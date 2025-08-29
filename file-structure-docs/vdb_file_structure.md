# VDB (Virtual Database) File Format Specification

## Overview

The Virtual Database (VDB) format is NCBI's columnar database system designed specifically for storing and accessing biological sequence data. VDB databases are stored within KAR archives and provide a structured, indexed, and compressed storage system optimized for genomic data access patterns. The format supports schema-driven data organization, multiple compression algorithms, and efficient random access to large datasets.

## Database Directory Structure

### Hierarchical Organization

VDB databases follow a standardized directory structure within their container (typically a KAR file):

```
Database Root/
├── md/                     # Metadata directory
│   ├── cur/               # Current metadata version
│   ├── vers/              # Version history
│   └── root               # Root metadata file
├── tbl/                   # Table directory
│   └── [TABLE_NAME]/      # Individual table directories
│       ├── md/            # Table-specific metadata
│       ├── col/           # Table column directory
│       │   └── [COLUMN]/  # Individual column directories
│       │       ├── data   # Column data file
│       │       ├── idx    # Primary index
│       │       ├── idx1   # Secondary index (optional)
│       │       └── idx2   # Tertiary index (optional)
│       └── idx/           # Table-level indexes
├── col/                   # Global column directory
│   └── [COLUMN_NAME]/     # Global column definitions
│       ├── data           # Column data blob
│       ├── idx            # Column index
│       └── vers           # Column version info
└── idx/                   # Database-level indexes
    └── [INDEX_NAME]       # Named index files
```

### Directory Naming Conventions

- **Metadata files**: Plain names (`root`, `cur`, `vers`)
- **Table names**: User-defined, typically uppercase (e.g., `SEQUENCE`, `QUALITY`)
- **Column names**: Schema-defined, often uppercase (e.g., `READ`, `SPOT_ID`)
- **Index files**: Numbered suffixes (`idx`, `idx1`, `idx2`)

## Column Storage Format

### Data Blob Structure

Column data is stored in compressed blobs with location descriptors and headers:

```c
typedef struct KColBlobLoc {
    uint64_t pg;                // Data page ID (file offset)
    union {
        struct {
            uint32_t size : 31; // Blob size in bytes
            uint32_t remove : 1; // Removal flag for journaling
        } blob;
        uint32_t gen;           // General 32-bit field
    } u;
    uint32_t id_range;          // Number of rows covered by this blob
    int64_t start_id;           // First row ID in this blob
} KColBlobLoc;

typedef struct VBlobHeaderData {
    int64_t *args;              // Argument array for transformations
    uint8_t *ops;               // Operation codes for transformations
    atomic32_t refcount;        // Reference count
    uint32_t op_count;          // Number of operations
    uint32_t arg_count;         // Number of arguments
    uint64_t osize;             // Original size before compression
    uint8_t flags;              // Control flags
    uint8_t version;            // Header format version
    uint32_t fmt;               // Format identifier
    bool read_only;             // Read-only flag
    bool args_alloc;            // Arguments allocated flag
    bool ops_alloc;             // Operations allocated flag
} VBlobHeaderData;
```

### Compression Algorithms

**Primary Encoding Types:**

1. **zip_encoding**: General-purpose compression using zlib/deflate
   - Used for: Text data, variable-length data
   - Compression ratio: High for repetitive data
   - Access pattern: Sequential decompression required

2. **izip_encoding**: Integer-specific compression
   - Used for: Numeric data, coordinates, IDs
   - Compression: Delta encoding + zlib
   - Access pattern: Optimized for numeric sequences

3. **pack_encoding**: Bit-packed encoding for small values
   - Used for: Boolean data, small integers, flags
   - Compression: Bit packing without additional compression
   - Access pattern: Direct bit-level access

4. **fzip_encoding**: Floating-point specific compression
   - Used for: Signal data, quality scores, measurements
   - Compression: Floating-point delta + quantization
   - Access pattern: Optimized for scientific data

### Index Structure

#### Primary Index (idx)

The primary index contains KColBlobLoc structures for direct blob access:

```c
typedef struct KDBHdr {
    uint32_t endian;         // Byte order: 0x05031988 (normal) or 0x88190305 (reversed)
    uint32_t version;        // Index format version
} KDBHdr;
```

**Index Organization:**
- **Level 0 Index (idx)**: Direct array of `KColBlobLoc` structures
- Each entry maps a row range to its blob location
- Sorted by `start_id` for binary search access
- Contains complete blob location information

#### Secondary Indexes (idx1, idx2)

**Level 1 Index (idx1):** Block-level organization
- Contains `KColBlockLoc` structures for grouped blob access
- Header includes `KDBHdr` plus column-specific metadata:
  - `data_eof`: End of data file
  - `idx2_eof`: End of level-2 index
  - `num_blocks`: Number of index blocks
  - `page_size`: Data page size (typically 1 byte for append mode)
  - `checksum`: Checksum type (none, CRC32, MD5)

**Level 2 Index (idx2):** Fine-grained access
- Contains compressed representations of multiple `KColBlobLoc` entries
- Four representation types:
  - **Type 0 (btypeRandom)**: Full specification for random access
  - **Type 1 (btypeUniform)**: Uniform sizes with different positions
  - **Type 2 (btypeMagnitude)**: Predictable sequence with varying sizes
  - **Type 3 (btypePredictable)**: Uniform contiguous sequence

## Data Type System

### Nucleotide Encodings

**2na_packed**: 2-bit nucleotide encoding
```c
typedef enum {
    NA_2na_A = 0,  // Adenine
    NA_2na_C = 1,  // Cytosine  
    NA_2na_G = 2,  // Guanine
    NA_2na_T = 3   // Thymine
} INSDC_2na_bin;
```
- Storage: 4 bases per byte
- Usage: Standard DNA sequences
- Compression: Excellent for uniform sequences

**4na_packed**: 4-bit nucleotide encoding
```c
typedef enum {
    NA_4na_A = 1,   // Adenine
    NA_4na_C = 2,   // Cytosine
    NA_4na_G = 4,   // Guanine
    NA_4na_T = 8,   // Thymine
    NA_4na_N = 15   // Ambiguous (any)
} INSDC_4na_bin;
```
- Storage: 2 bases per byte
- Usage: Sequences with ambiguous bases
- Compression: Good for mixed-quality sequences

**x2na_bin**: Extended 2-nucleotide format
- Supports colorspace encoding
- Used for SOLiD platform data
- Includes color-to-base conversion information

### Quality Score Types

**phred_33**: Standard Illumina quality encoding
```c
typedef struct {
    char quality_char;       // ASCII character
    uint8_t phred_value;    // quality_char - 33
} INSDC_quality_phred_33;
```
- Range: 0-93 (ASCII 33-126)
- Usage: Modern Illumina platforms
- Storage: 1 byte per base

**phred_64**: Legacy Illumina quality encoding  
```c
typedef struct {
    char quality_char;       // ASCII character  
    uint8_t phred_value;    // quality_char - 64
} INSDC_quality_phred_64;
```
- Range: 0-62 (ASCII 64-126)
- Usage: Older Illumina platforms
- Storage: 1 byte per base

### Platform-Specific Types

**Signal Intensities**: Floating-point arrays for raw signal data
- **fsamp4**: 4-channel floating-point samples (A, C, G, T)
- **swapped_fsamp4**: Channel-swapped version for called base optimization
- **rotated_fsamp4**: Rotated channels for colorspace data

**Coordinate Types**:
- **coord:zero**: Zero-based coordinates (uint32_t)
- **coord:one**: One-based coordinates (uint32_t) 
- **pos16**: 16-bit position values for memory optimization

## Schema System

### Schema Storage

Schema information is stored in the metadata directory as binary data:

```c
typedef struct SchemaData {
    uint32_t schema_version; // Schema format version
    uint32_t table_count;    // Number of tables defined
    TableDef tables[];       // Array of table definitions
} SchemaData;

typedef struct TableDef {
    uint16_t name_len;       // Length of table name
    char name[name_len];     // Table name string
    uint32_t version_maj;    // Major version number
    uint32_t version_min;    // Minor version number  
    uint32_t column_count;   // Number of columns
    ColumnDef columns[];     // Array of column definitions
} TableDef;

typedef struct ColumnDef {
    uint16_t name_len;       // Length of column name
    char name[name_len];     // Column name string
    uint16_t type_len;       // Length of type specification
    char type[type_len];     // VDB type specification
    uint8_t encoding;        // Physical encoding method
    uint8_t flags;           // Column attribute flags
} ColumnDef;
```

### Schema Inheritance

VDB supports table inheritance for schema reuse:

```
Base Table: INSDC:tbl:sequence
├── Inherited by: INSDC:SRA:tbl:sra  
│   ├── Inherited by: NCBI:SRA:tbl:sra
│   └── Mixed with: NCBI:SRA:tbl:spotdesc
└── Mixed with: INSDC:SRA:tbl:spotname
```

**Inheritance Rules:**
1. Child tables inherit all columns from parent tables
2. Child tables can override parent column definitions  
3. Multiple inheritance supported through mixins
4. Schema versioning maintains compatibility

### Type Specification Language

VDB uses a declarative type system:

```
// Basic types
column INSDC:dna:text READ;
column INSDC:quality:phred QUALITY;
column U64 SPOT_ID;

// Parameterized types  
column < INSDC:coord:zero > READ_START;
column < ascii > zip_encoding LABEL;

// Complex types
readonly column INSDC:SRA:platform_id PLATFORM = out_platform;
```

## Metadata System

### Database Metadata

**Root Metadata** (`md/root`):
```c
typedef struct DatabaseMeta {
    uint32_t format_version;     // VDB format version
    uint32_t schema_version;     // Schema format version
    uint64_t create_timestamp;   // Database creation time
    uint64_t modify_timestamp;   // Last modification time
    uint32_t table_count;        // Number of tables
    char software_version[];     // Creating software version
    char schema_text[];          // Embedded schema text
} DatabaseMeta;
```

**Version Metadata** (`md/vers`):
- Tracks schema evolution over time
- Maintains compatibility information
- Records software versions used

**Current Metadata** (`md/cur`):
- Points to current active version
- Used for version resolution
- Updated during schema migrations

### Table Metadata

Each table maintains its own metadata in `tbl/[TABLE]/md/`:

```c
typedef struct TableMeta {
    uint32_t table_version;      // Table schema version
    uint64_t row_count;          // Total number of rows
    uint64_t create_time;        // Table creation timestamp
    uint32_t column_count;       // Number of columns
    ColumnMeta columns[];        // Column metadata array
} TableMeta;

typedef struct ColumnMeta {
    char name[64];               // Column name
    uint32_t data_type;          // VDB data type ID
    uint32_t encoding_type;      // Physical encoding
    uint64_t element_count;      // Total elements stored
    uint64_t blob_count;         // Number of storage blobs
    uint32_t compression_ratio;  // Achieved compression ratio
} ColumnMeta;
```

## Access Patterns and Indexing

### Row-Based Access

**Sequential Row Reading:**
1. Read primary index to locate row range
2. Decompress relevant blobs
3. Extract row data from multiple columns
4. Reassemble complete rows

**Random Row Access:**
1. Use idx2 (element index) to locate exact row
2. Read minimal blob set containing target row
3. Decompress only necessary data
4. Extract specific row elements

### Column-Based Access

**Full Column Scan:**
1. Sequential read through all column blobs
2. Decompress blobs in storage order
3. Process data in compressed chunks
4. Optimal for analytical queries

**Column Range Queries:**
1. Use primary index to identify blob range
2. Read only blobs containing target range
3. Decompress minimal data set
4. Extract range from decompressed data

### Blob Management

**Blob Size Optimization:**
- Target blob size: 64KB - 1MB compressed
- Balance between compression ratio and access granularity
- Larger blobs: better compression, slower random access
- Smaller blobs: worse compression, faster random access

**Blob Caching:**
- LRU cache for recently accessed blobs
- Configurable cache size (default: 128MB)
- Cache warming for predictable access patterns
- Write-through caching for modified data

## Compression Implementation

### Algorithm Selection

**Data Type Mapping:**
```c
static const encoding_map[] = {
    {INSDC_dna_text,     zip_encoding},     // DNA sequences
    {INSDC_coord_zero,   izip_encoding},    // Coordinates  
    {INSDC_quality_phred, zip_encoding},    // Quality scores
    {U8,                 pack_encoding},    // Small integers
    {bool,               pack_encoding},    // Boolean flags
    {F32,                fzip_encoding},    // Floating point
    {ascii,              zip_encoding}      // Text data
};
```

**Compression Parameters:**
- **zip_encoding**: zlib level 6 (balanced speed/ratio)
- **izip_encoding**: Delta encoding + zlib level 9
- **pack_encoding**: Bit packing, no secondary compression
- **fzip_encoding**: 16-bit quantization + zlib level 6

### Performance Characteristics

**Compression Ratios (typical):**
- DNA sequences (2na_packed + zip): 95-98% reduction
- Quality scores (phred + zip): 60-80% reduction  
- Coordinates (izip): 80-95% reduction
- Signal data (fzip): 50-75% reduction

**Decompression Speed:**
- zip_encoding: ~100-200 MB/s
- izip_encoding: ~50-100 MB/s
- pack_encoding: ~500+ MB/s (minimal overhead)
- fzip_encoding: ~75-150 MB/s

## Error Handling and Validation

### Data Integrity

**Checksum Types:**
- **CRC32**: Fast checksum for blob headers
- **MD5**: Strong checksum for large data blocks  
- **Fletcher**: Fast checksum for indexes

**Validation Levels:**
1. **Header validation**: Magic numbers, versions, sizes
2. **Index validation**: Offset bounds, entry consistency
3. **Data validation**: Checksums, decompression success
4. **Schema validation**: Type compatibility, constraint checks

### Error Recovery

**Corruption Detection:**
- Checksum mismatches trigger error handling
- Invalid offsets detected during index traversal
- Decompression failures indicate data corruption

**Recovery Strategies:**
- Skip corrupted blobs when possible
- Use redundant indexes for cross-validation
- Graceful degradation for partial data loss
- Detailed error reporting with specific failure locations

### Consistency Checks

**Database Consistency:**
- Row counts match across all columns in table
- Index entries point to valid data locations
- Schema definitions match actual column data
- Version information consistent across metadata

## Implementation Constants

### Byte Order and Format Constants

**Byte Order Constants (from [`libs/kdb/kdbfmt.h`](https://github.com/ncbi/ncbi-vdb/blob/master/libs/kdb/kdbfmt.h)):**
```c
#define eByteOrderTag     0x05031988    // Normal byte order
#define eByteOrderReverse 0x88190305    // Reversed byte order
```

**Checksum Types (from [`interfaces/kdb/column.h`](https://github.com/ncbi/ncbi-vdb/blob/master/interfaces/kdb/column.h)):**
```c
typedef uint8_t KChecksum;
enum {
    kcsNone,     // No checksum
    kcsCRC32,    // CRC32 checksum
    kcsMD5       // MD5 checksum
};
```

### Size Limits

```c
#define VDB_MAX_COLUMNS_PER_TABLE    1024
#define VDB_MAX_TABLES_PER_DATABASE  256  
#define VDB_MAX_COLUMN_NAME_LEN      64
#define VDB_MAX_TABLE_NAME_LEN       64
#define VDB_MAX_BLOB_SIZE           (16*1024*1024)  // 16MB
#define VDB_MIN_BLOB_SIZE           (4*1024)        // 4KB
#define VDB_MAX_INDEX_ENTRIES       (1024*1024)     // 1M entries
#define VDB_MAX_ROW_ID              INT64_MAX       // Maximum row ID (signed)
#define VDB_DEFAULT_PAGE_SIZE       1               // 1 byte (append mode)
```

### Alignment Requirements

```c
#define VDB_BLOB_ALIGNMENT          8    // 8-byte alignment
#define VDB_INDEX_ALIGNMENT         4    // 4-byte alignment
#define VDB_HEADER_ALIGNMENT        8    // 8-byte alignment
```

### Default Values

```c
#define VDB_DEFAULT_BLOB_SIZE       (256*1024)     // 256KB
#define VDB_DEFAULT_CACHE_SIZE      (128*1024*1024) // 128MB  
#define VDB_DEFAULT_COMPRESSION     6               // zlib level 6
#define VDB_INDEX_VERSION           1               // Current index version
```

## Platform-Specific Optimizations

### Sequencing Platform Support

**Platform Constants** (from `insdc/sra.vschema:405-423`):
```c
const INSDC:SRA:platform_id SRA_PLATFORM_ILLUMINA = 2;
const INSDC:SRA:platform_id SRA_PLATFORM_454 = 1;
const INSDC:SRA:platform_id SRA_PLATFORM_PACBIO_SMRT = 6;
const INSDC:SRA:platform_id SRA_PLATFORM_OXFORD_NANOPORE = 9;
// ... (19 total platforms supported)
```

**Platform-Specific Optimizations:**
- **Illumina**: Optimized quality score compression, coordinate handling
- **454**: Variable-length read support, flow signal storage
- **PacBio**: Long-read optimizations, kinetic data support  
- **Nanopore**: Ultra-long read handling, event data compression
- **Ion Torrent**: Flow-based quality encoding, homopolymer handling

### Memory Management

**Cache Architecture:**
```c
typedef struct VDBCache {
    uint32_t max_size;           // Maximum cache size
    uint32_t current_size;       // Current cache usage
    uint32_t blob_count;         // Number of cached blobs
    CacheEntry *lru_head;        // LRU list head
    CacheEntry *lru_tail;        // LRU list tail
    HashTable *blob_hash;        // Hash table for fast lookup
} VDBCache;
```

**Memory Optimization:**
- Lazy loading of column data
- Reference counting for shared blobs
- Memory-mapped file access where supported
- Automatic garbage collection of unused blobs

## Integration with SRA Format

### SRA-VDB Relationship

SRA files are VDB databases packaged in KAR archives with specific schema constraints:

```
SRA File Structure:
└── KAR Archive (NCBI.sra header)
    └── VDB Database
        ├── Schema: NCBI:SRA:tbl:sra
        ├── Tables: SEQUENCE (primary)
        ├── Columns: READ, QUALITY, SPOT_ID, etc.
        └── Indexes: Optimized for sequence access
```

**Key SRA Tables:**
- **SEQUENCE**: Primary sequence data table
- **STATS**: Run-level statistics and metadata
- **SPOTCOORD**: Spatial coordinates for cluster-based platforms
- **SPOTNAME**: External spot identifiers and naming

**SRA-Specific Columns:**
- **READ**: DNA/RNA sequence data (2na_packed or 4na_packed)
- **QUALITY**: Base quality scores (phred_33 or phred_64)
- **SPOT_ID**: Unique spot/cluster identifier (uint64_t)
- **READ_TYPE**: Read classification (technical/biological)
- **READ_FILTER**: Pass/fail quality determination
- **PLATFORM**: Sequencing platform identifier

This specification provides complete technical details for implementing VDB-compatible software, including all data structures, algorithms, compression methods, and integration patterns used in the NCBI SRA ecosystem.
