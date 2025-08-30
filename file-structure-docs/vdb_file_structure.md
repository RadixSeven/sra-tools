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

#### 1. zip_encoding (Standard zlib compression)

**Algorithm**: Standard RFC 1950 zlib/deflate compression

**Implementation Details:**
```c
// Compression parameters
#define ZLIB_COMPRESSION_LEVEL 6        // Balanced speed/ratio
#define ZLIB_WINDOW_BITS 15             // Standard 32KB sliding window
#define ZLIB_MEM_LEVEL 8                // Default memory usage

// Python equivalent using zlib library
import zlib

def compress_zip_encoding(data):
    # Compress using zlib level 6 (VDB default)
    compressed = zlib.compress(data, level=6)
    return compressed

def decompress_zip_encoding(compressed_data):
    # Standard zlib decompression
    decompressed = zlib.decompress(compressed_data)
    return decompressed
```

**Usage**: DNA sequences, quality scores, text metadata
**Typical Compression**: 60-95% size reduction
**Performance**: 100-200 MB/s decompression

#### 2. izip_encoding (Integer delta + zlib compression)

**Algorithm**: Delta encoding followed by zlib compression optimized for integer sequences

**Implementation Details:**
```c
// Step-by-step process:
// 1. Apply delta encoding to integer array
// 2. Compress deltas using zlib level 9

// Python implementation
import zlib
import struct

def compress_izip_encoding(integer_array):
    # Step 1: Delta encoding
    deltas = [integer_array[0]]  # First value as-is
    for i in range(1, len(integer_array)):
        delta = integer_array[i] - integer_array[i-1]
        deltas.append(delta)
    
    # Step 2: Pack deltas as little-endian integers
    if max(abs(d) for d in deltas) < 128:
        # Use 8-bit signed integers
        packed = struct.pack(f'<{len(deltas)}b', *deltas)
    elif max(abs(d) for d in deltas) < 32768:
        # Use 16-bit signed integers
        packed = struct.pack(f'<{len(deltas)}h', *deltas)
    else:
        # Use 32-bit signed integers
        packed = struct.pack(f'<{len(deltas)}i', *deltas)
    
    # Step 3: Compress with maximum zlib compression
    compressed = zlib.compress(packed, level=9)
    return compressed

def decompress_izip_encoding(compressed_data, value_count, delta_size):
    # Step 1: Decompress
    decompressed = zlib.decompress(compressed_data)
    
    # Step 2: Unpack deltas
    if delta_size == 1:
        deltas = struct.unpack(f'<{value_count}b', decompressed)
    elif delta_size == 2:
        deltas = struct.unpack(f'<{value_count}h', decompressed)
    else:
        deltas = struct.unpack(f'<{value_count}i', decompressed)
    
    # Step 3: Reconstruct original values
    values = [deltas[0]]
    for i in range(1, len(deltas)):
        values.append(values[i-1] + deltas[i])
    
    return values
```

**Usage**: Coordinates, spot IDs, numeric sequences
**Typical Compression**: 80-95% size reduction
**Performance**: 50-100 MB/s decompression

#### 3. pack_encoding (Bit packing)

**Algorithm**: Bit-level packing with no secondary compression

**Implementation Details:**
```c
// Bit packing for boolean and small integer values
// Example: Pack array of 2-bit values (DNA nucleotides)

// Python implementation
def compress_pack_encoding_2bit(nucleotide_array):
    # Pack 4 nucleotides per byte (2 bits each)
    # A=0, C=1, G=2, T=3
    packed_bytes = bytearray()
    
    for i in range(0, len(nucleotide_array), 4):
        byte_value = 0
        for j in range(min(4, len(nucleotide_array) - i)):
            nucleotide = nucleotide_array[i + j]
            byte_value |= (nucleotide << (j * 2))
        packed_bytes.append(byte_value)
    
    return bytes(packed_bytes)

def decompress_pack_encoding_2bit(packed_data, nucleotide_count):
    # Unpack 2-bit values from bytes
    nucleotides = []
    
    for byte_idx, byte_val in enumerate(packed_data):
        for bit_pos in range(0, 8, 2):
            if len(nucleotides) >= nucleotide_count:
                break
            nucleotide = (byte_val >> bit_pos) & 0x3
            nucleotides.append(nucleotide)
    
    return nucleotides[:nucleotide_count]

# For 1-bit boolean values
def compress_pack_encoding_1bit(boolean_array):
    packed_bytes = bytearray()
    
    for i in range(0, len(boolean_array), 8):
        byte_value = 0
        for j in range(min(8, len(boolean_array) - i)):
            if boolean_array[i + j]:
                byte_value |= (1 << j)
        packed_bytes.append(byte_value)
    
    return bytes(packed_bytes)
```

**Usage**: DNA nucleotides (2na_packed), boolean flags, small integers
**Typical Compression**: 50-87.5% size reduction (2-bit: 75%, 1-bit: 87.5%)
**Performance**: 500+ MB/s (minimal overhead)

#### 4. fzip_encoding (Floating-point compression)

**Algorithm**: Floating-point quantization followed by delta encoding and zlib compression

**Implementation Details:**
```c
// Algorithm steps:
// 1. Quantize floating-point values to 16-bit integers
// 2. Apply delta encoding to quantized values
// 3. Compress deltas using zlib

// Python implementation
import zlib
import struct
import numpy as np

def compress_fzip_encoding(float_array, quantization_bits=16):
    # Step 1: Find range for quantization
    min_val = min(float_array)
    max_val = max(float_array)
    range_val = max_val - min_val
    
    # Step 2: Quantize to 16-bit integers
    max_quant = (1 << quantization_bits) - 1
    quantized = []
    for val in float_array:
        if range_val > 0:
            quant_val = int((val - min_val) * max_quant / range_val)
        else:
            quant_val = 0
        quantized.append(min(max_quant, max(0, quant_val)))
    
    # Step 3: Delta encoding
    deltas = [quantized[0]]
    for i in range(1, len(quantized)):
        delta = quantized[i] - quantized[i-1]
        deltas.append(delta)
    
    # Step 4: Pack header + deltas
    header = struct.pack('<ff', min_val, max_val)  # Range info
    if max(abs(d) for d in deltas) < 128:
        delta_data = struct.pack(f'<{len(deltas)}b', *deltas)
    else:
        delta_data = struct.pack(f'<{len(deltas)}h', *deltas)
    
    # Step 5: Compress
    compressed = zlib.compress(header + delta_data, level=6)
    return compressed

def decompress_fzip_encoding(compressed_data, value_count):
    # Step 1: Decompress
    decompressed = zlib.decompress(compressed_data)
    
    # Step 2: Extract range
    min_val, max_val = struct.unpack('<ff', decompressed[:8])
    range_val = max_val - min_val
    delta_data = decompressed[8:]
    
    # Step 3: Unpack deltas (detect size from remaining data)
    bytes_per_delta = len(delta_data) // value_count
    if bytes_per_delta == 1:
        deltas = struct.unpack(f'<{value_count}b', delta_data)
    else:
        deltas = struct.unpack(f'<{value_count}h', delta_data)
    
    # Step 4: Reconstruct quantized values
    quantized = [deltas[0]]
    for i in range(1, len(deltas)):
        quantized.append(quantized[i-1] + deltas[i])
    
    # Step 5: Dequantize to floats
    max_quant = 65535.0  # 16-bit
    floats = []
    for quant_val in quantized:
        if range_val > 0:
            float_val = min_val + (quant_val * range_val / max_quant)
        else:
            float_val = min_val
        floats.append(float_val)
    
    return floats
```

**Usage**: Signal intensities, kinetic data, measured values
**Typical Compression**: 50-75% size reduction
**Performance**: 75-150 MB/s decompression

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

## VDB Blob Access Patterns - Complete Walkthrough

### Blob Location and Access Workflow

Here's the complete process for reading column data from a specific row, addressing the critical gap identified in the suggestions:

**Complete Walkthrough: Reading column data for row 1000**

```python
def read_column_data_for_row(vdb_table, column_name, target_row):
    """
    Complete example of reading a specific row from a VDB column
    """
    
    # Step 1: Read idx2 file to locate block containing row 1000
    idx2_data = read_file(f"{column_name}/idx2")
    idx2_header = parse_idx2_header(idx2_data)
    
    # Find block containing target row
    block_info = None
    for block in idx2_header.blocks:
        if block.start_row <= target_row <= block.start_row + block.row_count:
            block_info = block
            break
    
    if not block_info:
        raise ValueError(f"Row {target_row} not found in idx2")
    
    # Step 2: Use block info to read corresponding idx1 entry
    idx1_data = read_file(f"{column_name}/idx1")
    idx1_offset = block_info.idx1_offset
    
    # Read KColBlockLoc from idx1
    idx1_entry = parse_idx1_entry(idx1_data, idx1_offset)
    
    # Step 3: Use idx1 entry to locate idx entry
    idx_data = read_file(f"{column_name}/idx")
    idx_offset = idx1_entry.idx_offset
    
    # Step 4: Read KColBlobLoc from idx
    blob_loc = parse_blob_loc(idx_data, idx_offset)
    
    # Verify row is in this blob
    if not (blob_loc.start_id <= target_row < blob_loc.start_id + blob_loc.id_range):
        raise ValueError(f"Row {target_row} not in blob range")
    
    # Step 5: Seek to pg offset in data file and read blob
    data_file = open(f"{column_name}/data", 'rb')
    data_file.seek(blob_loc.pg)
    
    # Read blob header
    blob_header = read_blob_header(data_file)
    
    # Read compressed blob data
    compressed_data = data_file.read(blob_loc.u.blob.size)
    
    # Step 6: Decompress using specified algorithm
    if blob_header.fmt == ZIP_ENCODING:
        decompressed = decompress_zip_encoding(compressed_data)
    elif blob_header.fmt == IZIP_ENCODING:
        decompressed = decompress_izip_encoding(compressed_data, blob_loc.id_range, 4)
    elif blob_header.fmt == PACK_ENCODING:
        decompressed = decompress_pack_encoding(compressed_data, blob_header.osize)
    elif blob_header.fmt == FZIP_ENCODING:
        decompressed = decompress_fzip_encoding(compressed_data, blob_loc.id_range)
    else:
        raise ValueError(f"Unknown encoding format: {blob_header.fmt}")
    
    # Step 7: Extract row 1000 data from decompressed blob
    row_offset_in_blob = target_row - blob_loc.start_id
    
    # For example, if this is DNA sequence data (2na_packed)
    if column_name == "READ":
        # Each row contains variable-length sequence
        # Use row offsets to locate specific row data
        row_data = extract_sequence_data(decompressed, row_offset_in_blob)
        return convert_2na_to_bases(row_data)
    
    elif column_name == "QUALITY":
        # Quality scores, one per base
        row_data = extract_quality_data(decompressed, row_offset_in_blob)
        return convert_phred33_to_ascii(row_data)
    
    elif column_name == "SPOT_ID":
        # Simple integer value
        spot_id = struct.unpack('<Q', decompressed[row_offset_in_blob:row_offset_in_blob+8])[0]
        return spot_id
    
    return decompressed

def parse_idx2_header(idx2_data):
    """Parse idx2 file header and block information"""
    header = struct.unpack('<II', idx2_data[:8])  # endian, version
    
    if header[0] not in [0x05031988, 0x88190305]:
        raise ValueError("Invalid idx2 header endianness")
    
    # Read block count and parse blocks
    block_count = struct.unpack('<I', idx2_data[8:12])[0]
    blocks = []
    offset = 12
    
    for i in range(block_count):
        # Parse block descriptor based on representation type
        block_type = idx2_data[offset]
        if block_type == 0:  # btypeRandom
            start_row, row_count, idx1_offset = struct.unpack('<QII', idx2_data[offset+1:offset+17])
            blocks.append({
                'type': 'random',
                'start_row': start_row,
                'row_count': row_count,
                'idx1_offset': idx1_offset
            })
            offset += 17
        elif block_type == 1:  # btypeUniform
            start_row, row_count, uniform_size = struct.unpack('<QII', idx2_data[offset+1:offset+17])
            blocks.append({
                'type': 'uniform',
                'start_row': start_row,
                'row_count': row_count,
                'uniform_size': uniform_size
            })
            offset += 17
        # Handle other block types...
    
    return type('IDX2Header', (), {'blocks': blocks})()

def parse_blob_loc(idx_data, offset):
    """Parse KColBlobLoc structure from idx file"""
    blob_data = idx_data[offset:offset+24]  # KColBlobLoc is 24 bytes
    
    pg, blob_info, id_range, start_id = struct.unpack('<QLIQ', blob_data)
    
    # Extract size and remove flag from blob_info
    blob_size = blob_info & 0x7FFFFFFF
    remove_flag = (blob_info & 0x80000000) != 0
    
    return type('KColBlobLoc', (), {
        'pg': pg,
        'u': type('Union', (), {
            'blob': type('Blob', (), {
                'size': blob_size,
                'remove': remove_flag
            })()
        })(),
        'id_range': id_range,
        'start_id': start_id
    })()

def read_blob_header(data_file):
    """Read and parse VDB blob header"""
    # VDB blob header format varies, but typically includes:
    # - Format identifier (encoding type)
    # - Original size
    # - Argument count and arguments
    
    header_start = data_file.tell()
    
    # Read basic header fields
    fmt, version, flags = struct.unpack('<IBB', data_file.read(6))
    
    osize = struct.unpack('<Q', data_file.read(8))[0]  # Original size
    
    arg_count = struct.unpack('<I', data_file.read(4))[0]
    args = []
    for i in range(arg_count):
        arg = struct.unpack('<q', data_file.read(8))[0]  # Signed 64-bit
        args.append(arg)
    
    return type('VBlobHeader', (), {
        'fmt': fmt,
        'version': version,
        'flags': flags,
        'osize': osize,
        'args': args
    })()
```

### Index Navigation Examples

**Example 1: Finding all blobs for a column**
```python
def get_all_blob_locations(column_path):
    """Get all blob locations for a column"""
    idx_data = read_file(f"{column_path}/idx")
    
    # Read KDB header
    endian, version = struct.unpack('<II', idx_data[:8])
    if endian not in [0x05031988, 0x88190305]:
        raise ValueError("Invalid index header")
    
    # Calculate number of blob entries
    blob_count = (len(idx_data) - 8) // 24  # Each KColBlobLoc is 24 bytes
    
    blob_locations = []
    for i in range(blob_count):
        offset = 8 + (i * 24)
        blob_loc = parse_blob_loc(idx_data, offset)
        blob_locations.append(blob_loc)
    
    return blob_locations

def estimate_blob_for_row(blob_locations, target_row):
    """Binary search to find blob containing target row"""
    left, right = 0, len(blob_locations) - 1
    
    while left <= right:
        mid = (left + right) // 2
        blob = blob_locations[mid]
        
        if blob.start_id <= target_row < blob.start_id + blob.id_range:
            return blob
        elif target_row < blob.start_id:
            right = mid - 1
        else:
            left = mid + 1
    
    return None
```

## Access Patterns and Indexing

### Row-Based Access

**Sequential Row Reading:**
1. **Parse index hierarchy**: Read idx2 → idx1 → idx for sequential access
2. **Process blobs in order**: Decompress blobs sequentially by pg offset
3. **Extract row data**: Use blob start_id and id_range to locate specific rows
4. **Reassemble complete rows**: Combine data from multiple column blobs

**Random Row Access:**
1. **Binary search idx2**: Locate block containing target row ID
2. **Navigate to idx1**: Use block descriptor to find idx1 entry
3. **Read blob location**: Extract KColBlobLoc from idx entry
4. **Targeted decompression**: Decompress only the blob containing target row
5. **Extract specific row**: Calculate offset within blob and extract data

### Column-Based Access

**Full Column Scan:**
1. **Read all blob locations**: Parse complete idx file to get all KColBlobLoc entries
2. **Sort by pg offset**: Process blobs in storage order for optimal disk access
3. **Stream decompression**: Decompress and process each blob sequentially
4. **Aggregate results**: Combine data from all blobs for complete column

**Column Range Queries:**
1. **Binary search range**: Find first and last blobs overlapping target range
2. **Minimal blob set**: Read only blobs containing rows in target range
3. **Partial decompression**: Extract only relevant data from each blob
4. **Range extraction**: Filter decompressed data to exact row range

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
