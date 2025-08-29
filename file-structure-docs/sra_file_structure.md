# SRA File Format Documentation

## Overview

The Sequence Read Archive (SRA) format is a KAR archive containing a VDB (Virtual/Vertical Database) with a specific directory tree structure optimized for genomic sequence data storage. The VDB architecture provides a columnar storage system that organizes data vertically by columns rather than horizontally by rows, enabling efficient analytical queries and superior compression ratios for biological data. The format supports two variants: **SRA Normalized** (full format) and **SRA Lite** (compressed format with simplified quality scores).

**Related Documentation:**
- [KAR File Format Specification](kar_file_structure.md) - The underlying archive format
- [VDB File Format Specification](vdb_file_structure.md) - The database storage system

## Layered Architecture

The SRA format is built on a three-layer architecture where each layer provides specific functionality:

```
┌─────────────────────────────────────┐
│    Layer 3: SRA Schema              │
│    - Biological data types          │
│    - Platform optimizations         │
│    - SRA Norm vs Lite variants      │
└─────────────────────────────────────┘
┌─────────────────────────────────────┐
│    Layer 2: VDB Database System     │
│    - Columnar storage               │
│    - Compression & indexing         │
│    - Schema-driven organization     │
└─────────────────────────────────────┘
┌─────────────────────────────────────┐
│    Layer 1: KAR Archive Container   │
│    - File packaging                 │
│    - Magic signature validation     │
│    - Directory structure            │
└─────────────────────────────────────┘
```

### Layer Integration

1. **KAR Archive Layer** provides file packaging and validation
2. **VDB Database Layer** implements columnar storage and indexing
3. **SRA Schema Layer** defines biological data semantics and optimizations

Each layer can be understood and implemented independently, making the format modular and extensible.

## Layer 1: KAR Archive Container

### Magic Signature and Identification

Every SRA file begins with the 8-byte magic signature: **`NCBI.sra`**
- Bytes 0-3: `"NCBI"` (0x4E434249)
- Bytes 4-7: `".sra"` (0x2E737261)

This signature identifies the file as an NCBI SRA format file and triggers KAR archive processing.

### Archive Structure

SRA files follow the standard KAR archive format:

```
[KAR Header: "NCBI.sra" + metadata] → [Table of Contents] → [File Data]
```

**Key KAR Features for SRA:**
- **4-byte aligned storage** for optimal access performance
- **Binary search tree TOC** for O(log n) file lookup
- **Sorted file storage** (by size) for access optimization
- **Chunked file support** for large data files

### Archive Contents

The KAR archive contains the complete VDB database as a directory tree:

```
Archive Contents:
├── md/                     # Database metadata
├── tbl/SEQUENCE/          # Primary sequence table
├── col/                   # Global column definitions
└── idx/                   # Database-level indexes
```

For complete KAR format details, see [KAR File Format Specification](kar_file_structure.md).

## Layer 2: VDB Database System

### Directory Tree Structure

The VDB database within the KAR archive follows a standardized hierarchy:

```
VDB Database Structure:
├── md/                     # Metadata directory
│   ├── cur                # Current version pointer
│   ├── vers               # Version history
│   └── root               # Root metadata (schema, timestamps)
├── tbl/                   # Table directory
│   └── SEQUENCE/          # Primary sequence table
│       ├── md/            # Table metadata
│       ├── col/           # Table columns
│       │   ├── READ/      # DNA sequence data
│       │   ├── QUALITY/   # Quality scores
│       │   ├── SPOT_ID/   # Spot identifiers
│       │   └── [others]/  # Additional columns
│       └── idx/           # Table indexes
├── col/                   # Global column definitions
│   └── [COLUMN_NAME]/     # Shared column specs
└── idx/                   # Database-level indexes
```

### Columnar Storage Implementation

**Column Organization:**
- Each column stored in separate directory (`col/[COLUMN_NAME]/`)
- Data compressed and stored in blobs (`data` file)
- Multi-level indexing (`idx`, `idx1`, `idx2`) for efficient access
- Column metadata tracks compression and statistics

**Blob Structure:**
```c
typedef struct KColBlobLoc {
    uint64_t pg;                // File offset to blob data
    uint32_t size : 31;         // Blob size in bytes
    uint32_t remove : 1;        // Deletion flag
    uint32_t id_range;          // Number of rows in blob
    int64_t start_id;           // Starting row ID
} KColBlobLoc;
```

### Compression and Encoding

**Algorithm Selection by Data Type:**
- **DNA Sequences**: `zip_encoding` (zlib compression)
- **Quality Scores**: `zip_encoding` or `pack_encoding`
- **Coordinates**: `izip_encoding` (integer delta + zlib)
- **Boolean Flags**: `pack_encoding` (bit packing)
- **Signal Data**: `fzip_encoding` (floating-point optimized)

**Compression Performance:**
- DNA sequences: 95-98% size reduction
- Quality scores: 60-80% reduction
- Coordinates: 80-95% reduction
- Decompression: 50-500 MB/s depending on algorithm

For complete VDB format details, see [VDB File Format Specification](vdb_file_structure.md).

## Layer 3: SRA Schema Implementation

### Biological Data Types

**Nucleotide Encodings:**
- **2na_packed**: 2-bit encoding (A=0, C=1, G=2, T=3), 4 bases per byte
- **4na_packed**: 4-bit encoding including ambiguous bases (N=15), 2 bases per byte
- **x2na_bin**: Extended format for colorspace data (SOLiD platform)

**Quality Score Formats:**
- **phred_33**: Standard Illumina (ASCII 33-126, quality = char - 33)
- **phred_64**: Legacy Illumina (ASCII 64-126, quality = char - 64)
- **Simplified**: SRA Lite format (30 = pass, 3 = fail)

### SRA-Specific Tables and Columns

**Primary Table: SEQUENCE**
- `READ`: DNA/RNA sequence data (2na_packed or 4na_packed)
- `QUALITY`: Base quality scores (format varies by variant)
- `SPOT_ID`: Unique spot/cluster identifier (uint64_t)
- `READ_TYPE`: Read classification (technical/biological)
- `READ_FILTER`: Pass/fail quality determination
- `READ_START`: Starting positions of reads within spots
- `READ_LEN`: Lengths of individual reads

**Supporting Tables:**
- **STATS**: Run-level statistics (base_count, spot_count, platform info)
- **SPOTCOORD**: Spatial coordinates (X_COORD, Y_COORD)
- **SPOTNAME**: External identifiers (NAME_FMT, SPOT_NAME)

### Platform Optimizations

**Supported Platforms** (19 total, platform IDs 0-18):
- **Illumina** (ID 2): Optimized quality compression, coordinate handling
- **454** (ID 1): Variable-length reads, signal intensity support
- **PacBio SMRT** (ID 6): Long-read optimizations, kinetic data
- **Oxford Nanopore** (ID 9): Ultra-long reads, signal compression
- **Ion Torrent** (ID 7): Flow-based quality encoding
- **ABI SOLiD** (ID 3): Colorspace encoding support
- [Additional platforms 4, 5, 8, 10-18]

### Format Variants

**SRA Normalized (.sra):**
- Full per-base quality scores as produced by sequencers
- Platform-specific quality encodings
- Complete data fidelity
- Larger file sizes

**SRA Lite (.sralite):**
- Simplified quality scores (30 = pass, 3 = fail)
- ~60-80% size reduction
- Full tool compatibility
- Faster transfer and processing

## Integration and Implementation

### Format Layer Integration

The three layers work together in a specific workflow:

```
1. KAR Archive Processing:
   ├── Validate magic signature: NCBI.sra
   ├── Check byte order: 0x05031988 or 0x88190305  
   ├── Parse Table of Contents (TOC)
   └── Extract VDB database files

2. VDB Database Access:
   ├── Read metadata from md/ directory
   ├── Open tables from tbl/ structure
   ├── Access columns via multi-level indexes
   └── Decompress blobs using schema algorithms

3. SRA Schema Interpretation:
   ├── Apply biological data type conversions
   ├── Handle platform-specific encodings
   ├── Manage SRA Norm vs Lite differences
   └── Provide genomic data access
```

### Access Patterns

**Sequential Access (Analytical Queries):**
1. Read all blobs for target columns
2. Decompress in storage order
3. Process data in chunks
4. Optimal for genome-wide analysis

**Random Access (Targeted Queries):**
1. Use idx2 to locate specific rows
2. Read minimal blob set
3. Decompress only necessary data
4. Extract specific elements

### Implementation Guide

**For Read-Only SRA Access:**

1. **Implement KAR Parser**
   - Validate magic signature and headers
   - Parse binary TOC using BSTree navigation
   - Extract files with 4-byte alignment handling

2. **Implement VDB Reader**  
   - Parse directory structure (md/, tbl/, col/)
   - Handle multi-level indexes (idx, idx1, idx2)
   - Implement blob decompression (zip, izip, pack, fzip)

3. **Implement SRA Schema**
   - Support 19 sequencing platforms
   - Handle nucleotide encodings (2na, 4na)
   - Manage quality score formats
   - Distinguish SRA Norm vs Lite variants

**Required Components:**
- Binary data parsing (little-endian integers, structures)
- Compression libraries (zlib, custom algorithms)
- Index traversal (binary search, range queries)
- Type conversion (nucleotide, quality, coordinate systems)

### Error Handling and Validation

**Validation Hierarchy:**
1. **KAR Level**: Magic signature, TOC consistency, file bounds
2. **VDB Level**: Index integrity, blob checksums, schema compatibility
3. **SRA Level**: Platform constraints, data type validation, biological consistency

**Recovery Strategies:**
- Graceful degradation for partial corruption
- Skip corrupted blobs when possible
- Use redundant indexes for cross-validation
- Detailed error reporting with layer-specific context

### Performance Characteristics

**Storage Efficiency:**
- **SRA Normalized**: Full fidelity, optimal compression for platform
- **SRA Lite**: 60-80% size reduction with quality simplification
- **Columnar benefits**: Superior compression through homogeneous data

**Access Performance:**
- **Index lookups**: O(log n) via binary search trees
- **Blob caching**: LRU cache (default 128MB)
- **Sequential reads**: 100-500 MB/s depending on compression
- **Random access**: Sub-linear with proper indexing

## Conversion Between Variants

### SRA Norm to SRA Lite (Delite Process)

**Conversion Steps:**
1. **Schema Translation**: Update to Lite-compatible schema version
2. **Quality Processing**: 
   - Preserve original as ORIGINAL_QUALITY
   - Generate simplified scores using READ_FILTER
   - Quality = 30 for pass reads, 3 for reject reads
3. **Column Optimization**: Remove platform-specific columns
4. **Validation**: Ensure compatibility and integrity

**Platform Restrictions:**
- Colorspace platforms (ABI SOLiD, ID 3) cannot be converted
- TRACE type archives not supported
- Validation occurs before conversion begins

**Conversion Tools:**
- `sra_delite.sh`: Primary conversion script
- `vdb-dump`: Data extraction and validation
- `sra-stat`: Statistics and verification

## Version and Compatibility

**Current Versions:**
- **SRA Format Version**: 1 (`FS_SRA_CUR_VERSION = 1`)
- **Schema Version**: 1.1.1 (for delite format)
- **Toolkit Compatibility**: SRA Lite requires 2.11.2+

**Compatibility Features:**
- Forward/backward compatibility maintained
- Automatic format detection
- Schema versioning system
- Graceful handling of unknown versions

## Complete Example: Minimal SRA Reader

```c
// Pseudo-code for minimal SRA file reader
int read_sra_file(const char* filename) {
    // Layer 1: KAR Archive
    if (!validate_kar_magic(filename)) return -1;
    KARArchive* kar = open_kar_archive(filename);
    
    // Layer 2: VDB Database  
    VDBDatabase* vdb = open_vdb_from_kar(kar, "/");
    VDBTable* seq_table = open_vdb_table(vdb, "SEQUENCE");
    
    // Layer 3: SRA Schema
    VDBColumn* read_col = open_column(seq_table, "READ");
    VDBColumn* qual_col = open_column(seq_table, "QUALITY");
    
    // Access data
    for (int64_t row = 1; row <= get_row_count(seq_table); row++) {
        read_sequence_data(read_col, row);
        read_quality_data(qual_col, row);
    }
    
    return 0;
}
```

This layered architecture documentation provides a clear path for understanding and implementing SRA format support, from the foundational KAR archive through the VDB database system to the biological data semantics of the SRA schema.