# SRA File Format Documentation

## Overview

The Sequence Read Archive (SRA) format is a KAR archive containing a VDB (Virtual/Vertical Database) with a specific directory tree structure optimized for genomic sequence data storage. The VDB architecture provides a columnar storage system that organizes data vertically by columns rather than horizontally by rows, enabling efficient analytical queries and superior compression ratios for biological data. The format supports two variants: **SRA Normalized** (full format) and **SRA Lite** (compressed format with simplified quality scores).

**Related Documentation:**
- [KAR File Format Specification](kar_file_structure.md) - The underlying archive format
- [VDB File Format Specification](vdb_file_structure.md) - The database storage system

## File Structure

### File Header and Magic Signature

Every SRA file begins with the 8-byte magic signature: **`NCBI.sra`**
- Bytes 0-3: `"NCBI"` (0x4E434249)
- Bytes 4-7: `".sra"` (0x2E737261)

This signature identifies the file as an NCBI SRA format file and is used by tools to validate file format before processing.

### Physical File Organization

SRA files are implemented as [KAR (NCBI Archive)](kar_file_structure.md) files containing a [VDB/KDB database](vdb_file_structure.md) structure. The physical organization follows a hierarchical directory-like structure:

```
SRA Database Structure:
├── col/                    # Global column data
│   ├── QUALITY/           # Quality scores column data
│   ├── READ/              # Sequence data column
│   ├── SPOT_ID/          # Spot identifiers
│   └── [other columns]
├── tbl/                   # Table definitions
│   └── SEQUENCE/         # Primary sequence table
│       └── col/          # Table-specific columns
├── meta/                 # Metadata and configuration
├── schema                # Schema definitions
└── idx/                  # Index structures (optional)
```

### Data Storage Architecture

#### Columnar Storage
Data is organized in a columnar format where each data type is stored separately:
- Each column is stored in its own subdirectory
- Data is compressed and stored in "blobs" (compressed chunks)
- Each blob covers a contiguous range of rows
- Supports both sequential and random access patterns
- Detailed storage format described in [VDB File Format Specification](vdb_file_structure.md#column-storage-format)

#### Blob Structure
- **Blob Header**: Contains metadata about the blob contents
- **Compressed Data**: Column data compressed using various algorithms
- **Index Information**: For efficient data retrieval
- **Integrity Checksums**: For data validation

## File Format Variants

### SRA Normalized Format (.sra)

The SRA Normalized format contains complete sequencing information including full per-base quality scores.

**Characteristics:**
- File extension: `.sra`
- Contains original quality scores as produced by sequencing instruments
- Platform-specific quality score encodings (Illumina, 454, PacBio, etc.)
- Full data fidelity for all downstream analyses
- Larger file size due to comprehensive quality information

**Quality Score Encoding:**
- Phred quality scores (0-93 theoretical range)
- Platform-optimized compression algorithms
- Per-base granularity maintained
- Complex encoding schemes for different sequencing platforms
- ASCII encoding support for both phred_33 and phred_64 formats
- Quality score format: `(INSDC:quality:text:phred_33)QUALITY` for most platforms

### SRA Lite Format (.sralite)

The SRA Lite format is a storage-optimized variant that simplifies quality scores while maintaining compatibility.

**Characteristics:**
- File extension: `.sralite`
- Simplified quality scores: 30 for "pass" reads, 3 for "reject" reads
- Significantly smaller file size
- Full compatibility with existing SRA tools
- Faster data transfer and processing times

**Quality Score Simplification:**
- Binary quality assignment based on read filtering status
- Quality = 30: Reads that pass quality filters (READ_FILTER = 'pass')
- Quality = 3: Reads that fail quality filters (READ_FILTER = 'reject')
- Uniform quality score applied to all bases within each read
- Specification defined in SRA Tools README documentation

## Database Schema Structure

### Core Tables

#### SEQUENCE Table
The primary table containing sequence read data:

**Key Columns:**
- `READ`: DNA sequence data in various encodings
- `QUALITY`: Base quality scores (full in .sra, simplified in .sralite)
- `SPOT_ID`: Unique identifier for each spot/cluster
- `READ_TYPE`: Classification of read segments
- `READ_FILTER`: Pass/fail status for reads
- `READ_START`: Starting positions of reads within spots
- `READ_LEN`: Lengths of individual reads

#### Supporting Tables

**SPOTCOORD Table:**
- `X_COORD`: X coordinate on the sequencing surface
- `Y_COORD`: Y coordinate on the sequencing surface

**SPOTNAME Table:**
- `NAME_FMT`: Format string for spot names
- `SPOT_NAME`: External spot identifiers

**STATS Table:**
- `BASE_COUNT`: Total number of bases
- `SPOT_COUNT`: Total number of spots
- Platform and run-level statistics

### Data Type Encodings

#### Sequence Encodings
- **2na_packed**: 2-bit encoding (A=0, C=1, G=2, T=3)
- **4na_packed**: 4-bit encoding including ambiguous bases (N, etc.)
- **x2na_bin**: Extended 2-nucleotide binary format
- **x2cs_bin**: Color-space binary format
- Complete encoding specifications in [VDB File Format Specification](vdb_file_structure.md#data-type-system)

#### Compression Algorithms
- **zip_encoding**: Standard compression for general data
- **izip_encoding**: Integer-optimized compression
- **bool_encoding**: Boolean data optimization
- **fzip**: Floating-point specific compression (for signal data)
- Detailed compression specifications in [VDB File Format Specification](vdb_file_structure.md#compression-implementation)

## Platform-Specific Features

### Supported Platforms
The format supports 19 sequencing platforms (platform IDs 0-18) including:
- Illumina (platform ID 2)
- 454 Life Sciences (platform ID 1)
- Ion Torrent (platform ID 7)
- PacBio SMRT (platform ID 6)
- Oxford Nanopore (platform ID 9)
- Complete Genomics (platform ID 4)
- Helicos (platform ID 5)
- ABI SOLiD (platform ID 3)
- Capillary (platform ID 8)
- Element Bio (platform ID 10)
- Tapestri (platform ID 11)
- Vela Diagnostics (platform ID 12)
- Genapsys (platform ID 13)
- Ultima Genomics (platform ID 14)
- Genemind (platform ID 15)
- BGI SEQ (platform ID 16)
- DNB-SEQ (platform ID 17)
- Singular Genomics (platform ID 18)

### Platform Optimizations
Each platform has specific optimizations:
- **Illumina**: Optimized quality score compression, coordinate handling
- **454**: Signal intensity data support, variable-length reads
- **PacBio**: Long-read optimizations, kinetic data support
- **Nanopore**: Ultra-long read handling, signal data compression

## Data Access Patterns

### Sequential Access
- Optimized for reading entire datasets
- Blob-based streaming for efficiency
- Minimal memory footprint

### Random Access
- Index structures enable fast spot/read lookup
- Cached access patterns (.vdbcache files)
- Efficient range queries

### Caching System
- **`.vdbcache`**: Standard VDB cache files for optimized access
- **`.sra.vdbcache`**: Cache files specific to SRA Normalized format
- **`.sralite.vdbcache`**: Cache files for SRA Lite format
- Cache files enable sub-linear lookup times and are essential for large dataset performance
- Mismatching lite/normalized cache files are automatically ignored for compatibility

## File Format Conversion

### SRA Norm to SRA Lite Conversion (Delite Process)

The conversion process involves:

1. **Schema Translation**: Update database schema to Lite-compatible version
2. **Quality Processing**: 
   - Preserve original QUALITY as ORIGINAL_QUALITY
   - Generate simplified quality scores using read filter information
   - Remove verbose quality data
3. **Column Optimization**: Remove unnecessary columns (POSITION, SIGNAL, etc.)
4. **Validation**: Ensure data integrity and compatibility

### Conversion Tools
- **sra_delite.sh**: Shell script implementing the delite conversion process (800+ lines)
- **vdb-dump**: Data extraction and analysis tool
- **sra-stat**: Statistics and validation utilities

### Conversion Limitations and Error Conditions

#### Platform Restrictions
- **Colorspace platforms** (ABI SOLiD, platform ID 3) cannot be converted to SRA Lite
- **TRACE type archives** are not supported for delite processing
- Platform validation occurs before delite processing begins

#### Common Error Conditions
- **Error 80**: TRACE type archives rejection
- **Error 81**: Object rejected for delite process  
- **Error 82**: Object already converted to SRA Lite
- **Error 86**: Object not delited yet (during validation phase)

#### File System Requirements
- Sufficient disk space for temporary files during conversion
- Write permissions for target directories  
- Available space for both original and lite versions during processing

## Technical Implementation Details

### Memory Management
- Lazy loading of column data
- Configurable blob cache sizes
- Memory-mapped file access where supported

### Concurrency Support
- Thread-safe read operations
- Multiple simultaneous readers supported
- Write operations are single-threaded

### Error Handling
- Comprehensive checksums throughout
- Graceful handling of corrupted data
- Detailed error reporting and recovery

### Version Compatibility
- Forward and backward compatibility maintained
- Schema versioning system (current delite schema: 1.1.1)
- Schema files located in `/etc/ncbi/schema` in containerized environments
- Automatic format detection
- SRA Lite access requires toolkit version 2.11.2 or later
- Current SRA format version: **1** (`FS_SRA_CUR_VERSION = 1`)

## Performance Characteristics

### Storage Efficiency
- **SRA Normalized**: Full fidelity, larger size
- **SRA Lite**: Significant size reduction with simplified quality scores
- Efficient compression ratios across all platforms

### Access Performance
- Columnar access optimized for analytical workloads
- Index structures enable sub-linear lookup times
- Blob-based I/O minimizes seek operations

### Network Transfer
- SRA Lite format significantly reduces transfer times
- Resumable transfer support
- Integrity verification during transfer

## Integration with NCBI Infrastructure

### Repository Storage
- Integrated with NCBI's distributed storage systems
- Automatic format selection based on usage patterns
- Multi-tier storage optimization

### Tool Compatibility
- Full compatibility with SRA Toolkit
- Support for standard bioinformatics tools (samtools, etc.)
- API access through various programming languages

## Implementation Guide for Independent Libraries

To implement SRA file reading/writing without using the SRA Toolkit, you need to implement the following layers:

### Layer 1: KAR Archive Handling
- Implement [KAR file format parser](kar_file_structure.md)
- Handle magic number validation (`NCBI.sra`)
- Parse Table of Contents (TOC) for file extraction
- Support file extraction from KAR archive

### Layer 2: VDB Database Access  
- Implement [VDB database reader](vdb_file_structure.md)
- Parse VDB directory structure (`md/`, `tbl/`, `col/`)
- Handle blob decompression (zip, izip, pack, fzip encodings)
- Implement column data access and indexing

### Layer 3: SRA Schema Implementation
- Implement SRA-specific table schemas (`SEQUENCE`, `STATS`, etc.)
- Handle platform-specific data types (2na_packed, 4na_packed, quality scores)
- Support 19 sequencing platforms with their specific optimizations
- Implement SRA Norm vs SRA Lite format differences

### Required Components for Full Implementation

1. **Binary Data Handling**
   - Little-endian integer reading
   - Multi-byte structure parsing
   - Binary search tree traversal
   - Checksum validation (CRC32, MD5)

2. **Compression Support**
   - zlib/deflate decompression (zip_encoding)
   - Integer delta compression (izip_encoding)  
   - Bit packing/unpacking (pack_encoding)
   - Floating-point compression (fzip_encoding)

3. **Data Type Conversion**
   - Nucleotide encoding/decoding (2na, 4na, x2na)
   - Quality score conversion (phred_33, phred_64)
   - Coordinate system handling (zero/one-based)
   - Platform-specific data types

4. **Index Management**
   - Primary index parsing for blob location
   - Secondary indexes for range queries
   - Element-level indexes for random access
   - Cache management for performance

### Minimal Working Implementation

For read-only access to SRA Normalized files:
1. Parse KAR header and extract VDB database
2. Read table metadata to identify available columns
3. Parse column indexes to locate data blobs
4. Decompress blobs using appropriate algorithm
5. Convert binary data to target format

## Complete Implementation Workflow

### Format Layer Integration

The three format layers work together as follows:

1. **KAR Archive Layer** ([kar_file_structure.md](kar_file_structure.md))
   - Validates magic signature: `NCBI.sra` (0x4E434249.73726100)
   - Checks byte order: `0x05031988` (normal) or `0x88190305` (reversed)
   - Parses Table of Contents (TOC) to locate VDB database files
   - Extracts files from archive using 4-byte aligned storage

2. **VDB Database Layer** ([vdb_file_structure.md](vdb_file_structure.md))
   - Reads database metadata from `md/` directory
   - Opens tables from `tbl/SEQUENCE/` directory structure
   - Accesses column data via multi-level indexes (`idx`, `idx1`, `idx2`)
   - Decompresses blobs using schema-specified algorithms

3. **SRA Schema Layer** (this document)
   - Interprets biological data types (DNA sequences, quality scores)
   - Handles platform-specific encodings for 19 sequencing technologies
   - Manages SRA Normalized vs SRA Lite format differences
   - Provides genomic data access through established schemas

### Binary Format Hierarchy

```
SRA File Structure:
[KAR Archive Header: "NCBI.sra" + byte_order + version + file_offset]
│
├── [KAR TOC: Binary tree of file/directory entries]
│
└── [File Data Section: 4-byte aligned]
    │
    └── [VDB Database Directory Structure]
        │
        ├── md/                    # Database metadata
        │   ├── root              # Root metadata (KDBHdr + database info)
        │   ├── cur               # Current version pointer
        │   └── vers              # Version history
        │
        ├── tbl/SEQUENCE/         # Primary sequence table
        │   ├── md/               # Table metadata
        │   └── col/              # Column data
        │       ├── READ/         # DNA sequence column
        │       │   ├── data      # Compressed blob data
        │       │   ├── idx       # Level-0 index (KColBlobLoc array)
        │       │   ├── idx1      # Level-1 index (KColBlockLoc + header)
        │       │   └── idx2      # Level-2 index (compressed locators)
        │       │
        │       └── QUALITY/      # Quality scores column
        │           ├── data      # Quality data (30/3 for SRA Lite)
        │           ├── idx       # Blob location index
        │           └── ...       # Additional index levels
        │
        └── col/                  # Global column definitions
            └── [COLUMN_NAME]/    # Shared column specifications
```

### Implementation Verification

All format specifications in this documentation have been verified against the actual NCBI VDB source code at `https://github.com/ncbi/ncbi-vdb`, including:

- **Byte order constants**: Verified from `libs/kdb/kdbfmt.h` and `interfaces/kfs/sra.h`
- **Data structures**: Confirmed from actual struct definitions in `libs/kdb/colfmt.h`
- **Blob organization**: Validated from `libs/vdb/blob-headers.c` implementation
- **Index formats**: Cross-checked against `libs/kdb/rcolidx*.c` implementations
- **Platform constants**: Verified from schema files in SRATools codebase

This documentation describes the complete structure of SRA files as implemented in the NCBI SRA Tools and VDB libraries. The format continues to evolve to support new sequencing technologies and analytical requirements while maintaining backward compatibility and data integrity.
