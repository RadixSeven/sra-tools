# KAR (NCBI Archive) File Format Specification

## Overview

KAR (NCBI Archive) is a custom archive format developed by NCBI for packaging and storing biological sequence data and metadata. Unlike standard archive formats like TAR or ZIP, KAR is specifically optimized for the storage patterns used by NCBI's Virtual Database (VDB) system. The format combines file archiving capabilities with structured metadata organization and supports efficient random access to contained files.

## File Header Structure

### Magic Signature and Header Layout

Every KAR file begins with a fixed-size header containing magic numbers and format information:

**Note**: The `KSraHeader` structure is defined in external VDB libraries (`kfs/sra.h`) rather than in the KAR loader implementation itself.

```c
typedef struct KSraHeader {
    char ncbi[4];       // Magic: "NCBI" (0x4E434249)
    char sra[4];        // Type: ".sra" (0x2E737261)
    uint32_t byte_order; // Byte order indicator
    uint32_t version;    // Format version (currently 1)
    union {
        struct {
            uint64_t file_offset;  // Offset where file data begins
        } v1;
    } u;
} KSraHeader;
```

**Header Components:**
- **Magic Numbers**: `"NCBI.sra"` - 8-byte signature identifying KAR format
- **Byte Order**: Determines endianness for multi-byte values
- **Version**: Format version number (currently only version 1 is supported)
- **File Offset**: 64-bit offset pointing to the start of file data section

### Byte Order Handling

**Byte Order Tags:**
- **Normal**: `eSraByteOrderTag = 0x05031988` - Native byte order
- **Reversed**: `eSraByteOrderReverse = 0x88190305` - Swapped byte order

**Note**: These constants are defined as `#define` macros rather than enum values due to integer size requirements.

All multi-byte integers in the file follow the byte order specified in the header.

### Header Size Calculation

For version 1 files, the header size is calculated as:
```c
header_size = base_header_size + sizeof(v1_union_data)
```
This allows for future format versions with different header sizes while maintaining backward compatibility.

## File Organization

### Overall File Layout

```
[Header] [TOC Section] [Alignment Padding] [File Data Section]
```

1. **Header Section**: Fixed-size header with format identification
2. **TOC Section**: Variable-size Table of Contents with file metadata
3. **Padding**: Zero-fill padding to align file data to 4-byte boundaries  
4. **File Data**: Actual file contents, sorted by size for optimal access

### File Sorting and Alignment

**File Sorting Strategy:**
- Files are sorted by size before storage (smallest to largest)
- This optimization improves access patterns for typical usage
- Sorting function: `kar_entry_sort_size()`

**Alignment Requirements:**
- All file data is aligned to **4-byte boundaries**
- Alignment function: `align_offset(offset, 4)`
- Padding filled with zero bytes between TOC and file data

## Table of Contents (TOC) Format

### Entry Types

```c
enum TOCEntryTypes {
    ktocentrytype_unknown = -1,
    ktocentrytype_notfound,     // Entry not found
    ktocentrytype_dir,          // Directory
    ktocentrytype_file,         // Regular file
    ktocentrytype_chunked,      // File stored in chunks
    ktocentrytype_softlink,     // Symbolic link
    ktocentrytype_hardlink,     // Hard link
    ktocentrytype_emptyfile,    // Zero-byte file
    ktocentrytype_zombiefile    // Deleted/invalid file
};
```

### Directory Structure

The TOC uses a hierarchical structure based on binary search trees for efficient access:

```c
typedef struct KAREntry {
    BSTNode node;             // Binary search tree node
    KTime_t mod_time;         // 64-bit modification timestamp
    const char *name;         // Entry name string
    KARDir *parentDir;        // Parent directory reference
    uint32_t access_mode;     // Unix-style permissions (rwxrwxrwx)
    uint8_t type;            // Entry type from enum above
    uint8_t eff_type;        // Effective type after resolution
    bool the_flag;           // General purpose flag
} KAREntry;

typedef struct KARDir {
    KAREntry base;           // Inherits from KAREntry
    BSTree contents;         // Binary search tree of child entries
} KARDir;

typedef struct KARFile {
    KAREntry base;           // Inherits from KAREntry
    uint64_t byte_size;      // File size in bytes
    uint64_t byte_offset;    // Offset in archive file
    KTocChunk *chunks;       // Array of chunks (for chunked files)
    uint32_t nun_chunks;     // Number of chunks (note: contains typo in implementation)
} KARFile;

typedef struct KARAlias {
    KAREntry base;           // Inherits from KAREntry
    KAREntry *resolved;      // Target entry for links
    const char *link;        // Link target path string
} KARAlias;
```

### Binary TOC Serialization

Each entry in the TOC is serialized in binary format as follows:

```c
struct SerializedTOCEntry {
    uint16_t name_len;       // Length of entry name
    char name[name_len];     // Entry name (not null-terminated)
    uint64_t mod_time;       // Modification timestamp (64-bit)
    uint32_t access_mode;    // Unix permissions (32-bit)
    uint8_t type_code;       // Entry type code
    
    // Type-specific data follows:
    
    // For files (ktocentrytype_file):
    uint64_t byte_offset;    // Offset to file data
    uint64_t byte_size;      // Size of file data
    
    // For directories (ktocentrytype_dir):
    // Child entries follow recursively
    
    // For symbolic links (ktocentrytype_softlink):
    uint16_t link_len;       // Length of link target
    char link[link_len];     // Link target path
    
    // For empty files (ktocentrytype_emptyfile):
    // No additional data needed
    
    // For chunked files (ktocentrytype_chunked):
    uint32_t num_chunks;     // Number of chunks
    // Array of chunk descriptors follows
};
```

### TOC Tree Organization

The TOC is organized as a **Persistent Binary Search Tree (PBSTree)** that provides:
- **Efficient Lookup**: O(log n) access to any entry
- **Ordered Traversal**: In-order traversal yields alphabetically sorted entries  
- **Compact Storage**: Binary tree structure minimizes metadata overhead
- **Random Access**: Direct access to any subtree without full parsing

## Chunked File Support

### Chunk Structure

For large files or files with non-contiguous storage requirements:

```c
typedef struct KTocChunk {
    uint64_t logical_offset;  // Offset within the logical file
    uint64_t chunk_offset;    // Offset within the archive
    uint64_t chunk_size;      // Size of this chunk
} KTocChunk;
```

**Note**: The exact KTocChunk structure definition varies across implementations. The above represents the logical structure based on usage patterns in the codebase.

### Chunked File Benefits
- **Non-contiguous Storage**: Files can be stored in multiple segments
- **Large File Support**: Files larger than available contiguous space
- **Future Compression**: Individual chunks can be compressed independently
- **Flexible Layout**: Allows for complex file organization strategies

## Data Storage and Access

### File Data Organization

**Storage Principles:**
1. Files sorted by size (ascending) before storage
2. All file data aligned to 4-byte boundaries
3. No compression at archive level (files stored as-is)
4. Zero-padding used for alignment

**Offset Calculation:**
```c
uint64_t calculate_aligned_offset(uint64_t current_offset) {
    return (current_offset + 3) & ~3ULL;  // Round up to next 4-byte boundary
}
```

### Access Patterns

**Sequential Access:**
- Optimized for reading entire archives
- TOC parsed once at open time
- File data accessed in storage order

**Random Access:**
- Direct file access via TOC lookup
- O(log n) entry location using BSTree
- Minimal seeks for file data retrieval

## Metadata and Permissions

### File Metadata

Each file entry preserves:
- **Name**: Variable-length string (max 65,535 characters)
- **Modification Time**: 64-bit timestamp (`KTime_t`)
- **Access Permissions**: 32-bit Unix-style permissions
- **Size**: 64-bit file size
- **Storage Offset**: 64-bit offset within archive

### Permission Handling

**Permission Format**: Standard Unix octal permissions (e.g., 0644, 0755)
**Permission Preservation**: Original file permissions maintained in archive
**Cross-platform**: Permissions preserved across different operating systems

### Timestamp Format

**Timestamp Type**: `KTime_t` (64-bit signed integer)
**Resolution**: Platform-dependent (typically nanosecond or microsecond)
**Epoch**: Unix epoch (January 1, 1970, 00:00:00 UTC)

## Error Handling and Validation

### Header Validation

1. **Magic Number Check**: Verify `"NCBI.sra"` signature
2. **Version Validation**: Ensure version <= current supported version
3. **Byte Order Check**: Validate byte order indicator
4. **Offset Validation**: Verify file_offset points within file bounds

### TOC Validation

1. **Name Length Check**: Ensure name_len > 0 and < 65536
2. **Type Code Validation**: Verify type_code is valid enum value
3. **Tree Structure**: Validate BSTree invariants and structure
4. **Offset Bounds**: Ensure all offsets point within file bounds
5. **Size Consistency**: Verify file sizes match available data

### Data Integrity

1. **Alignment Verification**: Check all file data starts at 4-byte boundaries
2. **Size Verification**: Confirm file sizes match TOC entries
3. **EOF Handling**: Graceful handling of truncated files
4. **Checksum Support**: Optional MD5/CRC32 validation (implementation-dependent)

## Implementation Constants

### Size Limits
- **Maximum file name length**: 65,535 characters (`UINT16_MAX`)
- **Maximum link target length**: 65,535 characters
- **Maximum file size**: 2^64 - 1 bytes (`UINT64_MAX`)
- **Maximum number of files**: Limited by available memory
- **Maximum number of chunks per file**: 2^32 - 1 (`UINT32_MAX`)

### Alignment Requirements
- **File data alignment**: 4 bytes
- **Header alignment**: No specific requirement
- **TOC alignment**: No specific requirement
- **String alignment**: 1 byte (no alignment required)

### Buffer Sizes
- **Default I/O buffer**: 128MB for file operations
- **TOC read buffer**: Platform-dependent (typically 64KB)
- **Maximum path length**: Platform PATH_MAX or 65,535 characters

## Version Management

### Current Version Support
- **Supported Version**: 1 only
- **Version Validation**: Strict - rejects any version > 1
- **Backward Compatibility**: Version 0 treated as invalid

### Future Extensibility

The header union structure allows for future format versions:
```c
union version_data {
    struct v1_data { uint64_t file_offset; } v1;
    struct v2_data { /* future version data */ } v2;
    // Additional versions can be added here
} u;
```

This design enables:
- **Backward Compatibility**: Older readers can detect newer versions
- **Forward Compatibility**: Newer readers can handle older versions  
- **Graceful Degradation**: Unknown versions can be rejected cleanly

## Binary Data Types

### Integer Types
- **uint8_t**: 8-bit unsigned integer
- **uint16_t**: 16-bit unsigned integer  
- **uint32_t**: 32-bit unsigned integer
- **uint64_t**: 64-bit unsigned integer
- **KTime_t**: 64-bit signed timestamp

### String Types
- **Length-prefixed strings**: uint16_t length + data bytes
- **No null termination**: Strings stored without trailing null byte in TOC
- **UTF-8 encoding**: All strings assumed to be UTF-8 encoded

### Endianness
- **Byte order**: Specified in header byte_order field
- **Multi-byte values**: All follow header byte order specification
- **String data**: Byte order not applicable (single-byte characters)

## Implementation Files Reference

### Core Implementation
- **`kar.c:795-914`**: Header structure and TOC serialization
- **`kar+.h:40-111`**: Data structure definitions
- **`kar+util.c`**: Utility functions for archive operations

### Supporting Components
- **`kar-args.c`**: Command line argument processing
- **`kar+meta.c`**: Metadata handling and validation
- **`ccsra.c`**: SRA-specific KAR integration

This specification provides complete technical details for implementing KAR file format readers and writers, including all binary layouts, data structures, algorithms, and validation requirements found in the NCBI SRATools codebase.