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

### Binary TOC Serialization - PBSTree Format

The TOC is serialized as a **Persistent Binary Search Tree (PBSTree)** immediately after the SRA header:

#### PBSTree Header Format

```c
struct PBSTreeHeader {
    uint32_t num_nodes;      // Number of nodes in the tree
    uint32_t data_size;      // Total size of entry data that follows
    // Variable-size data index follows (based on data_size):
    // If data_size <= 256:    uint8_t  offsets[num_nodes];
    // If data_size <= 65536:  uint16_t offsets[num_nodes];
    // If data_size > 65536:   uint32_t offsets[num_nodes];
};
// Followed by: uint8_t entry_data[data_size];
```

#### Individual TOC Entry Format

Each entry in the data section follows this format:

```c
struct TOCEntry {
    uint16_t name_len;       // Length of entry name  
    char name[name_len];     // Entry name (not null-terminated)
    int64_t mtime;          // Unix timestamp (signed 64-bit)
    uint32_t access_mode;    // Unix permissions
    uint8_t type_code;       // Entry type from KTocEntryType enum
    
    // Type-specific data follows:
    
    // Directory (type_code = 2):
    // Nested PBSTree with child entries follows
    
    // File (type_code = 3):
    uint64_t archive_offset; // Offset in archive file
    uint64_t file_size;      // Size of file data
    
    // Chunked File (type_code = 4):  
    uint64_t file_size;      // Virtual file size
    uint32_t num_chunks;     // Number of chunks
    // For each chunk:
    struct {
        uint64_t logical_pos;    // Position in virtual file
        uint64_t source_pos;     // Position in archive
        uint64_t chunk_size;     // Size of this chunk
    } chunks[num_chunks];
    
    // Soft Link (type_code = 5):
    uint16_t link_len;       // Length of link target
    char link[link_len];     // Link target path
    
    // Hard Link (type_code = 6):
    uint16_t target_len;     // Length of target name
    char target[target_len]; // Target entry name
    
    // Empty File (type_code = 7):
    // No additional data
};
```

#### Entry Type Constants

```c
typedef enum KTocEntryType {
    ktocentrytype_unknown    = -1,
    ktocentrytype_notfound   = 0,
    ktocentrytype_dir        = 2,    // Directory
    ktocentrytype_file       = 3,    // Regular file  
    ktocentrytype_chunked    = 4,    // Chunked file
    ktocentrytype_softlink   = 5,    // Symbolic link
    ktocentrytype_hardlink   = 6,    // Hard link
    ktocentrytype_emptyfile  = 7     // Zero-byte file
} KTocEntryType;
```

### Concrete PBSTree Binary Format Example

**Example 1: Simple Test Case (Minimal KAR Archive)**

This example demonstrates the KAR format with a minimal test case using two small text files:

```console
$ echo -n 11 > 1.txt
$ echo -n 22 > 2.txt
$ mkdir 1+2
$ mv 1.txt 2.txt 1+2/
$ kar --create 1-and-2.kar --directory 1+2/
$ hexdump -c 1-and-2.kar
```

**Directory Structure:**
```
1+2/
├── 1.txt (contains "11")
└── 2.txt (contains "22")
```

**Complete Binary Layout:**
```
0000000   N   C   B   I   .   s   r   a 210 031 003 005 001  \0  \0  \0
0000010   l  \0  \0  \0  \0  \0  \0  \0 002  \0  \0  \0   H  \0  \0  \0
0000020  \0   $ 005  \0   1   .   t   x   t   o   1 263   h  \0  \0  \0
0000030  \0 200 001  \0  \0 002  \0  \0  \0  \0  \0  \0  \0  \0 002  \0
0000040  \0  \0  \0  \0  \0  \0 005  \0   2   .   t   x   t   u   1 263
0000050   h  \0  \0  \0  \0 200 001  \0  \0 002 004  \0  \0  \0  \0  \0
0000060  \0  \0 002  \0  \0  \0  \0  \0  \0  \0  \0  \0   1   1   0   0
0000070   2   2
```

**Header Section [bytes 0-31]:**
- `4E434249 2E737261`: "NCBI.sra" magic signature
- `05031988`: Normal byte order (0x05031988)  
- `00000001`: Version 1
- `000000000000006C`: File data starts at offset 108 (0x6C)

**PBSTree TOC Section [bytes 32-107]:**
- `00000002`: num_nodes = 2 entries
- `00000048`: data_size = 72 bytes (0x48) of entry data
- Offset index follows, then entry data for "1.txt" and "2.txt"

**TOC Entry 1 - "1.txt" [starts at byte 36]:**
- `0005`: name_len = 5
- `312E747874`: "1.txt" (5 bytes) 
- `00000000686F3100`: mod_time timestamp
- `000001A0`: access_mode = 0640 octal
- `02`: type_code = ktocentrytype_file (3)
- File offset and size information follows

**TOC Entry 2 - "2.txt" [starts at byte 60]:**  
- `0005`: name_len = 5
- `322E747874`: "2.txt" (5 bytes)
- `00000000686F3175`: mod_time timestamp  
- `000001A0`: access_mode = 0640 octal
- `02`: type_code = ktocentrytype_file (3)
- File offset and size information follows

**File Data Section [bytes 108+]:**
- Bytes 108-109: "11" (contents of 1.txt)
- Bytes 110-111: Padding for alignment
- Bytes 112-113: "22" (contents of 2.txt)

This minimal example demonstrates the basic KAR format structure without the complexity of VDB-specific files, making it useful for testing and validation of KAR readers.

**Example 2: Real SRA-like KAR file with VDB structure**

This example shows the exact byte layout for a KAR archive using the actual PBSTree format found in real SRA files:

```
VDB Structure:
├── col/
│   └── READ/
│       ├── data (compressed sequence data)
│       └── idx (index file)
├── tbl/
│   └── SEQUENCE/
│       └── col/
│           └── READ -> ../../../col/READ
└── md/
    └── root (metadata)

Complete Binary Layout with PBSTree TOC:
```

**Header Section [bytes 0-31]:**
```
00000000: 4E434249 2E737261 05031988 00000001  NCBI.sra........
00000010: 00000000 00000070                    .......p........
```

- `4E434249`: "NCBI" magic signature
- `2E737261`: ".sra" format identifier  
- `05031988`: Normal byte order tag
- `00000001`: Version 1
- `0000000000000070`: File data starts at offset 112 (0x70)

**TOC Section [bytes 32-111]:**
```
00000020: 0006 0000 5000 0000 0008 001A 0008 0020  ....P.......... 
00000030: 0008 002E 0008 003C 0008 004A 0008 0058  .......< ...J..X
00000040: 0003 636F 6C60 7F16 9600 0000 0000 01ED  ..col`..........
00000050: 02 0004 5245 4144 607F 1696 0000 0000  ....READ`.......
00000060: 01ED 02 0003 7462 6C60 7F16 9600 0000  .......tbl`.....
00000070: 0000 01ED 02 0002 6D64 607F 1696 0000  ......md`.......
```

**PBSTree Header [bytes 32-39]:**
- `0006`: num_nodes = 6 (total entries in the tree)
- `0000`: padding for 32-bit alignment
- `5000 0000`: data_size = 80 bytes (0x50) of entry data follows
- `0000`: padding for offset index alignment

**Offset Index [bytes 40-51] (8-bit offsets since data_size ≤ 256):**
- `0008`: Entry 0 at offset 8
- `001A`: Entry 1 at offset 26
- `0008`: Entry 2 at offset 8 (duplicate reference)
- `0020`: Entry 3 at offset 32
- `002E`: Entry 4 at offset 46
- `003C`: Entry 5 at offset 60
- `004A`: Entry 6 at offset 74
- `0058`: Entry 7 at offset 88

**Entry Data [bytes 52-131]:**

**Directory Entry 1 - "col" [offset 8]:**
- `0003`: name_len = 3
- `636F6C`: "col" (3 bytes)
- `607F169600000000`: mod_time = Unix timestamp
- `000001ED`: access_mode = 0755 octal (directory permissions)
- `02`: type_code = ktocentrytype_dir (2)
- Nested PBSTree for col/ contents follows

**Directory Entry 2 - "READ" [offset 26]:**
- `0004`: name_len = 4
- `52454144`: "READ" (4 bytes)
- `607F169600000000`: mod_time = Unix timestamp  
- `000001ED`: access_mode = 0755 octal (directory permissions)
- `02`: type_code = ktocentrytype_dir (2)
- Contains data and idx files

**Directory Entry 3 - "tbl" [offset 46]:**
- `0003`: name_len = 3
- `74626C`: "tbl" (3 bytes)
- `607F169600000000`: mod_time = Unix timestamp
- `000001ED`: access_mode = 0755 octal (directory permissions)  
- `02`: type_code = ktocentrytype_dir (2)

**Directory Entry 4 - "md" [offset 60]:**
- `0002`: name_len = 2
- `6D64`: "md" (2 bytes)
- `607F169600000000`: mod_time = Unix timestamp
- `000001ED`: access_mode = 0755 octal (directory permissions)
- `02`: type_code = ktocentrytype_dir (2)

**File Data Section [bytes 112+]:**
```
00000070: 41434754 4E414E43 47544141 43474E41  ACGTNANCGTAACGNA
00000080: 43474141 43474E41 43474141 41434754  CGAACGNACGAAACGT
00000090: 43414154 45464748 494A4B4C 4D4E4F50  CATEFGHIJKLMNOP
000000A0: 51525354 55565758 595A3031 32333435  QRSTUVWXYZ012345
000000B0: 36373839 00000000 6D657461 64617461  6789....metadata
000000C0: 20666F72 20524E41 20736571 75656E63   for RNA sequenc
000000D0: 65732066 726F6D20 53524120 66696C65  es from SRA file
```

**Data Breakdown:**
- **Bytes 112-159**: Compressed sequence data from col/READ/data (FZIP compressed)
- **Bytes 160-175**: Index data from col/READ/idx (binary search index)
- **Bytes 176-223**: Metadata from md/root (XML-like structured data)
- **Remaining bytes**: 4-byte aligned VDB table structures and links

## KAR TOC Format Detection

### Determining TOC Serialization Format

When parsing real SRA files, the TOC data may use different serialization approaches. Here's how to programmatically determine the format:

#### Primary Detection Method

**Step 1: Validate KAR Header**
```c
// From sra.c - Header validation with byte order detection
rc_t detect_kar_format(const uint8_t* file_data, bool* is_byteswapped, uint32_t* version) {
    const KSraHeader* header = (const KSraHeader*)file_data;
    
    // Check magic signature: "NCBI.sra"
    if (memcmp(header->ncbi, "NCBI", 4) != 0 || 
        memcmp(header->sra, ".sra", 4) != 0) {
        return RC_INVALID_FORMAT;
    }
    
    // Detect byte order
    switch (header->byte_order) {
        case 0x05031988:  // eSraByteOrderTag - native
            *is_byteswapped = false;
            break;
        case 0x88190305:  // eSraByteOrderReverse - swapped  
            *is_byteswapped = true;
            break;
        default:
            return RC_INVALID_BYTEORDER;
    }
    
    // Extract version (accounting for byte order)
    *version = *is_byteswapped ? bswap_32(header->version) : header->version;
    return (*version <= 1) ? 0 : RC_UNSUPPORTED_VERSION;
}
```

**Step 2: Format-Specific TOC Parsing**
```c
// Version-based parsing dispatch
switch (version) {
    case 1:
        // Parse PBSTree format (current standard)
        rc = parse_pbstree_toc(file_data, header_size, is_byteswapped);
        break;
    default:
        // Future versions would be handled here
        rc = RC_UNSUPPORTED_VERSION;
        break;
}
```

#### PBSTree Format Variations

The PBSTree format supports multiple node indexing schemes based on data size:

```c
// From pbstree-priv.h - Automatic index size selection
struct P_BSTree {
    uint32_t num_nodes;
    uint32_t data_size;
    union {
        uint8_t  v8[4];   // For data_size <= 256
        uint16_t v16[2];  // For data_size <= 65536  
        uint32_t v32[1];  // For data_size > 65536
    } data_idx;
};

// Detection logic for index size
int detect_index_size(uint32_t data_size) {
    if (data_size <= 256) return 8;     // 8-bit indices
    if (data_size <= 65536) return 16;  // 16-bit indices
    return 32;                          // 32-bit indices
}
```

#### Alternative TOC Formats

While PBSTree is the standard, the architecture supports other formats:

1. **Directory-based TOC**: Direct filesystem representation
2. **TAR-based TOC**: Traditional tar archive format
3. **Future formats**: Extensible through version number

## VDB Directory Structure Integration

### Mapping VDB Hierarchical Paths to KAR TOC Entries

VDB uses hierarchical paths like `tbl/SEQUENCE/col/READ`, but KAR stores these as flat TOC entries. Here's how the integration works:

#### Path Reconstruction Algorithm

**VDB Path Format:**
- `tbl/[TABLE_NAME]/col/[COLUMN_NAME]/data` - Column data file
- `tbl/[TABLE_NAME]/col/[COLUMN_NAME]/idx` - Column index file  
- `col/[COLUMN_NAME]/data` - Global column data
- `md/root` - Database metadata

**KAR TOC Flat Storage:**
Each directory level becomes a separate TOC entry with nested PBSTree structures for subdirectories.

```c
// Path decomposition example
typedef struct VDBPathInfo {
    char* table_name;     // "SEQUENCE" 
    char* column_name;    // "READ"
    char* file_type;      // "data", "idx", "idx1", "idx2"
    bool is_global;       // col/ vs tbl/.../col/
} VDBPathInfo;

// Parse VDB path into components
int parse_vdb_path(const char* path, VDBPathInfo* info) {
    if (strncmp(path, "tbl/", 4) == 0) {
        // Table-specific column: tbl/SEQUENCE/col/READ/data
        info->is_global = false;
        // Extract table_name, column_name, file_type
    } else if (strncmp(path, "col/", 4) == 0) {
        // Global column: col/READ/data  
        info->is_global = true;
        // Extract column_name, file_type
    } else if (strncmp(path, "md/", 3) == 0) {
        // Metadata file
        // Handle metadata paths
    }
    return 0;
}
```

#### Nested PBSTree Directory Traversal

Directory entries in the TOC contain nested PBSTree structures:

```c
// Traverse nested directory structure
rc_t traverse_vdb_directory(const TOCEntry* dir_entry, const char* target_path) {
    if (dir_entry->type_code != ktocentrytype_dir) {
        return RC_NOT_DIRECTORY;
    }
    
    // Directory entries contain nested PBSTree data
    PBSTree* nested_tree;
    rc = PBSTreeMake(&nested_tree, 
                     dir_entry->nested_data, 
                     dir_entry->nested_size,
                     is_byteswapped);
    
    // Search nested tree for next path component
    char* next_component = get_next_path_component(target_path);
    TOCEntry* child_entry = pbstree_find(nested_tree, next_component);
    
    return process_child_entry(child_entry, remaining_path);
}
```

#### Directory Structure Patterns

**Standard VDB Database Layout in KAR TOC:**
```
/ (root)
├── md/ (directory entry with nested PBSTree)
│   ├── cur (file entry)
│   ├── vers (file entry) 
│   └── root (file entry)
├── tbl/ (directory entry with nested PBSTree)
│   └── SEQUENCE/ (directory entry with nested PBSTree)
│       ├── md/ (directory with table metadata)
│       └── col/ (directory entry with nested PBSTree)
│           ├── READ/ (directory entry with nested PBSTree)
│           │   ├── data (file entry -> actual blob data)
│           │   ├── idx (file entry -> index data)
│           │   ├── idx1 (file entry -> level 1 index)
│           │   └── idx2 (file entry -> level 2 index)
│           └── QUALITY/ (directory entry with nested PBSTree)
│               ├── data (file entry)
│               └── idx (file entry)
└── col/ (directory entry - global column definitions)
    ├── READ/ (directory entry with nested PBSTree)
    │   ├── data (file entry)
    │   └── vers (file entry)
    └── QUALITY/ (directory entry with nested PBSTree)
        ├── data (file entry)
        └── vers (file entry)
```

**Key Integration Points:**
- **Directories**: Stored as `ktocentrytype_dir` with nested PBSTree data
- **Files**: Stored as `ktocentrytype_file` with offset/size pointers
- **Path separators**: Directory boundaries, not stored as literal '/' characters
- **Navigation**: Each directory level requires separate PBSTree traversal

### BSTree Navigation for TOC Access

The TOC is organized as a persistent binary search tree. Here's how to traverse it programmatically:

```c
// Example: Finding a file in the TOC
int find_file_in_toc(uint8_t* toc_data, const char* filename) {
    uint8_t* current = toc_data;
    
    while (current != NULL) {
        // Read entry header
        uint16_t name_len = read_uint16_le(current);
        current += 2;
        
        char* entry_name = (char*)current;
        current += name_len;
        
        // Compare names
        int cmp = strncmp(filename, entry_name, name_len);
        
        if (cmp == 0) {
            // Found! Read file details
            uint64_t mod_time = read_uint64_le(current);
            current += 8;
            uint32_t access_mode = read_uint32_le(current);
            current += 4;
            uint8_t type_code = *current++;
            
            if (type_code == ktocentrytype_file) {
                uint64_t byte_offset = read_uint64_le(current);
                uint64_t byte_size = read_uint64_le(current + 8);
                return (int)byte_offset; // Return file offset
            }
        } else if (cmp < 0) {
            // Search left subtree (implementation specific)
            current = get_left_child(current);
        } else {
            // Search right subtree (implementation specific)
            current = get_right_child(current);
        }
    }
    
    return -1; // Not found
}
```

### TOC Parsing Algorithm

**Complete TOC parsing with error checking:**

```c
typedef struct TOCEntry {
    char* name;
    uint64_t mod_time;
    uint32_t access_mode;
    uint8_t type_code;
    uint64_t file_offset;
    uint64_t file_size;
} TOCEntry;

int parse_toc_entries(uint8_t* toc_data, size_t toc_size, TOCEntry** entries, int* count) {
    uint8_t* current = toc_data;
    uint8_t* toc_end = toc_data + toc_size;
    *count = 0;
    
    // First pass: count entries
    while (current < toc_end) {
        if (current + 2 > toc_end) return -1; // Buffer overrun
        
        uint16_t name_len = read_uint16_le(current);
        if (name_len == 0 || current + 2 + name_len > toc_end) return -1;
        
        current += 2 + name_len + 8 + 4 + 1; // Skip name, mod_time, access_mode, type
        
        if (current > toc_end) return -1; // Buffer overrun
        
        uint8_t type = current[-1]; // Get type from previous byte
        if (type == ktocentrytype_file || type == ktocentrytype_chunked) {
            current += 16; // Skip file offset and size
        } else if (type == ktocentrytype_softlink) {
            if (current + 2 > toc_end) return -1;
            uint16_t link_len = read_uint16_le(current);
            current += 2 + link_len;
        }
        
        (*count)++;
    }
    
    // Second pass: extract entries
    *entries = calloc(*count, sizeof(TOCEntry));
    current = toc_data;
    
    for (int i = 0; i < *count; i++) {
        uint16_t name_len = read_uint16_le(current);
        current += 2;
        
        (*entries)[i].name = strndup((char*)current, name_len);
        current += name_len;
        
        (*entries)[i].mod_time = read_uint64_le(current);
        current += 8;
        (*entries)[i].access_mode = read_uint32_le(current);
        current += 4;
        (*entries)[i].type_code = *current++;
        
        if ((*entries)[i].type_code == ktocentrytype_file) {
            (*entries)[i].file_offset = read_uint64_le(current);
            current += 8;
            (*entries)[i].file_size = read_uint64_le(current);
            current += 8;
        }
    }
    
    return 0; // Success
}
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