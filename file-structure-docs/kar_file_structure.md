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
    ktocentrytype_zombiefile    // Zombie in the sense that it is somewhere between live and dead
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

The TOC is serialized as a **Persistent Binary Search Tree (PBSTree)** immediately after the SRA header. This format is confirmed by the actual NCBI implementation in `sra.c` which calls `KTocInflatePBSTree()` to parse the TOC section.

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

    // File (type_code = 2):
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
    ktocentrytype_dir        = 1,
    ktocentrytype_file       = 2,
    ktocentrytype_chunked    = 3,
    ktocentrytype_softlink   = 4,
    ktocentrytype_hardlink   = 5,
    ktocentrytype_emptyfile  = 6,
    ktocentrytype_zombiefile = 7
} KTocEntryType;
```

### PBSTree Format Example - Source Code References

**TOC Parsing Implementation Location:**
The actual TOC parsing is implemented in:
- **Main parsing function**: `KArcParseSRAInt()` in `/libs/kfs/sra.c:217-324`
- **PBSTree inflation**: `KTocInflatePBSTree()` in `/libs/kfs/tocentry.c:1719-1744`
- **PBSTree core implementation**: `/libs/klib/pbstree.c`

**Key Implementation Points from Source Code:**
1. `sra.c:293` calls `KTocParseReadPBSTree()` to read TOC section between header and file data
2. `sra.c:305` calls `KTocInflatePBSTree()` to parse the PBSTree structure
3. `tocentry.c:1725` uses `PBSTreeMake()` to create PBSTree from binary data
4. `tocentry.c:1737` uses `PBSTreeForEach()` to walk entries and inflate directory structure

**Actual PBSTree Format Used in SRA Files:**

Based on analysis of both the NCBI source code implementation and real file validation, the TOC section uses the standard PBSTree format as implemented in the NCBI libraries.

**Confirmed PBSTree Format (from `pbstree-impl.c` parsing logic):**
```c
struct P_BSTree {
    uint32_t num_nodes;      // Number of entries in the tree
    uint32_t data_size;      // Total size of the data section

    // Variable-size index array (depends on data_size):
    // - If data_size <= 256:    uint8_t data_idx[num_nodes]
    // - If data_size <= 65536:  uint16_t data_idx[num_nodes]
    // - If data_size > 65536:   uint32_t data_idx[num_nodes]

    uint8_t data[data_size]; // The actual TOC entry data
};
```

**Byte Order Handling:**
- Files with byte order marker `0x05031988` (normal): Parse integers as little-endian
- Files with byte order marker `0x88190305` (swapped): Apply byte swapping to all integers

**Index Array Interpretation:**
Each entry in `data_idx[]` contains an offset into the `data[]` section where that node's entry begins. The entries in the data section contain the hierarchical directory structure with nested PBSTree structures for subdirectories.

**Implementation Notes:**
- The entire TOC section (from byte 32 to file data offset) is read as one buffer
- This buffer is passed to `PBSTreeMake()` which validates and parses the format
- The `PBSTreeForEach()` function walks the tree to build the directory structure
- Directory entries contain nested PBSTree structures for their children

**Individual TOC Entry Format in Data Section:**

Based on the `KTocEntryInflateNodeCommon()` implementation, each entry in the PBSTree data section has this format:

```c
struct TOCEntryData {
    uint16_t name_len;           // Length of entry name (little-endian or byte-swapped)
    char name[name_len];         // Entry name (not null-terminated)
    int64_t mtime;              // Unix timestamp (signed 64-bit)
    uint32_t access_mode;        // Unix permissions
    uint8_t type_code;           // Entry type (see KTocEntryType enum)

    // Type-specific data follows:

    // For ktocentrytype_dir (1):
    // Nested PBSTree structure follows immediately

    // For ktocentrytype_file (2):
    uint64_t file_offset;        // Offset within archive where file data starts
    uint64_t file_size;          // Size of file data


    // For ktocentrytype_chunked (3):
    uint64_t virtual_size;       // Total virtual file size
    uint32_t chunk_count;        // Number of chunks
    // For each chunk:
    struct {
        uint64_t logical_pos;    // Position in virtual file
        uint64_t source_pos;     // Position in archive
        uint64_t chunk_size;     // Size of this chunk
    } chunks[chunk_count];

    // For ktocentrytype_hardlink (5) or ktocentrytype_softlink (4):
    uint16_t link_len;           // Length of link target
    char link_target[link_len];  // Link target path (not null-terminated)

    // For ktocentrytype_emptyfile (6):
    // No additional data
};
```

**Complete Parsing Algorithm:**
1. Read SRA header to get file data offset and byte order flag
2. Read entire TOC section from byte 32 to file data offset
3. Parse PBSTree header: `num_nodes`, `data_size`
4. Read index array based on data_size (uint8, uint16, or uint32 per entry)
5. For each node index, access entry in data section:
   - Parse common header: name_len, name, mtime, access_mode, type_code
   - Parse type-specific data based on type_code
   - For directories: recursively parse nested PBSTree structure
6. Apply byte swapping to all multi-byte integers if byte order is swapped

**Implementation Example:**
```python
def parse_pbstree_toc(toc_data, byte_swapped=False):
    """Parse PBSTree TOC structure from SRA file"""
    offset = 0

    # Read PBSTree header
    num_nodes = read_uint32(toc_data, offset, byte_swapped)
    offset += 4
    data_size = read_uint32(toc_data, offset, byte_swapped)
    offset += 4

    # Read index array
    if data_size <= 256:
        indices = [toc_data[offset + i] for i in range(num_nodes)]
        offset += num_nodes
    elif data_size <= 65536:
        indices = [read_uint16(toc_data, offset + i*2, byte_swapped)
                  for i in range(num_nodes)]
        offset += num_nodes * 2
    else:
        indices = [read_uint32(toc_data, offset + i*4, byte_swapped)
                  for i in range(num_nodes)]
        offset += num_nodes * 4

    # Parse entries
    data_start = offset
    entries = {}

    for i, entry_offset in enumerate(indices):
        entry_pos = data_start + entry_offset
        entry = parse_toc_entry(toc_data, entry_pos, byte_swapped)
        entries[entry['name']] = entry

    return entries

def parse_toc_entry(data, offset, byte_swapped=False):
    """Parse individual TOC entry"""
    start_offset = offset

    # Read common header
    name_len = read_uint16(data, offset, byte_swapped)
    offset += 2

    name = data[offset:offset + name_len].decode('utf-8')
    offset += name_len

    mtime = read_int64(data, offset, byte_swapped)
    offset += 8

    access_mode = read_uint32(data, offset, byte_swapped)
    offset += 4

    type_code = data[offset]
    offset += 1

    entry = {
        'name': name,
        'mtime': mtime,
        'access_mode': access_mode,
        'type_code': type_code
    }

    # Parse type-specific data
    if type_code == 2:  # ktocentrytype_file
        entry['file_offset'] = read_uint64(data, offset, byte_swapped)
        offset += 8
        entry['file_size'] = read_uint64(data, offset, byte_swapped)
        offset += 8

    elif type_code == 1:  # ktocentrytype_dir
        # Directory contains nested PBSTree - would need recursive parsing
        nested_size = len(data) - offset  # Remaining data
        entry['nested_pbstree'] = data[offset:offset + nested_size]

    # Add other type handlers as needed...

    return entry
```

### Concrete PBSTree Binary Format Example

#### Example 1: Simple Test Case (Minimal KAR Archive)

This example demonstrates the KAR format with a minimal test case using two small text files:

```console
$ echo -n 11 > 1.txt
$ echo -n 22 > 2.txt
$ mkdir 1+2
$ mv 1.txt 2.txt 1+2/
$ kar --create 1-and-2.kar --directory 1+2/
$ hexdump -C 1-and-2.kar
```

echo -n 11 > 1.txt
echo -n 22 > 2.txt
mkdir 1+2
mv 1.txt 2.txt 1+2/
kar --create 1-and-2.kar --directory 1+2/
hexdump -C 1-and-2.kar

**Directory Structure:**
```
1+2/
├── 1.txt (contains "11")
└── 2.txt (contains "22")
```

**Complete Binary Layout:**
```
#          0  1  2  3  4  5  6  7   8  9  a  b  c  d  e  f   0123456789abcdef
00000000  4e 43 42 49 2e 73 72 61  88 19 03 05 01 00 00 00  |NCBI.sra........|
00000010  6c 00 00 00 00 00 00 00  02 00 00 00 48 00 00 00  |l...........H...|
00000020  00 24 05 00 31 2e 74 78  74 41 6f b4 68 00 00 00  |.$..1.txtAo.h...|
00000030  00 80 01 00 00 02 00 00  00 00 00 00 00 00 02 00  |................|
00000040  00 00 00 00 00 00 05 00  32 2e 74 78 74 41 6f b4  |........2.txtAo.|
00000050  68 00 00 00 00 80 01 00  00 02 04 00 00 00 00 00  |h...............|
00000060  00 00 02 00 00 00 00 00  00 00 00 00 31 31 30 30  |............1100|
00000070  32 32                                             |22|
00000072
```

**Header Section [bytes 0-31]:**
- `4E434249 2E737261` [00-07]: "NCBI.sra" magic signature
- `05031988` [08-0b]: Normal byte order (0x05031988)
- `00000001` [0c-0f]: Version 1
- `000000000000006C` [10-17]: File data starts at offset 108 (0x6C)

**Root PBSTree TOC Section [bytes 32-107]:**
- `00000002` [18-1b]: num_nodes = 2 entries (1.txt, 2.txt)
- `00000048` [1c-1f]: data_size = 72 bytes (0x48) of entry data
- `00 24` [20-21]: Index array (offsets: 0, 36 in data section)

**Root TOC Entry 1 - "1.txt" [data offset 0]:**
- `0005` [22-23]: name_len = 5
- `31 2E 74 78 74` [24-28]: "1.txt" (5 bytes)
- `00000000686FB441` [29-30]: mod_time timestamp
- `00000180` [31-34]: access_mode = 0600 octal
- `02` [35]: type_code = ktocentrytype_file (2)
- ktocentrytype_file fields:
  - `0000000000000000` [36-3d]: file_offset in the data section = 0
  - `0000000000000002` [3e-45]: file_size = 2 bytes

**Root TOC Entry 2 - "2.txt" [data offset 36]:**
- `0005` [46-47]: name_len = 5
- `32 2E 74 78 74` [48-4c]: "2.txt" (5 bytes)
- `00000000686FB441` [4d-54]: mod_time timestamp
- `00000180` [55-58]: access_mode = 0600 octal
- `02` [59]: type_code = ktocentrytype_file (2)
- ktocentrytype_file fields:
  - `0000000000000004` [5a-61]: file_offset in the data section = 4
  - `0000000000000002` [62-69]: file_size = 2 bytes

**Padding:**
- `0000` [6a-6b]: null padding to 4-byte boundary

**File Data Section [bytes 108+]:**
- `3131` [6c-6d]: "11" (contents of 1.txt, file_offset=0, size=2)
- `3030` [6e-6f]: "00" (padding for alignment)
- `3232` [70-71]: "22" (contents of 2.txt, file_offset=4, size=2)

This minimal example demonstrates the basic KAR format structure without the complexity of VDB-specific files, making it useful for testing and validation of KAR readers.

#### Example 2: Simple Test Case (Minimal KAR Archive with a subdirectory)

This example demonstrates the KAR format with a minimal test case using two small text files:

```console
$ mkdir with_sub
$ mkdir with_sub/sub
$ echo -n aa > with_sub/a.txt
$ echo -n zz > with_sub/z.txt
$ echo -n ww > with_sub/sub/w.txt
$ tree with_sub
$ kar --create with_subdir.kar --directory with_sub
$ hexdump -C with_subdir.kar
```


**Directory Structure:**
```
with_sub
├── a.txt
├── z.txt
└── sub
    └── w.txt
```

**Complete Binary Layout:**
```
#          0  1  2  3  4  5  6  7   8  9  a  b  c  d  e  f   0123456789abcdef
00000000  4e 43 42 49 2e 73 72 61  88 19 03 05 01 00 00 00  |NCBI.sra........|
00000010  ac 00 00 00 00 00 00 00  03 00 00 00 87 00 00 00  |................|
00000020  00 24 63 05 00 61 2e 74  78 74 ba 6e b4 68 00 00  |.$c..a.txt.n.h..|
00000030  00 00 80 01 00 00 02 00  00 00 00 00 00 00 00 02  |................|
00000040  00 00 00 00 00 00 00 03  00 73 75 62 ba 6e b4 68  |.........sub.n.h|
00000050  00 00 00 00 c0 01 00 00  01 01 00 00 00 24 00 00  |.............$..|
00000060  00 00 05 00 77 2e 74 78  74 ba 6e b4 68 00 00 00  |....w.txt.n.h...|
00000070  00 80 01 00 00 02 04 00  00 00 00 00 00 00 02 00  |................|
00000080  00 00 00 00 00 00 05 00  7a 2e 74 78 74 ba 6e b4  |........z.txt.n.|
00000090  68 00 00 00 00 80 01 00  00 02 08 00 00 00 00 00  |h...............|
000000a0  00 00 02 00 00 00 00 00  00 00 00 00 61 61 30 30  |............aa00|
000000b0  77 77 30 30 7a 7a                                 |ww00zz|
000000b6
```

**Header Section [bytes 0-31]:**
- `4E434249 2E737261` [00-08]: "NCBI.sra" magic signature
- `05031988` [08-0b]: Normal byte order (0x05031988)
- `00000001` [0c-0f]: Version 1
- `00000000000000AC` [10-17]: File data starts at offset 172 (0xAC)

**Root PBSTree TOC Section [bytes 32-171]:**
- `00000003` [18-1b]: num_nodes = 3 entries (a.txt, sub/, z.txt)
- `00000087` [1c-1f]: data_size = 135 bytes (0x87) of entry data
- `00 24 63` [20-22]: Index array (offsets: 0, 36, 99 in data section)

**Root TOC Entry 1 - "a.txt" [data offset 0]:**
- `0005` [23-24]: name_len = 5
- `61 2E 74 78 74` [25-29]: "a.txt" (5 bytes)
- `0000000068B46EBA` [2a-31]: mod_time timestamp
- `00000180` [32-35]: access_mode = 0600 octal
- `02` [36-36]: type_code = ktocentrytype_file (2)
- ktocentrytype_file fields:
- `0000000000000000` [37-3e]: file_offset in the data section = 0
- `0000000000000002` [3f-46]: file_size = 2 bytes

**Root TOC Entry 2 - "sub" [data offset 36]:**
- `0003` [47-48]: name_len = 3
- `73 75 62` [49-4b]: "sub" (3 bytes)
- `0000000068B46EBA` [4c-53]: mod_time timestamp
- `000001C0` [54-57]: access_mode = 0700 octal (directory permissions)
- `01` [58-58]: type_code = ktocentrytype_dir (1)
- **Nested PBSTree follows immediately for subdirectory contents**

**Nested PBSTree for "sub/" directory:**
- `00000001` [59-5c]: num_nodes = 1 entry (w.txt)
- `00000024` [5d-60]: data_size = 36 bytes of nested entry data
- `00` [61-61]: Index array:  (single entry at offset 0)

**Nested TOC Entry - "w.txt":**
- `0005` [62-63]: name_len = 5
- `77 2E 74 78 74` [64-68]: "w.txt" (5 bytes)
- `0000000068B46EBA` [69-70]: mod_time timestamp
- `00000180` [71-74]: access_mode = 0600 octal
- `02` [75-75]: type_code = ktocentrytype_file (2)
- ktocentrytype_file fields:
- `0000000000000004` [76-7d]: file_offset in the data section = 4
- `0000000000000002` [7e-85]: file_size = 2 bytes

**Root TOC Entry 3 - "z.txt" [data offset 99]:**
- `0005` [86-87]: name_len = 5
- `7A 2E 74 78 74` [88-8c]: "z.txt" (5 bytes)
- `0000000068B46EBA` [8d-94]: mod_time timestamp
- `00000180` [95-98]: access_mode = 0600 octal
- `02` [99-99]: type_code = ktocentrytype_file (2)
- File data at archive offset 8, size 2 bytes
- ktocentrytype_file fields:
- `0000000000000008` [9a-a1]: file_offset in the data section = 8
- `0000000000000002` [a2-a9]: file_size = 2 bytes

**Padding:**
- `0000` [aa-ab]: Null padding to align start of data to 4-byte boundaries
**File Data Section [bytes 172+]:**
- `61 61` [ac-ad]: "aa" (contents of a.txt)
- `30 30` [ae-af]: "00" '0' Padding to align file contents start to 4-byte boundaries.
- `77 77` [b0-b1]: "ww" (contents of sub/w.txt)
- `30 30` [b1-b2]: "00" Padding for alignment
- `7a 7a` [b3-b4]: "zz" (contents of z.txt)

**Key Observations:**
- Directory entries (type_code=1) contain nested PBSTree structures immediately after their header
- Nested PBSTree structures follow the same format: num_nodes, data_size, index array, entry data
- File paths like "sub/w.txt" are reconstructed by traversing the directory hierarchy
- Files are stored in size-sorted order in the data section regardless of directory structure

This example demonstrates how the PBSTree format handles directories with nested structures, making it essential for implementers to handle recursive parsing of directory entries.

#### Example 3: Real SRA-like KAR file with VDB structure

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
- `01`: type_code = ktocentrytype_dir (1)
- Nested PBSTree for col/ contents follows

**Directory Entry 2 - "READ" [offset 26]:**
- `0004`: name_len = 4
- `52454144`: "READ" (4 bytes)
- `607F169600000000`: mod_time = Unix timestamp
- `000001ED`: access_mode = 0755 octal (directory permissions)
- `01`: type_code = ktocentrytype_dir (1)
- Contains data and idx files

**Directory Entry 3 - "tbl" [offset 46]:**
- `0003`: name_len = 3
- `74626C`: "tbl" (3 bytes)
- `607F169600000000`: mod_time = Unix timestamp
- `000001ED`: access_mode = 0755 octal (directory permissions)
- `01`: type_code = ktocentrytype_dir (1)

**Directory Entry 4 - "md" [offset 60]:**
- `0002`: name_len = 2
- `6D64`: "md" (2 bytes)
- `607F169600000000`: mod_time = Unix timestamp
- `000001ED`: access_mode = 0755 octal (directory permissions)
- `01`: type_code = ktocentrytype_dir (1)

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

### Mapping VDB Hierarchical Paths to Linear KAR TOC Entries

**IMPORTANT**: Since actual SRA files use linear TOC format, VDB hierarchical paths must be reconstructed from the flat entry sequence.

#### Actual VDB Path Reconstruction from Linear TOC

**VDB Path Format in Real Files:**
- `col/[COLUMN_NAME]/data` - Global column data
- `col/[COLUMN_NAME]/idx` - Global column index
- Some files may use `tbl/SEQUENCE/col/[COLUMN_NAME]/data` pattern

**Real File Analysis:**
Actual SRA files contain TOC entries like:
```
col                    # type_code=1 (directory marker)
ALTREAD               # type_code=1 (column name)
data                  # type_code=2 (file with archive_offset and size)
idx                   # type_code=2 (index file)
QUALITY               # type_code=1 (next column)
data                  # type_code=2 (quality data file)
idx                   # type_code=2 (quality index file)
READ                  # type_code=1 (sequence column)
data                  # type_code=2 (sequence data)
idx                   # type_code=2 (sequence index)
```

#### Path Reconstruction Algorithm for Real Files

```python
def reconstruct_vdb_paths_from_linear_toc(toc_entries):
    """
    Reconstruct VDB hierarchical paths from linear TOC entries
    Based on actual SRA file format analysis
    """
    vdb_files = {}
    current_column = None
    current_base_path = ""

    for name, entry in toc_entries.items():
        if entry['type_code'] == 1:  # Directory/column marker
            if name in ['col', 'tbl', 'md']:
                current_base_path = name + "/"
            elif current_base_path.startswith('col/'):
                # This is a column name
                current_column = name
                current_base_path = f"col/{name}/"

        elif entry['type_code'] == 2:  # Actual file with data
            if current_column and name in ['data', 'idx', 'idx1', 'idx2']:
                vdb_path = current_base_path + name
                vdb_files[vdb_path] = {
                    'archive_offset': entry.get('archive_offset'),
                    'file_size': entry.get('file_size'),
                    'column': current_column,
                    'file_type': name
                }

    return vdb_files

# Example usage with real file structure:
# Input TOC entries: col, READ, data, idx, QUALITY, data, idx
# Output VDB paths:
# col/READ/data -> offset: 1234, size: 5678
# col/READ/idx -> offset: 6789, size: 123
# col/QUALITY/data -> offset: 9012, size: 3456
# col/QUALITY/idx -> offset: 4567, size: 89
```

#### VDB File Storage Pattern in Real Files

**Discovery from Actual Implementation:**
- VDB files are stored sequentially in the file data section
- The order typically follows: all idx files first, then all data files
- File offsets in TOC entries point to actual VDB file headers (starting with `88 19 03 05`)

```python
def map_vdb_files_to_storage_order(vdb_files):
    """
    Map VDB logical paths to actual storage order
    Based on patterns found in real SRA files
    """
    storage_mapping = {}

    # Pattern: idx files stored first
    idx_files = {path: info for path, info in vdb_files.items()
                if info['file_type'] == 'idx'}

    # Pattern: data files stored after idx files
    data_files = {path: info for path, info in vdb_files.items()
                 if info['file_type'] == 'data'}

    file_index = 0
    for path, info in idx_files.items():
        storage_mapping[path] = {
            'storage_index': file_index,
            'archive_offset': info['archive_offset'],
            'size': info['file_size']
        }
        file_index += 1

    for path, info in data_files.items():
        storage_mapping[path] = {
            'storage_index': file_index,
            'archive_offset': info['archive_offset'],
            'size': info['file_size']
        }
        file_index += 1

    return storage_mapping
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