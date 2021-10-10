#!/usr/bin/env python3

import os
import sys

import struct

# https://www.htslib.org/doc/bgzip.html

"""
The file contents are:
    uint64_t number_entries

followed by number_entries pairs of:
    uint64_t compressed_offset
    uint64_t uncompressed_offset
"""

def print_index(index_file:str) -> None:
    with open(index_file, 'rb') as fhd:
        tgtfile        = index_file[:-4]
        filesize       = os.path.getsize(tgtfile)
        number_entries = struct.unpack_from("Q", fhd.read(struct.calcsize('Q')))[0]
        
        print(f"number_entries: {number_entries:15,d}")
        print(f"filesize      : {filesize:15,d}")

        struct_len = struct.calcsize('QQ')
        for pos in range(number_entries):
            frag_bytes = fhd.read(struct_len)
            compressed_offset, uncompressed_offset = struct.unpack_from('QQ', frag_bytes)
            print(f"pos: {pos:15,d} compressed_offset {compressed_offset:15,d} uncompressed_offset {uncompressed_offset:15,d}")

        print(f"number_entries: {number_entries:15,d}")
        print(f"filesize      : {filesize:15,d}")

def main():
    index_file = sys.argv[1]
    print_index(index_file)

if __name__ == "__main__":
    main()