#!/usr/bin/env python3
"""
merge to gff3 files into one, resolve by Name attribute
"""
import argparse
import tempfile

import tc_utils as tc
# include version number
from version import __version__



def main():
    # parse arguments:
    parser = argparse.ArgumentParser(
            description='Merge two gff3 files into one, resolve by Name attribute.'
    )
    parser.add_argument('-1', '--gff3_1', required=True, help='gff3 file 1')
    parser.add_argument('-2', '--gff3_2', required=True, help='gff3 file 2')
    parser.add_argument('-o', '--output', required=True, help='output file gff3')

    args = parser.parse_args()

    # concatenate  gff3 1 and gff3 2 files into one to temp file
    tmp_gff3 = tempfile.NamedTemporaryFile(mode='w', delete=False)
    with open(tmp_gff3.name, 'w') as fout:
        with open(args.gff3_1) as f:
            for line in f:
                fout.write(line)
        with open(args.gff3_2) as f:
            for line in f:
                fout.write(line)

    tc.merge_overlapping_gff3_intervals(tmp_gff3.name, args.output)

if __name__ == "__main__":
    main()