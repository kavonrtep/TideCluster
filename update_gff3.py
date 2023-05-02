#!/usr/bin/env python3
"""
This script take gff3 and update attribute using conversion table.
"""

import argparse
import tc_utils as tc

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Update gff3 attributes based on conversion table.'
            )
    parser.add_argument('-g', '--gff3', required=True, help='gff3 file')
    parser.add_argument(
            '-t', '--table', required=True,
            help='Conversion table as tab-delimited file. First column if '
                 'original attribute value, second column is new attribute '
                 'value.'
            )
    parser.add_argument('-o', '--output', required=True, help='output file gff3')
    parser.add_argument(
            '-a', '--attribute_name', default="Name",
            help='attribute name to update, default attribute is "%(default)s"'
            )

    args = parser.parse_args()
    ori_attribute = args.attribute_name + "_old"
    # read table to dictionary
    table = {}
    with open(args.table) as f:
        for line in f:
            line = line.rstrip()
            if line:
                old, new = line.split("\t")
                table[old] = new

    # parse gff3 file and update attributes ans save to output file
    with open(args.gff3) as infile, open(args.output, "w") as outfile:
        for line in infile:
            line = line.rstrip()
            if line.startswith("#"):
                outfile.write(line + "\n")
            else:
                gff3 = tc.Gff3Feature(line)
                if args.attribute_name in gff3.attributes_dict:
                    try:
                        gff3.attributes_dict[ori_attribute] = gff3.attributes_dict[
                            args.attribute_name]
                        gff3.attributes_dict[args.attribute_name] = table[
                            gff3.attributes_dict[args.attribute_name]]
                    except KeyError:
                        # keep original value
                        pass
                else:
                    print(F"Attribute {args.attribute_name} not found in {str(gff3)}")
                outfile.write(str(gff3) + "\n")
