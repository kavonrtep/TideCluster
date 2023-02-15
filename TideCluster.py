#!/usr/bin/env python
"""Wrapper of TideHunter
input is a fasta file of DNA sequences
output is a table of detected tandem repeats by tidehunter
this wrapper split DNA sequence to chunks, run tidehunter on each chunk and merge results
coordinates in the table must be recalculated to the original sequence

"""
import os
import sys
import tc_utils
import argparse
import subprocess

# minimal python version is 3.6
assert sys.version_info >= (3, 6), "Python 3.6 or newer is required"


def run_tidehunter(fasta_file, num_threads=4, min_period=40, max_period=3000):
    """run tidehunter on fasta file
    require TideHunter to be in PATH
    version of TideHunter must be 1.4.3
    return path to tidehunter output file
    """
    # verify tidehunter version
    tidehunter_version = subprocess.check_output(
        "TideHunter -v", shell=True
        ).decode().strip()
    assert tidehunter_version == "1.4.3", "TideHunter version 1.4.3 is required"

    # run tidehunter
    tmp_file = fasta_file + ".out"
    tidehunter_cmd = (F"TideHunter -c 5 -p 40 -P 3000 -f 2 -o {tmp_file} -t {num_threads}"
                      F" -p {min_period} -P {max_period}"
                      F" {fasta_file}")

    subprocess.check_call(tidehunter_cmd, shell=True)
    return tmp_file


class TideHunterFeature:
    """class to store features from tidehunter output
    """

    def __init__(self, line):
        self.seq_name = line.split()[0]
        self.repeat_ID = line.split()[1] + "_" + line.split()[0]
        self.seq_length = int(line.split()[2])
        self.start = int(line.split()[3])
        self.end = int(line.split()[4])
        self.cons_length = int(line.split()[5])
        self.copy_numer = float(line.split()[6])
        self.aver_match = float(line.split()[7])
        self.full_length = int(line.split()[8])
        self.sub_position = line.split()[9]
        self.consensus = line.split()[10]

    def recalculate_coordinates(self, matching_table):
        """recalculate coordinates in the table using matching table
        """
        ori_name, ori_start, ori_end, _ = tc_utils.get_original_header_and_coordinates(
            self.seq_name, self.start, self.end, matching_table
            )
        self.seq_name = ori_name
        self.start = ori_start
        self.end = ori_end

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.seq_name, self.repeat_ID, self.seq_length, self.start, self.end,
            self.cons_length, self.copy_numer, self.aver_match, self.full_length,
            self.sub_position, self.consensus
            )

    def gff3(self):
        """return line in gff3 format
        """
        columns = [self.seq_name, "TideHunter", "repeat_region", self.start, self.end,
                   self.aver_match, ".", "."]
        attributes = {"ID": self.repeat_ID, "Consensus_sequence": self.consensus,
                      "Consensus_length": self.cons_length,
                      "Copy_number": self.copy_numer, "Average_match": self.aver_match}

        # put attributes in the right order and format, some are floats or ints - must be
        # converted to string
        columns_str = "\t".join([str(i) for i in columns])
        atr_str = ";".join(["{}={}".format(k, str(v)) for k, v in attributes.items()])
        gff_line = columns_str + "\t" + atr_str
        return gff_line


def main(args):
    """
    run tidehunter on fasta file
    :param args: object with command line arguments

    """
    # get size of input file
    chunk_size = 500000
    overlap = 50000
    # this fill split sequences to chunk and all is stored in single file
    fasta_file_chunked, matching_table = tc_utils.split_fasta_to_chunks(
        args.input, chunk_size, overlap
        )
    results = run_tidehunter(
        fasta_file_chunked, num_threads=args.threads, min_period=args.min_period,
        max_period=args.max_period
        )
    with open(args.output, "w") as out:
        # write GFF3 header
        out.write("##gff-version 3\n")

        with open(results) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                feature = TideHunterFeature(line)
                # TideHunte is also returning sequences of NNN - to not include them
                # in the output
                if feature.consensus == "N" * feature.cons_length:
                    continue
                feature.recalculate_coordinates(matching_table)
                out.write(feature.gff3() + "\n")
    # clean up
    os.remove(fasta_file_chunked)
    os.remove(results)

    with open(args.output + "_chunks.bed", "w") as out:
        for m in matching_table:
            out.write(F'{m[0]}\t{m[2]}\t{m[3]}\t{m[4]}\n')


if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input', type=str, required=True,
        help='path to input DNA sequence in fasta format'
        )
    parser.add_argument(
        '-o', '--output', type=str, required=True, help='path to output file'
        )

    # TideHunter parameters
    parser.add_argument(
        '-t', '--threads', type=int, default=4,
        help='number of threads to use, default 4'
        )
    parser.add_argument(
        "-p", '--min_period', type=int, default=40,
        help="minimum period of repeats to be detected, default 40"
        )
    parser.add_argument(
        "-P", '--max_period', type=int, default=3000,
        help="maximum period of repeats to be detected, default 3000"
        )

    # add program description
    parser.description = "Wrapper of TideHunter"
    cmd_args = parser.parse_args()
    main(cmd_args)
