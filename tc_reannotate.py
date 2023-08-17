#!/usr/bin/env python3
"""Reanotate tandem repeats using library of repeats generated by TideCluster/TAREAN
"""

import argparse
import os
import shutil
import tempfile
import tc_utils as tc

def get_seq_lengths(fasta_file):
    """Get length of sequences in fasta file
    Args:
        fasta_file (str): path to fasta file
    Returns:
        dict: dictionary with sequence names as keys and sequence lengths as values
        """
    # same sequences names occure multiple times in fasta file
    # read length of the sequence and keep is only if it was not read before or if it
    # is longer then the one read before
    seq_lengths = {}
    index = 0
    with open(fasta_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                index = index + 1
                seq_name = F" {line[1:]} {index}"
                seq_lengths[seq_name] = 0
            else:
                seq_lengths[seq_name] += len(line)

    seq_lengths_unique = {}
    for seq_name in seq_lengths:
        seq_name_unique = seq_name.split()[0].split("#")[0]
        if seq_name_unique not in seq_lengths_unique:
            seq_lengths_unique[seq_name_unique] = seq_lengths[seq_name]
        elif seq_lengths_unique[seq_name_unique] < seq_lengths[seq_name]:
            seq_lengths_unique[seq_name_unique] = seq_lengths[seq_name]

    return seq_lengths_unique


def repeatmasker_to_gff3(repeatmasker_file, gff3_file):
    """Parse RepeatMasker output file and create gff3 file.
    Args:
        repeatmasker_file (str): path to RepeatMasker output file
        gff3_file (str): path to gff3 file to create
    Returns:
        None
    """
    # parse RepeatMasker output file

    with open(repeatmasker_file) as f, open(gff3_file, "w") as gff3_out:
        gff3_out.write("##gff-version 3\n")
        # skip header 3 lines
        for i in range(3):
            next(f)
        for line in f:
            items = line.split()
            seq_name = items[4]
            start = items[5]
            end = items[6]
            strand = items[8]
            strand = "+" if strand == "C" else "-"
            Name = items[10]
            gff_line = (F"{seq_name}\tRepeatMasker\ttandem repeat\t{start}"
                        F"\t{end}\t.\t{strand}\t.\tCluster_ID={Name};Name={Name}")
            gff_record = tc.Gff3Feature(gff_line)
            gff3_out.write(gff_record.print_line())


def main():
    # parse arguments
    parser = argparse.ArgumentParser(
            description=(
                "Annotate tandem repeats in a FASTA file using the tandem repeat library "
                "from "
                "TideCluster. The reference sequence is annotated by RepeatMasker, "
                "and the TRC library serves as a custom repeat library. The output from "
                "RepeatMasker is filtered to retain only high-quality TRC hits. "
                "Overlapping TRC hits are consolidated. Regions shorter than twice the monomer "
                "length are excluded from the output. Ambiguous regions are also removed."
                )
            )
    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument(
        "-r", "--repeatmasker_file", help="RepeatMasker output file", type=str)
    group.add_argument(
        "-s", "--ref_seq", help="FASTA file to annotated by TRC library", type=str
        )

    parser.add_argument(
        "-f" , "--fasta_file", help="Fasta file wiht TRC  library used for RepeatMasker "
                                "search",
            type=str,
            required=True
        )

    parser.add_argument("-c", "--cpu", help="Number of CPUs to use", type=int, default=1)


    parser.add_argument(
        "-o", "--output", help="GFF3 output file", type=str, required=True
        )

    parser.add_argument("-d", "--debug", help="Keep temp files for debuging",
                        action="store_true", default=False)

    args = parser.parse_args()

    # create temp directory
    temp_dir = tempfile.mkdtemp()
    if args.ref_seq:
        # run RepeatMasker in temp directory
        cmds = (F"RepeatMasker -dir {temp_dir} {args.ref_seq} "
                F"-nolow -no_is -e ncbi  -lib {args.fasta_file} -pa {args.cpu}")
        print("Running RepeatMasker...")
        os.system(cmds)
        # list files in temp directory
        files = os.listdir(temp_dir)
        args.repeatmasker_file = F"{temp_dir}/{os.path.basename(args.ref_seq)}.out"


    # get length of sequences in fasta file - this is the length of sat dimers
    seq_lengths = get_seq_lengths(args.fasta_file)
    print(seq_lengths)
    # intermediate files - in temp directory:

    rm_gff3 = F"{temp_dir}/rm.gff3"
    rm_merged_gff3 = F"{temp_dir}/rm.merged.gff3"
    rm_merged_ext_gff3 = F"{temp_dir}/rm.merged.extended.gff3"
    rm_merged_ext_merged_gff3 = F"{temp_dir}/rm.merged.extended.merged.gff3"
    gr_filter_gff3 = F"{temp_dir}/gr_filter.gff3"

    # parse RepeatMasker output file and create gff3 file
    repeatmasker_to_gff3(args.repeatmasker_file, rm_gff3)
    # merge overlapping intervals of the same TRC
    tc.merge_overlapping_gff3_intervals(rm_gff3, rm_merged_gff3, use_strand=True)
    # extend intervals by offset - this will be used as filter for the original gff3 file
    tc.extend_gff3_intervals(rm_merged_gff3, rm_merged_ext_gff3, seq_lengths)
    # this extended and merged will be used as filter for the original gff3 file
    tc.merge_overlapping_gff3_intervals(rm_merged_ext_gff3,
                                        rm_merged_ext_merged_gff3, use_strand=True)
    # filter on min size

    gr_filter = []
    with open(rm_merged_ext_merged_gff3) as f:
        for line in f:
            if line.startswith("#"):
                continue
            gff_record = tc.Gff3Feature(line)
            if gff_record.end - gff_record.start < seq_lengths[
                gff_record.attributes_dict["Name"]]:
                continue
            else:
                gr_filter.append(
                    tc.GRange(
                        gff_record.seqid, gff_record.start,
                        gff_record.end,
                        gff_record.attributes_dict["Name"],
                        gff_record.strand
                        )
                    )

    # export gr_filter to gff3 file

    with open(gr_filter_gff3, "w") as f:
        f.write("##gff-version 3\n")
        for gr in gr_filter:
            f.write(gr.print_gff3())
            f.write("\n")

    # tmp_gff_extended will be used as filter for gff
    gr = tc.read_gff3_to_grange_list(rm_merged_gff3)
    granges_filtered = tc.filter_intervals(gr_filter,gr)

    granges_parsed = tc.split_intervals(granges_filtered)
    with open(args.output, "w") as f:
        f.write("##gff-version 3\n")
        for gr in granges_parsed:
            # check ambiguity - if there is ambiguity, skip
            if "|" in gr.name:
                continue
            f.write(gr.print_gff3())
            f.write("\n")

    if not args.debug:
        shutil.rmtree(temp_dir)
    else:
        print(F"Temp directory: {temp_dir}")
        # content of temp directory
        print("Temp directory content:")
        print(os.listdir(temp_dir))


if __name__ == "__main__":
    main()
