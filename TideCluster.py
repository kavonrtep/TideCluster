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
from tc_utils import run_tidehunter, TideHunterFeature
from tc_utils import gff3_to_fasta, find_cluster_by_mmseqs2
from tc_utils import find_clusters_by_blast_connected_component
from tc_utils import get_cluster_size, add_cluster_info_to_gff3
from tc_utils import merge_overlapping_gff3_intervals
from tc_utils import save_consensus_files

# minimal python version is 3.6
assert sys.version_info >= (3, 6), "Python 3.6 or newer is required"


def clustering(args):
    """
    Run clustering on sequences defined in gff3 file and fasta file
    :param args: command line arguments
    :return:
    """
    fasta = args.fasta
    gff3 = args.gff3
    gff3_out = args.output

    # get consensus sequences for clustering
    consensus = {}
    consensus_dimers = {}
    for seq_id, seq, cons in gff3_to_fasta(gff3, fasta):
        mult = round(1 + 10000 / len(cons))
        consensus[seq_id] = cons * mult
        consensus_dimers[seq_id] = cons * 4
        if len(consensus[seq_id]) > 10000:
            consensus[seq_id] = consensus[seq_id][0:10000]

    # first round of clustering by mmseqs2
    clusters1 = find_cluster_by_mmseqs2(consensus)

    representative_id = set(clusters1.values())
    consensus_representative = {k: consensus_dimers[k] for k in representative_id}
    # second round of clustering by blastn
    clusters2 = find_clusters_by_blast_connected_component(consensus_representative)
    # combine clusters
    clusters_final = clusters1.copy()

    for k, v in clusters1.items():
        if v in clusters2:
            clusters_final[k] = clusters2[v]
        else:
            clusters_final[k] = v

    # get total size of each cluster, store in dict
    cluster_size = get_cluster_size(gff3, clusters_final)

    # representative id sorted by cluster size
    representative_id = sorted(cluster_size, key=cluster_size.get, reverse=True)
    # rename values in clusters dictionary
    cluster_names = {}
    for i, v in enumerate(representative_id):
        cluster_names[v] = F"TRC_{i + 1}"

    for k, v in clusters_final.items():
        clusters_final[k] = cluster_names[v]

    cons_cls, cons_cls_dimer = add_cluster_info_to_gff3(gff3, gff3_out, clusters_final)
    merge_overlapping_gff3_intervals(gff3_out, gff3_out)

    # save also first round of clustering for debugging
    cons_cls1, cons_cls_dimer1_ = add_cluster_info_to_gff3(
        gff3, gff3_out + "_1.gff3", clusters1
        )
    merge_overlapping_gff3_intervals(gff3_out + "_1.gff3", gff3_out + "_1.gff3")

    # for debugging
    # write consensus sequences by clusters to directory
    #  used gff3_out as base name for directory
    consensus_dir = gff3_out + "_consensus"
    # create directory if it does not exist
    save_consensus_files(consensus_dir, cons_cls, cons_cls_dimer)
    save_consensus_files(consensus_dir + "_1", cons_cls1, cons_cls_dimer1_)


def tidehunter(args):
    """
    run tidehunter on fasta file
    :param args: object with command line arguments

    """
    # get size of input file
    chunk_size = 500000
    overlap = 50000
    # this fill split sequences to chunk and all is stored in single file
    fasta_file_chunked, matching_table = tc_utils.split_fasta_to_chunks(
        args.fasta, chunk_size, overlap
        )
    results = run_tidehunter(
        fasta_file_chunked, args.tidehunter_arguments
        )
    with open(args.output, "w") as out:
        # write GFF3 header
        out.write("##gff-version 3\n")

        with open(results) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                feature = TideHunterFeature(line)
                # TideHunter is also returning sequences of NNN - to not include them
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
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest='command', help='TideHunter wrapper')
    parser_tidehunter = subparsers.add_parser(
        'tidehunter', help='Run wrapper of TideHunter'
        )

    parser_tidehunter.add_argument(
        '-f', '--fasta', type=str, required=True,
        help='Path to reference sequence in fasta format'
        )
    parser_tidehunter.add_argument(
        '-o', '--output', type=str, required=True, help='path to output file'
        )
    # TideHunter parameters as single arguments in quotes
    parser_tidehunter.add_argument(
        '-T', '--tidehunter_arguments', type=str, nargs="?", required=False, default="",
        help='additional arguments for TideHunter in quotes'
        )
    # add program description

    parser_clustering = subparsers.add_parser(
        'clustering', help='Run clustering on TideHunter output'
        )
    parser_clustering.add_argument(
        "-f", "--fasta", help="Reference fasta", required=True
        )
    parser_clustering.add_argument(
        "-g", "--gff3", help="gff3 file from TideHunter wrapper", required=True
        )
    parser_clustering.add_argument(
        "-o", "--output", help="Output gff3 file with clusters", required=True
        )

    parser.description = """Wrapper of TideHunter
    This script enable to run TideHunter on large fasta files in parallel. It splits
    fasta file into chunks and run TideHunter on each chunk. Then it merges results
    and clusters similar tandem repeats.
    """

    # make epilog, in epilog keep line breaks as preformatted text

    parser.epilog = ('''
    Example of usage:
    
    # first run tidehunter on fasta file to generate raw gff3 output
    TideCluster.py tidehunter -f test.fasta -o test.gff3 -T "-p 40 -P 3000 -c 5 -e 0.25 
    -t 46"
    
    # then run clustering on the output to cluster similar tandem repeats
    TideCluster.py clustering -f test.fasta -g test.gff3 -o test_clustered.gff3
    
    
    Recommended parameters for TideHunter:
    short monomers: -T "-p 10 -P 39 -c 5 -e 0.25"
    long monomers: -T "-p 40 -P 3000 -c 5 -e 0.25"
    
    For parallel processing include TideHunter -t option. 
    
    For more information about TideHunter parameters see TideHunter manual.
    ''')

    cmd_args = parser.parse_args()
    if cmd_args.command == "tidehunter":
        tidehunter(cmd_args)
    elif cmd_args.command == "clustering":
        clustering(cmd_args)
    else:
        parser.print_help()
        sys.exit(1)
