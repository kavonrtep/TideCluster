#!/usr/bin/env python
"""Wrapper of TideHunter
input is a fasta file of DNA sequences
output is a table of detected tandem repeats by tidehunter
this wrapper split DNA sequence to chunks, run tidehunter on each chunk and merge results
coordinates in the table must be recalculated to the original sequence

"""
import glob
import os
import subprocess
import sys
from multiprocessing import Pool

import tc_utils as tc
import argparse
from tc_utils import run_tidehunter, TideHunterFeature
from tc_utils import gff3_to_fasta, find_cluster_by_mmseqs2
from tc_utils import find_clusters_by_blast_connected_component
from tc_utils import get_cluster_size, add_cluster_info_to_gff3
from tc_utils import merge_overlapping_gff3_intervals
from tc_utils import save_consensus_files, filter_gff_by_length
from tc_utils import read_fasta_sequence_size, get_repeatmasker_annotation
from tc_utils import add_attribute_to_gff, read_single_fasta_to_dictionary
from tc_utils import extract_sequences_from_gff3
from tc_utils import group_sequences_by_orientation, save_fasta_dict_to_file
from tc_utils import reverse_complement

# minimal python version is 3.6
assert sys.version_info >= (3, 6), "Python 3.6 or newer is required"


def tarean(prefix, gff, fasta=None, cpu=4):
    """
    Run tarean on genomic sequences specified in gff3 file
    from gff record extract sequence end analyse it with tarean algorithm
    :param gff: gff3 file with tidehunter results
    :param fasta: reference fasta
    :param cpu: number of cpu cores to use
    :return:
    """
    if gff is None:
        # use preferentially gff3 file with annotation
        if os.path.exists(prefix + "_annotation.gff3"):
            gff = prefix + "_annotation.gff3"
        elif os.path.exists(prefix + "_clustering.gff3"):
            gff = prefix + "_clustering.gff3"
        else:
            print(F"gff3 file annotation not found in {prefix}")
            return

    # directory for tarean results
    tarean_dir = prefix + "_tarean"
    # directory for fasta files extracted from gff3, it could be deleted after tarean run
    tarean_dir_fasta = tarean_dir + "/fasta"
    # create directory for tarean results
    if not os.path.exists(tarean_dir):
        # create directory for tarean results
        os.mkdir(tarean_dir)
    # create directory for fasta files
    if not os.path.exists(tarean_dir_fasta):
        os.mkdir(tarean_dir_fasta)

    fasta_dict = extract_sequences_from_gff3(gff, fasta, tarean_dir_fasta)
    l_debug = 0
    cmd_list = []
    print("preparaing sequence for TAREAN")
    for k, v in fasta_dict.items():
        l_debug += 1
        v_basename = os.path.basename(v)
        with open(v, "r") as f:
            seqs = read_single_fasta_to_dictionary(f)
        # does seqs contain only one sequence?
        if len(seqs) > 1:
            # more than one sequence in fasta file, set orientation
            seqs2 = group_sequences_by_orientation(seqs, k=8)
            for id in seqs:
                if id in seqs2['reverse']:
                    seqs[id] = reverse_complement(seqs[id])
        # save oriented sequences to file
        save_fasta_dict_to_file(seqs, v)
        # get current script path
        script_path = os.path.dirname(os.path.realpath(__file__))
        tarean_out = F"{tarean_dir}/{v_basename}_tarean"
        cmd = F"{script_path}/tarean/tarean.R -i {v} -s 0 -n {1} -o {tarean_out}"
        cmd_list.append(cmd)
    # run cmd tarean in parallel using multiprocessing module
    print("running TAREAN in parallel")
    with Pool(cpu) as p:
        total_jobs = len(cmd_list)
        completed_jobs = 0
        for _ in p.imap(tc.run_cmd, cmd_list):
            completed_jobs += 1
            print(F"completed {completed_jobs} of {total_jobs}")

    print("TAREAN finished")
    # get SSR info for tarean report from gff3 file
    SSR = {}
    with open(gff, "r") as f:
        for i in f:
            if i.startswith("#"):
                continue
            gff_record = tc.Gff3Feature(i)
            if "SSR" in gff_record.attributes:
                SSR[gff_record.attributes_dict['Name']] = gff_record.attributes_dict[
                "SSR"]
    # export SSR info to csv file
    with open(F"{tarean_dir}/SSRS_summary.csv", "w") as f:
        for k, v in SSR.items():
            f.write(F"{k}\t{v}\n")



    # final report:
    cmd = (F"{script_path}/tarean/tarean_report.R -i {tarean_dir} -o"
           F" {prefix}_tarean_report")
    print("Making final report")
    tc.run_cmd(cmd)



def annotation(prefix, library, gff=None, consensus_dir=None, cpu=1):
    """
    Run annotation on sequences defined in gff3 file based on coresponding
    consensus sequences in stored in consensu directory
    produce gff3 file with updated annotation  information
    :param prefix: prefix - base naame for input and output files
    :param library: library file for RepeatMasker
    :param gff: gff3 file with tidehunter results
    :param consensus_dir: directory with consensus sequences
    :param cpu: number of cpu cores to use
    :return:
    """
    if consensus_dir is None:
        consensus_dir = prefix + "_consensus"
    if gff is None:
        gff = prefix + "_clustering.gff3"
    gff_out = prefix + "_annotation.gff3"

    # get list consensus sequences from consensus directory
    # naming scheme is TRC_10_dimer.fasta
    # use glob to get all files in directory
    consensus_files = glob.glob(consensus_dir + "/TRC*dimers.fasta")
    # conncatenate all consensus sequences to one file
    consensus_files_concat = consensus_dir + "/consensus_sequences_all.fasta"
    if not os.path.exists(consensus_files_concat):
        with open(consensus_files_concat, "w") as f:
            for consensus_file in consensus_files:
                with open(consensus_file, "r") as f_in:
                    for line in f_in:
                        f.write(line)
    seq_lengths = read_fasta_sequence_size(consensus_files_concat)
    # run RepeatMasker on concatenated consensus sequences
    cmd = (F"RepeatMasker -pa {cpu} -lib {library} -e ncbi -s -no_is -norna "
           F"-nolow "
           F"-dir {consensus_dir} {consensus_files_concat}")
    subprocess.check_call(cmd, shell=True)

    rm_file = F"{consensus_files_concat}.out"
    # parse RepeatMasker output
    rm_annotation = get_repeatmasker_annotation(rm_file, seq_lengths, prefix)
    # add annotation to gff3 file, only if gff3 file exists
    if os.path.exists(gff):
        add_attribute_to_gff(gff, gff_out, "Name", "Annotation", rm_annotation)
    else:
        print(F"gff3 file {gff} does not exist, no annotation added to gff3 file")


def clustering(fasta, prefix, gff3=None, min_length=None, dust=False, cpu=4):
    """
    Run clustering on sequences defined in gff3 file and fasta file
    produce gff3 file with cluster information
    :param fasta: fasta file with sequences
    :param prefix: prefix - base naame for input and output files
    :param gff3: gff3 file with tidehunter results
    :param min_length: minimal length of repeat to be included in clustering
    :param cpu: number of cpu cores to use
    :return:

    """
    fasta = fasta
    if gff3 is None:
        gff3 = prefix + "_tidehunter.gff3"
    if min_length is not None:
        gff3 = filter_gff_by_length(gff3, min_length=min_length)
    gff3_out = prefix + "_clustering.gff3"

    # get consensus sequences for clustering
    consensus = {}
    consensus_dimers = {}
    for seq_id, seq, cons in gff3_to_fasta(gff3, fasta, "Consensus_sequence"):
        mult = round(1 + 10000 / len(cons))
        consensus[seq_id] = cons * mult
        consensus_dimers[seq_id] = cons * 4
        if len(consensus[seq_id]) > 10000:
            consensus[seq_id] = consensus[seq_id][0:10000]

    # run dustmasker first, sequences which are completely masked
    # will not be used in clastering.
    mask_prop = tc.get_ssrs_proportions(consensus_dimers)

    # if count mask_prop above 0.9
    ssrs_id = [k for k, v in mask_prop.items() if v > 0.9]
    # remove ssrs from consensus sequences, they will be added back later
    # but not used for clustering
    ssrs_description = {}
    for k in ssrs_id:
        consensus.pop(k, None)
        fdimer = consensus_dimers.pop(k, None)
        # for each ssrs caclulate ssr proportions
        ssrs_description[k] = tc.get_ssrs_description(fdimer)

    # first round of clustering by mmseqs2
    clusters1 = find_cluster_by_mmseqs2(consensus, cpu=cpu)

    representative_id = set(clusters1.values())
    consensus_representative = {k: consensus_dimers[k] for k in representative_id}
    # second round of clustering by blastn
    clusters2 = find_clusters_by_blast_connected_component(
        consensus_representative, dust=dust, cpu=cpu
        )
    # combine clusters
    clusters_final = clusters1.copy()

    for k, v in clusters1.items():
        if v in clusters2:
            clusters_final[k] = clusters2[v]
        else:
            clusters_final[k] = v

    # add ssrs_id back to clusters_final
    # and also to clusters1 - these are saved as well
    for k in ssrs_id:
        clusters_final[k] = k
        clusters1[k] = k

    # get total size of each cluster, store in dict
    cluster_size = get_cluster_size(gff3, clusters_final)

    # representative id sorted by cluster size
    representative_id = sorted(cluster_size, key=cluster_size.get, reverse=True)
    # rename values in clusters dictionary
    cluster_names = {}
    for i, v in enumerate(representative_id):
        cluster_names[v] = F"TRC_{i + 1}"

    ssrs_info = {}  # store ssrs info for gff3 file
    for k, v in clusters_final.items():
        clusters_final[k] = cluster_names[v]
        if k in ssrs_id:
            ssrs_info[cluster_names[v]] = "SSR"
        else:
            ssrs_info[cluster_names[v]] = "TR"
    #
    ssrs_description_final = {} # with TRC names
    for k, v in ssrs_description.items():
        ssrs_description_final[cluster_names[k]] = v

    cons_cls, cons_cls_dimer = add_cluster_info_to_gff3(gff3, gff3_out, clusters_final)

    merge_overlapping_gff3_intervals(gff3_out, gff3_out)

    gff_tmp = gff3_out + "_tmp"

    add_attribute_to_gff(gff3_out, gff_tmp, "Name", "repeat_type", ssrs_info)
    os.rename(gff_tmp, gff3_out)
    add_attribute_to_gff(gff3_out, gff_tmp, "Name", "SSR", ssrs_description_final)
    os.rename(gff_tmp, gff3_out)


    # save also first round of clustering for debugging
    cons_cls1, cons_cls_dimer1_ = add_cluster_info_to_gff3(
        gff3, gff3_out + "_1.gff3", clusters1
        )
    merge_overlapping_gff3_intervals(gff3_out + "_1.gff3", gff3_out + "_1.gff3")

    # for debugging
    # write consensus sequences by clusters to directory
    #  used gff3_out as base name for directory
    consensus_dir = prefix + "_consensus"
    save_consensus_files(consensus_dir, cons_cls, cons_cls_dimer)
    consensus_dir = prefix + "_consensus_1"
    save_consensus_files(consensus_dir, cons_cls1, cons_cls_dimer1_)


def tidehunter(fasta, tidehunter_arguments, prefix, cpu=4):
    """
    run tidehunter on fasta file
    :param fasta: fasta file with sequences
    :param tidehunter_arguments: tidehunter arguments
    :param prefix: prefix - base name for input and output files
    :param cpu: number of cpu cores to use
    :return:

    """
    # get size of input file
    chunk_size = 500000
    overlap = 50000
    output = prefix + "_tidehunter.gff3"
    output_chunks = prefix + "_chunks.bed"
    # check is tidehunter _arguments contain specification for
    # number of threads to use - format is -t <number>
    # if not add it from cpu variable
    if " -t" not in tidehunter_arguments:
        tidehunter_arguments += F" -t {cpu}"

    # this fill split sequences to chunk and all is stored in single file
    fasta_file_chunked, matching_table = tc.split_fasta_to_chunks(
        fasta, chunk_size, overlap
        )
    results = run_tidehunter(
        fasta_file_chunked, tidehunter_arguments
        )
    with open(output, "w") as out:
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

    with open(output_chunks, "w") as out:
        for m in matching_table:
            out.write(F'{m[0]}\t{m[2]}\t{m[3]}\t{m[4]}\n')


if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-c", "--cpu", type=int, default=4, help="Number of CPUs to use")

    subparsers = parser.add_subparsers(dest='command', help='TideHunter wrapper')

    # TideHunter
    parser_tidehunter = subparsers.add_parser(
        'tidehunter', help='Run wrapper of TideHunter'
        )

    parser_tidehunter.add_argument(
        '-f', '--fasta', type=str, required=True,
        help='Path to reference sequence in fasta format'
        )
    parser_tidehunter.add_argument(
        '-pr', '--prefix', type=str, required=True, help='Base name for output files'
        )
    # TideHunter parameters as single arguments in quotes
    parser_tidehunter.add_argument(
        '-T', '--tidehunter_arguments', type=str, nargs="?", required=False,
        default="-p 40 -P 3000 -c 5 -e 0.25",
        help=('additional arguments for TideHunter in quotes'
              ' default: -p 40 -P 3000 -c 5 -e 0.25')
        )

    # Clustering
    parser_clustering = subparsers.add_parser(
        'clustering', help='Run clustering on TideHunter output'
        )
    parser_clustering.add_argument(
        "-f", "--fasta", help="Reference fasta", required=True
        )

    parser_clustering.add_argument(
        "-m", "--min_length", help="Minimum length of tandem repeat", required=False,
        default=None, type=int
        )

    parser_clustering.add_argument(
        "-pr", "--prefix", help="Base name used for input and output files from ",
        required=True
        )
    parser_clustering.add_argument(
        "-g", "--gff", help=("gff output file from tidehunter ster. If not provided "
                             "the file named 'prefix_tidehunter.gff3' will be used"),
        required=False, default=None
        )
    parser_clustering.add_argument(
        '-d', '--dust', required=False, default=False, action='store_true',
        help='Use dust filter in blastn when clustering'
        )

    # Annotation
    parser_annotation = subparsers.add_parser(
        'annotation', help=('Run annotation on TideHunter output using reference '
                            'library of tandem repeats')
        )
    parser_annotation.add_argument(
        "-pr", "--prefix", help="Base name used for input and output files from ",
        required=True
        )

    parser_annotation.add_argument(
        "-g", "--gff", help=("gff output file from clustering step. If not provided "
                             "the file named 'prefix_clustering.gff3' will be used"),
        required=False, default=None
        )
    parser_annotation.add_argument(
        "-cd", "--consensus_directory",
        help=("Directory with consensus sequences which are to be "
              "annotated. If not provided the directory named 'prefix_consensus' "
              "will be used"), required=False, default=None
        )

    parser_annotation.add_argument(
        "-l", "--library", help="Path to library of tandem repeats", required=True, )

    # tarean
    parser_tarean = subparsers.add_parser(
        'tarean', help='Run TAREAN on custers to extract representative sequences'
        )
    parser_tarean.add_argument(
        "-g", "--gff", help="Base name used for input and output files from ",
        required=False, default=None
        )
    parser_tarean.add_argument(
        "-f", "--fasta", help="Reference fasta", required=True
        )
    parser_tarean.add_argument(
        "-pr", "--prefix", help="Base name used for input and output files from ",
        required=True
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
        tidehunter(
            cmd_args.fasta, cmd_args.tidehunter_arguments, cmd_args.prefix, cmd_args.cpu
            )
    elif cmd_args.command == "clustering":
        clustering(
            cmd_args.fasta, cmd_args.prefix, cmd_args.gff, cmd_args.min_length,
            cmd_args.dust, cmd_args.cpu
            )
    elif cmd_args.command == "annotation":
        annotation(
            cmd_args.prefix, cmd_args.library, cmd_args.gff, cmd_args.consensus_directory,
            cmd_args.cpu
            )
    elif cmd_args.command == "tarean":
        tarean(
            cmd_args.prefix, cmd_args.gff, cmd_args.fasta, cmd_args.cpu
            )
    else:
        parser.print_help()
        sys.exit(1)
