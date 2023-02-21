#!/usr/bin/env python
"""Wrapper of TideHunter
input is a fasta file of DNA sequences
output is a table of detected tandem repeats by tidehunter
this wrapper split DNA sequence to chunks, run tidehunter on each chunk and merge results
coordinates in the table must be recalculated to the original sequence

"""
import os
import sys
import textwrap

import tc_utils
import argparse
import subprocess
import tempfile
from shutil import rmtree
from tc_utils import Gff3Feature

# minimal python version is 3.6
assert sys.version_info >= (3, 6), "Python 3.6 or newer is required"


# clustering related functions:
def fasta_to_dict(fasta_file):
    """
    convert fasta file to dictionary
    :param fasta_file: path to fasta file
    :return: dictionary with fasta sequences
    """
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        for line in f:

            if line.startswith(">"):
                seq_name = line.split()[0].replace(">", "")
                fasta_dict[seq_name] = [""]  # initialize list
            else:
                fasta_dict[seq_name].append(line.strip())
    # concatenate list in dictionary
    fasta_dict = {k: "".join(v) for k, v in fasta_dict.items()}
    return fasta_dict


def get_seq_from_fasta(fasta_dict, seq_id, start=None, end=None):
    """
    get sequence from fasta dictionary
    :param fasta_dict: dictionary with fasta sequences
    :param seq_id: sequence ID
    :param start: start position
    :param end: end position
    :return: sequence
    """
    if start and end:
        return fasta_dict[seq_id][start:end]
    return fasta_dict[seq_id]


def gff3_to_fasta(gff3_file, fasta_file):
    """
    extract fasta sequences from gff3 file
    it is generator, returns one sequence at time and seq ID
    :param gff3_file: path to gff3 file
    :param fasta_file: path to fasta file
    :return:
    """

    fasta_dict = fasta_to_dict(fasta_file)
    with open(gff3_file, 'r') as f1:
        for line in f1:
            if line.startswith("#"):
                continue
            gff3_feature: Gff3Feature = Gff3Feature(line)
            s = get_seq_from_fasta(
                fasta_dict, gff3_feature.seqid, gff3_feature.start, gff3_feature.end
                )
            yield [gff3_feature.attributes_dict['ID'], s,
                   gff3_feature.attributes_dict['Consensus_sequence']]


def save_fasta_dict_to_file(fasta_dict, fasta_file):
    """
    save fasta dictionary to file
    :param fasta_dict: dictionary with fasta sequences
    :param fasta_file: path to fasta file
    :return:
    """
    with open(fasta_file, 'w') as f:
        for k, v in fasta_dict.items():
            f.write(">{}\n{}\n".format(k, v))


# run mmseqs2 on consensus sequences
def find_cluster_by_mmseqs2(sequences):
    """
    run mmseqs2 on consensus sequences
    :param sequences:
    :return: clusters
    """
    tmp_dir = tempfile.mkdtemp()
    input_fasta_file = tempfile.NamedTemporaryFile(mode='w', delete=False)
    save_fasta_dict_to_file(sequences, input_fasta_file.name)

    cmd = (F'mmseqs easy-cluster {input_fasta_file.name} {input_fasta_file.name}.clu'
           F' {tmp_dir} --cluster-mode 1 '
           F'--mask 0  ')

    subprocess.check_call(cmd, shell=True)

    # read clusters to dictionary
    with open(F"{input_fasta_file.name}.clu_cluster.tsv", 'r') as f:
        clusters = {}
        for line in f:
            cluster, seq_id = line.split()
            clusters[seq_id] = cluster
    # remove temporary files
    os.remove(F"{input_fasta_file.name}.clu_all_seqs.fasta")
    os.remove(F"{input_fasta_file.name}.clu_cluster.tsv")
    os.remove(F"{input_fasta_file.name}.clu_rep_seq.fasta")
    rmtree(tmp_dir)
    os.remove(input_fasta_file.name)
    return clusters


def make_graph(pairs):
    """
    make graph from list of pairs
    :param pairs:
    :return:
    """
    graph = {}
    for pair in pairs:
        v1, v2 = pair
        if v1 not in graph:
            graph[v1] = set()
        if v2 not in graph:
            graph[v2] = set()
        graph[v1].add(v2)
        graph[v2].add(v1)
    return graph


def find_connected_components(graph):
    """
    find connected components in graph
    :param graph:
    :return: components
    """
    visited = set()  # create a set to keep track of visited vertices
    components = []  # create an empty list to keep track of components

    # define a recursive function to perform DFS
    def dfs(vertex, component):
        """
        perform depth-first search on the graph to find all connected components
        :param vertex:
        :param component:
        """
        visited.add(vertex)
        component.append(vertex)
        for neighbor in graph[vertex]:
            if neighbor not in visited:
                dfs(neighbor, component)

    # iterate over all vertices to perform DFS
    for vertex in graph:
        if vertex not in visited:
            component = []
            dfs(vertex, component)
            components.append(component)

    return components


def run_blastn(sequences):
    """
    run blastn on fasta file
    :param sequences : dictionary with sequences
    :return: dictionary with clusters
    """
    tmp_dir = tempfile.mkdtemp()
    fasta_file = F'{tmp_dir}/seqs.fasta'
    blast_out = F'{tmp_dir}/blastn.out'
    save_fasta_dict_to_file(sequences, fasta_file)
    # make blast database
    cmd = F"makeblastdb -in {tmp_dir}/seqs.fasta -dbtype nucl"
    subprocess.check_call(cmd, shell=True)
    # run blastn

    outfmt = "'6 qseqid sseqid pident length evalue bitscore qlen slen'"

    cmd = (F"blastn -query {fasta_file} -db {fasta_file} -outfmt {outfmt}"
           F" -out {blast_out} -num_threads 5 -evalue 1e-5 -perc_identity 80 "
           F"-word_size 21")

    subprocess.check_call(cmd, shell=True)
    # read pairs to list, exclude self hits and duplicates
    pairs = set()
    with open(blast_out, 'r') as f:
        for line in f:
            qseqid, sseqid, pident, length, evalue, bitscore, qlen, slen = line.split()
            if qseqid != sseqid:
                # overlap should be at least 50% over the shorter sequence
                if min([int(qlen), int(slen)]) / int(length) > 0.5:
                    pairs.add(tuple(sorted([qseqid, sseqid])))
    rmtree(tmp_dir)
    return pairs


def find_clusters_by_blast(consensus_representative):
    """
    find clusters by blastn, return dictionary with clusters
    cluaste are connected components in graph
    :param consensus_representative:
    :return: clusters
    """
    pairs = run_blastn(consensus_representative)
    graph = make_graph(pairs)
    components = find_connected_components(graph)

    clusters = {}
    for vertices in components:
        v_representative = sorted(vertices)[0]
        for v in vertices:
            clusters[v] = v_representative

    return clusters


def merge_overlapping_gff3_intervals(gff3_file, gff3_out_file):
    """
    merge overlapping intervals in gff3 file
    merge only if they have same Cluster_ID
    :param gff3_file: path to gff3 file
    :param gff3_out_file: path to output gff3 file
    :return: nothing
    """
    # read gff3 file, split to lists by cluster ID and seqname
    gff3_dict = {}
    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            gff3_feature = Gff3Feature(line)
            cluster_id = gff3_feature.attributes_dict['Cluster_ID']
            seqname = gff3_feature.seqid
            seq_cluster_id = (seqname, cluster_id)
            if seq_cluster_id not in gff3_dict:
                gff3_dict[seq_cluster_id] = []
            gff3_dict[seq_cluster_id].append(gff3_feature)
    with open(gff3_out_file, 'w') as f:
        f.write("##gff-version 3\n")
        for i in gff3_dict:
            intervals = []
            for j in gff3_dict[i]:
                intervals.append((j.start, j.end))
            intervals = merge_intervals(intervals)
            for start, end in intervals:
                gff_line = (F'{i[0]}\tTideCluster\trepeat_region'
                            F'\t{start}\t{end}\t{1}\t.\t.\tName={i[1]}\n')

                f.write(gff_line)


def merge_intervals(intervals):
    """
    :param intervals:list of (start, end) tuples
    :return:
    list of merged (start, end) tuples
    """
    # Sort intervals by start point
    intervals.sort()
    # Initialize the merged intervals list with the first interval
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        # Check if the current interval overlaps with the last merged interval
        if start <= merged[-1][1] + 1:
            # Merge the intervals by extending the end point of the last interval
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            # Add the current interval to the list of merged intervals
            merged.append((start, end))
    return merged


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
        consensus_dimers[seq_id] = cons * 2
        if len(consensus[seq_id]) > 10000:
            consensus[seq_id] = consensus[seq_id][0:10000]

    # first round of clustering by mmseqs2
    clusters1 = find_cluster_by_mmseqs2(consensus)
    representative_id = set(clusters1.values())
    consensus_representative = {k: consensus_dimers[k] for k in representative_id}

    # second round of clustering by blastn
    clusters2 = find_clusters_by_blast(consensus_representative)
    # combine clusters
    clusters_final = clusters1.copy()
    for k, v in clusters1.items():
        if v in clusters2:
            clusters_final[k] = clusters2[v]
        else:
            clusters_final[k] = v

    # rename values in clusters dictionary
    representative_id = sorted(set(clusters_final.values()))
    cluster_names = {}
    for i, v in enumerate(representative_id):
        cluster_names[v] = F"TRC_{i}"

    for k, v in clusters_final.items():
        clusters_final[k] = cluster_names[v]

    # read gff3 and append cluster ID
    with open(gff3, 'r') as f1, open(gff3_out, 'w') as f2:
        for line in f1:
            if line.startswith("#"):
                f2.write(line)
                continue
            gff3_feature = Gff3Feature(line)
            cluster_id = clusters_final[gff3_feature.attributes_dict['ID']]
            gff3_feature.attributes_dict['Cluster_ID'] = cluster_id
            gff3_feature.attributes_dict['Name'] = cluster_id
            unique_id = cluster_id + "_" + gff3_feature.attributes_dict['ID']
            gff3_feature.attributes_dict['ID'] = unique_id
            f2.write(gff3_feature.print_line())

    merge_overlapping_gff3_intervals(gff3_out, gff3_out)


# TideHunter related functions:
def run_tidehunter(fasta_file, tidehunter_arguments):
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
    tidehunter_cmd = (F"TideHunter -f 2 -o {tmp_file} {tidehunter_arguments}"
                      F" {fasta_file}")

    print("running TideHunter")
    print(tidehunter_cmd)
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
    TideCluster.py tidehunter -f test.fasta -o test.gff3 -T "-p 40 -P 3000 -c 5 -e 0.25 -t 46"
    
    # then run clustering on the output to cluster similar tantem repeats
    TideCluster.py clustering -f test.fasta -g test.gff3 -o test_clustered.gff3
    
    
    Recommended parameters for TideHunter:
    short monomers: -T "-p 10 -P 39 -c 5 -e 0.25"
    long monomers: -T "-p 40 -P 3000 -c 5 -e 0.25"
    
    For parallel processing include -t option. 
    
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
