""" Utilities for the project. """
import os
import subprocess
import tempfile
from itertools import cycle, permutations
from shutil import rmtree
import networkx as nx

def get_kmers(dna, k):
    """ Return all kmers of length k from dna.
    param dna: dict of dna sequences
    param k: length of kmers
    return: dictionary kmer:count
    """
    kmers = {}
    for seq in dna.values():
        for i in range(len(seq)-k+1):
            kmer = seq[i:i+k]
            if kmer in kmers:
                kmers[kmer] += 1
            else:
                kmers[kmer] = 1
    return kmers

def kmers2debruijn(kmers):
    """ Return debruijn graph from kmers.
    param kmers: list of kmers
    return: dictionary of nodes
    """
    graph = nx.DiGraph()
    for kmer in kmers:
        edge_weight = kmers[kmer]
        graph.add_edge(kmer[:-1], kmer[1:], weight=edge_weight)
    return graph


def smallest_circular_paths(G):

    cycles = nx.simple_cycles(G)
    cycles = sorted(
        cycles, key=lambda cycle: sum(
            G[u][v]['weight'] for u, v in zip(cycle, cycle[1:] + cycle[:1])
            ), reverse=True
        )
    paths = []
    while cycles:
        cycle = cycles.pop(0)
        H = G.copy()
        H.remove_edges_from(zip(cycle, cycle[1:] + cycle[:1]))
        H.remove_nodes_from([v for v in cycle if H.degree[v] == 0])
        while True:
            subcycles = nx.simple_cycles(H)
            print(len(subcycles))
            if not subcycles:
                break
            subcycles = sorted(
                subcycles, key=lambda subcycle: sum(
                    H[u][v]['weight'] for u, v in
                    zip(subcycle, subcycle[1:] + subcycle[:1])
                    ), reverse=True
                )
            subcycle = subcycles[0]
            H.remove_edges_from(zip(subcycle, subcycle[1:] + subcycle[:1]))
            H.remove_nodes_from([v for v in subcycle if H.degree[v] == 0])
            cycles.insert(0, subcycle)
        paths.append(cycle)
    if len(paths) == 1:
        return paths
    else:
        assigned = set()
        result = []
        for path in paths:
            unassigned = set(path) - assigned
            result.append([v for v in path if v in unassigned])
            assigned |= unassigned
        return result

def debruijn2contigs(graph):
    """ Return contigs from debruijn graph.
    param graph: debruijn graph
    return: list of contigs
    """
    contigs_kmers = {}
    # get list of nodes as set:
    nodes = set(list(graph.nodes()))
    while nodes:
        node = nodes.pop()
        neighbor = next(graph.neighbors(node))
        path = nx.shortest_path(graph, neighbor, node, weight='weight')
        frozenset_path = frozenset(sorted(path))
        contigs_kmers[frozenset_path] = path
        nodes -= frozenset_path
    return contigs_kmers

def contigs_kmers2strings(contigs_kmers):
    """ Return contigs from contigs_kmers.
    param contigs_kmers: dictionary of contigs
    return: list of contigs
    """
    if type(contigs_kmers) is dict:
        contigs = []
        for contig in contigs_kmers.values():
            contig_string = contig[0]
            for kmer in contig[1:]:
                contig_string += kmer[-1]
            contigs.append(contig_string)
    else:
        contigs = []
        for contig in contigs_kmers:
            contig_string = contig[0]
            for kmer in contig[1:]:
                contig_string += kmer[-1]
            contigs.append(contig_string)

    return contigs





class Gff3Feature:
    """
    Class for gff3 feature
    """

    def __init__(self, line):
        self.line = line
        self.items = line.strip().split('\t')
        self.seqid = self.items[0]
        self.source = self.items[1]
        self.type = self.items[2]
        self.start = int(self.items[3])
        self.end = int(self.items[4])
        self.score = self.items[5]
        self.strand = self.items[6]
        self.frame = self.items[7]
        self.attributes = self.items[8]
        self._attributes_dict = {}
        for item in self.attributes.split(';'):
            if item != '':
                key, value = item.split('=')
                self._attributes_dict[key] = value

        self._attributes_str = ';'.join(
            ['{}={}'.format(key, value) for key, value in self._attributes_dict.items()]
            )

    @property
    def attributes_dict(self):
        """
        store attributes as dictionary
        :return:
        """
        return self._attributes_dict

    @property
    def attributes_str(self):
        """
        store attributes as string
        :return:
        """
        return self._attributes_str

    @attributes_str.getter
    def attributes_str(self):
        """
        store attributes as string
         :return:
        """
        self._attributes_str = ';'.join(
            ['{}={}'.format(key, value) for key, value in self._attributes_dict.items()]
            )
        return self._attributes_str

    @attributes_dict.setter
    def attributes_dict(self, value):
        self._attributes_dict = value
        self._attributes_str = ';'.join(
            ['{}={}'.format(key, value) for key, value in self._attributes_dict.items()]
            )

    def __str__(self):
        return '\t'.join(
            [self.seqid, self.source, self.type, str(self.start), str(self.end),
             self.score, self.strand, self.frame, self.attributes_str]
            ) + '\n'

    def __repr__(self):
        return '\t'.join(
            [self.seqid, self.source, self.type, str(self.start), str(self.end),
             self.score, self.strand, self.frame, self.attributes_str]
            ) + '\n'

    def __eq__(self, other):
        return self.line_recalculated() == other.line_recalculated()

    def __hash__(self):
        return hash(self.line_recalculated())

    def get_line(self):
        """returns original line"""
        return self.line

    def overlap(self, other):
        """
        Check if two features overlap
        :param other:
        :return:
        """
        if self.start <= other.end and self.end >= other.start:
            return True
        else:
            return False

    def line_recalculated(self):
        """
        :return:
        string with recalculated line
        """
        return '\t'.join(
            [self.seqid, self.source, self.type, str(self.start), str(self.end),
             self.score, self.strand, self.frame, self.attributes_str]
            ) + '\n'

    def __lt__(self, other):
        width = self.end - self.start
        other_width = other.end - other.start
        return width < other_width

    def __gt__(self, other):
        width = self.end - self.start
        other_width = other.end - other.start
        return width > other_width

    def identical_region(self, other):
        """
        Check if two features are identical
        :param other:
        :return:
        """
        if self.start == other.start and self.end == other.end and self.seqid == \
                other.seqid:
            return True
        else:
            return False

    def print_line(self):
        """
        :return:
        string with recalculated line
        """
        columns = [self.seqid, self.source, self.type, str(self.start), str(self.end),
                   self.score, self.strand, self.frame]
        attributes_list = ['{}={}'.format(key, value) for key, value in
                           self.attributes_dict.items()]
        attributes = [';'.join(attributes_list)]
        return '\t'.join(columns + attributes) + '\n'


class RepeatMaskerFeature:
    """
    class for parsing repeatmasker ouput from .out file
    """

    def __init__(self, line):
        items = line.split()
        if len(items) < 12:
            raise ValueError('Line does not contain enough columns')
        self.seqid = items[4]
        self.start = int(items[5])
        self.end = int(items[6])
        self.length = self.end - self.start
        self.strand = items[8]
        self.annot = items[10]
        self.refid = items[9]
        self.family = items[11]


def get_repeatmasker_annotation(rm_file, seq_lengths, prefix):
    """
    :parse repeatmasker output and calculate proportion of each annotation
    :param rm_file:
    :param seq_lengths:
    :param prefix:
    :return:
    """

    seq_rm_info = {}
    with open(rm_file, 'r') as f:
        # parse repeatmasker output, first three lines are header
        for i in range(3):
            next(f)
        for line in f:
            rm_feature = RepeatMaskerFeature(line)
            if rm_feature.seqid not in seq_rm_info:
                seq_rm_info[rm_feature.seqid] = {}
                if rm_feature.annot not in seq_rm_info[rm_feature.seqid]:
                    seq_rm_info[rm_feature.seqid][rm_feature.annot] = 0
                prop = rm_feature.length / seq_lengths[rm_feature.seqid]
                seq_rm_info[rm_feature.seqid][rm_feature.annot] += prop
    # for each TRC calculate mean value for each annotation
    # not all sequences are necessarily annotated
    annot_summary = {}
    trc_consensus_count = {}  # how many consensus sequences are in each TRC
    # scan all consensus sequences id
    for seqid in seq_lengths:
        trc_id = "TRC_" + seqid.split("_")[1]
        if trc_id not in trc_consensus_count:
            trc_consensus_count[trc_id] = 0
        trc_consensus_count[trc_id] += 1
        if trc_id not in annot_summary:
            annot_summary[trc_id] = {}
        # does seqid have any annotations
        if seqid in seq_rm_info:
            for annot in seq_rm_info[seqid]:
                if annot not in annot_summary[trc_id]:
                    annot_summary[trc_id][annot] = 0
                annot_summary[trc_id][annot] += seq_rm_info[seqid][annot]

    # divide by number of consensus sequences in each TRC
    for trc_id in annot_summary:
        for annot in annot_summary[trc_id]:
            annot_summary[trc_id][annot] /= trc_consensus_count[trc_id]

    # export annotation to table:
    annot_table = prefix + "_annotation.tsv"
    annot_description = {}
    with open(annot_table, "w") as f:
        f.write("TRC\tannotation\tproportion_coverage\n")
        for trc_id in annot_summary:
            annot_description[trc_id] = "NA"
            if len(annot_summary[trc_id]) == 0:
                continue
            # sort keys by values
            key_sorted = sorted(
                annot_summary[trc_id], key=annot_summary[trc_id].get, reverse=True
                )
            for i, v in enumerate(key_sorted):
                p = round(annot_summary[trc_id][v], 3)
                if i > 0:
                    # round to 2 decimal places
                    annot_description[trc_id] += F"; {v} ({p})"
                    f.write(F"\t{v}\t{annot_summary[trc_id][v]}")
                else:
                    annot_description[trc_id] = F"{v} ({p})"
                    f.write(F"{trc_id}\t{v}\t{annot_summary[trc_id][v]}")
            f.write("\n")
    print("Annotation exported to: " + annot_table)
    return annot_description


def read_fasta_sequence_size(fasta_file):
    """Read size of sequence into dictionary"""
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line[0] == '>':
                header = line.strip().split(' ')[0][1:]  # remove part of name after space
                fasta_dict[header] = 0
            else:
                fasta_dict[header] += len(line.strip())
    return fasta_dict


def read_single_fasta_to_dictionary(fh):
    """
    Read fasta file into dictionary
    :param fh: - file handle
    :return:
    fasta_dict
    """
    fasta_dict = {}
    for line in fh:
        if line[0] == '>':
            header = line.strip().split(' ')[0][1:]  # remove part of name after space
            fasta_dict[header] = []
        else:
            fasta_dict[header] += [line.strip()]
    fasta_dict = {k: ''.join(v) for k, v in fasta_dict.items()}
    return fasta_dict


def split_fasta_to_chunks(fasta_file, chunk_size=100000000, overlap=100000):
    """
    Split fasta file to chunks, sequences longer than chunk size are split to overlaping
    pieces. If sequences are shorter, chunk with multiple sequences are created.
    :param fasta_file:

    :param fasta_file:
    :param chunk_size:
    :param overlap:
    :return:
    fasta_file_split:
    matching_table
    """
    min_chunk_size = chunk_size * 2
    fasta_dict = read_fasta_sequence_size(fasta_file)
    # calculates ranges for splitting of fasta files and store them in list
    matching_table = []
    fasta_file_split = tempfile.NamedTemporaryFile(delete=False).name
    for header, size in fasta_dict.items():
        if size > min_chunk_size:
            number_of_chunks = int(size / chunk_size)
            adjusted_chunk_size = int(size / number_of_chunks)
            for i in range(number_of_chunks):
                start = i * adjusted_chunk_size
                end1 = (i + 1) * adjusted_chunk_size + overlap
                end = end1 if i + 1 < number_of_chunks else size
                new_header = header + '_' + str(i)
                matching_table.append([header, i, start, end, new_header])
        else:
            new_header = header + '_0'
            matching_table.append([header, 0, 0, size, new_header])
    # read sequences from fasta files and split them to chunks according to matching table
    # open output and input files, use with statement to close files
    fasta_dict = read_single_fasta_to_dictionary(open(fasta_file, 'r'))
    with open(fasta_file_split, 'w') as fh_out:
        for header in fasta_dict:
            matching_table_part = [x for x in matching_table if x[0] == header]
            for header2, i, start, end, new_header in matching_table_part:
                fh_out.write('>' + new_header + '\n')
                fh_out.write(fasta_dict[header][start:end] + '\n')
    return fasta_file_split, matching_table


def make_temp_files(number_of_files):
    """
    Make named temporary files, file will not be deleted upon exit!
    :param number_of_files:
    :return:
    filepaths
    """
    temp_files = []
    for i in range(number_of_files):
        temp_files.append(tempfile.NamedTemporaryFile(delete=False).name)
        os.remove(temp_files[-1])
    return temp_files


def get_new_header_and_coordinates(header, start, end, matching_table):
    """
    Get new header and coordinates for sequence
    :param header:
    :param start:
    :param end:
    :param matching_table:
    :return:
    new_header
    new_start
    new_end
    """
    matching_table_part = [x for x in matching_table if x[0] == header]
    new_coords = []
    for chunk in matching_table_part:
        if chunk[2] <= start < chunk[3]:
            new_header = chunk[4]
            new_start = start - chunk[2]
            new_end = end - chunk[2]
            new_sequence_length = chunk[3] - chunk[2]
            new_coords.append([new_header, new_start, new_end, new_sequence_length])
    return new_coords


def get_original_header_and_coordinates(new_header, new_start, new_end, matching_table):
    """
    Get original header and coordinates for sequence
    :param new_header:
    :param new_start:
    :param new_end:
    :param matching_table:
    :return:
    original_header
    original_start
    original_end
    """
    matching_table_part = [x for x in matching_table if x[4] == new_header]
    real_chunk_size = matching_table_part[0][3] - matching_table_part[0][2]
    ori_header = matching_table_part[0][0]
    start = matching_table_part[0][2]
    ori_start = new_start + start
    ori_end = new_end + start
    return ori_header, ori_start, ori_end, real_chunk_size


# recalculate gff3 coordinates, use gff3_feature class
def recalculate_gff3_coordinates(gff3_file, matching_table):
    """
    Recalculate gff3 coordinates, use gff3_feature class
    :param gff3_file:
    :param matching_table:
    :return:
    gff3_file_recalculated
    """
    gff3_file_recalculated = tempfile.NamedTemporaryFile(delete=False).name

    with open(gff3_file, 'r') as fh_in:
        with open(gff3_file_recalculated, 'w') as fh_out:
            for line in fh_in:
                if line[0] == '#':
                    fh_out.write(line)
                else:
                    feature = Gff3Feature(line)
                    new_coords = get_new_header_and_coordinates(
                        feature.seqid, feature.start, feature.end, matching_table
                        )
                    for new_header, new_start, new_end, sequence_length in new_coords:
                        if new_start >= 1 and new_end <= sequence_length:
                            feature.seqid = new_header
                            feature.start = new_start
                            feature.end = new_end
                            fh_out.write(str(feature))
    return gff3_file_recalculated


def write_temp_fasta_chunks(fasta_seq_size, fasta_file, chunk_size):
    """
    Write temporary fasta chunks
    :param fasta_seq_size: dictionary of fasta sequences and their sizes
    :param fasta_file: path to fasta file
    :param chunk_size: int size of chunk
    :return:
    temp_files_fasta: list of temporary fasta files

    input fasta file is deleted after this function!
    """
    number_of_chunks = len(fasta_seq_size)
    seq_id_size_sorted = [i[0] for i in sorted(
        fasta_seq_size.items(), key=lambda x: int(x[1]), reverse=True
        )]
    number_of_temp_files = int(os.path.getsize(fasta_file) / chunk_size) + 1
    if number_of_temp_files > number_of_chunks:
        number_of_temp_files = number_of_chunks
    temp_files_fasta = make_temp_files(number_of_temp_files)
    seq_id_file_path_dict = dict(zip(seq_id_size_sorted, cycle(temp_files_fasta)))
    # write sequences to temporary files
    with open(fasta_file, 'r') as f:
        for line in f:
            if line[0] == '>':
                header = line.strip().split(' ')[0][1:]
                # append to file
                with open(seq_id_file_path_dict[header], 'a') as fout:
                    fout.write(line)
            else:
                with open(seq_id_file_path_dict[header], 'a') as fout:
                    fout.write(line)
    os.remove(fasta_file)
    return temp_files_fasta


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
        ori_name, ori_start, ori_end, _ = get_original_header_and_coordinates(
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
def find_cluster_by_mmseqs2(sequences, cpu=4):
    """
    run mmseqs2 on consensus sequences
    :param cpu:
    :param sequences:
    :return: clusters
    """
    tmp_dir = tempfile.mkdtemp()
    input_fasta_file = tempfile.NamedTemporaryFile(mode='w', delete=False)
    save_fasta_dict_to_file(sequences, input_fasta_file.name)

    cmd = (F'mmseqs easy-cluster {input_fasta_file.name} {input_fasta_file.name}.clu'
           F' {tmp_dir} --cluster-mode 0 '
           F'--mask 0  -s 1 --threads {cpu}')
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


def get_connected_component_clusters(pairs):
    """
    find connected components in graph
    :param pairs: list of pairs
    :return: clusters as list of lists
    """
    # Create a NetworkX graph from the list of edges
    g = nx.Graph()
    g.add_edges_from(pairs)

    # Find the connected components in the graph
    components = list(nx.connected_components(g))

    # Convert the components to a list of clusters
    clusters = [list(component) for component in components]

    return clusters


def run_blastn(sequences, dust = False, cpu=4):
    """
    run blastn on fasta file
    :param cpu:
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
    # blast parameters:
    dust = "yes" if dust else "no"
    outfmt = "'6 qseqid sseqid pident length evalue bitscore qlen slen'"
    # run blastn
    cmd = (F"blastn -task blastn -query {fasta_file} -db {fasta_file} -outfmt {outfmt}"
           F" -out {blast_out} -num_threads {cpu} -evalue 1e-20 -perc_identity 75"
           F" -word_size 9 -max_target_seqs 1000000 -dust {dust}"
           F" -gapextend 1 -gapopen 2 -reward 1 -penalty -1")
    subprocess.check_call(cmd, shell=True)
    # read pairs to list, exclude self hits and duplicates
    pairs = set()
    with open(blast_out, 'r') as f:
        for line in f:
            qseqid, sseqid, pident, length, evalue, bitscore, qlen, slen = line.split()
            if qseqid != sseqid:
                # overlap should be at least 50% over the shorter sequence
                if int(length) / min([int(qlen), int(slen)]) > 0.8:
                    pairs.add(tuple(sorted([qseqid, sseqid])))
    # add self hits separately, so they are in the graph later
    for seq_id in sequences:
        pairs.add((seq_id, seq_id))

    rmtree(tmp_dir)
    return pairs


def find_clusters_by_blast_connected_component(consensus_representative, dust=False,
                                               cpu=4):
    """
    find clusters by blastn, return dictionary with clusters
    cluaste are connected components in graph
    :param cpu:
    :param consensus_representative:
    :return: clusters
    """
    pairs = run_blastn(consensus_representative, dust=dust, cpu=cpu)
    # graph = make_graph(pairs)
    # components = find_connected_components(graph)
    components = get_connected_component_clusters(pairs)
    clusters = {}
    for vertices in components:
        v_representative = sorted(vertices)[0]
        for v in vertices:
            clusters[v] = v_representative

    return clusters


def get_cluster_size(fin, clusters):
    """
    get cluster size
    :param fin: input gff
    :param clusters: clusters dictionary id:cluster_id
    :return: cluster size as total length of intervals
    """
    cluster_size = {}
    with open(fin, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            gff3_feature = Gff3Feature(line)
            cluster_id = clusters[gff3_feature.attributes_dict['ID']]
            if cluster_id not in cluster_size:
                cluster_size[cluster_id] = 0
            cluster_size[cluster_id] += gff3_feature.end - gff3_feature.start
    return cluster_size


def add_attribute_to_gff(fin, fout, attr2match, new_attr, attr_dict):
    """
    add attribute to gff file
    :param fin: input gff
    :param fout: output gff
    :param attr2match: which attribute to use to match with attribute_dict
    :param new_attr: name of new attribute to be added to gff
    :param attr_dict: dictionary with attributes
    :return:
    """
    with open(fin, 'r') as f1, open(fout, 'w') as f2:
        for line in f1:
            if line.startswith("#"):
                f2.write(line)
                continue
            gff3_feature = Gff3Feature(line)
            attr_value = gff3_feature.attributes_dict[attr2match]
            if attr_value in attr_dict:
                gff3_feature.attributes_dict[new_attr] = attr_dict[attr_value]
            f2.write(str(gff3_feature))


def add_cluster_info_to_gff3(fin, fout, clusters):
    """
    add cluster info to gff3 file
    :param fin: input gff
    :param fout: output gff
    :param clusters: clusters dictionary id:cluster_id
    :return: consensus and consensus dimers of representative sequences
    """
    consensus_clusters = {}
    consensus_clusters_dimers = {}
    with open(fin, 'r') as f1, open(fout, 'w') as f2:
        for line in f1:
            if line.startswith("#"):
                f2.write(line)
                continue
            gff3_feature = Gff3Feature(line)
            cluster_id = clusters[gff3_feature.attributes_dict['ID']]
            gff3_feature.attributes_dict['Cluster_ID'] = cluster_id
            gff3_feature.attributes_dict['Name'] = cluster_id
            unique_id = cluster_id + "_" + gff3_feature.attributes_dict['ID']
            gff3_feature.attributes_dict['ID'] = unique_id
            # add consensus sequence to dictionary of dictionaries
            consensus = gff3_feature.attributes_dict['Consensus_sequence']
            if cluster_id not in consensus_clusters:
                consensus_clusters[cluster_id] = {}
                consensus_clusters_dimers[cluster_id] = {}
            consensus_clusters[cluster_id][unique_id] = consensus
            consensus_clusters_dimers[cluster_id][unique_id] = consensus * 2

            f2.write(gff3_feature.print_line())
    return consensus_clusters, consensus_clusters_dimers


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


def filter_gff_by_length(gff3_file, min_length=1000):
    """
    Filter gff3 file by size of intervals
    :param gff3_file:
    :param min_length:
    :return: filtered gff3 file path
    """
    gff_out = tempfile.NamedTemporaryFile(delete=False).name
    with open(gff3_file, 'r') as f1, open(gff_out, 'w') as f2:
        for line in f1:
            if line.startswith("#"):
                f2.write(line)
                continue
            gff3_feature = Gff3Feature(line)
            if gff3_feature.end - gff3_feature.start > min_length:
                f2.write(line)
    return gff_out


def save_consensus_files(consensus_dir, cons_cls, cons_cls_dimer):
    """
    Save consensus sequences to fasta files
    :param consensus_dir: directory to save consensus sequences
    :param cons_cls: consensus sequences for each cluster
    :param cons_cls_dimer: consensus sequences for each cluster, dimers
    """
    if not os.path.exists(consensus_dir):
        os.makedirs(consensus_dir)
    for cluster_id in cons_cls:
        f = os.path.join(consensus_dir, cluster_id + ".fasta")
        save_fasta_dict_to_file(cons_cls[cluster_id], f)
        f = os.path.join(consensus_dir, cluster_id + "_dimers.fasta")
        save_fasta_dict_to_file(cons_cls_dimer[cluster_id], f)
