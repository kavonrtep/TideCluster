""" Utilities for the project. """
import glob
import gzip
import itertools
import os
import shutil
import subprocess
import tempfile
from collections import OrderedDict
from itertools import cycle
from shutil import rmtree
import re
import networkx as nx

class GRange:
    def __init__(self, seqid, start, end, name, strand):
        self.seqid = seqid
        self.start = start
        self.end = end
        self.name = name
        self.strand = strand
    def __repr__(self):
        return f"GRange(seqid={self.seqid}, start={self.start}, end={self.end}, name={self.name}, strand={self.strand})"
    def print_bed(self):
        return f"{self.seqid}\t{self.start}\t{self.end}\t{self.name}\t0\t{self.strand}"

    def print_gff3(self):
        return f"{self.seqid}\tRepeatMasker\trepeat\t{self.start}\t{self.end}\t.\t{self.strand}\t.\tName={self.name}"

    def overlap(self, other):
        return self.seqid == other.seqid and self.start <= other.end and self.end >= other.start

    def within(self, other):
        return self.seqid == other.seqid and self.start >= other.start and self.end <= other.end


def read_gff3_to_grange_list(gff3_file):
    granges = []
    with open(gff3_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            gff_record = Gff3Feature(line)
            granges.append(GRange(gff_record.seqid,
                                  gff_record.start,
                                  gff_record.end,
                                  gff_record.attributes_dict["Name"],
                                  gff_record.strand))
    return granges

def split_intervals(granges):
    events = []
    for grange in granges:
        events.append((grange.seqid, grange.start, 1, grange))  # 1 for start
        events.append((grange.seqid, grange.end, -1, grange))  # -1 for end
    # Sort events by seqid, position, and type (start before end)
    events.sort(key=lambda x: (x[0], x[1], -x[2]))
    split = []
    active_granges = set()
    prev_pos = None
    prev_seqid = None
    count = 0
    last_end = None
    for seqid, pos, event_type, grange in events:
        # Reset active granges and counters for a new seqid
        if prev_seqid and seqid != prev_seqid:
            active_granges.clear()
            count = 0
            last_end = None
        if count > 0 and prev_pos is not None:
            end_pos = pos if event_type == -1 else pos - 1
            if last_end is not None and prev_pos <= last_end:
                prev_pos = last_end + 1
            new_grange = GRange(seqid=prev_seqid, start=prev_pos, end=end_pos, name="", strand="")
            ori_names = [g.name for g in active_granges]
            names = "|".join(set(ori_names))
            new_grange.name = names
            ori_strands = [g.strand for g in active_granges]
            strand = set(ori_strands)
            if len(strand) == 1:
                new_grange.strand = strand.pop()
            else:
                new_grange.strand = "."
            split.append(new_grange)
            last_end = end_pos
        if event_type == 1:
            active_granges.add(grange)
        else:
            active_granges.remove(grange)
        count += event_type
        prev_pos = pos
        prev_seqid = seqid
    return split


def filter_intervals(gr_filter, gr):
    """
    Filters intervals in listB based on listA.

    Args:
    - gr_filter (list of GRange): List of GRange intervals to use as a filter.
    - gr (list of GRange): List of GRange intervals to be filtered.

    Returns:
    - List of Grange intervals from listB that are within intervals in listA.
    """


    # Split and sort lists by seqid and start position
    dict_filter = {}
    dict_gr = {}
    for grange in gr_filter:
        dict_filter.setdefault(grange.seqid, []).append(grange)
    for grange in gr:
        dict_gr.setdefault(grange.seqid, []).append(grange)

    for seqid in dict_filter:
        dict_filter[seqid].sort(key=lambda x: x.start)
    for seqid in dict_gr:
        dict_gr[seqid].sort(key=lambda x: x.start)

    filtered = []

    for seqid in dict_filter:
        if seqid not in dict_gr:
            continue
        for grange_out in dict_gr[seqid]:
            for grange_filter in dict_filter[seqid]:
                if grange_out.name == grange_filter.name:
                    if grange_out.within(grange_filter):
                        filtered.append(grange_out)
                        break

    return filtered


def get_kmers(dna, k):
    """ Return all kmers of length k from dna.
    param dna: dict of dna sequences
    param k: length of kmers
    return: dictionary kmer:count
    """
    kmers = {}
    for seq in dna.values():
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if kmer in kmers:
                kmers[kmer] += 1
            else:
                kmers[kmer] = 1
    return kmers


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
                try:
                    key, value = item.split('=')
                    self._attributes_dict[key] = value
                except ValueError:
                    print(item)
                    print(self.attributes)
                    raise

        self._attributes_str = ';'.join(
                ['{}={}'.format(key, value) for key, value in
                 self._attributes_dict.items()]
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
                ['{}={}'.format(key, value) for key, value in
                 self._attributes_dict.items()]
                )
        return self._attributes_str

    @attributes_dict.setter
    def attributes_dict(self, value):
        self._attributes_dict = value
        self._attributes_str = ';'.join(
                ['{}={}'.format(key, value) for key, value in
                 self._attributes_dict.items()]
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


def split_gff3_by_cluster_name(gff3_file, dirname):
    """
    Split gff3 file by cluster name
    :param gff3_file:
    :param dirname: output dir where gff3 files will exported

    :return:
    """
    # create directory if not exists
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    else:
        # if not empty, remove all files
        # this is necessary as we use append mode
        for file in os.listdir(dirname):
            os.remove(os.path.join(dirname, file))

    with open(gff3_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                feature = Gff3Feature(line)
                cluster_name = feature.attributes_dict['Name']
                with open(
                        os.path.join(dirname, F'{cluster_name}.gff3'),
                        'a'
                        ) as out:
                    out.write(feature.print_line())


def annotate_gff(gff, gff_out, library, cpu=1):
    """
    Run annotation on sequences included in gff3, sequences are part of gff attribute
    "consensus_sequence". consensus sequences are extracted from gff3 file and saved
    as dimers and analyzed by RepeatMasker against custom library
    :param gff: gff3 file with tidehunter results containing consensus sequences
    :param gff_out: gff3 file with updated annotation  information
    :param library: custom library for RepeatMasker
    :param cpu: number of cpu cores to use
    :return:
    """
    # read gff3 file and make dictionary with dimer of consensus sequences
    consensus_sequences = {}
    with open(gff, "r") as f:
        for i in f:
            if i.startswith("#"):
                continue
            gff_record = Gff3Feature(i)
            try:
                id_name = gff_record.attributes_dict['ID']
            except KeyError:
                print(F"ID attribute not found in gff, \ngff line:\n{i}")
                exit(1)
            try:
                consensus_sequences[id_name] = gff_record.attributes_dict[
                                                   "consensus_sequence"] * 2
            except KeyError:
                print(F"consensus_sequence attribute not found in gff, \ngff line:\n{i}")
                exit(1)
    # save consensus sequences to fasta file
    tmp_consensus_file = tempfile.NamedTemporaryFile(delete=False).name
    save_fasta_dict_to_file(consensus_sequences, tmp_consensus_file)
    # run repeatmasker with automatic sequence name renaming
    rm_file = run_repeatmasker_with_renaming(tmp_consensus_file, library, cpu)
    # all rm generated files:
    all_rm_files = glob.glob(tmp_consensus_file + ".*")

    seq_lengths = {k: len(v) for k, v in consensus_sequences.items()}
    prefix = gff_out.replace("_annotation.gff3", "")
    rm_annotation = get_repeatmasker_annotation(
            rm_file, seq_lengths, prefix,
            parse_id=False
            )
    add_attribute_to_gff(gff, gff_out, "ID", "annotation", rm_annotation)
    # remove tmp files
    os.remove(tmp_consensus_file)
    for file in all_rm_files:
        os.remove(file)


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


def get_repeatmasker_annotation(rm_file, seq_lengths, prefix, parse_id=True):
    """
    :parse repeatmasker output and calculate proportion of each annotation
    :param rm_file:
    :param seq_lengths:
    :param prefix:
    :param parse_id: if True, parse ID from repeatmasker output, if False, use as it is
    :return:
    """
    seq_rm_info = {}
    with open(rm_file, 'r') as f:
        # parse repeatmasker output, first three lines are header
        for i in range(3):
            # file could be empty!
            try:
                next(f)
            except StopIteration:
                return {}
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
        if parse_id:
            trc_id = "TRC_" + seqid.split("_")[1]
        else:
            trc_id = seqid
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
                    annot_description[trc_id] += F", {v} ({round(p * 100, 1)}%25)"
                    f.write(F"\t{v}\t{annot_summary[trc_id][v]}")
                else:
                    annot_description[trc_id] = F"{v} ({round(p * 100, 1)}%25)"
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
                header = line.strip().split()[0][1:]  # remove part of name after space
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
    # verify that fh is a file handle
    if not hasattr(fh, 'read'):
        raise TypeError("fh must be a file handle")
    fasta_dict = {}
    header = None
    for line in fh:
        if line[0] == '>':
            header = line.strip().split()[0][1:]  # remove part of name after space
            fasta_dict[header] = []
        else:
            if header is None:
                raise ValueError('Fasta file does not start with header')
            fasta_dict[header] += [line.strip()]
    fasta_dict = {k: ''.join(v) for k, v in fasta_dict.items()}
    return fasta_dict


def read_single_fasta_as_generator(fh):
    """
    Read fasta file as a generator
    :param fh: - file handle
    :yield:
    header, sequence
    """
    # Verify that fh is a file handle
    if not hasattr(fh, 'read'):
        raise TypeError("fh must be a file handle")

    header = None
    sequence = []
    for line in fh:
        if line[0] == '>':
            if header is not None:
                yield header, ''.join(sequence)
                sequence = []
            header = line.strip().split()[0][1:]  # Remove part of name after space
        else:
            if header is None:
                raise ValueError('Fasta file does not start with header')
            sequence.append(line.strip())

    if header is not None:
        yield header, ''.join(sequence)



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
    # fasta_dict = read_single_fasta_to_dictionary(open(fasta_file, 'r'))
    with open(fasta_file_split, 'w') as fh_out:
        with open(fasta_file) as fh:
            for header, sequence in read_single_fasta_as_generator(fh):
                matching_table_part = [x for x in matching_table if x[0] == header]
                for header2, i, start, end, new_header in matching_table_part:
                    fh_out.write('>' + new_header + '\n')
                    fh_out.write(sequence[start:end] + '\n')
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
                header = line.strip().split()[0][1:]
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
    For large files it splits fasta file into chunks and run tidehunter on each chunk
    to limit memory usage, on large files, tidehunter is consuming excessive amount of memory

    """
    # verify tidehunter version
    tidehunter_version = subprocess.check_output(
            "TideHunter -v", shell=True
            ).decode().strip()
    assert tidehunter_version == "1.4.3", "TideHunter version 1.4.3 is required"
    # check size of fasta file, if > 1GB, split into chunks to limit memory usage
    # each sequence in fasta ~ 500000nt
    file_size = os.path.getsize(fasta_file)
    file_size_limit = 50000000
    if file_size > file_size_limit:
        # split fasta file into parts
        # run tidehunter on each part
        # merge results
        #
        print("splitting fasta file into parts")
        number_of_parts = int(file_size / file_size_limit) + 1
        print("Number of parts:", number_of_parts)
        fasta_file_parts = split_fasta_to_parts(fasta_file, number_of_parts)
        print("Number of parts:", len(fasta_file_parts))
        tidehunter_parts = []
        for f in fasta_file_parts:
            tmp_file_out = f + ".out"
            tidehunter_cmd = (F"TideHunter -f 2 -o {tmp_file_out} {tidehunter_arguments}"
                              F" {f}")
            print("running TideHunter")
            print(tidehunter_cmd)
            subprocess.check_call(tidehunter_cmd, shell=True)
            tidehunter_parts.append(tmp_file_out)
        # merge results
        print("merging TideHunter results")
        tidehunter_out = fasta_file + ".out"
        with open(fasta_file + ".out", 'w') as fout:
            for p in tidehunter_parts:
                with open(p, 'r') as fin:
                    for line in fin:
                        fout.write(line)
        # remove directory with temporary files (all files are in the same directory)
        tmp_dir_name = os.path.dirname(fasta_file_parts[0])
        shutil.rmtree(tmp_dir_name)
    else:
        # run tidehunter on all
        tidehunter_out = fasta_file + ".out"
        tidehunter_cmd = (F"TideHunter -f 2 -o {tidehunter_out} {tidehunter_arguments}"
                          F" {fasta_file}")

        print("running TideHunter")
        print(tidehunter_cmd)
        subprocess.check_call(tidehunter_cmd, shell=True)

    return tidehunter_out

def split_fasta_to_parts(fasta_file, number_of_parts):
    """split fasta file to parts
    :param fasta_file:
    :param number_of_parts:
    :return: List of paths to temporary fasta files

    """
    total_sequences = len(read_fasta_sequence_size(fasta_file))
    tmp_dir = tempfile.mkdtemp()
    file_list = []
    chunk_size = int(total_sequences / number_of_parts)
    if total_sequences % number_of_parts != 0:
        chunk_size += 1

    print("Number of sequences:", total_sequences)
    print("Number of parts:", number_of_parts)
    print("Chunk size:", chunk_size)
    with open(fasta_file, 'r') as f:
        fasta_gen = read_single_fasta_as_generator(f)
        for n in range(number_of_parts):
            f_out = F"{tmp_dir}/part_{n}.fasta"
            file_list.append(f_out)
            with open(f_out, 'w') as out_file:
                for _ in range(chunk_size):
                    try:
                        header, sequence = next(fasta_gen)
                        out_file.write(">" + header + "\n")
                        out_file.write(sequence + "\n")
                    except StopIteration:
                        break
    print(file_list)
    return file_list


def split_fasta_to_parts2(fasta_file, number_of_parts):
    """split fasta file to parts
    :param fasta_file:
    :param number_of_parts:
    :return: List of paths to temporary fasta files

    """
    with open(fasta_file, 'r') as f:
        s = read_single_fasta_to_dictionary(f)
    l = len(s)
    chunk_size = int(l / number_of_parts) + 1
    # check sanity of chunk size
    if chunk_size < 1:
        chunk_size = 1
    tmp_dir = tempfile.mkdtemp()
    file_list = []
    n = 0
    while True:
        k = list(itertools.islice(s.keys(), chunk_size))
        s_part = {i: s.pop(i) for i in k}
        f_out = F"{tmp_dir}/part_{n}.fasta"
        save_fasta_dict_to_file(s_part, f_out)
        file_list.append(f_out)
        if len(s) == 0:
            break
        n += 1
    return file_list



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
        columns = [self.seq_name, "TideHunter", "tandem_repeat", self.start, self.end,
                   self.aver_match, ".", "."]
        attributes = {
            "ID": self.repeat_ID, "consensus_sequence": self.consensus,
            "consensus_length": self.cons_length,
            "copy_number": self.copy_numer, "cverage_match": self.aver_match
            }

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

def fasta_to_list(fasta_file):
    """
    convert fasta file to list
    :param fasta_file:
    :return: list of tuples with name and sequences
    use this instead of fasta_todit if there are duplicated names in fasta file

    """
    fasta_list = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                seq_name = line.split()[0].replace(">", "")
                fasta_list.append([seq_name, ""])
            else:
                fasta_list[-1][1] += line.strip()
    return fasta_list

def save_fasta_list_to_file(fasta_list, fasta_file):
    """
    save fasta list to file
    :param fasta_list: list
    :param fasta_file: str
    :return:
    """
    with open(fasta_file, 'w') as f:
        for name, seq in fasta_list:
            f.write(F">{name}\n{seq}\n")


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


def gff3_to_fasta(gff3_file, fasta_file, additonal_attribute=None):
    """
    extract fasta sequences from gff3 file
    it is generator, returns one sequence at time and seq ID plus additional attribute
    if provided
    :param additonal_attribute: yield additional attribute from gff3 file
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
            if "ID" not in gff3_feature.attributes_dict:
                gff3_feature.attributes_dict["ID"] = (gff3_feature.seqid + "_" +
                                                      str(gff3_feature.start) + "_" +
                                                      str(gff3_feature.end))
            if additonal_attribute:
                yield [gff3_feature.attributes_dict['ID'], s,
                       gff3_feature.attributes_dict[additonal_attribute]]
            else:
                yield [gff3_feature.attributes_dict['ID'], s]


def save_fasta_dict_to_file(fasta_dict, fasta_file):
    """
    save fasta dictionary to file
    :param fasta_dict: dictionary with fasta sequences
    :param fasta_file: path to fasta file
    :return:
    """
    with open(fasta_file, 'w') as f:
        for k, v in fasta_dict.items():
            f.write(">{}\n{}\n".format(k, v.upper()))


# run mmseqs2 on consensus sequences
def find_cluster_by_mmseqs2(fasta_file, cpu=4, memory_limit=64000):
    """
    run mmseqs2 on consensus sequences
    :param cpu:
    :param sequences:
    :param memory_limit: in megabytes
    :return: clusters
    """
    # check fasta file size
    if os.path.getsize(fasta_file) == 0:
        return {}

    print("Clustering by mmseqs2")
    tmp_dir = tempfile.mkdtemp()

    cmd = (F'mmseqs easy-cluster {fasta_file} {fasta_file}.clu'
           F' {tmp_dir} --cluster-mode 0 -v 1 '
           F'--mask 0  -s 1 --threads {cpu} '
           F'--split-memory-limit {memory_limit}'
           )
    subprocess.check_call(cmd, shell=True)

    # read clusters to dictionary
    with open(F"{fasta_file}.clu_cluster.tsv", 'r') as f:
        clusters = {}
        for line in f:
            cluster, seq_id = line.split()
            clusters[seq_id] = cluster
    # remove temporary files
    os.remove(F"{fasta_file}.clu_all_seqs.fasta")
    os.remove(F"{fasta_file}.clu_cluster.tsv")
    os.remove(F"{fasta_file}.clu_rep_seq.fasta")
    rmtree(tmp_dir)
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


def get_connected_component_clusters(pairs, representative=False):
    """
    find connected components in graph
    :param pairs: list of pairs
    :param representative: return clusters as dictionary with representative
    :return: clusters as list of lists or dictionary with representative

    """
    # Create a NetworkX graph from the list of edges
    g = nx.Graph()
    g.add_edges_from(pairs)

    # Find the connected components in the graph
    components = list(nx.connected_components(g))

    # Convert the components to a list of clusters
    clusters = [list(component) for component in components]
    if representative:
        clusters_representative = {}
        for vertices in components:
            v_representative = sorted(vertices)[0]
            for v in vertices:
                clusters_representative[v] = v_representative
        return clusters_representative
    else:
        return clusters


def extract_sequences_from_gff3(gff3_file, fasta_file, output_dir):
    """
    extract sequences from gff3 file
    :param gff3_file: path to gff3 file
    :param fasta_file: path to fasta file
    :param output_dir: path to output directory
    :return: dictionary with path to fasta files
    """
    # make dir:
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    else:
        # remove old files
        for f in os.listdir(output_dir):
            try:
                os.remove(F"{output_dir}/{f}")
            except OSError:
                pass
    fasta_files = {}
    for seq_id, seq, name in gff3_to_fasta(gff3_file, fasta_file, "Name"):
        f = F"{output_dir}/{name}.fasta"
        fasta_files[name] = f
        out_file = F"{output_dir}/{name}.fasta"
        with open(out_file, 'a') as f:
            f.write(F">{seq_id}\n{seq}")
            f.write(seq)
            f.write("\n")
    return fasta_files


def group_sequences_by_orientation(sequences, k):
    """
    group sequences by orientation
    :param sequences: dictionary with sequences
    :param k: kmers size
    :return:
    """
    # some sequences may have lower case letters, change them to upper case
    for i in sequences:
        sequences[i] = sequences[i].upper()
    groups = {'forward': [], 'reverse': []}
    reference = max(list(sequences.values()), key=len)
    ref_kmers = generate_kmers(reference, k)

    for i in sequences:
        seq = sequences[i]
        rc_seq = reverse_complement(seq)

        kmers = generate_kmers(seq, k)
        kmers_sorted = sorted(kmers, key=lambda key: kmers[key], reverse=True)
        kmers = set(list(kmers_sorted[0:int(len(kmers_sorted) * 0.5)]))

        rc_kmers = generate_kmers(rc_seq, k)
        rc_kmers_sorted = sorted(rc_kmers, key=lambda key: rc_kmers[key], reverse=True)
        rc_kmers = set(rc_kmers_sorted[0:int(len(rc_kmers_sorted) * 0.5)])

        overlap = len(kmers.intersection(ref_kmers))
        rc_overlap = len(rc_kmers.intersection(ref_kmers))
        if overlap > rc_overlap:
            groups['forward'].append(i)
        else:
            groups['reverse'].append(i)
    return groups


def reverse_complement(dna):
    """
    reverse complement of dna sequence, including ambiguous bases
    :param dna:
    :return: reverse complement of dna sequence
    """
    # complement including all IUPAC codes
    complement = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R': 'Y', 'Y': 'R', 'S': 'S',
        'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B',
        'N': 'N'
        }
    return ''.join(complement[base] for base in reversed(dna.upper()))


def generate_kmers(dna, k):
    """
    generate kmers from dna sequence
    :param dna: string
    :param k: length of kmers
    :return: dictionary with kmers, value is frequency of each kmer
    """
    kmers = dict()
    for i in range(len(dna) - k + 1):
        kmer = dna[i:i + k]
        if kmer not in kmers:
            kmers[kmer] = 0
        kmers[kmer] += 1
    return kmers


def run_blastn(fasta_file, dust=False, cpu=4):
    """
    run blastn on fasta file
    :param cpu:
    :param dust: dust filter for blastn
    :param sequences : dictionary with sequences
    :return: dictionary with clusters
    """
    tmp_dir = tempfile.mkdtemp()
    # make symbolic link to fasta file

    fasta_file_symlink = F'{tmp_dir}/seqs.fasta'
    os.symlink(fasta_file, fasta_file_symlink)
    blast_out = F'{tmp_dir}/blastn.out'
    # make blast database
    cmd = F"makeblastdb -in {tmp_dir}/seqs.fasta -dbtype nucl"
    subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL)
    # blast parameters:
    dust = "yes" if dust else "no"
    outfmt = "'6 qseqid sseqid pident length evalue bitscore qlen slen'"
    # run blastn
    cmd = (F"blastn -task blastn -query {fasta_file_symlink} -db {fasta_file_symlink} -outfmt"
           F" {outfmt}"
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
    with open(fasta_file, 'r') as f:
        for id, s in read_single_fasta_as_generator(f):
            pairs.add((id, id))
    rmtree(tmp_dir)
    return pairs


def get_ssrs_proportions(fasta_file):
    """
    Run dustmasker on fasta file to identifie simple repeats from masking. if
    repeats is masked by more than 80% of the sequence, then ths sequnece is
    scanned for SSRs and the proportion of SSRs is calculated.
    it is indended to be used for masking tandem repeats with monomer up to 6 nt
    :param sequences: dictionary with sequences
    :return: dictionary with dustmasker results
    """
    tmp_dir = tempfile.mkdtemp()
    dust_out = F'{tmp_dir}/dust.out'
    # run dustmasker
    cmd = (F"dustmasker -in {fasta_file} -out {dust_out} -outfmt fasta -window 64 "
           F"-level 20")
    subprocess.check_call(cmd, shell=True)
    # read dustmasker results to dictionary
    dust = {}
    lengths = {}
    # read fasta and calucutate proportion of lower case letters
    with open(dust_out, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line[1:].strip()
                dust[seq_id] = 0
                lengths[seq_id] = 0
            else:
                lengths[seq_id] += len(line.strip())
                dust[seq_id] += sum(1 for char in line.strip() if char.islower())
    ssrs_prop = {}
    dust08_seqs = F'{tmp_dir}/dust08_seqs.fasta'
    dust08_id = []
    for id_name in dust:
        dust[id_name] /= lengths[id_name]
        ssrs_prop[id_name] = 0
        if dust[id_name] > 0.8:
            dust08_id.append(id_name)
    # filter fasta file in dust08_id
    with open(fasta_file, 'r') as f, open(dust08_seqs, 'w') as fout:
        for id, s in read_single_fasta_as_generator(f):
            if id in dust08_id:
                ssrs = find_ssrs(s)
                ssrs_prop[id] = single_sequence_ssrs_proportion(ssrs)
    rmtree(tmp_dir)
    return ssrs_prop


def homopolymer(motif, motiflength):
    """
    return true if motif is repeat of single nucleotide
    :param motif:
    :param motiflength:
    :return:
    """
    reps = motiflength - 1
    return True if re.search(r'([gatc])\1{%d}' % reps, motif) else False


def novel(position, locations):
    """
    check if position is novel
    :param position:
    :param locations:
    :return:
    """
    if position in locations:
        return False
    else:
        locations[position] = True
        return True


def find_ssrs(sequence, specs=None):
    """
    :param sequence: string
    :param specs: list of tuples with length of motif and minimum number of repeats
    :return: list of tuples with motif and number of repeats
    """
    if specs is None:
        specs = [(1,20),
                 (2, 9),  # dinucl. with >= 9 repeats
                 (3, 6),  # trinucl. with >= 6 repeats
                 (4, 5),  # tetranucl. with >= 5 repeats
                 (5, 5),  # pentanucl. with >= 5 repeats
                 (6, 5),  # hexanucl. with >= 5 repeats
                 (7, 5),  # heptanucl. with >= 5 repeats
                 (8, 5),  # octanucl. with >= 5 repeats
                 (9, 5),  # nonanucl. with >= 5 repeats
                 (10, 5),  # decanucl. with >= 5 repeats
                 (11, 5),  # undecanucl. with >= 5 repeats
                 (12, 5),  # dodecanucl. with >= 5 repeats
                 (13, 5),  # tridecanucl. with >= 5 repeats
                 (14, 5)]  # tetradecanucl. with >= 5 repeats
    results = []
    homopolymers = []
    locations = {}
    for i in range(len(specs)):
        motiflength, minreps = specs[i]
        regexp = r'(([gatc]{%d})\2{%d,})' % (motiflength, minreps - 1)
        for match in re.finditer(regexp, sequence.lower()):
            ssr = match.group(1)
            motif = match.group(2).lower()
            if homopolymer(motif, motiflength):
                if motiflength == 1:
                    ssrlength = len(ssr)
                    repeats = ssrlength // motiflength
                    end = match.end()
                    start = end - ssrlength + 1
                    homopolymers.append(
                            {
                                'ssr_number': len(results) + 1,
                                'motif_length': motiflength,
                                'motif_sequence': motif,
                                'repeats': repeats,
                                'start': start,
                                'end': end,
                                'seq_length': len(sequence) - sequence.count('N'),
                                }
                            )
                continue
            ssrlength = len(ssr)
            repeats = ssrlength // motiflength
            end = match.end()
            start = end - ssrlength + 1
            if novel(start, locations):
                results.append(
                        {
                            'ssr_number': len(results) + 1,
                            'motif_length': motiflength,
                            'motif_sequence': motif,
                            'repeats': repeats,
                            'start': start,
                            'end': end,
                            'seq_length': len(sequence) - sequence.count('N'),
                            }
                        )
    if len(homopolymers) > 0 and len(results) == 0:
        results += homopolymers
    return results


def normalize_motif(motif):
    """
    find the lexicographically smallest string in the set of all possible cyclic
    rotations. This will be the normalized version of the motif.
    :param motif:
    :return: normalized motif
    """
    doubled_motif = (motif + motif).upper()
    rc_motif: str = reverse_complement(doubled_motif)
    fm = min(doubled_motif[i:i + len(motif)] for i in range(len(motif)))
    rcm = min(rc_motif[i:i + len(motif)] for i in range(len(motif)))
    return min(fm, rcm)



def get_unique_motifs(motifs):
    """
    return unique motifs
    :param motifs: list
    :return: motifs: list
    exclude motifs which are multiple of shorter motif
    """
    motifs.sort(key=len)  # Sort motifs by length to ensure shorter motifs are checked first
    unique_motifs = []

    for i, motif in enumerate(motifs):
        is_unique = True
        for shorter_motif in unique_motifs:
            if len(motif) % len(shorter_motif) == 0 and shorter_motif * (
                    len(motif) // len(shorter_motif)) == motif:
                is_unique = False
                break
        if is_unique:
            unique_motifs.append(motif)
    return unique_motifs


def sum_up_ssrs_motifs(results):
    """
    sum up motifs
    :param results: list of dictionaries with ssr results
    """
    motifs = {}
    for i in results:
        motif = i['motif_sequence']
        motif = normalize_motif(motif)
        if motif in motifs:
            motifs[motif] += i['repeats']
        else:
            motifs[motif] = i['repeats']
    unique_motifs = get_unique_motifs(list(motifs.keys()))
    # keep only unique motifs
    motifs = {k: motifs[k] for k in unique_motifs}

    return motifs


def get_ssrs_description_multiple(seq_strs):
    """
    return description of ssrs in multiple sequences
    :param seq_strs:
    :return:
    """
    ssrs = []
    for seq_str in seq_strs:
        ssrs += find_ssrs(seq_str)
    motif_count = sum_up_ssrs_motifs(ssrs)
    length = sum([len(seq_str) for seq_str in seq_strs])
    motif_percent = {k: 100 * v * len(k) / length for k, v in motif_count.items()}
    desc = ""
    # iterate over sorted motifs by percent
    for motif, percent in sorted(
            motif_percent.items(),
            key=lambda x: x[1],
            reverse=True
            ):
        desc += F"{motif} ({percent:.1f}%25), "
    return desc[:-2]


def get_ssrs_description(seq_str):
    """
    return description of ssrs in single sequence
    :param seq_str:
    :return:
    """
    ssrs = find_ssrs(seq_str)

    motif_count = sum_up_ssrs_motifs(ssrs)

    motif_percent = {k: 100 * v * len(k) / len(seq_str) for k, v in motif_count.items()}
    desc = ""
    # iterate over sorted motifs by percent
    for motif, percent in sorted(
            motif_percent.items(),
            key=lambda x: x[1], reverse=True
            ):
        desc += F"{motif} ({percent:.2f}%25), "
    # remove_last_comma
    return desc[:-2]


def single_sequence_ssrs_proportion(ssrs):
    """
    calculate proportion of sequence covered by ssrs
    :param ssrs: ssrs is list of dictionaries with ssr information
    :return: ssrs proportion
    """
    if len(ssrs) == 0:
        return 0
    # boolean list with length of sequence
    is_ssr = [False] * len(ssrs[0])
    # set True for each position covered by ssr
    for ssr in ssrs:
        is_ssr[ssr['start'] - 1:ssr['end']] = [True] * ssr['seq_length']
    # calculate proportion of True values
    return sum(is_ssr) / len(is_ssr)


def find_clusters_by_blast_connected_component(
        fasta_file, dust=False,
        cpu=4
        ):
    """
    find clusters by blastn, return dictionary with clusters
    cluaste are connected components in graph
    :param cpu:
    :param dust: use dust filter
    :param consensus_representative:
    :return: clusters
    """
    # check fasta file size
    if os.path.getsize(fasta_file) == 0:
        return {}

    print("Clustering by BLASTN")
    pairs = run_blastn(fasta_file, dust=dust, cpu=cpu)
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


def get_cluster_size2(fin, clusters):
    """
    get cluster size
    :param fin: input gff
    :param clusters: clusters dictionary id:cluster_id
    :return: cluster size as total length of intervals
    Calculate total size of region for particular cluster, note that regions can overlap
    but each region will be counted only once.
    """
    # first read all regions and sort them by chromosome and start position for each
    # cluster
    cluster_regions = {}
    with open(fin, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            gff3_feature = Gff3Feature(line)
            cluster_id = clusters[gff3_feature.attributes_dict['ID']]
            if cluster_id not in cluster_regions:
                cluster_regions[cluster_id] = []
            cluster_regions[cluster_id].append(
                    (gff3_feature.seqid,
                     gff3_feature.start,
                     gff3_feature.end)
                    )
    # merge overlapping regions before calculating size
    cluster_size = {}
    for cluster_id in cluster_regions:
        cluster_regions[cluster_id] = merge_genomic_intervals(cluster_regions[cluster_id])
        # calculate size
        cluster_size[cluster_id] = sum(
                [end - start for seqid, start, end in cluster_regions[cluster_id]]
                )
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
            consensus = gff3_feature.attributes_dict['consensus_sequence']
            if cluster_id not in consensus_clusters:
                consensus_clusters[cluster_id] = {}
                consensus_clusters_dimers[cluster_id] = {}
            consensus_clusters[cluster_id][unique_id] = consensus
            consensus_clusters_dimers[cluster_id][unique_id] = consensus * 2

            f2.write(gff3_feature.print_line())
    return consensus_clusters, consensus_clusters_dimers

def extend_gff3_intervals(gff3_file, gff3_out_file, rep_size, extend_prop=0.05):
    """
    extend gff3 intervals by rep_size * extend_prop
    :param gff3_file:
    :param gff3_out_file:
    :param rep_size: dict with repeat sizes
    :param extend_prop:
    :return: None

    This is helper function to for evaluation which TRC region to keep
    """
    with open(gff3_file, 'r') as f1, open(gff3_out_file, 'w') as f2:
        for line in f1:
            if line.startswith("#"):
                f2.write(line)
                continue
            gff3_feature = Gff3Feature(line)
            size1 = rep_size[gff3_feature.attributes_dict['Name']]
            size2 = (gff3_feature.end - gff3_feature.start) * 2 * extend_prop
            size_ext = int(min(size1, size2))
            gff3_feature.start = max(1, gff3_feature.start - size_ext)
            gff3_feature.end += size_ext
            f2.write(gff3_feature.print_line())

def merge_overlapping_gff3_intervals(gff3_file, gff3_out_file, use_strand=False):
    """
    merge overlapping intervals in gff3 file
    merge only if they have same Cluster_ID
    :param gff3_file: path to gff3 file
    :param gff3_out_file: path to output gff3 file
    :param use_strand: if True, merge only intervals on the same strand
    :return: nothing
    """
    # read gff3 file, split to lists by cluster ID and seqname
    gff3_dict = {}
    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            gff3_feature = Gff3Feature(line)
            # Cluster_ID could be missing in some gff3 files, use Name instead
            if 'Cluster_ID' not in gff3_feature.attributes_dict:
                cluster_id = gff3_feature.attributes_dict['Name']
            else:
                cluster_id = gff3_feature.attributes_dict['Cluster_ID']
            seqname = gff3_feature.seqid
            strand = gff3_feature.strand
            if use_strand:
                seq_cluster_id = (seqname, cluster_id,strand)
            else:
                seq_cluster_id = (seqname, cluster_id)
            if seq_cluster_id not in gff3_dict:
                gff3_dict[seq_cluster_id] = []
            gff3_dict[seq_cluster_id].append(gff3_feature)
    with open(gff3_out_file, 'w') as f:
        f.write("##gff-version 3\n")
        for i in gff3_dict:
            if use_strand:
                strand = i[2]
            else:
                strand = "."
            intervals = []
            for j in gff3_dict[i]:
                intervals.append((j.start, j.end))
            intervals = merge_intervals(intervals)
            for start, end in intervals:
                gff_line = (F'{i[0]}\tTideCluster\ttandem_repeat'
                            F'\t{start}\t{end}\t{1}\t{strand}\t.\tName={i[1]}\n')

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


def merge_genomic_intervals(intervals):
    """
    :param intervals:list of (seqid, start, end) tuples
    :return:
    list of merged (seqid, start, end) tuples
    """
    # same as merge_intervals but uses also seqid
    # split it by seqid and merge each seqid separately using merge_intervals
    intervals_dict = {}
    for seqid, start, end in intervals:
        if seqid not in intervals_dict:
            intervals_dict[seqid] = []
        intervals_dict[seqid].append((start, end))
    merged = []
    for seqid in intervals_dict:
        merged += [(seqid, start, end) for start, end in
                   merge_intervals(intervals_dict[seqid])]
    return merged


def filter_gff_remove_duplicates(gff3_file):
    """
    Filter gff3 file by removing duplicates, diplicates can have different ID but same
    exactly genomic region
    :param gff3_file:
    :return: filtered gff3 file path
    """
    gff_out = tempfile.NamedTemporaryFile(delete=False).name
    gff_data = {}
    with open(gff3_file, 'r') as f1:
        for line in f1:
            if line.startswith("#"):
                continue
            gff3_feature = Gff3Feature(line)
            gff_data[gff3_feature.attributes_dict['ID']] = gff3_feature
    gff_data = OrderedDict(
            sorted(gff_data.items(), key=lambda t: (t[1].seqid, t[1].start))
            )
    duplicated_ids = set()
    for i, (k1, v1) in enumerate(gff_data.items()):
        if i == len(gff_data) - 1:
            break
        k2, v2 = list(gff_data.items())[i + 1]
        if v1.seqid == v2.seqid and v1.start == v2.start and v1.end == v2.end:
            duplicated_ids.add(k1)  # add just one, second will be kept
    if len(duplicated_ids) > 0:
        with open(gff3_file, 'r') as f1, open(gff_out, 'w') as f2:
            for line in f1:
                if line.startswith("#"):
                    f2.write(line)
                    continue
                gff3_feature = Gff3Feature(line)
                if gff3_feature.attributes_dict['ID'] not in duplicated_ids:
                    f2.write(line)
        return gff_out
    else:
        return gff3_file


def filter_gff_by_length(gff3_file, gff_short, min_length=1000):
    """
    Filter gff3 file by size of intervals
    :param gff3_file:
    :param min_length:
    :param gff_short: path to output gff3 file with short intervals
    :return: filtered gff3 file path
    """
    gff_out = tempfile.NamedTemporaryFile(delete=False).name
    with open(gff3_file, 'r') as f1, open(gff_out, 'w') as f2, open(gff_short, 'w') as f3:
        for line in f1:
            if line.startswith("#"):
                f2.write(line)
                f3.write(line)
                continue
            gff3_feature = Gff3Feature(line)
            if gff3_feature.end - gff3_feature.start > min_length:
                f2.write(line)
            else:
                f3.write(line)
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



def run_cmd(cmd):
    """
    runs shell command and returns status
    :param cmd:
    :return: list
    """
    try:
        # run command and capture warning and error messages
        result = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        # if command fails, print error message and return error
        if result.returncode != 0:
            print(result.stderr.decode('utf-8'))
            return [cmd, 'error']
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode('utf-8'))
        return [cmd, 'error']
    return [cmd, 'ok']

def filter_fasta_file(input_fasta, output_fasta, ids):
    """
    filter fasta file by ids
    :param input_fasta:
    :param output_fasta:
    :param ids:
    :return: None
    """
    with open(input_fasta, 'r') as f1, open(output_fasta, 'w') as f2:
        for id, seq in read_single_fasta_as_generator(f1):
            if id in ids:
                f2.write(F">{id}\n{seq}\n")


def rename_fasta_sequences_to_indices(input_fasta, output_fasta):
    """
    Rename FASTA sequences to numerical indices to handle RepeatMasker 50-character limit.
    :param input_fasta: path to input FASTA file
    :param output_fasta: path to output FASTA file with renamed sequences
    :return: dictionary mapping new names (indices) to original names
    """
    name_mapping = {}
    with open(input_fasta, 'r') as f1, open(output_fasta, 'w') as f2:
        for i, (original_id, seq) in enumerate(read_single_fasta_as_generator(f1)):
            new_id = str(i + 1)  # Start from 1
            name_mapping[new_id] = original_id
            f2.write(F">{new_id}\n{seq}\n")
    return name_mapping


def restore_sequence_names_in_repeatmasker_output(rm_file, output_file, name_mapping):
    """
    Restore original sequence names in RepeatMasker output file.
    :param rm_file: path to RepeatMasker .out file with renamed sequences
    :param output_file: path to output file with restored names
    :param name_mapping: dictionary mapping indices to original names
    :return: None
    """
    with open(rm_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line_num, line in enumerate(f_in):
            if line_num < 3:  # Skip RepeatMasker header (first 3 lines)
                f_out.write(line)
                continue
            
            # Skip empty lines
            if not line.strip():
                f_out.write(line)
                continue
            
            # Split the line and replace the query sequence name (5th column, index 4)
            parts = line.split()
            if len(parts) >= 5:
                query_name = parts[4]
                if query_name in name_mapping:
                    parts[4] = name_mapping[query_name]
                    f_out.write(' '.join(parts) + '\n')
                else:
                    f_out.write(line)  # Keep original line if mapping not found
            else:
                f_out.write(line)  # Keep malformed lines as-is


def run_repeatmasker_with_renaming(input_fasta, library, cpu=1, additional_params=""):
    """
    Run RepeatMasker with automatic sequence name renaming to handle 50-character limit.

    Args:
        input_fasta: Path to input FASTA file
        library: Path to RepeatMasker library
        cpu: Number of CPU cores to use
        additional_params: Additional RepeatMasker parameters

    Returns:
        Path to RepeatMasker output file (.out) with restored sequence names
    """
    # Create renamed version to handle RepeatMasker 50-character sequence name limit
    input_dir = os.path.dirname(input_fasta)
    input_basename = os.path.basename(input_fasta)  # do not assume input_fasta has .fasta extension!
    input_renamed = F"{input_dir}/{input_basename}_renamed.fasta"


    print("input fasta file:", input_fasta)
    print("renamed fasta file:", input_renamed)
    # Rename sequences to numerical indices
    name_mapping = rename_fasta_sequences_to_indices(input_fasta, input_renamed)

    # Build RepeatMasker command
    cmd = (F"RepeatMasker -pa {cpu} -lib {library} -e ncbi -s -no_is -norna "
           F"-nolow {additional_params} -dir {input_dir} {input_renamed}")

    print(cmd)
    print("---------")
    subprocess.run(cmd, shell=True)

    # Get RepeatMasker output files
    rm_file_renamed = F"{input_renamed}.out"
    rm_file_original = F"{input_fasta}.out"

    # Restore original sequence names in RepeatMasker output
    restore_sequence_names_in_repeatmasker_output(rm_file_renamed, rm_file_original, name_mapping)

    # Clean up renamed files
    os.remove(input_renamed)
    if os.path.exists(rm_file_renamed):
        os.remove(rm_file_renamed)

    return rm_file_original


def mask_fasta_with_gff3(fasta_file, gff3_file):
    """
    Mask FASTA sequences using bedtools maskfasta based on GFF3 coordinates.
    Regions in GFF3 file will be masked with 'N'.

    :param fasta_file: Path to input FASTA file
    :param gff3_file: Path to GFF3 file with regions to mask
    :return: Path to masked FASTA file
    """
    # Create temporary BED file from GFF3
    bed_file = tempfile.NamedTemporaryFile(delete=False, suffix=".bed").name

    with open(gff3_file, 'r') as f_in, open(bed_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                continue
            gff_feature = Gff3Feature(line)
            # BED format: seqid, start (0-based), end, name
            # GFF3 is 1-based, BED is 0-based
            bed_line = F"{gff_feature.seqid}\t{gff_feature.start - 1}\t{gff_feature.end}\n"
            f_out.write(bed_line)

    # Create output masked FASTA file
    masked_fasta = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta").name

    # Run bedtools maskfasta
    cmd = F"bedtools maskfasta -fi {fasta_file} -bed {bed_file} -fo {masked_fasta}"
    print(f"Masking FASTA: {cmd}")
    subprocess.check_call(cmd, shell=True)

    # Clean up BED file
    os.remove(bed_file)

    return masked_fasta


def prepare_fasta_input(fasta_path):
    """
    Check if fasta is gzipped and decompress to temp directory if needed.
    Returns path to uncompressed fasta file and cleanup function.

    :param fasta_path: Path to input fasta (can be .gz or .fasta.gz)
    :return: (path_to_uncompressed_fasta, cleanup_function)
    """
    # Check if file is gzipped (by extension or magic bytes)
    is_gzipped = fasta_path.endswith('.gz')

    if not is_gzipped:
        # Check magic bytes as fallback
        try:
            with open(fasta_path, 'rb') as f:
                is_gzipped = f.read(2) == b'\x1f\x8b'
        except:
            pass

    if not is_gzipped:
        # Not gzipped, return original path with no-op cleanup
        return fasta_path, lambda: None

    # Create temp uncompressed file
    # Use basename without .gz for the temp file
    base_name = os.path.basename(fasta_path)
    if base_name.endswith('.gz'):
        base_name = base_name[:-3]

    temp_dir = tempfile.gettempdir()
    temp_fasta = os.path.join(temp_dir, f"tidecluster_{os.getpid()}_{base_name}")

    print(f"Decompressing gzipped FASTA: {fasta_path} -> {temp_fasta}")

    with gzip.open(fasta_path, 'rb') as f_in:
        with open(temp_fasta, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Create cleanup function
    def cleanup():
        if os.path.exists(temp_fasta):
            print(f"Cleaning up temporary FASTA: {temp_fasta}")
            os.remove(temp_fasta)

    return temp_fasta, cleanup
