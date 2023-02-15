""" Utilities for the project. """
import os
import tempfile
from itertools import cycle

class Gff3Feature:
    """
    Class for gff3 feature
    """

    def __init__(self, line):
        self.line = line
        self.items = line.strip().split('\t')
        self.header = self.items[0]
        self.source = self.items[1]
        self.type = self.items[2]
        self.start = int(self.items[3])
        self.end = int(self.items[4])
        self.score = self.items[5]
        self.strand = self.items[6]
        self.frame = self.items[7]
        self.attributes = self.items[8]
        self.attributes_dict = {}
        for item in self.attributes.split(';'):
            if item != '':
                key, value = item.split('=')
                self.attributes_dict[key] = value

        self.attributes_str = ';'.join(
            ['{}={}'.format(key, value) for key, value in self.attributes_dict.items()]
            )

    def __str__(self):
        return '\t'.join(
            [self.header, self.source, self.type, str(self.start), str(self.end),
             self.score, self.strand, self.frame, self.attributes_str]
            ) + '\n'

    def __repr__(self):
        return '\t'.join(
            [self.header, self.source, self.type, str(self.start), str(self.end),
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
            [self.header, self.source, self.type, str(self.start), str(self.end),
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
        if self.start == other.start and self.end == other.end and self.header == \
                other.header:
            return True
        else:
            return False


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
    :param fh:
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
                end = ((
                                   i + 1) * adjusted_chunk_size + overlap) if i + 1 < \
                                                                              number_of_chunks else size
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
                        feature.header, feature.start, feature.end, matching_table
                        )
                    for new_header, new_start, new_end, sequence_length in new_coords:
                        if new_start >= 1 and new_end <= sequence_length:
                            feature.header = new_header
                            feature.start = new_start
                            feature.end = new_end
                            fh_out.write(str(feature))
    return gff3_file_recalculated


# recalculate gff3 back to original coordinates, use gff3_feature class
def recalculate_gff3_back_to_original_coordinates(
        gff3_file, matching_table, chunk_size, overlap
        ):
    """
    Recalculate gff3 back to original coordinates, use gff3_feature class
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
                    ori_header, ori_start, ori_end, real_chunk_size = \
                        get_original_header_and_coordinates(
                        feature.header, feature.start, feature.end, matching_table
                        )
                    # if feature is too close (less than 100 nt) to the end ends of
                    # chunk,skip it
                    if feature.start < 100 or feature.end > real_chunk_size - 100:
                        continue
                    feature.header = ori_header
                    feature.start = ori_start
                    feature.end = ori_end
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
