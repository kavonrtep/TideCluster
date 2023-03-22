#!/usr/bin/env python3
import logging
logger = logging.getLogger(__name__)
import sys
import operator

REQUIRED_VERSION = (3, 4)
MAX_PRINT = 10
MEGABLAST = "-task megablast"
HITSORT = "-task megablast"

if sys.version_info < REQUIRED_VERSION:
    raise Exception("\n\npython 3.4 or higher is required!\n")

# additional functions
class Sequence:

    def __init__(self, seq, name="", paired=False):
        # the mode os seq storage can be changed later to make it more
        # memory efficient
        self._seq = bytes(str(seq), "ascii")
        self.name = str(name)

    @property
    def seq(self):
        return self._seq.decode("utf-8")

    @seq.setter
    def seq(self, value):
        self._seq = bytes(str(value), "ascii")

    def __str__(self):
        return "{0} : {1}".format(self.name, self.seq)

    @staticmethod
    def read_fasta(fasta_file_name):
        '''
        generator - reads sequences from fasta file
        return sequence one by one
        '''
        with open(fasta_file_name, 'r') as f:
            header = None
            seqstr = None
            for rawline in f:
                line = rawline.strip()
                if line == "":
                    continue
                if line[0] == ">":
                    if header and seqstr:
                        yield Sequence(seqstr, header)
                        # reset
                        seqstr = None
                        header = line[1:]
                    elif seqstr:
                        Warning("sequence was not preceeded by header")
                    else:
                        header = line[1:]
                else:
                    seqstr = line if not seqstr else seqstr + line
        # skip empty lines:
        if header and seqstr:
            yield Sequence(seqstr, header)
        return

    def write2fasta(self, file_object):
        file_object.write(">{0}\n{1}\n".format(self.name, self.seq))


def get_kmers(string, width=11):
    L = len(string)
    parts = [string[i:i + width] for i in range(L - width + 0)]
    return parts


def count_kmers_from_file(f, width=11):
    counts = {}
    for i in Sequence.read_fasta(f):
        a = get_kmers(i.seq, width)
        for km in a:
            if "N" in km:
                continue
            if km in counts:
                counts[km] += 1
            else:
                counts[km] = 1
    sorted_counts = sorted(counts.items(),
                           key=operator.itemgetter(1),
                           reverse=True)
    return sorted_counts


if __name__ == "__main__":
    L = len(sys.argv) - 1
    kmer_length = int(sys.argv[-1])
    files = sys.argv[1:-1]
    for fin in files:
        counts = count_kmers_from_file(fin, kmer_length)
        fout = "{}_{}.kmers".format(fin, kmer_length)
        with open(fout, "w") as f:
            for i in counts:
                f.write("{}\t{}\n".format(*i))
        print(fout)
