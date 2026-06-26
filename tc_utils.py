""" Utilities for the project. """
import bisect
import collections
import csv
import glob
import gzip
import itertools
import os
import shutil
import statistics
import subprocess
import tempfile
from collections import OrderedDict
from itertools import cycle
from multiprocessing import Pool
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

    A GRange ``g`` from ``gr`` is kept iff some interval ``f`` in ``gr_filter``
    on the same seqid has ``f.name == g.name`` and ``g`` is contained in ``f``
    (``f.start <= g.start`` and ``g.end <= f.end``), i.e. ``g.within(f)``.

    Implementation: a per-(seqid, name) sweep. For each name the filter
    intervals are sorted by start and a running maximum of their ends is
    precomputed. Every filter with ``start <= g.start`` is a containment
    candidate, so containment reduces to "is the largest end among those
    candidates >= g.end" — a single binary search answers it. This is
    O((N+M) log M) and returns exactly the same objects, in the same order, as
    the previous O(N*M) nested loop (the outer iteration order over ``gr`` per
    seqid, sorted by start, is preserved).
    """
    # Split and sort lists by seqid and start position
    dict_filter = {}
    dict_gr = {}
    for grange in gr_filter:
        dict_filter.setdefault(grange.seqid, []).append(grange)
    for grange in gr:
        dict_gr.setdefault(grange.seqid, []).append(grange)

    for seqid in dict_gr:
        dict_gr[seqid].sort(key=lambda x: x.start)

    filtered = []

    for seqid in dict_filter:
        if seqid not in dict_gr:
            continue
        # Index this seqid's filter intervals by name: starts sorted ascending
        # plus the running maximum of ends (prefix max), so a binary search on
        # start gives the best containment candidate.
        by_name = {}
        for f in dict_filter[seqid]:
            by_name.setdefault(f.name, []).append(f)
        name_index = {}
        for name, flist in by_name.items():
            flist.sort(key=lambda x: x.start)
            starts = [f.start for f in flist]
            prefix_max_end = []
            running = None
            for f in flist:
                running = f.end if running is None else max(running, f.end)
                prefix_max_end.append(running)
            name_index[name] = (starts, prefix_max_end)

        for grange_out in dict_gr[seqid]:
            idx = name_index.get(grange_out.name)
            if idx is None:
                continue
            starts, prefix_max_end = idx
            pos = bisect.bisect_right(starts, grange_out.start) - 1
            if pos >= 0 and prefix_max_end[pos] >= grange_out.end:
                filtered.append(grange_out)

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
        self.length = self.end - self.start + 1
        self.strand = items[8]
        self.annot = items[10]
        self.refid = items[9]
        self.family = items[11]


def merged_interval_length(intervals):
    """Length of the union of (start, end) intervals (inclusive, 1-based).

    Overlapping / adjacent intervals are merged so the total never
    double-counts a position. Used for coverage fractions (covered length
    divided by a sequence length stays in [0, 1]).
    """
    total = 0
    cur_start = cur_end = None
    for start, end in sorted(intervals):
        if cur_end is None or start > cur_end + 1:
            if cur_end is not None:
                total += cur_end - cur_start + 1
            cur_start, cur_end = start, end
        else:
            cur_end = max(cur_end, end)
    if cur_end is not None:
        total += cur_end - cur_start + 1
    return total


def get_repeatmasker_annotation(rm_file, seq_lengths, prefix, parse_id=True):
    """
    :parse repeatmasker output and calculate proportion of each annotation
    :param rm_file:
    :param seq_lengths:
    :param prefix:
    :param parse_id: if True, parse ID from repeatmasker output, if False, use as it is
    :return:
    """
    # Collect, for every consensus sequence and every reference annotation, the
    # set of query intervals covered by RepeatMasker hits. The coverage of an
    # annotation is the length of the UNION of its hit intervals divided by the
    # sequence length, so multiple monomer hits along a tandem consensus
    # accumulate (rather than only the first hit being counted) while
    # overlapping hits are not double-counted (each annotation stays in [0, 1]).
    seq_annot_intervals = {}
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
            interval = (min(rm_feature.start, rm_feature.end),
                        max(rm_feature.start, rm_feature.end))
            seq_annot_intervals.setdefault(
                rm_feature.seqid, {}).setdefault(
                rm_feature.annot, []).append(interval)
    seq_rm_info = {}
    for seqid, annot_intervals in seq_annot_intervals.items():
        seq_rm_info[seqid] = {
            annot: merged_interval_length(intervals) / seq_lengths[seqid]
            for annot, intervals in annot_intervals.items()
        }
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


# ---------------------------------------------------------------------------
# rDNA identification
#
# Similarity-search TRC consensus (or, as a fallback, genomic arrays) against a
# reference library of rRNA genes in RepeatMasker `name#class` format (classes
# `rDNA_45S/{18S,5.8S,25S}` and `rDNA_5S/5S`). A TRC is called rDNA from the
# *best single-subunit reference coverage*: it is rDNA_45S/5S if one library
# reference (e.g. an 18S or 25S gene) is matched near-full-length at high
# identity. This is robust to the ITS/IGS spacer (absent from the library)
# diluting whole-monomer coverage and to partial consensus. The search engine
# is blastn because the chosen metric is reference coverage, which blastn yields
# directly (query = reference) — the same blast tooling the superfamily step
# uses.
# ---------------------------------------------------------------------------

RDNA_TOP_TYPES = ("rDNA_45S", "rDNA_5S")


def blastn_rdna_reference_coverage(rdna_library, subject_fasta, cpu,
                                   min_identity=85.0):
    """
    blastn the rDNA library (as query) against ``subject_fasta`` and return, per
    subject sequence, the best single-reference coverage for each rDNA top type.

    Reference coverage = union of the aligned query (reference) intervals for one
    library entry, divided by that entry's length. Per subject and per top type
    (rDNA_45S / rDNA_5S) we keep the maximum reference coverage over all
    references of that type whose best hit meets ``min_identity``.

    :param rdna_library: rDNA reference FASTA (RepeatMasker `name#class` headers)
    :param subject_fasta: sequences to test (TRC consensus or genomic arrays)
    :param cpu: blastn threads
    :param min_identity: minimum percent identity for a reference to count
    :return: {subject_seqid: {top_type: (coverage_frac, identity)}}
    """
    ref_len = read_fasta_sequence_size(rdna_library)
    ref_top = {}
    for h in ref_len:
        cls = h.split("#", 1)[1] if "#" in h else ""
        ref_top[h] = cls.split("/")[0]
    work = tempfile.mkdtemp(prefix="tc_rdna_blast_")
    db = os.path.join(work, "subjdb")
    out = os.path.join(work, "hits.tsv")
    subprocess.run(F"makeblastdb -in {subject_fasta} -dbtype nucl -out {db}",
                   shell=True, check=True, stdout=subprocess.DEVNULL)
    subprocess.run(
        F"blastn -query {rdna_library} -db {db} -evalue 1e-5 -dust no "
        F"-num_threads {cpu} -outfmt '6 qseqid sseqid pident qstart qend' > {out}",
        shell=True, check=True)
    # per (subject, reference): aligned reference intervals + best identity
    intervals = {}
    identity = {}
    with open(out) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                continue
            qseqid, sseqid, pid = c[0], c[1], float(c[2])
            qs, qe = int(c[3]), int(c[4])
            key = (sseqid, qseqid)
            intervals.setdefault(key, []).append((min(qs, qe), max(qs, qe)))
            identity[key] = max(identity.get(key, 0.0), pid)
    shutil.rmtree(work, ignore_errors=True)
    best = {}
    for (sseqid, qseqid), ivs in intervals.items():
        if identity[(sseqid, qseqid)] < min_identity:
            continue
        top = ref_top.get(qseqid)
        if top not in RDNA_TOP_TYPES:
            continue
        cov = merged_interval_length(ivs) / ref_len[qseqid]
        cur = best.setdefault(sseqid, {})
        if cov > cur.get(top, (0.0, 0.0))[0]:
            cur[top] = (cov, identity[(sseqid, qseqid)])
    return best


def assign_rdna_to_trcs(best_by_seqid, min_coverage):
    """
    Aggregate per-sequence best reference coverage (from
    :func:`blastn_rdna_reference_coverage`) to per-TRC rDNA calls.

    Subject seqids are ``TRC_<n>_...`` so the TRC id is
    ``"TRC_" + seqid.split("_")[1]``. The maximum coverage per top type across a
    TRC's sequences is taken, then the type with the higher coverage is called
    if it clears ``min_coverage``.

    :return: {TRC_id: (label, coverage)} e.g. {'TRC_1': ('45S', 0.92)}
    """
    per_trc = {}
    for sseqid, tops in best_by_seqid.items():
        trc = "TRC_" + sseqid.split("_")[1]
        agg = per_trc.setdefault(trc, {})
        for top, (cov, _ident) in tops.items():
            if cov > agg.get(top, 0.0):
                agg[top] = cov
    calls = {}
    for trc, agg in per_trc.items():
        top = max(agg, key=agg.get)
        if agg[top] >= min_coverage:
            calls[trc] = (top.replace("rDNA_", ""), round(agg[top], 3))
    return calls


def _write_trc_consensus_subject(prefix, consensus_dir, out_fasta):
    """Write a blast subject of per-TRC consensus dimers with clean
    ``TRC_<n>_c<i>`` headers. Prefers the TAREAN dimer library
    (``<prefix>_consensus_dimer_library.fasta``, headers ``TRC_x#TRC_x``) and
    falls back to the raw clustering consensus (``TRC_x_dimers.fasta``) per TRC.
    Returns the set of TRC ids that had a usable consensus."""
    tarean = {}
    dimer_lib = prefix + "_consensus_dimer_library.fasta"
    if os.path.exists(dimer_lib):
        for name, seq in fasta_to_list(dimer_lib):
            tarean.setdefault(name.split("#")[0], []).append(seq)
    clustering = {}
    for fpath in glob.glob(os.path.join(consensus_dir, "TRC_*_dimers.fasta")):
        trc = os.path.basename(fpath)[:-len("_dimers.fasta")]
        clustering[trc] = [s for _n, s in fasta_to_list(fpath)]
    covered = set()
    with open(out_fasta, "w") as out:
        for trc in sorted(set(tarean) | set(clustering)):
            seqs = tarean.get(trc) or clustering.get(trc) or []
            for i, seq in enumerate(seqs):
                out.write(F">{trc}_c{i}\n{seq}\n")
            if seqs:
                covered.add(trc)
    return covered


def _write_trc_region_subject(gff, fasta, trc_set, out_fasta):
    """Write a blast subject of the genomic arrays of the given TRCs, with
    ``TRC_<n>_r<i>`` headers. Used as the last-resort fallback for TRCs without
    a usable consensus."""
    tmp_gff = out_fasta + ".gff3"
    with open(gff) as fin, open(tmp_gff, "w") as fo:
        for line in fin:
            if line.startswith("#"):
                fo.write(line)
                continue
            if Gff3Feature(line).attributes_dict.get("Name") in trc_set:
                fo.write(line)
    n = 0
    counts = {}
    with open(out_fasta, "w") as out:
        for _seqid, seq, name in gff3_to_fasta(tmp_gff, fasta, "Name"):
            i = counts.get(name, 0)
            counts[name] = i + 1
            out.write(F">{name}_r{i}\n{seq}\n")
            n += 1
    os.remove(tmp_gff)
    return n


def _write_rdna_attributes(gff, calls):
    """Add ``rDNA_type`` and ``rDNA_coverage`` attributes to the matched TRCs in
    ``gff`` (rewritten in place). TRCs absent from ``calls`` are left as-is."""
    if not calls:
        return
    type_dict = {trc: lbl for trc, (lbl, _cov) in calls.items()}
    cov_dict = {trc: str(cov) for trc, (_lbl, cov) in calls.items()}
    tmp = gff + ".rdna.tmp"
    add_attribute_to_gff(gff, tmp, "Name", "rDNA_type", type_dict)
    add_attribute_to_gff(tmp, gff, "Name", "rDNA_coverage", cov_dict)
    os.remove(tmp)


def identify_rdna(prefix, fasta, rdna_library, cpu, min_coverage=0.7,
                  min_identity=85.0, gff=None, consensus_dir=None):
    """
    Identify which TRCs are rDNA (45S / 5S) and write ``rDNA_type`` +
    ``rDNA_coverage`` attributes into the clustering GFF3 (and the annotation
    GFF3 if present, so the report picks them up).

    Strategy (consensus-first, genome fallback):
      1. blastn the rDNA library against per-TRC consensus dimers (TAREAN dimer
         library where available, else the raw TideHunter clustering consensus).
      2. For TRCs with no usable consensus, blastn against their genomic arrays
         extracted from the GFF3.
    A TRC is labelled from the best single-subunit reference coverage
    (>= ``min_coverage`` at >= ``min_identity``); see
    :func:`blastn_rdna_reference_coverage`.

    :return: {TRC_id: (label, coverage)}
    """
    gff = gff or (prefix + "_clustering.gff3")
    consensus_dir = consensus_dir or (prefix + "_consensus")
    work = tempfile.mkdtemp(prefix="tc_rdna_")
    all_trcs = set()
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            all_trcs.add(Gff3Feature(line).attributes_dict.get("Name"))
    all_trcs.discard(None)

    calls = {}
    # 1) consensus-based search
    subject = os.path.join(work, "consensus.fasta")
    covered = _write_trc_consensus_subject(prefix, consensus_dir, subject)
    if os.path.getsize(subject) > 0:
        best = blastn_rdna_reference_coverage(rdna_library, subject, cpu,
                                              min_identity)
        calls.update(assign_rdna_to_trcs(best, min_coverage))

    # 2) genomic fallback for TRCs that had no usable consensus
    no_consensus = all_trcs - covered
    if no_consensus:
        regions = os.path.join(work, "regions.fasta")
        if _write_trc_region_subject(gff, fasta, no_consensus, regions) > 0:
            best_g = blastn_rdna_reference_coverage(rdna_library, regions, cpu,
                                                    min_identity)
            calls.update(assign_rdna_to_trcs(best_g, min_coverage))

    # 3) write attributes into the GFF3(s) the report reads
    _write_rdna_attributes(gff, calls)
    annot_gff = prefix + "_annotation.gff3"
    if os.path.exists(annot_gff):
        _write_rdna_attributes(annot_gff, calls)

    # 4) side table
    with open(prefix + "_rdna.tsv", "w") as f:
        f.write("TRC\trDNA_type\tcoverage\n")
        for trc in sorted(calls):
            lbl, cov = calls[trc]
            f.write(F"{trc}\t{lbl}\t{cov}\n")

    shutil.rmtree(work, ignore_errors=True)
    n45 = sum(1 for v in calls.values() if v[0] == "45S")
    n5 = sum(1 for v in calls.values() if v[0] == "5S")
    print(F"rDNA identification: {len(calls)} TRC(s) labelled "
          F"({n45} x 45S, {n5} x 5S)")
    return calls


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


def split_fasta_to_chunk_files(fasta_file, out_dir, chunk_size=50000000,
                               overlap=100000):
    """
    Split a (genome) FASTA into multiple chunk FASTA files for parallel,
    per-chunk RepeatMasker.

    Sequences longer than ``2 * chunk_size`` are cut into overlapping pieces of
    about ``chunk_size`` bases (the overlap is added to the right edge of each
    interior piece so a repeat straddling a cut is fully contained in at least
    one piece); shorter sequences are emitted whole. Every piece is given a
    short numeric header to avoid RepeatMasker's 50-character sequence-name
    limit, and recorded in a ``matching_table`` compatible with
    :func:`get_original_header_and_coordinates` so chunk-local hit coordinates
    can be mapped back to genome coordinates.

    Pieces are packed greedily into files of roughly ``chunk_size`` bases, so
    the number of RepeatMasker processes stays bounded while each carries a
    similar amount of work.

    :param fasta_file: path to input FASTA
    :param out_dir: directory to write chunk FASTA files into
    :param chunk_size: target piece / file size in bases
    :param overlap: bases of right-edge overlap added to each interior piece
    :return: (list_of_chunk_fasta_paths, matching_table) where matching_table
             rows are [orig_header, piece_index, start, end, token]
    """
    fasta_dict = read_fasta_sequence_size(fasta_file)
    min_chunk_size = chunk_size * 2
    # 1) plan pieces; each piece gets a short numeric token as its header
    matching_table = []
    token = 0
    for header, size in fasta_dict.items():
        if size > min_chunk_size:
            number_of_chunks = int(size / chunk_size)
            adjusted_chunk_size = int(size / number_of_chunks)
            for i in range(number_of_chunks):
                start = i * adjusted_chunk_size
                end1 = (i + 1) * adjusted_chunk_size + overlap
                end = end1 if i + 1 < number_of_chunks else size
                matching_table.append([header, i, start, end, str(token)])
                token += 1
        else:
            matching_table.append([header, 0, 0, size, str(token)])
            token += 1
    # 2) bin-pack pieces into files of ~chunk_size bases
    token_to_file = {}
    file_paths = []
    current_load = 0
    for row in matching_table:
        piece_len = row[3] - row[2]
        if not file_paths or current_load + piece_len > chunk_size:
            file_paths.append(os.path.join(out_dir, F"chunk_{len(file_paths)}.fasta"))
            current_load = 0
        token_to_file[row[4]] = file_paths[-1]
        current_load += piece_len
    # group pieces per original header for a single streaming pass
    pieces_by_header = {}
    for row in matching_table:
        pieces_by_header.setdefault(row[0], []).append(row)
    # 3) stream the FASTA once, writing each piece to its assigned file
    open_handles = {p: open(p, "w") for p in file_paths}
    try:
        with open(fasta_file) as fh:
            for header, sequence in read_single_fasta_as_generator(fh):
                for _orig, _i, start, end, tok in pieces_by_header.get(header, []):
                    out_fh = open_handles[token_to_file[tok]]
                    out_fh.write(F">{tok}\n")
                    out_fh.write(sequence[start:end] + "\n")
    finally:
        for fh in open_handles.values():
            fh.close()
    return file_paths, matching_table


def _repeatmasker_chunk_worker(task):
    """Run a single-threaded RepeatMasker on one chunk FASTA.

    Module-level (picklable) worker for the process pool used by
    :func:`run_repeatmasker_genome_chunked`. Each chunk is masked in its own
    output directory so concurrent jobs cannot clobber each other's files.

    :param task: (chunk_fasta, rm_library, sensitivity_flag)
    :return: path to the chunk's ``.out`` file, or None if RepeatMasker wrote
             none (an empty / hit-less / all-N chunk).
    """
    chunk_fasta, rm_library, sensitivity_flag = task
    out_dir = chunk_fasta + "_rmdir"
    os.makedirs(out_dir, exist_ok=True)
    cmds = ["RepeatMasker", "-dir", out_dir, "-nolow", "-no_is", "-e", "ncbi"]
    if sensitivity_flag:
        cmds.append(sensitivity_flag)
    cmds.extend(["-lib", rm_library, "-pa", "1", chunk_fasta])
    # Tolerate the no-hit case the way run_repeatmasker_with_renaming does:
    # RepeatMasker writes no .out for an empty / all-N / hit-less input. Do not
    # use check=True, so such chunks return None instead of aborting the pool.
    subprocess.run(cmds)
    out_file = os.path.join(out_dir, os.path.basename(chunk_fasta) + ".out")
    return out_file if os.path.exists(out_file) else None


def run_repeatmasker_genome_chunked(ref_seq, rm_library, cpu, sensitivity_flag,
                                    out_gff3, chunk_size=50000000,
                                    overlap=100000, debug=False):
    """
    Parallel, output-preserving replacement for a single whole-genome
    ``RepeatMasker -pa {cpu}`` call.

    RepeatMasker's ``-pa`` does not parallelise effectively with a custom
    ``-lib`` on the RMBlast/NCBI engine (Dfam #274), and its single-threaded
    ``ProcessRepeats`` tail is not parallelised by ``-pa`` at all. This splits
    ``ref_seq`` into overlapping chunks, runs one single-threaded RepeatMasker
    per chunk in a pool of ``cpu`` workers (so both the search *and* each
    chunk's ProcessRepeats run concurrently), then maps every hit back to
    genome coordinates and writes them to ``out_gff3``.

    Output equivalence: RepeatMasker treats each sequence independently, the
    chunk overlap exceeds any library-entry (hence any hit) length so every hit
    is fully contained in at least one chunk, and duplicate / partial hits in
    the overlap zones are collapsed by the downstream
    :func:`merge_overlapping_gff3_intervals`. The GFF3 written here is sorted,
    so the result does not depend on pool scheduling order.

    :param ref_seq: genome FASTA to annotate
    :param rm_library: RepeatMasker custom library (already length-adjusted)
    :param cpu: number of concurrent single-threaded RepeatMasker processes
    :param sensitivity_flag: RepeatMasker sensitivity flag ("", "-q" or "-qq")
    :param out_gff3: path to write the (genome-coordinate) RepeatMasker GFF3
    :param chunk_size: target chunk size in bases
    :param overlap: right-edge overlap per interior chunk in bases
    :param debug: if True, keep the temporary chunk workspace
    :return: out_gff3
    """
    work_dir = tempfile.mkdtemp(prefix="tc_rm_chunks_")
    chunk_files, matching_table = split_fasta_to_chunk_files(
        ref_seq, work_dir, chunk_size=chunk_size, overlap=overlap
    )
    print(F"Chunked RepeatMasker: {len(chunk_files)} chunk file(s), "
          F"{len(matching_table)} piece(s), {cpu} worker(s)")

    tasks = [(c, rm_library, sensitivity_flag) for c in chunk_files]
    pool_size = max(1, min(cpu, len(tasks)))
    with Pool(pool_size) as pool:
        out_files = pool.map(_repeatmasker_chunk_worker, tasks)

    # O(1) token -> matching_table row lookup (avoid the O(M) linear scan in
    # get_original_header_and_coordinates per hit; there can be millions).
    token_row = {row[4]: row for row in matching_table}

    records = []
    for out_file in out_files:
        if out_file is None or not os.path.exists(out_file):
            continue
        with open(out_file) as f:
            for _ in range(3):  # skip the 3 RepeatMasker header lines
                next(f, None)
            for line in f:
                items = line.split()
                if len(items) < 11:  # blank / malformed line
                    continue
                row = token_row.get(items[4])
                if row is None:
                    continue
                offset = row[2]
                ori_header = row[0]
                start = int(items[5]) + offset
                end = int(items[6]) + offset
                # mirror repeatmasker_to_gff3 strand mapping verbatim
                strand = "+" if items[8] == "C" else "-"
                name = items[10]
                records.append((ori_header, start, end, strand, name))

    # sort for a deterministic intermediate (final output is sorted downstream
    # regardless, but this keeps the GFF3 reproducible run-to-run)
    records.sort(key=lambda r: (r[0], r[1], r[2], r[3], r[4]))
    with open(out_gff3, "w") as gff3_out:
        gff3_out.write("##gff-version 3\n")
        for seqid, start, end, strand, name in records:
            gff_line = (F"{seqid}\tRepeatMasker\ttandem repeat\t{start}"
                        F"\t{end}\t.\t{strand}\t.\tCluster_ID={name};Name={name}")
            gff3_out.write(Gff3Feature(gff_line).print_line())

    if not debug:
        shutil.rmtree(work_dir, ignore_errors=True)
    else:
        print(F"Chunked RepeatMasker workspace kept: {work_dir}")
    return out_gff3


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


def run_blastn(fasta_file, dust=False, cpu=4, perc_identity=75, min_coverage=0.8):
    """
    run blastn on fasta file
    :param cpu:
    :param dust: dust filter for blastn
    :param perc_identity: minimum BLASTN percent identity for a pair to be kept
    :param min_coverage: minimum alignment coverage over the shorter sequence
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
           F" -out {blast_out} -num_threads {cpu} -evalue 1e-20 -perc_identity {perc_identity}"
           F" -word_size 9 -max_target_seqs 1000000 -dust {dust}"
           F" -gapextend 1 -gapopen 2 -reward 1 -penalty -1")
    subprocess.check_call(cmd, shell=True)
    # read pairs to list, exclude self hits and duplicates
    pairs = set()
    with open(blast_out, 'r') as f:
        for line in f:
            qseqid, sseqid, pident, length, evalue, bitscore, qlen, slen = line.split()
            if qseqid != sseqid:
                # overlap should be at least min_coverage over the shorter sequence
                if int(length) / min([int(qlen), int(slen)]) > min_coverage:
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
        cpu=4, perc_identity=75, min_coverage=0.8
        ):
    """
    find clusters by blastn, return dictionary with clusters
    cluaste are connected components in graph
    :param cpu:
    :param dust: use dust filter
    :param perc_identity: minimum BLASTN percent identity for a clustering edge
    :param min_coverage: minimum alignment coverage over the shorter sequence
    :param consensus_representative:
    :return: clusters
    """
    # check fasta file size
    if os.path.getsize(fasta_file) == 0:
        return {}

    print("Clustering by BLASTN")
    pairs = run_blastn(fasta_file, dust=dust, cpu=cpu,
                       perc_identity=perc_identity, min_coverage=min_coverage)
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


def resolve_trc_overlaps(gff3_file, gff3_out_file):
    """
    Make a clustering GFF3 non-overlapping across TRCs ("annotate each region
    once").

    Different TRCs (different ``Name``) may cover the same genomic bp — e.g.
    variant arrays of a satellite (such as rDNA) clustered separately that
    interleave and overlap at their boundaries. This rewrites the GFF3 so every
    base belongs to at most one TRC, assigning each contested span to the
    **dominant** TRC: the one with the largest total array length over the whole
    genome (tie-break: ``Name``). No base in the union is lost — only
    reassigned — and per-TRC attributes (repeat_type, ssr, rDNA_type, …) are
    carried through unchanged.

    Implementation: a breakpoint sweep per seqid. Every elementary segment
    between consecutive feature boundaries is covered by a fixed set of TRCs;
    it goes to the dominant one. Adjacent segments kept by the same (Name,
    strand) are then merged back into single features.

    :param gff3_file: input clustering GFF3 (features carry Name=TRC_x)
    :param gff3_out_file: output GFF3 (may equal the input — written via a temp)
    """
    feats = {}        # seqid -> list of (start, end, name, strand)
    attrs = {}        # name -> attributes dict (from its first feature)
    total = {}        # name -> total array length (dominance weight)
    header = "##gff-version 3\n"
    with open(gff3_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            feat = Gff3Feature(line)
            name = feat.attributes_dict.get("Name")
            if name is None:
                continue
            feats.setdefault(feat.seqid, []).append(
                (feat.start, feat.end, name, feat.strand))
            attrs.setdefault(name, dict(feat.attributes_dict))
            total[name] = total.get(name, 0) + (feat.end - feat.start)

    kept = set()
    out_rows = []
    for seqid, fl in feats.items():
        # half-open breakpoints: each feature [s, e] -> [s, e+1)
        points = sorted({p for s, e, _n, _st in fl for p in (s, e + 1)})
        segments = []
        for k in range(len(points) - 1):
            a, b = points[k], points[k + 1]
            covering = [(n, st) for s, e, n, st in fl if s <= a and (e + 1) >= b]
            if not covering:
                continue
            # dominant TRC: largest total array length, then Name for determinism
            win_name, win_strand = max(covering, key=lambda x: (total[x[0]], x[0]))
            segments.append((a, b - 1, win_name, win_strand))
        # merge adjacent segments kept by the same TRC (and strand)
        for s, e, name, strand in segments:
            if (out_rows and out_rows[-1][0] == seqid
                    and out_rows[-1][3] == name and out_rows[-1][4] == strand
                    and s == out_rows[-1][2] + 1):
                last = out_rows[-1]
                out_rows[-1] = (seqid, last[1], e, name, strand)
            else:
                out_rows.append((seqid, s, e, name, strand))
            kept.add(name)

    lost = sorted(set(total) - kept)
    if lost:
        print(F"resolve_trc_overlaps: {len(lost)} TRC(s) fully absorbed by "
              F"dominant overlapping TRCs and dropped: {', '.join(lost)}")

    out_rows.sort(key=lambda r: (r[0], r[1], r[2]))
    tmp = gff3_out_file + ".resolve.tmp"
    with open(tmp, "w") as out:
        out.write(header)
        for seqid, s, e, name, strand in out_rows:
            ad = dict(attrs[name])
            ad["Name"] = name  # keep Name as the first attribute
            ad = {"Name": name, **ad}
            attr_str = ";".join(F"{k}={v}" for k, v in ad.items())
            out.write(F"{seqid}\tTideCluster\ttandem_repeat\t{s}\t{e}\t1\t"
                      F"{strand}\t.\t{attr_str}\n")
    os.replace(tmp, gff3_out_file)


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

    # RepeatMasker writes no .out file when the input FASTA is empty or no
    # repeats are detected. Treat that as "zero annotations" (an empty .out
    # that get_repeatmasker_annotation parses to {}) instead of crashing on the
    # missing file.
    if os.path.exists(rm_file_renamed):
        # Restore original sequence names in RepeatMasker output
        restore_sequence_names_in_repeatmasker_output(rm_file_renamed, rm_file_original, name_mapping)
        os.remove(rm_file_renamed)
    else:
        open(rm_file_original, 'w').close()

    # Clean up renamed files
    if os.path.exists(input_renamed):
        os.remove(input_renamed)

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


# ----------------------------------------------------------------------
# kitehor IO helpers (replacement for tarean/kite.R)
# ----------------------------------------------------------------------

def build_kite_multifasta(fasta_dir, out_fa):
    """Concatenate every TRC_*.fasta in `fasta_dir` into one multi-FASTA
    suitable for `kitehor analyze`.

    The per-TRC fastas under `{prefix}_tarean/fasta/` store sequences as
    dimers (sequence concatenated with itself). kitehor expects one TR
    array per record, so each sequence is halved back to its single-array
    length. Headers are rewritten to `<TRC>:<original_header>` so the
    TRC ID survives kitehor's per-record output and `tc_rerender_report`
    can split it back."""
    n_records = 0
    with open(out_fa, "w") as fh_out:
        for path in sorted(glob.glob(os.path.join(fasta_dir, "TRC_*.fasta"))):
            trc = os.path.basename(path)[:-len(".fasta")]
            with open(path) as fh:
                for hdr, seq in read_single_fasta_as_generator(fh):
                    half = seq[:len(seq) // 2] if len(seq) >= 2 else seq
                    if not half:
                        continue
                    fh_out.write(f">{trc}:{hdr}\n{half}\n")
                    n_records += 1
    return n_records


_LONG_RESCORE_ID_GATE = 0.7   # below-cap best id_med under which we re-search


def extend_long_period_rescore(kite_dir, multi_fa, max_period, ext_max_period,
                               top_n, threads):
    """Selectively re-rescore arrays whose dominant monomer exceeds `max_period`.

    rescore emits `identity_med = NA` for periods above its `--max-period` cap,
    so an array whose true monomer is longer than the cap (e.g. a ~16 kb
    satellite under a 15 kb cap) falls back to a spurious short peak. This finds
    arrays that have a peak above `max_period` AND no confident founder below the
    cap (best below-cap `identity_med` < _LONG_RESCORE_ID_GATE), re-runs
    `kitehor rescore` for *only those* records at `ext_max_period`, and merges the
    new rows back into `kitehor.rescored.peaks.tsv`. Returns the re-rescored
    case_ids. rescore is O(period²) so confining the high cap to the few flagged
    arrays keeps the cost negligible (measured ~3 s per 112 kb array)."""
    rescored = F"{kite_dir}/kitehor.rescored.peaks.tsv"
    if not os.path.exists(rescored):
        return []
    by_case = collections.OrderedDict()
    with open(rescored, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        header = rd.fieldnames
        for r in rd:
            by_case.setdefault(r["case_id"], []).append(r)

    def _n(v):
        try: return float(v)
        except (TypeError, ValueError): return None

    flagged = set()
    for cid, peaks in by_case.items():
        best_id = max([(_n(p.get("identity_med")) or -1.0) for p in peaks],
                      default=-1.0)
        above_cap = any((_n(p.get("period")) or 0.0) > max_period for p in peaks)
        if above_cap and best_id < _LONG_RESCORE_ID_GATE:
            flagged.add(cid)
    if not flagged:
        return []

    sub_fa = F"{kite_dir}/_kite_input_longext.fasta"
    with open(multi_fa) as fin, open(sub_fa, "w") as fout:
        keep = False
        for line in fin:
            if line.startswith(">"):
                keep = line[1:].strip().split()[0] in flagged
            if keep:
                fout.write(line)
    kite_peaks = F"{kite_dir}/kitehor.kite.peaks.tsv"
    sub_peaks = F"{kite_dir}/_kite_peaks_longext.tsv"
    with open(kite_peaks, newline="") as fin, open(sub_peaks, "w", newline="") as fout:
        rd = csv.reader(fin, delimiter="\t")
        hrow = next(rd)
        ci = hrow.index("case_id")
        fout.write("\t".join(hrow) + "\n")
        for row in rd:
            if row and row[ci] in flagged:
                fout.write("\t".join(row) + "\n")

    out_prefix = F"{kite_dir}/_rescored_longext"
    run_cmd(F"kitehor rescore --peaks {sub_peaks} --out {out_prefix}"
            F" --max-period {ext_max_period} --top-n {top_n}"
            F" --threads {threads} {sub_fa}")

    ext_rows = collections.defaultdict(list)
    with open(F"{out_prefix}.peaks.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ext_rows[r["case_id"]].append(r)
    with open(rescored, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for cid, peaks in by_case.items():
            w.writerows(ext_rows.get(cid, peaks) if cid in flagged else peaks)
    return sorted(flagged)


def parse_trc_ssr_motif_len(gff3_path):
    """Per-TRC SSR motif length from the clustering GFF3's `repeat_type=SSR`
    classification. SSR families are decided at clustering — a TRA *consensus*
    that is > 90 % simple-sequence is pulled out, grouped by motif, and tagged
    `repeat_type=SSR` + `ssr=MOTIF (pct%)` (e.g. `ssr=ATC (98.6%)`). Returns
    {TRC_ID: len(MOTIF)} for SSR TRCs only (ATC -> 3, AAACCCT -> 7); TR TRCs and
    older runs without the attributes yield an empty/partial map. This is the
    authoritative SSR founder signal — kitehor's per-array ssr-scan is secondary
    (it annotates SSR content inside any TRC, e.g. a 500 bp founder that is also
    SSR-rich)."""
    out = {}
    if not gff3_path or not os.path.exists(gff3_path):
        return out
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or "repeat_type=SSR" not in line:
                continue
            attrs = line.rstrip("\n").split("\t")[-1]
            m_name = re.search(r"Name=([^;]+)", attrs)
            m_ssr  = re.search(r"ssr=\s*([A-Za-z]+)", attrs)
            if m_name and m_ssr:
                out[m_name.group(1)] = len(m_ssr.group(1))
    return out


def _hor_order_confidence(multiplicity, multiplicity_raw, irregular,
                          founder_method, fallback):
    """Confidence that the array's ×k *higher-order order* is real, as distinct
    from "founder recovered" (the basic monomer is reliable but the ×k may be an
    estimate). See docs/hor_order_confidence_design.md. Pure function of columns
    the row already carries — it never changes founder selection.

      none      multiplicity <= 1 — founder == strongest, no HOR claim.
      strict    clean integer divisor, low-k, strict path → confident HOR order.
      supported family/anchor/ladder-rescued (pass2 / kh_deeper / ladder),
                clean, low-k — order plausible, leans on corroborating evidence.
      weak      founder solid but order not confidently integer: irregular
                multiplicity, OR k >= _HOR_ORDER_HIGH_K (integer test vacuous),
                OR a relaxed rescue path (pass3 / cluster / fallback).

    For reporting, "HOR" = strict ∪ supported."""
    try:
        m = int(multiplicity)
    except (TypeError, ValueError):
        return "none"
    if m <= 1:
        return "none"
    if fallback or irregular:
        return "weak"
    k = multiplicity_raw if multiplicity_raw is not None else m
    if k >= _HOR_ORDER_HIGH_K:
        return "weak"
    if founder_method == "strict":
        return "strict"
    if founder_method in ("pass2", "kh_deeper", "ladder"):
        return "supported"
    return "weak"   # pass3 / cluster / fallback / other relaxed rescue


_FOUNDER_ID_MIN  = 0.7    # rescore --subrepeat-founder-id-min default
_KMAX            = 30     # rescore --rule-k-max default
_RATIO_TOL       = 0.05   # how close to integer multiplicity must be
# Prevalent-founder anchor (Lever 4) — k-scaled integer tolerance for deepening a
# long-HOR founder to the TRC consensus when the array's dominant peak IS the
# consensus. The bp residual |S - m*P| ~ period-estimation error, so the
# tolerance on |k - round(k)| grows ~linearly with k (measured on confident
# TRC_4 divisors: dev/k ≈ 0.3 % constant across k=2..150). Above k ≈ 0.5/frac the
# integer test is vacuous (every k is within 0.5 of an integer) -> trust the
# family consensus instead. frac ≈ 3x the measured scatter.
_RATIO_TOL_FRAC  = 0.01   # per-k integer tolerance for the consensus anchor
_CONS_ANCHOR_BP  = 5      # rank-1 must be within max(bp, pct) of the consensus
_CONS_ANCHOR_PCT = 0.05
# HOR-order confidence (see docs/hor_order_confidence_design.md). Above this k the
# integer-multiplicity test is vacuous (every k is within 0.5 of an integer), so a
# clean-looking ×k is a *scaling estimate*, not a verified higher-order order ->
# downgrade to "weak". Tied to _RATIO_TOL_FRAC so it tracks the same constant the
# prevalent-founder anchor uses (0.5 / 0.01 = 50).
_HOR_ORDER_HIGH_K = round(0.5 / _RATIO_TOL_FRAC)
_SUBREP_RATIO    = 0.33   # period <= founder/3 ⇒ qualifies as candidate band

# Harmonic-ladder founder (Pass 7) — recover the basic monomer of a divergent HOR
# satellite when it sits below the strict id gate but a clean harmonic ladder
# (peaks at integer multiples m·P0) corroborates it. See
# docs/harmonic_ladder_founder_plan.md. Gates are deliberately conservative;
# tools/founder_diff.py is the blast-radius check.
_LADDER_KMIN      = 2     # strongest must be >= 2x the fundamental (k=2 only via the exceptional path)
_LADDER_MIN_RUNGS = 3     # >=3 distinct integer multiples present (incl. P0 and S) for the k>=3 path
_LADDER_ID_FLOOR  = 0.60  # P0 / rung id_med floor — a diverged real monomer, not noise
_LADDER_IQR_MAX   = 0.10  # P0 id_iqr ceiling — consistent identity (k>=3 path)
_LADDER_OCC_MIN   = 0.50  # P0 scan_occupancy floor — tiles a real fraction (k>=3 path)
# Exceptionally-clean ×2 path (only 2 rungs ⇒ stricter, plus a HOR-conservation test):
_LADDER_X2_IQR_MAX   = 0.05   # P0 id_iqr ceiling (stricter)
_LADDER_X2_OCC_MIN   = 0.85   # P0 scan_occupancy floor (stricter)
_LADDER_X2_HOR_DELTA = 0.10   # id_med(2*P0) - id_med(P0): the double must be a genuine HOR unit

# --- Short-founder review aids (kitehor >= 0.13.2) ---------------------------
# kitehor 0.13.1's "shorter-monomer scanning" computed a pairwise identity_med
# for very short periods (P≈5-20) that was high *by chance* on low-complexity
# arrays (skewed composition ⇒ frequent matches) with no real periodicity
# (kite `score` ≈ 0.005); these cleared the identity_med >= 0.7 founder gate and
# collapsed real satellite monomers to microsatellite harmonics. kitehor 0.13.2
# fixed this at source — those phantom short periods now get a realistic LOW
# identity_med (≈0.54, below the gate), so TideCluster's founder logic no longer
# adopts them and needs no suppression gate of its own. What remains useful is
# *visibility*: a short founder kept on weak kite support is worth a manual look,
# so we flag it and surface its dominant longer-period alternative. These are
# pure diagnostics — they never change the founder call.
_FOUNDER_SCORE_MIN      = 0.20  # kite-score floor — at/above ⇒ a credible peak
_FOUNDER_SCAN_ID_STRONG = 0.80  # scan_identity_med at/above ⇒ a credible peak
_SHORT_FOUNDER_MAX      = 30    # founder period (bp) at/below which a founder is
                                # "short" (flagged + given a longer alternative)


def _peak_credible(peak):
    """True when *peak* shows real periodicity support — kite ``score`` >=
    _FOUNDER_SCORE_MIN OR ``scan_identity_med`` >= _FOUNDER_SCAN_ID_STRONG (or
    either value is absent). Used only to pick a *credible* longer-period
    alternative for the short-founder review column; never gates the founder."""
    sc  = _num_or_none(peak.get("score"))
    sid = _num_or_none(peak.get("scan_identity_med"))
    if sc is None or sid is None:
        return True
    return sc >= _FOUNDER_SCORE_MIN or sid >= _FOUNDER_SCAN_ID_STRONG

# --- kitehor >= 0.13.0 "go-deeper" founder adoption (Pass 1b) ---
# kitehor's rescore now emits its own per-array HOR decomposition
# (`hor_basic_period`, constant per case_id) using a strict + scan-relaxed
# divisor search. When that basic is *deeper* than TideCluster's strict
# Pass-1 pick AND the rescored peak at that period is a clean, near-full
# array-wide tandem, TideCluster adopts it. The coverage/occupancy gate is
# the sharp discriminator (validated on drapa): real deeper monomers tile
# the whole array (cov >= 0.75, scan_occ >= 0.90), whereas kitehor's
# spurious near-2x over-splits have ~0 coverage and are rejected.
_KH_DEEPER_COV_MIN   = 0.75   # coverage_frac of the deeper-basic peak
_KH_DEEPER_OCC_MIN   = 0.90   # scan_occupancy_frac of the deeper-basic peak
_KH_DEEPER_RATIO_TOL = 0.10   # |k - round(k)| tol (relaxed; high-k HORs)
# kitehor's basic must be a MEANINGFULLY smaller monomer, not a ~few-%
# lateral re-estimate of the same period. The confirmed deep-HOR cases
# drop the founder by >=25 % (mult rises >=1.33x); spurious lateral
# shifts move it only ~5 %. Require at least a 20 % drop.
_KH_DEEPER_MAX_RATIO = 0.80   # kh_basic <= 0.80 * current founder

# --- tandem-validate (kitehor >= 0.13.0) subrepeat adoption gate ---
# We surface a localized_subrepeat only when it is a genuine partial-
# occupancy nested motif that does NOT coincide with TideCluster's founder.
_TV_DENSITY_MIN   = 0.10   # min occupancy (density) to count as "strong"
_TV_PRESENCE_MIN  = 0.10   # min n_windows_present / n_windows_total
_TV_FOUNDER_TOL   = 0.05   # cand ~= founder ⇒ it IS the founder, suppress
# Edge diagnostic: a candidate that is a clean integer divisor of the
# founder at high occupancy may indicate the founder itself is miscalled.
_TV_EDGE_DENSITY  = 0.50


def _classify_subrepeat_tier(peak, founder_period):
    """Per-peak subrepeat tier, ported from kitehor's decision tree
    (docs/rule_proto.md). Reads only rescored peak columns:
       HIGH               period<=founder/4, scan_occ>=0.15,
                          AND (subrepeat=true OR phaseC>=0.10 OR autoF>=0.4)
       LIKELY             period<=founder/4, scan_occ>=0.20, scan_n>=10
       AMBIGUOUS          founder/4 < period <= founder/3 with mild support
       OBSERVATIONAL      founder is NA but scan_occ >= 0.05
       REJECT_*           any of the disqualifying rules (phantom,
                          ratio > 1/3, founder-itself, no tandem run, ...)
    The report's main 'Subrepeat' cell surfaces HIGH+LIKELY only; the
    Details panel shows everything that wasn't rejected outright."""
    def _num(v):
        if v in (None, "", "NA"): return None
        try: return float(v)
        except (TypeError, ValueError): return None

    if peak.get("phantom", "").strip() == "true":
        return "REJECT_PHANTOM"
    period = _num(peak.get("period"))
    if period is None or period <= 0:
        return "REJECT_BAD"
    occ    = _num(peak.get("scan_occupancy_frac"))
    if founder_period is None or founder_period <= 0:
        # rule 4 — rescore over-flagged, no founder in record
        if occ is None or occ < 0.05:
            return "REJECT_NO_FOUNDER"
        return "OBSERVATIONAL"
    ratio = period / founder_period
    idm    = _num(peak.get("identity_med"))
    covf   = _num(peak.get("coverage_frac"))
    scan_n = _num(peak.get("scan_n_intervals"))
    phaseC = _num(peak.get("kmer_phase_contrast"))
    autoF  = _num(peak.get("kmer_autocorr_founder"))
    srflag = peak.get("subrepeat", "").strip() == "true"
    # rule 1d — this peak IS a clean array-wide tandem at its own period
    if (idm is not None and idm >= 0.85 and
        covf is not None and covf >= 0.95 and
        occ  is not None and occ  >= 0.95):
        return "REJECT_IS_FOUNDER"
    if ratio > _SUBREP_RATIO:
        return "REJECT_TOO_LARGE"
    if ratio <= 0.25:
        # HIGH: per-base scan confirms tandem AND an independent signal
        # (rescore subrepeat flag, phaseC, or autoF) agrees.
        if (occ is not None and occ >= 0.15
            and (srflag
                 or (phaseC is not None and phaseC >= 0.10)
                 or (autoF  is not None and autoF  >= 0.4))):
            return "HIGH"
        # LIKELY: per-base scan finds it even without alignment support.
        if (occ is not None and occ >= 0.20
            and scan_n is not None and scan_n >= 10):
            return "LIKELY"
        # KMER_SUPPORT: scan didn't catch it (e.g. interspersed not
        # contiguous), but k-mer autocorrelation / phase contrast fires.
        # Cheat-sheet calls this 'probably real'; we surface it in the
        # Details panel rather than the main Subrepeat cell.
        if ((phaseC is not None and phaseC >= 0.10)
            or (autoF is not None and autoF  >= 0.4)):
            return "KMER_SUPPORT"
        return "WEAK"
    # 0.25 < ratio <= 0.33 — ambiguous band
    if ((occ is not None and occ >= 0.20)
        or (phaseC is not None and phaseC >= 0.10)
        or (autoF  is not None and autoF  >= 0.4)):
        return "AMBIGUOUS"
    return "REJECT_AMBIGUOUS_LOW"


# Two-pass founder reassignment thresholds (see build_monomer_size_csv).
_TRC_CONSENSUS_MIN_ARRAYS = 3      # min Pass-1 founders to compute consensus
_TRC_CONSENSUS_MIN_FRAC   = 0.25   # dominant cluster must hold >= this fraction
_TRC_CONSENSUS_TOL        = 0.10   # ±10% to belong to the dominant period cluster
_RESCUE_ID_MIN            = 0.5    # relaxed id_med gate when TRC consensus exists
_RESCUE_WINDOW_PCT        = 0.10   # ± this fraction of trc_consensus_founder
_RESCUE_WINDOW_BP_MIN     = 20     # ... or this many bp, whichever is larger
_RELAXED_RATIO_TOL        = 0.20   # solo-TRA fallback: relaxed k tolerance

# (The former coverage-based SSR-founder override was removed: SSR founders now
# come from the clustering repeat_type classification, not per-array kitehor SSR
# coverage. The threshold constant that drove it is gone.)
_RELAXED_ID_MIN           = 0.7    # solo-TRA fallback: keeps strict id_med gate

# Pass 5 — peak-cluster fallback rescue (kitehor.rescored.peaks.tsv).
# When rescore returns founder_period=NA the legacy fallback set
# founder = rank-1 peak by single-peak score. On noisy centromeric
# arrays this can pick a spurious short-period peak (an SSR-like
# sub-monomer or off-target alignment hit) over the real ~10 kb HOR
# signal that is fragmented across many adjacent rescored periods
# (e.g. 9383, 9386, 10461, 10466, 10468, 10689). Pass 5 clusters
# nearby periods, sums their kite scores, and accepts the dominant
# cluster as founder when it also passes coverage / spatial-contrast
# gates. Fires only on rows where rescore returned NA (fallback=True)
# and no SSR override applied — Pass-1 successes / Pass-2/3 rescues
# are not touched.
_CLUSTER_WINDOW_PCT       = 0.05   # ±5 % single-link relative window
_CLUSTER_WINDOW_BP_MIN    = 100    # ... or ±100 bp, whichever is larger
_CLUSTER_COV_FRAC_MIN     = 0.20   # max(coverage_frac) over cluster members
_CLUSTER_SPATIAL_MIN      = 0.50   # max(spatial_contrast) over cluster members
_CLUSTER_SCORE_MARGIN     = 1.0    # cluster.score_sum >= rank1.score × this
_CLUSTER_MIN_PEAKS        = 2      # singleton clusters do not qualify
_ALT_CLUSTERS_REPORTED    = 2      # alt_cluster_1, alt_cluster_2 in CSV


def _num_or_none(v):
    if v in (None, "", "NA"):
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _parse_tv_candidates(cell):
    """Parse a tandem-validate `candidates` cell into a list of dicts.

    Each `;`-separated entry is `kind/period:d=..:sc=..:pc=..:label`,
    e.g. `other/305:d=0.084:sc=0.500:pc=0.111:localized`."""
    out = []
    if not cell or cell in ("", "NA"):
        return out
    for entry in cell.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        parts = entry.split(":")
        head = parts[0]
        if "/" not in head:
            continue
        kind, per = head.split("/", 1)
        d = sc = pc = None
        label = ""
        for f in parts[1:]:
            if f.startswith("d="):    d  = _num_or_none(f[2:])
            elif f.startswith("sc="): sc = _num_or_none(f[3:])
            elif f.startswith("pc="): pc = _num_or_none(f[3:])
            else:                     label = f
        out.append({"kind": kind, "period": _num_or_none(per),
                    "d": d, "sc": sc, "pc": pc, "label": label})
    return out


def _tv_subrepeats(tv_row, founder_period):
    """Select adopted subrepeats from a tandem-validate row, gated against
    TideCluster's founder.

    Returns `(sr_list, edge_flag, edge_k)` where `sr_list` is up to two
    `(period, density, "HIGH")` tuples that pass the strong partial-
    occupancy gate and do NOT coincide with the founder, and
    `edge_flag`/`edge_k` flag a candidate that is a clean integer divisor
    of the founder at high occupancy (the possible-founder-miscall
    watch-list). Returns empty / False when no tandem-validate row."""
    sr = []
    edge_flag = False
    edge_k = None
    if not tv_row:
        return sr, edge_flag, edge_k
    hint  = tv_row.get("decision_hint", "")
    ntot  = _num_or_none(tv_row.get("n_windows_total"))
    npres = _num_or_none(tv_row.get("n_windows_present"))
    presfrac = (npres / ntot) if (npres and ntot and ntot > 0) else None
    cands = _parse_tv_candidates(tv_row.get("candidates"))
    if not cands and _num_or_none(tv_row.get("best_candidate_period")) is not None:
        cands = [{"kind":   tv_row.get("best_candidate_kind", ""),
                  "period": _num_or_none(tv_row.get("best_candidate_period")),
                  "d":      _num_or_none(tv_row.get("density")),
                  "sc":     _num_or_none(tv_row.get("spatial_contrast")),
                  "pc":     _num_or_none(tv_row.get("phase_contrast")),
                  "label":  hint}]
    for c in cands:
        per  = c["period"]
        dens = c["d"]
        if per is None or per <= 0:
            continue
        # Coincides with founder, or kitehor labelled it the founder/host
        # ⇒ it IS the HOR subunit, not a subrepeat. Suppress.
        if founder_period and abs(per - founder_period) <= max(
                12.0, _TV_FOUNDER_TOL * max(per, founder_period)):
            continue
        if c["kind"] == "founder":
            continue
        # Edge diagnostic: a clean integer divisor of the founder at high
        # occupancy is not a partial-occupancy subrepeat — it hints the
        # founder may be miscalled. Flag it (watch-list); don't surface as
        # a subrepeat.
        if founder_period and per < founder_period:
            kk = founder_period / per
            kr = round(kk)
            if (2 <= kr <= _KMAX and abs(kk - kr) <= 0.10
                    and dens is not None and dens >= _TV_EDGE_DENSITY):
                edge_flag = True
                if edge_k is None:
                    edge_k = kr
                continue
        # Strong partial-occupancy gate.
        if "localized" not in (c["label"] or "") and hint != "localized_subrepeat":
            continue
        if dens is None or dens < _TV_DENSITY_MIN:
            continue
        if presfrac is not None and presfrac < _TV_PRESENCE_MIN:
            continue
        if founder_period and per >= founder_period:
            continue
        sr.append((per, dens, "HIGH"))
    sr.sort(key=lambda x: -(x[1] or 0.0))
    return sr[:2], edge_flag, edge_k


def _trc_consensus_founder(periods):
    """Median of the largest ±_TRC_CONSENSUS_TOL period cluster across the
    strict (Pass-1) per-TRA founder assignments in one TRC. Returns
    (consensus_period, n_supporting) or (None, 0) when the cluster is
    too small (< _TRC_CONSENSUS_MIN_ARRAYS) or doesn't hold enough of
    the TRC's reassigned arrays (< _TRC_CONSENSUS_MIN_FRAC)."""
    if len(periods) < _TRC_CONSENSUS_MIN_ARRAYS:
        return None, 0
    best_n, best_members = 0, None
    for anchor in periods:
        members = [p for p in periods
                   if abs(p - anchor) / max(anchor, 1.0) <= _TRC_CONSENSUS_TOL]
        if len(members) > best_n:
            best_n, best_members = len(members), members
    if (best_n < _TRC_CONSENSUS_MIN_ARRAYS
        or best_n / len(periods) < _TRC_CONSENSUS_MIN_FRAC):
        return None, 0
    return statistics.median(best_members), best_n


def _rescue_founder_from_trc(peaks, strongest, trc_consensus):
    """Pass 2 — TRC-aware rescue. Look for a rescored peak whose period
    is within max(±10%, ±20bp) of the TRC consensus founder, with the
    relaxed id_med gate (>= 0.5). Returns (peak, period) or None."""
    if trc_consensus is None or strongest is None:
        return None
    window = max(trc_consensus * _RESCUE_WINDOW_PCT, _RESCUE_WINDOW_BP_MIN)
    candidates = []
    for p in peaks:
        P = _num_or_none(p.get("period"))
        if P is None or P >= strongest:
            continue
        if abs(P - trc_consensus) > window:
            continue
        idm = _num_or_none(p.get("identity_med"))
        if idm is None or idm < _RESCUE_ID_MIN:
            continue
        candidates.append((abs(P - trc_consensus), P, p))
    if not candidates:
        return None
    candidates.sort()              # closest to consensus wins
    _, period, peak = candidates[0]
    return peak, period


def _relaxed_individual_founder(peaks, strongest, kmax=_KMAX):
    """Pass 3 — solo / no-consensus fallback. Relaxed ±0.20 multiplicity
    tolerance, strict id_med >= 0.7 gate. Returns (peak, P, k_round,
    k_raw) or None."""
    if strongest is None:
        return None
    cands = []
    for p in peaks:
        P = _num_or_none(p.get("period"))
        if P is None or P >= strongest:
            continue
        idm = _num_or_none(p.get("identity_med"))
        if idm is None or idm < _RELAXED_ID_MIN:
            continue
        k = strongest / P
        kr = round(k)
        if not (2 <= kr <= kmax):
            continue
        if abs(k - kr) > _RELAXED_RATIO_TOL:
            continue
        cands.append((P, p, kr, k))
    if not cands:
        return None
    cands.sort(key=lambda c: c[0])
    P, peak, kr, k = cands[0]
    return peak, P, kr, k


def _cluster_peaks_by_period(peaks):
    """Greedy single-link clustering of rescored peaks by period.

    A peak joins an existing cluster iff its period is within
    ``max(member_period * _CLUSTER_WINDOW_PCT, _CLUSTER_WINDOW_BP_MIN)``
    of *any* member already in the cluster. Peaks are processed in
    descending kite ``score`` order so the highest-scoring peak in
    each cluster anchors it (used as `highest_score_peak` for the
    Pass-5 rescue to point founder_peak/strongest_peak at).

    Returns clusters sorted by ``score_sum`` desc. Each cluster dict::

        {
          "median_period":    float (median of member periods),
          "score_sum":        float (Σ of member kite scores),
          "n_peaks":          int,
          "cov_frac_max":     float | None,
          "spatial_max":      float | None,
          "id_med_max":       float | None,
          "members":          [peak dict, ...],
          "highest_score":    float,
          "highest_score_peak": peak dict,
        }
    """
    enriched = []
    for p in peaks:
        period = _num_or_none(p.get("period"))
        score  = _num_or_none(p.get("score"))
        if period is None or period <= 0 or score is None:
            continue
        enriched.append((score, period, p))
    # Score descending; ties broken by period ascending to keep deterministic.
    enriched.sort(key=lambda t: (-t[0], t[1]))
    clusters = []
    for score, period, p in enriched:
        joined = False
        for c in clusters:
            for q_period in c["_periods"]:
                window = max(q_period * _CLUSTER_WINDOW_PCT,
                             _CLUSTER_WINDOW_BP_MIN)
                if abs(period - q_period) <= window:
                    c["members"].append(p)
                    c["_periods"].append(period)
                    c["score_sum"] += score
                    joined = True
                    break
            if joined:
                break
        if not joined:
            clusters.append({
                "members":            [p],
                "_periods":           [period],
                "score_sum":          score,
                "highest_score":      score,
                "highest_score_peak": p,
            })
    # Finalize per-cluster aggregates
    for c in clusters:
        c["n_peaks"]       = len(c["members"])
        c["median_period"] = statistics.median(c["_periods"])
        cov  = [_num_or_none(m.get("coverage_frac"))    for m in c["members"]]
        spat = [_num_or_none(m.get("spatial_contrast")) for m in c["members"]]
        idm  = [_num_or_none(m.get("identity_med"))     for m in c["members"]]
        c["cov_frac_max"] = max((v for v in cov  if v is not None), default=None)
        c["spatial_max"]  = max((v for v in spat if v is not None), default=None)
        c["id_med_max"]   = max((v for v in idm  if v is not None), default=None)
        del c["_periods"]
    clusters.sort(key=lambda c: -c["score_sum"])
    return clusters


def _cluster_rescue_founder(peaks):
    """Pass 5 — peak-cluster fallback rescue.

    Aggregates the rescored peaks into period clusters (single-link,
    mixed ±5 % / ±100 bp window) and accepts the highest-``score_sum``
    cluster as founder iff *all* of:

    - ``n_peaks >= _CLUSTER_MIN_PEAKS`` (singleton clusters skipped)
    - ``max(coverage_frac)    >= _CLUSTER_COV_FRAC_MIN``
    - ``max(spatial_contrast) >= _CLUSTER_SPATIAL_MIN``
    - ``score_sum >= rank-1 peak score * _CLUSTER_SCORE_MARGIN``

    Returns ``(cluster_dict, founder_period)`` on success, ``None``
    otherwise. ``founder_period`` is the rounded median period of the
    chosen cluster.
    """
    if not peaks:
        return None
    rank1_score = _num_or_none(peaks[0].get("score"))
    if rank1_score is None:
        return None
    clusters = _cluster_peaks_by_period(peaks)
    if not clusters:
        return None
    best = clusters[0]   # already sorted by score_sum desc
    if best["n_peaks"] < _CLUSTER_MIN_PEAKS:
        return None
    if (best["cov_frac_max"] is None
        or best["cov_frac_max"] < _CLUSTER_COV_FRAC_MIN):
        return None
    if (best["spatial_max"] is None
        or best["spatial_max"] < _CLUSTER_SPATIAL_MIN):
        return None
    if best["score_sum"] < rank1_score * _CLUSTER_SCORE_MARGIN:
        return None
    founder_period = int(round(best["median_period"]))
    return best, founder_period


def _harmonic_ladder_founder(peaks, strongest_period, strongest_id):
    """Pass 7 helper. Recover the basic monomer P0 of a divergent HOR satellite
    when the strongest period S sits atop a clean harmonic ladder (peaks at
    integer multiples m·P0). Returns (founder_peak, founder_period, multiplicity)
    or None. P0 may be below the strict id gate — that is the point; the ladder
    + tight-IQR + occupancy gates substitute for the missing identity. See
    docs/harmonic_ladder_founder_plan.md.

    Two regimes:
      k ≥ 3  → require ≥ _LADDER_MIN_RUNGS distinct integer multiples present
               (incl. P0 and S), tight IQR, decent occupancy.
      k == 2 → "exceptionally clean" only: stricter IQR + occupancy, AND the
               double must be a genuine HOR unit (id_med(2·P0) − id_med(P0) ≥
               _LADDER_X2_HOR_DELTA) — the conservation gap that tells the half
               is the diverged fundamental, not a sidelobe of a real 2·P0 monomer.
    """
    S = _num_or_none(strongest_period)
    if S is None or S <= 0:
        return None
    cand = []
    for p in peaks:
        P = _num_or_none(p.get("period"))
        idm = _num_or_none(p.get("identity_med"))
        if (P is None or P <= 0 or P >= S
                or idm is None or idm < _LADDER_ID_FLOOR):
            continue
        cand.append((P, p, idm))
    cand.sort(key=lambda t: t[0])               # smallest period (true fundamental) first
    for P0, p0, id0 in cand:
        k = int(round(S / P0))
        if k < _LADDER_KMIN or abs(S / P0 - k) > _RATIO_TOL:
            continue                            # S must be a clean rung of P0
        iqr0 = _num_or_none(p0.get("identity_iqr"))
        occ0 = _num_or_none(p0.get("scan_occupancy_frac"))
        if iqr0 is None or occ0 is None:
            continue
        if k >= _LADDER_MIN_RUNGS:              # --- k >= 3: rung-count path ---
            if iqr0 > _LADDER_IQR_MAX or occ0 < _LADDER_OCC_MIN:
                continue
            rungs = set()
            for m in range(1, k + 1):
                target = m * P0
                tol = max(_RATIO_TOL * target, 1.0)
                for q in peaks:
                    qp = _num_or_none(q.get("period"))
                    qid = _num_or_none(q.get("identity_med"))
                    if (qp is not None and abs(qp - target) <= tol
                            and qid is not None and qid >= _LADDER_ID_FLOOR):
                        rungs.add(m)
                        break
            if len(rungs) >= _LADDER_MIN_RUNGS and 1 in rungs and k in rungs:
                return p0, P0, k
        elif k == 2:                            # --- exceptionally-clean ×2 ---
            sid = _num_or_none(strongest_id)
            if (iqr0 <= _LADDER_X2_IQR_MAX and occ0 >= _LADDER_X2_OCC_MIN
                    and sid is not None and (sid - id0) >= _LADDER_X2_HOR_DELTA):
                return p0, P0, 2
    return None


# Dominant-score harmonic-ladder founder (Pass 7b) — the off-ladder counterpart
# to Pass 7. Pass 7 anchors the ladder on `strongest` (argmax id_med) and only
# helps when strongest = k·P0 is a CLEAN integer multiple of the fundamental.
# For a divergent satellite whose highest-identity peak is a long, OFF-ladder
# period (e.g. FabTR-53: strongest=3150, monomer=188, 3150/188=16.76), that
# decomposition either collapses (founder=strongest, mult=1) or lands on a
# non-monomer period (165, 168, 347, 1808). Pass 7b instead anchors on the
# rank-1-by-SCORE peak P0 — the true monomer of these satellites — and adopts it
# when P0 carries a self-standing harmonic ladder, regardless of where strongest
# sits. See docs/dominant_ladder_founder_plan.md.
_DOM_LADDER_MIN_RUNGS = 3      # distinct integer multiples m·P0 present (incl. m=1)
_DOM_LADDER_KH_TOL    = 0.05   # |kitehor_founder − P0| / max ≤ this ⇒ kitehor corroborates P0


def _commensurate(a, b, tol=_RATIO_TOL):
    """True when a and b are integer multiples of one another within an
    *absolute* |k − round(k)| ≤ tol (not k-scaled): a ≈ b (k=1), a ≈ k·b, or
    b ≈ k·a. Used as Pass 7b's safety guard — an array whose founder is already
    commensurate with its dominant-score peak is left untouched."""
    if not a or not b or a <= 0 or b <= 0:
        return False
    for x, y in ((a, b), (b, a)):
        k = x / y
        if round(k) >= 1 and abs(k - round(k)) <= tol:
            return True
    return False


def _dominant_ladder_founder(peaks, strongest_period, dominant_peak):
    """Pass 7b helper. Adopt the rank-1-by-score peak P0 as the founder when it
    sits atop a clean harmonic ladder (peaks at integer multiples m·P0),
    independent of whether `strongest` (= argmax id_med) is a clean multiple of
    P0. Returns (P0_peak, P0, multiplicity) or None. multiplicity = round(S/P0)
    — irregular when non-integer, which is the honest state for an off-ladder
    strongest. Gates mirror the k≥3 path of _harmonic_ladder_founder: P0 id_med
    above the ladder floor, tight IQR, decent occupancy, and
    ≥ _DOM_LADDER_MIN_RUNGS distinct rungs present (incl. m=1)."""
    S = _num_or_none(strongest_period)
    if dominant_peak is None or S is None or S <= 0:
        return None
    P0 = _num_or_none(dominant_peak.get("period"))
    if P0 is None or P0 <= 0 or P0 >= S:
        return None
    id0  = _num_or_none(dominant_peak.get("identity_med"))
    iqr0 = _num_or_none(dominant_peak.get("identity_iqr"))
    occ0 = _num_or_none(dominant_peak.get("scan_occupancy_frac"))
    if (id0 is None or id0 < _LADDER_ID_FLOOR
            or iqr0 is None or iqr0 > _LADDER_IQR_MAX
            or occ0 is None or occ0 < _LADDER_OCC_MIN):
        return None
    kmax = int(round(S / P0))
    # Assign each peak to its single nearest integer multiple of P0 (so one far
    # peak cannot fill several adjacent rungs through a width-growing tolerance).
    rungs = set()
    for q in peaks:
        qp  = _num_or_none(q.get("period"))
        qid = _num_or_none(q.get("identity_med"))
        if qp is None or qp <= 0 or qid is None or qid < _LADDER_ID_FLOOR:
            continue
        m = int(round(qp / P0))
        if 1 <= m <= kmax and abs(qp - m * P0) <= max(_RATIO_TOL * m * P0, 2.0):
            rungs.add(m)
    # Require a genuine LOW ladder: P0 (m=1), its dimer (m=2), and >= 3 distinct
    # rungs overall — not P0 plus a couple of scattered high multiples.
    if {1, 2} <= rungs and len(rungs) >= _DOM_LADDER_MIN_RUNGS:
        return dominant_peak, P0, max(1, kmax)
    return None


def build_monomer_size_csv(kite_tsv, ssr_tsv, rescored_peaks_tsv, out_csv,
                           tandem_validate_tsv=None, trc_repeat_type=None):
    """Build the per-array summary CSV (one row per array) from kitehor
    >=0.12.0 outputs: the rescored peaks TSV (`<prefix>.rescored.peaks.tsv`,
    24 cols = kite + 15 rescore) and the ssr-scan summary
    (`<prefix>.ssr.tsv`). Filename stays `monomer_size_top3_estimats.csv`
    so tc_per_tra_consensus.py + the R consensus prototype still find it.

    Per array we compute, in five passes:

    - **Pass 1 (per-TRA strict)**: the *strongest* period is the peak
      with the highest `identity_med` among peaks that pass the
      founder identity gate (≥ 0.7). On well-behaved arrays this is
      the same peak kitehor's rescore stamps as `founder_period`, but
      on the subset of arrays where kitehor picked by single-peak
      score rather than by identity (e.g. a 178 bp monomer that
      out-scores the 1786 bp 10× HOR period it tiles), this
      argmax-based pick recovers the higher-identity HOR period as
      strongest. Kitehor's original pick is retained as the
      `kitehor_founder_period` column for audit. The *founder* is
      reassigned by searching for divisor peaks P with
      `strongest = k · P`, `2 ≤ round(k) ≤ 30`, `|k − round(k)| ≤
      0.05`, and `identity_med(P) ≥ 0.7`. Candidates are then
      grouped by `kr = round(k)`; within each kr group the best
      representative is picked (highest id_med, then cleanest
      `|k − kr|`, then smallest P), and across kr groups the
      smallest period wins — so a deeper-decomposition (basic-
      monomer) candidate still outranks shallower ones, but at the
      same kr level a clean k=2.00 / id_med=0.995 divisor beats a
      fuzzy k=2.04 / id_med=0.972 one. If no candidate qualifies,
      founder = strongest, multiplicity = 1.
    - **Pass 2 (TRC consensus rescue)**: for each TRC, compute the
      consensus founder from Pass-1 successes (median of the largest
      ±10% cluster; gated by ≥ 3 arrays AND ≥ 25 % of the TRC's
      reassigned arrays). Arrays where Pass 1 left founder=strongest
      look for a rescored peak within max(±10%, ±20bp) of the
      consensus with id_med ≥ 0.5; if found, set founder = that peak's
      period, multiplicity = round(strongest/founder), and mark
      **irregular_multiplicity = true** (the raw fractional k is
      stored in `multiplicity_raw`).
    - **Pass 3 (solo / no-consensus fallback)**: when Pass 2 can't help
      (no TRC consensus, or no peak near it), relax the per-array
      multiplicity tolerance to ±0.20 (still id_med ≥ 0.7), flagging
      success as irregular_multiplicity = true.
    - **Pass 4 (SSR-founder override)**: see inline comment — forces
      founder = rank-1 kite peak on SSR-pure (≥ 95 %) arrays.
    - **Pass 5 (peak-cluster fallback rescue)**: when rescore returned
      NA *and* no SSR override applied, cluster the rescored peaks by
      period (single-link, ±5 % / ±100 bp) and accept the highest-
      `score_sum` cluster as founder when it has ≥ 2 peaks AND
      `max(coverage_frac) ≥ 0.20` AND `max(spatial_contrast) ≥ 0.5`
      AND `score_sum ≥ rank-1 peak score`. This recovers the real
      ~10 kb HOR signal on arrays where rescore's identity threshold
      filtered out every individual long-period peak but the
      fragmented evidence is collectively dominant. See
      `_cluster_rescue_founder()` for the gate logic.

    - **Pass 7 (harmonic-ladder founder)**: final rescue for arrays still
      at `founder == strongest` (mult == 1) — the basic monomer is below
      the strict id gate, so no deeper founder was found. When the
      strongest sits atop a clean harmonic ladder (peaks at integer
      multiples m·P0; ≥ 3 distinct rungs, or an *exceptionally-clean* ×2
      whose double is a genuine HOR unit), adopt the fundamental P0 with
      ×k. The ladder + tight-IQR + occupancy substitute for the missing
      identity. Recovers divergent-HOR satellites (e.g. 179 ×6, 111 ×3).
      Marks `founder_method = "ladder"`. See `_harmonic_ladder_founder()`
      and docs/harmonic_ladder_founder_plan.md.

    Other columns:
    - **kitehor_founder_period**: kitehor's original rescore pick for
      this case_id, preserved unchanged. When it differs from
      `strongest_period` the argmax(identity_med) override in Pass 1
      has fired (kitehor picked a lower-id_med peak by score).
    - **founder_fallback**: true when rescore returned NA *and* Pass 5
      did not rescue (i.e. truly no signal). Cleared on cluster rescue.
    - **Pass 1b (kitehor go-deeper)**: before Pass 2, consult kitehor's
      own per-array basic monomer (`hor_basic_period`). Adopt it as
      founder only when it is *deeper* than the strict pick AND its
      rescored peak is a near-full-coverage array-wide tandem
      (coverage_frac ≥ 0.75, scan_occupancy_frac ≥ 0.90). Recovers deep
      HORs the strict ±0.05 gate rejects at large k; coverage-gated so
      kitehor's spurious near-2× over-splits are not adopted. Marks
      `founder_method = "kh_deeper"`.
    - **founder_method**: enum recording which pass produced the final
      founder — `strict` / `none` / `kh_deeper` / `ladder` / `pass2` /
      `pass3` / `ssr` / `cluster` / `fallback`.
    - **cluster_rescue**: true when Pass 5 promoted a cluster median
      to founder. Diagnostic counterpart of `founder_method = cluster`.
    - **delta_id_pp**: identity_med(strongest) − identity_med(founder)
      in percentage points; NA when founder = strongest.
    - **top-5 monomer estimates by score**.
    - **top-2 subrepeat candidates**: when `tandem_validate_tsv` is
      given (kitehor >= 0.13.0), drawn from `tandem-validate`'s
      per-array `localized_subrepeat` candidates, gated against
      TideCluster's founder — adopted only when the candidate period
      does NOT coincide with the founder (not the HOR subunit) and shows
      genuine partial occupancy (density ≥ 0.10, presence ≥ 10 %). The
      candidate's occupancy lands in `subrepeat_*_occ` and its tier is
      `HIGH`. When no tandem-validate file is supplied, falls back to the
      legacy per-peak `_classify_subrepeat_tier` rescore-column heuristic.
    - **subrepeat_founder_divisor_flag / _k**: edge diagnostic — set when
      a tandem-validate candidate is a clean integer divisor of the
      founder at high occupancy (density ≥ 0.50). Surfaces arrays where
      the founder itself may be miscalled (a smaller unit tiles it), for
      review rather than an automatic founder change.
    - **alt_cluster_1..2_***: top-2 period clusters that are *not* the
      founder's, surfaced so the report can show secondary signals
      (e.g. a spurious 310-bp peak that was correctly demoted).
    - **SSR fields** lifted from `<prefix>.ssr.tsv`.

    `kite_tsv` is currently unused but kept on the signature to keep
    the call site self-documenting (it pairs with the rescored peaks
    file). Empty `hor_status` / `hor_confidence` cells are emitted so
    the per-TRA consensus R script's HOR_* alias shim doesn't choke
    (the spectral HOR classifier they fed is retired — see
    `hor_order_confidence` below and docs/hor_order_confidence_design.md).

    `hor_order_confidence` (none | strict | supported | weak) is the
    founder-arithmetic HOR tier that replaces the spectral classifier:
    it distinguishes a recovered founder from a confidently-ordered
    HOR. "HOR" for the report = strict ∪ supported."""
    import csv
    import collections

    _num = _num_or_none

    def _fmt_bp(v):
        if v is None: return ""
        try:
            iv = int(round(float(v)))
            return str(iv) if abs(float(v) - iv) < 1e-6 else f"{float(v):g}"
        except (TypeError, ValueError):
            return ""

    def _fmt_sc(v):
        if v is None: return ""
        return f"{float(v):.10f}".rstrip("0").rstrip(".")

    def _fmt_id(v):
        if v is None: return ""
        return f"{float(v):.4f}".rstrip("0").rstrip(".")

    def _fmt_pp(v):
        if v is None: return ""
        return f"{float(v):.2f}"

    # ---- Read all rescored peaks, group by record_id, sort by rank ----
    peaks_by_rec = collections.defaultdict(list)
    with open(rescored_peaks_tsv, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            peaks_by_rec[r["case_id"]].append(r)
    for rs in peaks_by_rec.values():
        rs.sort(key=lambda r: int(r.get("rank") or 0))

    # ---- Read SSR scan summary keyed by record_id ----
    ssr_by_rec = {}
    if ssr_tsv and os.path.exists(ssr_tsv):
        with open(ssr_tsv, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                ssr_by_rec[r["record_id"]] = r

    # ---- Read tandem-validate subrepeat detector keyed by record_id ----
    # kitehor >= 0.13.0's unified nested-TR detector (spec v5). Each row
    # carries the dominant candidate plus a `candidates` cell listing all
    # of them, `;`-separated, each entry `kind/period:d=..:sc=..:pc=..:label`
    # (e.g. `other/305:d=0.084:sc=0.500:pc=0.111:localized`).
    tv_by_rec = {}
    if tandem_validate_tsv and os.path.exists(tandem_validate_tsv):
        with open(tandem_validate_tsv, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                tv_by_rec[r["record_id"]] = r

    # ---- Pass 1: per-array strict reassignment ----
    # Stash one dict per record with the strict result + the raw peaks
    # list (used again in Pass 2/3 and for the subrepeat tier classifier).
    pass1 = []
    for record_id, peaks in peaks_by_rec.items():
        if ":" not in record_id:
            continue
        trc, rest = record_id.split(":", 1)
        parts = rest.rsplit("_", 2)
        if len(parts) != 3:
            continue
        seqid, start, end = parts
        rank1 = peaks[0]
        # kitehor's rescore writes the same `founder_period` on every row
        # of a case_id (its own pick of the best peak). We keep it as a
        # diagnostic only — see `kitehor_founder_period` in the CSV —
        # because on some arrays kitehor's pick is the rank-1-by-score
        # peak rather than the peak with the highest `identity_med`
        # (e.g. S. pimpinellifolium TRC_8 chr12:1314680 picks the 178bp
        # monomer at id_med=0.736 over the 1786bp HOR period at
        # id_med=0.818, which then suppresses the 10×178 HOR call).
        # An NA value still signals "no rescore signal at all", which
        # we keep as the Pass 5 cluster-rescue trigger.
        kitehor_fp = _num(rank1.get("founder_period"))
        fallback = False
        if kitehor_fp is None:
            fallback        = True
            strongest_peak  = rank1
            strongest_period = _num(rank1.get("period"))
            founder_peak    = rank1
            founder_period  = strongest_period
            multiplicity    = 1
            multiplicity_raw = 1.0 if strongest_period else None
        else:
            # Strongest = peak with the highest identity_med among peaks
            # that pass the founder identity gate (>= _FOUNDER_ID_MIN).
            # On well-behaved arrays kitehor's founder_period IS already
            # the argmax(id_med), so this is a no-op; on the subset where
            # kitehor picked by score and missed a higher-id_med long-
            # period peak, this gives Pass-1's divisor search the right
            # strongest target and recovers the HOR call.
            strongest_peak   = None
            strongest_period = None
            best_idm         = -1.0
            for p in peaks:
                idm    = _num(p.get("identity_med"))
                period = _num(p.get("period"))
                if idm is None or period is None or period <= 0:
                    continue
                if idm < _FOUNDER_ID_MIN:
                    continue
                if idm > best_idm:
                    best_idm         = idm
                    strongest_peak   = p
                    strongest_period = period
            if strongest_period is None:
                # No peak cleared the id gate even though kitehor reported
                # a founder. Defensive: keep kitehor's pick rather than
                # falling through (changes are then exactly the same as
                # 1.12.0 for these rows).
                strongest_period = kitehor_fp
                strongest_peak   = next(
                    (p for p in peaks if _num(p.get("period")) == kitehor_fp),
                    rank1)
            # Lever 2 — cluster-mean basic-monomer search. Cluster the
            # rescored peaks by period (single-link, mixed ±5% / ±100 bp
            # window, the same _cluster_peaks_by_period helper Pass 5
            # uses) and test each cluster's score-weighted mean period
            # as a divisor of `strongest_period` under the strict
            # ±_RATIO_TOL gate. This recovers cases where the basic
            # monomer is supported by multiple nearby peaks but no
            # single peak forms a clean integer divisor — e.g. drapa
            # TRC_26 chr9:1816989 strongest=1880 with peaks
            # 154/158/161/165 each at id_med ≈ 0.95 but |k − kr| up to
            # 0.21 individually; the score-weighted mean ≈ 156.8 gives
            # k = 11.99 (clean). A singleton cluster's mean equals the
            # peak's own period, so the prior single-peak behaviour is
            # preserved exactly when no near-neighbour exists.
            #
            # Each cluster yields one candidate tuple:
            #   (P_mean, rep_id_med, kr, k, rep_peak, rep_period)
            # where `rep_peak` = the cluster member with the highest
            # id_med (drives the founder_id_med and delta_id_pp cells)
            # and `rep_period` = that member's actual reported period
            # (used as the deterministic tiebreaker in the per-kr sort).
            clusters = _cluster_peaks_by_period(peaks)
            candidates = []
            for cluster in clusters:
                # Only members passing the founder id gate contribute to
                # the cluster's representation — otherwise a noise peak
                # with low id_med that happens to fall inside the cluster
                # window (max(±5%, ±100 bp)) can both drag the mean and
                # pretend to be a representative.
                qualified = [m for m in cluster["members"]
                             if (_num(m.get("identity_med")) or 0.0)
                                 >= _FOUNDER_ID_MIN]
                if not qualified:
                    continue
                # Strategy: prefer the cluster's score-weighted mean if
                # it forms a clean integer divisor of strongest (the
                # Lever 2 use case: multiple nearby peaks averaging to
                # the true basic). If the mean fails the strict gate
                # but an individual cluster member passes (the TRC_45
                # chr12:19629243 case: cluster {334,337,340,344} mean
                # fails k-gap 0.32 but 344 alone is clean), fall back
                # to per-member candidate generation within the cluster
                # — this preserves the pre-Lever-2 single-peak behaviour
                # for arrays where one peak is structurally cleaner than
                # its noisy neighbours.
                num = 0.0
                den = 0.0
                for m in qualified:
                    mp = _num(m.get("period"))
                    sc = _num(m.get("score"))
                    if mp is None or sc is None or sc <= 0:
                        continue
                    num += mp * sc
                    den += sc
                mean_added = False
                if den > 0:
                    P_mean = num / den
                    if P_mean < strongest_period:
                        k_mean = strongest_period / P_mean
                        kr_mean = round(k_mean)
                        if (2 <= kr_mean <= _KMAX
                            and abs(k_mean - kr_mean) <= _RATIO_TOL):
                            rep_peak = max(
                                qualified,
                                key=lambda m: _num(m.get("identity_med")) or 0.0)
                            rep_idm = _num(rep_peak.get("identity_med"))
                            rep_period = _num(rep_peak.get("period")) or P_mean
                            candidates.append(
                                (P_mean, rep_idm, kr_mean, k_mean,
                                 rep_peak, rep_period))
                            mean_added = True
                if mean_added:
                    continue
                # Fallback: per-member single-peak candidates inside
                # this cluster (same gate the pre-Lever-2 code used).
                for m in qualified:
                    P = _num(m.get("period"))
                    if P is None or P >= strongest_period:
                        continue
                    idm = _num(m.get("identity_med"))
                    k = strongest_period / P
                    kr = round(k)
                    if not (2 <= kr <= _KMAX):
                        continue
                    if abs(k - kr) > _RATIO_TOL:
                        continue
                    candidates.append((P, idm, kr, k, m, P))
            if candidates:
                # Two-level pick (unchanged from 9254c93, applied to
                # cluster-mean candidates). Across kr groups: smallest
                # mean wins (basic-monomer rule, so a kr=10 cluster at
                # mean 178 still beats a kr=2 cluster at mean 892).
                # Within a kr group: highest id_med, then cleanest
                # |k - kr|, then smallest mean.
                by_kr = collections.defaultdict(list)
                for c in candidates:
                    by_kr[c[2]].append(c)
                per_kr_winners = []
                for group in by_kr.values():
                    group.sort(key=lambda c: (-c[1], abs(c[3] - c[2]), c[0]))
                    per_kr_winners.append(group[0])
                per_kr_winners.sort(key=lambda c: c[0])
                P_winner, _id, kr, k_raw, fp_peak, _rep_period = per_kr_winners[0]
                # P_winner is either the cluster mean (Lever 2 path) or a
                # single member period (per-member fallback) — round in
                # both cases to give an integer founder_period.
                founder_period   = int(round(P_winner))
                founder_peak     = fp_peak
                multiplicity     = kr
                multiplicity_raw = k_raw
            else:
                founder_period   = strongest_period
                founder_peak     = strongest_peak
                multiplicity     = 1
                multiplicity_raw = 1.0
        # Harmonic-basis deepening (case 2). The dominant kite-score periodicity
        # (rank-1 by score) is the basic monomer of a harmonic series; when the
        # founder chosen above is a clean integer multiple of it, deepen the
        # founder to it. Recovers TRC_4-type cases where `strongest =
        # argmax(id_med)` landed on a long noisy period (7617) and the k<=_KMAX
        # divisor search settled on an intermediate multiple (318 = 6*53) instead
        # of the real basic (53) — even though 53 is the top kite peak and passes
        # the id gate. Guards: rank-1 must pass the id gate, be SHORTER than the
        # founder, and the founder/rank1 ratio must be a clean integer in
        # [2,_KMAX] (so a real 178 bp satellite, whose rank-1 peak IS 178, never
        # collapses to an internal sub-period).
        irregular_h = False
        if not fallback:
            p0_period = _num(rank1.get("period"))
            p0_idm    = _num(rank1.get("identity_med"))
            if (p0_period is not None and p0_period > 0
                    and founder_period is not None
                    and p0_period < founder_period
                    and p0_idm is not None and p0_idm >= _FOUNDER_ID_MIN):
                k_basic  = founder_period / p0_period
                kr_basic = round(k_basic)
                if (2 <= kr_basic <= _KMAX
                        and abs(k_basic - kr_basic) <= _RATIO_TOL):
                    founder_peak   = rank1
                    founder_period = int(round(p0_period))
                    if strongest_period:
                        multiplicity_raw = strongest_period / founder_period
                        multiplicity     = int(round(multiplicity_raw))
                        # flag irregular when strongest isn't a clean multiple
                        # of the basic (the basic still wins — the HOR period
                        # just isn't an exact integer multiple).
                        irregular_h = (abs(multiplicity_raw - multiplicity)
                                       > _RATIO_TOL)
        pass1.append({
            "record_id":         record_id,
            "trc":               trc,
            "seqid":             seqid,
            "start":             start,
            "end":               end,
            "peaks":             peaks,
            "rank1":             rank1,
            "strongest_peak":    strongest_peak,
            "strongest_period":  strongest_period,
            # kitehor's original rescore pick (may differ from
            # strongest_period when TideCluster's argmax(id_med)
            # override fires).
            "kitehor_founder_period": kitehor_fp,
            "founder_peak":      founder_peak,
            "founder_period":    founder_period,
            "multiplicity":      multiplicity,
            "multiplicity_raw":  multiplicity_raw,
            "fallback":          fallback,
            "irregular":         irregular_h,
            # rescue_method records which Pass produced the final founder.
            # Default reflects Pass 1's outcome: "strict" when multiplicity
            # >= 2 (real HOR call), "none" when mult=1 (rescore returned a
            # founder but no smaller-period subrepeat passed), "fallback"
            # when rescore returned NA. Passes 2/3/4/5 overwrite as needed.
            "rescue_method":     ("strict" if (not fallback and multiplicity > 1)
                                  else ("fallback" if fallback else "none")),
        })

    # ---- Pass 2 prep: TRC consensus founder from Pass-1 successes ----
    records_by_trc = collections.defaultdict(list)
    for e in pass1:
        records_by_trc[e["trc"]].append(e)
    consensus_by_trc = {}
    for trc, recs in records_by_trc.items():
        strict_founders = [e["founder_period"] for e in recs
                           if e["multiplicity"] > 1 and not e["fallback"]
                           and e["founder_period"] is not None]
        consensus_period, _n = _trc_consensus_founder(strict_founders)
        consensus_by_trc[trc] = consensus_period

    # ---- Pass 2 + Pass 3: rescue arrays where Pass 1 left mult=1
    #                       (or where Pass-1's founder sits at an
    #                       integer multiple of the TRC consensus —
    #                       the Lever 3 expansion) ----
    # Original gate: only Pass-1 mult=1 (or fallback) rows. Lever 3
    # adds a second escape clause covering rows where Pass 1 succeeded
    # with mult>1 but the founder ended up at 2..5× the TRC consensus
    # basic — the signature of "no peak near the basic monomer formed
    # a clean integer divisor of strongest individually AND no nearby
    # cluster existed either, so Pass 1 settled on a multiple". We
    # re-evaluate those rows against the TRC consensus oracle.
    # Examples on drapa: TRC_138 chr4:26085301 (founder=262=2*131,
    # consensus=131); the rescue from TRC consensus yields
    # founder=131, mult=round(strongest/131)=26, irregular=true.
    for e in pass1:
        if e["fallback"]:
            continue
        if e["strongest_period"] is None:
            continue
        if e["multiplicity"] > 1:
            # Lever 3 gate — only rescue mult>1 rows that look like an
            # integer-multiple-of-the-basic pick. Requires: a TRC
            # consensus exists, founder/consensus is approximately
            # integer in [2..5], founder is significantly off the
            # consensus (else there's nothing to fix), AND the
            # consensus is plausibly nested inside founder (founder
            # divides into strongest in proportion to consensus).
            tc_cons = consensus_by_trc.get(e["trc"])
            if tc_cons is None or tc_cons <= 0:
                continue
            fp = e["founder_period"]
            if fp is None or fp <= 0:
                continue
            ratio = fp / tc_cons
            rr = round(ratio)
            if not (2 <= rr <= 5):
                continue
            if abs(ratio - rr) > 0.10:
                continue
            if abs(fp - tc_cons) / tc_cons <= 0.20:
                continue
        # Pass 2 — TRC consensus rescue
        trc_consensus = consensus_by_trc.get(e["trc"])
        rescued = _rescue_founder_from_trc(
            e["peaks"], e["strongest_period"], trc_consensus)
        if rescued:
            founder_peak, founder_period = rescued
            raw = e["strongest_period"] / founder_period
            kr  = int(round(raw))
            # Only count as a real rescue when the rounded multiplicity
            # is at least 2 — strongest ≈ founder rescues add no HOR
            # signal and would inflate the irregular count.
            if kr >= 2:
                e["founder_peak"]     = founder_peak
                e["founder_period"]   = founder_period
                e["multiplicity_raw"] = raw
                e["multiplicity"]     = kr
                e["irregular"]        = True
                e["rescue_method"]    = "pass2"
                continue
        # Pass 3 — solo / no-consensus relaxed individual fallback
        relaxed = _relaxed_individual_founder(e["peaks"], e["strongest_period"])
        if relaxed:
            founder_peak, P, kr, k_raw = relaxed
            e["founder_peak"]     = founder_peak
            e["founder_period"]   = P
            e["multiplicity"]     = kr
            e["multiplicity_raw"] = k_raw
            e["irregular"]        = True
            e["rescue_method"]    = "pass3"

    # ---- Lever 4: prevalent-founder anchor (consensus-anchored deepening) ----
    # Recover long-HOR arrays whose founder landed on a non-consensus intermediate
    # because the strict k<=_KMAX / |k-round|<=_RATIO_TOL divisor gate rejected a
    # high or non-clean multiplicity (e.g. TRC_4 chr8:13133296: founder 150,
    # strongest 2106, while 53 — the family's prevalent founder — is the rank-1
    # kite peak; 2106/53 = 39.7 fails the strict gate). When the array's DOMINANT
    # (rank-1-by-score) peak IS the TRC consensus and passes the id gate, deepen
    # the founder to it, tolerating high/non-clean k via the k-scaled integer test
    # (the family vouches for the basic; above k≈0.5/frac the test is vacuous and
    # consensus alone carries it). Only deepens (founder > consensus), so Group-1
    # near-neighbour founders below the consensus are left alone and satellites
    # (already at consensus) are untouched. Marks `~` irregular when not clean.
    _rt = trc_repeat_type or {}
    for e in pass1:
        if e["fallback"] or e.get("ssr_override") or e["trc"] in _rt:
            continue
        tc_cons = consensus_by_trc.get(e["trc"])
        sp = e["strongest_period"]
        fp = e["founder_period"]
        if (tc_cons is None or tc_cons <= 0 or sp is None or fp is None
                or fp <= tc_cons * (1.0 + _CONS_ANCHOR_PCT)):
            continue
        p0     = _num(e["rank1"].get("period"))
        p0_idm = _num(e["rank1"].get("identity_med"))
        if (p0 is None or p0 <= 0 or p0_idm is None
                or p0_idm < _FOUNDER_ID_MIN
                or abs(p0 - tc_cons) > max(_CONS_ANCHOR_BP,
                                           tc_cons * _CONS_ANCHOR_PCT)):
            continue
        k  = sp / p0
        kr = int(round(k))
        if kr < 2:
            continue
        tol = _RATIO_TOL_FRAC * k
        if abs(k - kr) > tol and tol < 0.5:   # integer test still meaningful + fails
            continue
        e["founder_peak"]     = e["rank1"]
        e["founder_period"]   = int(round(p0))
        e["multiplicity"]     = kr
        e["multiplicity_raw"] = k
        e["irregular"]        = abs(k - kr) > _RATIO_TOL
        e["rescue_method"]    = "pass2"

    # ---- Pass 4: SSR founder from the clustering repeat_type classification ----
    # SSR families are identified up front at clustering: a TRA *consensus* that
    # is > 90 % simple-sequence is pulled out of clustering, grouped by motif, and
    # tagged `repeat_type=SSR` + `ssr=MOTIF` in the clustering GFF3 (see
    # TideCluster.py). `trc_repeat_type` maps such TRCs -> the fundamental motif
    # length (ATC -> 3, AAACCCT -> 7). For EVERY array of an SSR-family TRC the
    # founder IS that motif length, regardless of the per-array kite peaks /
    # coverage / rescore — the consensus call already settled the monomer, and the
    # per-array kite signal is noisy (one ATC family spans 16-100 % coverage on
    # the full sequence because the consensus hides per-array noise). This both
    # (a) propagates the motif to low-coverage members the old per-array 95 %
    # coverage override missed (-> founder 6/9), and (b) stops a TR satellite's
    # lone SSR-rich array from being flipped to the motif: TR families are absent
    # from `trc_repeat_type`, so e.g. TRC_2's 96 %-ATC array keeps its real 9490
    # founder (the SSR content stays in the row as annotation). No HOR
    # decomposition is attempted on an SSR motif. kitehor's per-array ssr-scan is
    # now purely secondary (annotates SSR-rich content inside any TRC).
    trc_repeat_type = trc_repeat_type or {}
    for e in pass1:
        e["ssr_override"] = False
        motif_len = trc_repeat_type.get(e["trc"])
        if not motif_len:
            continue
        peak_at = next((p for p in e["peaks"]
                        if _num(p.get("period")) == motif_len), e["rank1"])
        e["founder_period"]   = motif_len
        e["strongest_period"] = motif_len
        e["multiplicity"]     = 1
        e["multiplicity_raw"] = 1.0
        e["irregular"]        = False
        e["fallback"]         = False
        e["founder_peak"]     = peak_at
        e["strongest_peak"]   = peak_at
        e["ssr_override"]     = True
        e["rescue_method"]    = "ssr"

    # ---- Pass 5: peak-cluster fallback rescue ----
    # Fires only when rescore returned NA (fallback=True) AND no SSR
    # override applied. Pass-1 successes, Pass-2/3 rescues, and SSR
    # override rows are not entered. See _cluster_rescue_founder() for
    # the cluster + gating logic. On success: clear `fallback`, point
    # founder/strongest at the cluster's highest-score peak, mark
    # rescue_method="cluster". `multiplicity=1` because the rescue
    # makes no claim about a smaller-period subrepeat (and Pass 2/3
    # already had the opportunity to find one if rescore's founder
    # had been usable in the first place).
    for e in pass1:
        e["cluster_rescue"] = False
        if not e["fallback"] or e.get("ssr_override"):
            continue
        rescued = _cluster_rescue_founder(e["peaks"])
        if rescued is None:
            continue
        cluster, founder_period = rescued
        e["founder_period"]     = founder_period
        e["strongest_period"]   = founder_period
        e["founder_peak"]       = cluster["highest_score_peak"]
        e["strongest_peak"]     = cluster["highest_score_peak"]
        e["multiplicity"]       = 1
        e["multiplicity_raw"]   = 1.0
        e["fallback"]           = False     # cluster rescue clears the badge
        e["cluster_rescue"]     = True
        e["rescue_method"]      = "cluster"

    # ---- Pass 6: kitehor "go-deeper" founder adoption (Bucket 1) ----
    # Consult kitehor's own per-array basic monomer (hor_basic_period,
    # constant per case_id, from rescore >= 0.13.0). Adopt it as founder
    # ONLY when it is deeper than TideCluster's FINAL founder (after all
    # prior passes) AND the rescored peak at that period is a clean,
    # near-full-coverage array-wide tandem. This recovers real deep HORs
    # that TC's strict ±0.05 ratio gate rejects at large k and that the
    # TRC-consensus passes couldn't reach (the 4 confirmed drapa cases:
    # TRC_173 496->165 k=27.07; TRC_276 205->103 k=21.9; TRC_178 64->48;
    # TRC_203 3148->1522). Running LAST and gating "deeper than the final
    # founder" is deliberate: it never walks a correct Pass-2 deep founder
    # back to kitehor's shallower intermediate divisor (TRC_26 final 154
    # rejects kitehor's 462=3*154), and the hard coverage/occupancy gate
    # keeps kitehor's spurious near-2x over-splits out (TRC_1 1692->890 at
    # cov≈0.0). Deeper-only — never overrides with a shallower founder.
    for e in pass1:
        if e["fallback"] or e.get("ssr_override"):
            continue
        if e["strongest_period"] is None:
            continue
        kh_basic = _num(e["rank1"].get("hor_basic_period"))
        if kh_basic is None or kh_basic <= 0:
            continue
        cur_fp = e["founder_period"]
        # must be a meaningfully smaller monomer than TC's final founder
        # (not a ~few-% lateral re-estimate of the same period)
        if cur_fp is None or kh_basic > cur_fp * _KH_DEEPER_MAX_RATIO:
            continue
        # kitehor's basic must form a clean integer multiple of strongest
        k  = e["strongest_period"] / kh_basic
        kr = round(k)
        if not (2 <= kr <= _KMAX) or abs(k - kr) > _KH_DEEPER_RATIO_TOL:
            continue
        # locate the rescored peak at kh_basic; require a (near-)full
        # coverage array-wide tandem there (Bucket-1 vs Bucket-3 split).
        win  = max(12.0, 0.05 * kh_basic)
        near = [p for p in e["peaks"]
                if (_num(p.get("period")) is not None
                    and abs(_num(p.get("period")) - kh_basic) <= win)]
        if not near:
            continue
        bp  = max(near, key=lambda p: (_num(p.get("coverage_frac")) or 0.0))
        cov = _num(bp.get("coverage_frac"))
        occ = _num(bp.get("scan_occupancy_frac"))
        if cov is None or occ is None:
            continue
        if cov < _KH_DEEPER_COV_MIN or occ < _KH_DEEPER_OCC_MIN:
            continue
        e["founder_peak"]     = bp
        e["founder_period"]   = int(round(kh_basic))
        e["multiplicity"]     = kr
        e["multiplicity_raw"] = k
        # flag irregular only when the integer ratio needed the relaxed gate
        e["irregular"]        = abs(k - kr) > _RATIO_TOL
        e["rescue_method"]    = "kh_deeper"

    # ---- Pass 7: harmonic-ladder founder (divergent HOR basic monomer) ----
    # Final rescue for arrays the earlier passes left at founder == strongest
    # (multiplicity == 1): the basic monomer was below the strict id gate, so no
    # deeper founder was found. When the strongest sits atop a clean harmonic
    # ladder (peaks at integer multiples m·P0), adopt the fundamental P0 as the
    # founder with ×k. See docs/harmonic_ladder_founder_plan.md.
    _rt7 = trc_repeat_type or {}
    for e in pass1:
        if e["fallback"] or e.get("ssr_override") or e["trc"] in _rt7:
            continue
        if e["multiplicity"] != 1 or e["strongest_period"] is None:
            continue
        sp = e.get("strongest_peak") or {}
        res = _harmonic_ladder_founder(
            e["peaks"], e["strongest_period"], _num(sp.get("identity_med")))
        if res is None:
            continue
        p0, P0, k = res
        e["founder_peak"]     = p0
        e["founder_period"]   = P0
        e["multiplicity"]     = k
        e["multiplicity_raw"] = e["strongest_period"] / P0
        e["irregular"]        = abs(e["multiplicity_raw"] - k) > _RATIO_TOL
        e["rescue_method"]    = "ladder"

    # ---- Pass 7b: dominant-score harmonic-ladder founder ----
    # Complementary to Pass 7. Pass 7 flips founder==strongest only when
    # strongest = k·P0 is a CLEAN multiple of a ladder fundamental. Pass 7b
    # handles the off-ladder case: the argmax(id_med) strongest is NOT a clean
    # multiple of the dominant monomer (e.g. FabTR-53 strongest=3150,
    # monomer=188, 3150/188=16.76), so the divisor search either collapsed
    # (mult=1) or landed on a non-monomer period (165, 168, 347, 1808). When the
    # rank-1-by-SCORE peak P0 carries a strong, self-standing harmonic ladder,
    # the current founder is incommensurate with P0 (neither founder=k·P0 nor
    # P0=k·founder within _RATIO_TOL), AND kitehor's own k-mer-autocorrelation
    # founder_period corroborates P0 (same monomer), adopt P0. Two guards keep it
    # tight: the incommensurate test leaves every array already consistent with
    # its dominant monomer untouched (incl. the clean k·P0 cases Pass 7 just
    # fixed — 179×6, 111×3 — and ordinary HOR calls); the kitehor-corroboration
    # test requires an independent method to agree the dominant-score peak is the
    # founder, so a lone high-score peak that kitehor itself does not call a
    # founder (e.g. an array with a genuine ~12 kb HOR unit) is not overridden.
    # See docs/dominant_ladder_founder_plan.md.
    for e in pass1:
        if e["fallback"] or e.get("ssr_override") or e["trc"] in _rt7:
            continue
        if e["strongest_period"] is None:
            continue
        P0v = _num(e["rank1"].get("period"))
        if P0v is None or _commensurate(e["founder_period"], P0v):
            continue
        # kitehor's own founder pick must agree that P0 is the monomer.
        kh = _num(e["kitehor_founder_period"])
        if kh is None or abs(kh - P0v) > _DOM_LADDER_KH_TOL * max(kh, P0v):
            continue
        res = _dominant_ladder_founder(
            e["peaks"], e["strongest_period"], e["rank1"])
        if res is None:
            continue
        p0, P0, k = res
        e["founder_peak"]     = p0
        e["founder_period"]   = P0
        e["multiplicity"]     = k
        e["multiplicity_raw"] = e["strongest_period"] / P0
        e["irregular"]        = abs(e["multiplicity_raw"] - k) > _RATIO_TOL
        e["rescue_method"]    = "dominant_ladder"

    # ---- Build output rows ----
    rows = []
    for e in pass1:
        founder_id   = _num(e["founder_peak"].get("identity_med"))
        strongest_id = _num(e["strongest_peak"].get("identity_med"))
        delta_id_pp = (
            (strongest_id - founder_id) * 100.0
            if (founder_id is not None and strongest_id is not None
                and e["multiplicity"] > 1)
            else None)

        # Top-5 peaks by score (rescored peaks are emitted in rank order).
        peaks = e["peaks"]
        top5 = peaks[:5] + [None] * max(0, 5 - len(peaks))
        ms_scores = [
            (None, None) if p is None
            else (_num(p.get("period")), _num(p.get("score")))
            for p in top5]

        # Subrepeat candidates. When kitehor >= 0.13.0's tandem-validate
        # output is available, draw the top-2 from its per-array
        # localized_subrepeat candidates, gated against TideCluster's
        # founder (suppress candidates that coincide with / are the HOR
        # subunit; flag clean high-occupancy founder divisors). Otherwise
        # fall back to the legacy per-peak rescore-column tiering. Tuples
        # are stored as (occupancy, period, tier) for the row writer.
        sr_edge_flag = False
        sr_edge_k    = None
        if tv_by_rec:
            tv_sr, sr_edge_flag, sr_edge_k = _tv_subrepeats(
                tv_by_rec.get(e["record_id"]), e["founder_period"])
            sr_keep = [(dens, per, tier) for (per, dens, tier) in tv_sr]
        else:
            sr_keep = []
            for p in peaks:
                t = _classify_subrepeat_tier(p, e["founder_period"])
                if t in ("HIGH", "LIKELY"):
                    sr_keep.append((
                        _num(p.get("scan_occupancy_frac")) or 0.0,
                        _num(p.get("period")),
                        t))
            sr_keep.sort(key=lambda c: -c[0])
        sr1 = sr_keep[0] if len(sr_keep) > 0 else None
        sr2 = sr_keep[1] if len(sr_keep) > 1 else None

        # Alternative periodicities: cluster every row's rescored peaks
        # and report the top-N clusters that are *not* the founder's.
        # This surfaces secondary monomer signals (e.g. the spurious
        # short-period peak that fallback rescue rejected, or genuine
        # sub-monomer evidence on Pass-1 successes). The same clustering
        # is used by Pass 5's rescue but is recomputed here cheaply.
        fp_pick = e["founder_period"]
        clusters_all = _cluster_peaks_by_period(peaks)
        if fp_pick is None:
            alt_clusters = clusters_all
        else:
            fp_window = max(fp_pick * _CLUSTER_WINDOW_PCT,
                            _CLUSTER_WINDOW_BP_MIN)
            def _is_founder_cluster(c):
                for m in c["members"]:
                    mp = _num(m.get("period"))
                    if mp is not None and abs(mp - fp_pick) <= fp_window:
                        return True
                return False
            alt_clusters = [c for c in clusters_all if not _is_founder_cluster(c)]
        alt1 = alt_clusters[0] if len(alt_clusters) > 0 else None
        alt2 = alt_clusters[1] if len(alt_clusters) > 1 else None

        # Short-founder review aids (pure diagnostics — they never change the
        # founder). For any short founder (<= _SHORT_FOUNDER_MAX) surface the
        # dominant *longer*-period alternative so a reviewer can sanity-check the
        # short call at a glance, and flag the subset whose founder peak has weak
        # kite support (score < _FOUNDER_SCORE_MIN) — kitehor 0.13.2 already
        # suppresses the by-chance short phantoms, but a surviving weak short
        # founder still warrants a manual look. Strongly supported short founders
        # (e.g. a clean 15 bp monomer at score 0.87) surface the alternative but
        # are not flagged.
        fp_val      = e["founder_period"]
        weak_short  = False
        alt_longer  = None
        if fp_val is not None and fp_val <= _SHORT_FOUNDER_MAX:
            # Surface the most credible longer alternative: the highest
            # kite-score peak above the founder that itself shows real
            # periodicity support (_peak_credible). Blank ⇒ no credible longer
            # monomer (the short call stands on its own), which is itself
            # informative; requiring support avoids surfacing a noise peak.
            longer = [p for p in peaks
                      if (_num(p.get("period")) or 0) > fp_val
                      and _num(p.get("score")) is not None
                      and _peak_credible(p)]
            if longer:
                alt_longer = _num(
                    max(longer, key=lambda p: _num(p.get("score")) or 0.0)
                    .get("period"))
            fsc = _num(e["founder_peak"].get("score"))
            if fsc is not None and fsc < _FOUNDER_SCORE_MIN:
                weak_short = True

        ssr = ssr_by_rec.get(e["record_id"], {})
        rows.append({
            "TRC_ID":            e["trc"],
            "seqid":             e["seqid"],
            "start":             e["start"],
            "end":               e["end"],
            "array_length":      e["rank1"].get("array_length", ""),
            "monomer_size":      _fmt_bp(ms_scores[0][0]),
            "score":             _fmt_sc(ms_scores[0][1]),
            "monomer_size_2":    _fmt_bp(ms_scores[1][0]),
            "score_2":           _fmt_sc(ms_scores[1][1]),
            "monomer_size_3":    _fmt_bp(ms_scores[2][0]),
            "score_3":           _fmt_sc(ms_scores[2][1]),
            "monomer_size_4":    _fmt_bp(ms_scores[3][0]),
            "score_4":           _fmt_sc(ms_scores[3][1]),
            "monomer_size_5":    _fmt_bp(ms_scores[4][0]),
            "score_5":           _fmt_sc(ms_scores[4][1]),
            # hor_status / hor_confidence: the old spectral HOR classifier is
            # retired from reporting (docs/hor_order_confidence_design.md); these
            # stay as empty vestigial columns only because the per-TRA consensus
            # R alias shim (consensus_ensemble.R) selects HOR_status/_confidence.
            "hor_status":        "",
            "hor_confidence":    "",
            # hor_order_confidence: founder-arithmetic HOR tier that replaces the
            # spectral classifier — none | strict | supported | weak.
            "hor_order_confidence": _hor_order_confidence(
                e["multiplicity"], e["multiplicity_raw"], e["irregular"],
                e.get("rescue_method"), e["fallback"]),
            "founder_period":    _fmt_bp(e["founder_period"]),
            "strongest_period":  _fmt_bp(e["strongest_period"]),
            # kitehor's original rescore pick for this case_id. When it
            # differs from `strongest_period` above, TideCluster's
            # argmax(identity_med) override has fired (kitehor picked a
            # lower-id_med peak by score). Useful for auditing the
            # override; not required by downstream consumers.
            "kitehor_founder_period": _fmt_bp(e.get("kitehor_founder_period")),
            "multiplicity":      str(e["multiplicity"]),
            "multiplicity_raw":  (f"{e['multiplicity_raw']:.4f}"
                                  if e["multiplicity_raw"] is not None else ""),
            "irregular_multiplicity": "true" if e["irregular"] else "false",
            "delta_id_pp":       _fmt_pp(delta_id_pp),
            "founder_id_med":    _fmt_id(founder_id),
            "strongest_id_med":  _fmt_id(strongest_id),
            "founder_fallback":  "true" if e["fallback"] else "false",
            "ssr_founder_override": "true" if e.get("ssr_override") else "false",
            # Cluster-rescue diagnostics + founder_method enum
            # (strict | none | kh_deeper | ladder | pass2 | pass3 | ssr |
            #  cluster | fallback).
            "founder_method":     e.get("rescue_method") or "",
            "cluster_rescue":     "true" if e.get("cluster_rescue") else "false",
            "subrepeat_1_period": _fmt_bp(sr1[1]) if sr1 else "",
            "subrepeat_1_occ":    _fmt_id(sr1[0]) if sr1 else "",
            "subrepeat_1_tier":   sr1[2] if sr1 else "",
            "subrepeat_2_period": _fmt_bp(sr2[1]) if sr2 else "",
            "subrepeat_2_occ":    _fmt_id(sr2[0]) if sr2 else "",
            "subrepeat_2_tier":   sr2[2] if sr2 else "",
            "subrepeat_founder_divisor_flag": "true" if sr_edge_flag else "false",
            "subrepeat_founder_divisor_k":    str(sr_edge_k) if sr_edge_k else "",
            # Short-founder review aids (kitehor >= 0.13.2 short-monomer scan):
            # weak_short_founder_flag marks a short founder (<= 30 bp) whose
            # founder peak has weak kite support (score < 0.20); alt_longer_period
            # is the dominant longer-period alternative to eyeball against it.
            "weak_short_founder_flag": "true" if weak_short else "false",
            "alt_longer_period":       _fmt_bp(alt_longer) if alt_longer else "",
            # Alternative periodicities (top-N clusters that aren't the founder).
            "alt_cluster_1_period":       _fmt_bp(alt1["median_period"])  if alt1 else "",
            "alt_cluster_1_score_sum":    _fmt_sc(alt1["score_sum"])      if alt1 else "",
            "alt_cluster_1_n_peaks":      str(alt1["n_peaks"])            if alt1 else "",
            "alt_cluster_1_cov_frac_max": _fmt_id(alt1["cov_frac_max"])   if alt1 else "",
            "alt_cluster_2_period":       _fmt_bp(alt2["median_period"])  if alt2 else "",
            "alt_cluster_2_score_sum":    _fmt_sc(alt2["score_sum"])      if alt2 else "",
            "alt_cluster_2_n_peaks":      str(alt2["n_peaks"])            if alt2 else "",
            "alt_cluster_2_cov_frac_max": _fmt_id(alt2["cov_frac_max"])   if alt2 else "",
            "ssr_flag":           ssr.get("ssr_flag", ""),
            # Per-array CONSENSUS SSR (kitehor consensus_single): coverage of the
            # dominant motif measured against a motif-periodic *consensus* — can
            # be inflated (a consensus artifact on motif-rich satellites, e.g.
            # 96 % on a 6 %-ATC array). Kept for transparency, named so it cannot
            # be mistaken for the real per-array content.
            "ssr_consensus_dominant_motif": ssr.get("dominant_motif", ""),
            "ssr_consensus_dominant_motif_coverage_pct":
                                  ssr.get("dominant_motif_coverage_pct", ""),
            "ssr_consensus_total_coverage_pct": ssr.get("total_ssr_coverage_pct", ""),
            "ssr_consensus_top_motifs": ssr.get("top_motifs", ""),
            # Per-array RAW SSR: actual per-sequence composition (varies between
            # arrays, usually much lower than the consensus). The honest per-array
            # SSR content — this is what the report's SSR column shows.
            "ssr_raw_dominant_motif": ssr.get("ssr_raw_dominant_motif", ""),
            "ssr_raw_dominant_motif_coverage_pct":
                                  ssr.get("ssr_raw_dominant_motif_coverage_pct", ""),
            "ssr_raw_total_coverage_pct": ssr.get("ssr_raw_total_coverage_pct", ""),
            "ssr_raw_n_regions":  ssr.get("ssr_raw_n_regions", ""),
            "ssr_raw_top_motifs": ssr.get("ssr_raw_top_motifs", ""),
            "consensus_period_bp":ssr.get("consensus_period_bp", ""),
        })

    def _sort_key(r):
        try: return (0, int(r["TRC_ID"].split("_")[1]))
        except (IndexError, ValueError): return (1, r["TRC_ID"])
    rows.sort(key=_sort_key)

    header = [
        "TRC_ID", "seqid", "start", "end", "array_length",
        "monomer_size",   "score",
        "monomer_size_2", "score_2",
        "monomer_size_3", "score_3",
        "monomer_size_4", "score_4",
        "monomer_size_5", "score_5",
        "hor_status", "hor_confidence", "hor_order_confidence",
        "founder_period", "strongest_period", "kitehor_founder_period",
        "multiplicity",
        "multiplicity_raw", "irregular_multiplicity",
        "delta_id_pp", "founder_id_med", "strongest_id_med",
        "founder_fallback", "ssr_founder_override",
        "founder_method", "cluster_rescue",
        "subrepeat_1_period", "subrepeat_1_occ", "subrepeat_1_tier",
        "subrepeat_2_period", "subrepeat_2_occ", "subrepeat_2_tier",
        "subrepeat_founder_divisor_flag", "subrepeat_founder_divisor_k",
        "weak_short_founder_flag", "alt_longer_period",
        "alt_cluster_1_period", "alt_cluster_1_score_sum",
        "alt_cluster_1_n_peaks", "alt_cluster_1_cov_frac_max",
        "alt_cluster_2_period", "alt_cluster_2_score_sum",
        "alt_cluster_2_n_peaks", "alt_cluster_2_cov_frac_max",
        "ssr_flag",
        "ssr_consensus_dominant_motif",
        "ssr_consensus_dominant_motif_coverage_pct",
        "ssr_consensus_total_coverage_pct", "ssr_consensus_top_motifs",
        "ssr_raw_dominant_motif", "ssr_raw_dominant_motif_coverage_pct",
        "ssr_raw_total_coverage_pct", "ssr_raw_n_regions", "ssr_raw_top_motifs",
        "consensus_period_bp",
    ]
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return len(rows)
