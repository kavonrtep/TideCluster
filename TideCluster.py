#!/usr/bin/env python
"""Wrapper of TideHunter
Input is a fasta file of DNA sequences
Output is a table of detected tandem repeats by tidehunter
this wrapper split DNA sequence to chunks, run tidehunter on each chunk and merge results
coordinates in the table must be recalculated to the original sequence

"""
import glob
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool
import tc_utils as tc
import argparse
from version import __version__

# minimal python version is 3.6
assert sys.version_info >= (3, 6), "Python 3.6 or newer is required"


def tarean(prefix, gff, fasta=None, cpu=4, min_total_length=50000, args=None,
           version=__version__):
    """
    Run tarean on genomic sequences specified in gff3 file
    from gff record extract sequence end analyse it with tarean algorithm
    :param prefix: prefix for output files, if gff is not specied it is used to find gff3
    file
    :param gff: gff3 file with tidehunter results
    :param fasta: reference fasta
    :param cpu: number of cpu cores to use
    :param min_total_length: minimal total length of sequences to run tarean
    :return:
    """
    script_path = os.path.dirname(os.path.realpath(__file__))
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

    fasta_dict = tc.extract_sequences_from_gff3(gff, fasta, tarean_dir_fasta)
    # get total length of all sequences in fasta file
    input_fasta_length = sum([i for i in tc.read_fasta_sequence_size(fasta).values()])
    l_debug = 0
    cmd_list = []

    ssr = {}
    with open(gff, "r") as f:
        for i in f:
            if i.startswith("#"):
                continue
            gff_record = tc.Gff3Feature(i)
            if "ssr" in gff_record.attributes:
                ssr[gff_record.attributes_dict['Name']] = gff_record.attributes_dict[
                    "ssr"]

    print("preparaing sequences for TAREAN")
    omitted_clusters = []
    trc_total_length = int(0)
    trc_total_length_omitted = 0
    for k, v in fasta_dict.items():
        l_debug += 1
        v_basename = os.path.basename(v)
        with open(v, "r") as f:
            seqs = tc.read_single_fasta_to_dictionary(f)
        total_length = sum([len(i) for i in seqs.values()]) / 2  # it is dimer!!
        trc_total_length += total_length
        if total_length < min_total_length:
            trc_total_length_omitted += total_length
            print(
                    F"total length of sequences in {k} is less than {min_total_length} "
                    F"nt, "
                    F"skipping"
                    )
            omitted_clusters.append((k, total_length, len(seqs)))
            continue
        # does seqs contain only one sequence?
        if len(seqs) > 1:
            # more than one sequence in fasta file, set orientation
            seqs2 = tc.group_sequences_by_orientation(seqs, k=8)
            for n in seqs:
                if n in seqs2['reverse']:
                    seqs[n] = tc.reverse_complement(seqs[n])
        # save oriented sequences to file
        tc.save_fasta_dict_to_file(seqs, v)
        # get current script path
        # do not run if ssrs
        if k in ssr:
            print(F"{k} is SSR, skipping")
            continue
        tarean_out = F"{tarean_dir}/{v_basename}_tarean"
        cmd = F"{script_path}/tarean/tarean.R -i {v} -s 0 -n {1} -o {tarean_out}"
        cmd_list.append(cmd)
    # run cmd tarean in parallel using multiprocessing module
    if len(omitted_clusters) > 0:
        with open(F"{tarean_dir}/omitted_clusters.txt", "w") as f:
            f.write("cluster_id\ttotal_length\tnumber_of_arrays\n")
            # sort by total length
            omitted_clusters.sort(key=lambda x: x[1], reverse=True)
            for i in omitted_clusters:
                f.write(F"{i[0]}\t{i[1]}\t{i[2]}\n")

        # RUN kite
    cmd = (F"{script_path}/tarean/kite.R -d {prefix}_tarean/fasta"
           F" -p {prefix} -c {cpu}")
    print("Making tandem repeat array profiles.")
    tc.run_cmd(cmd)


    # copy index.html as prefix_index.html
    html_src = script_path + "/tarean/index.html"
    html_dst = F"{prefix}_index.html"
    # Format string to report analysis settings (from args and version)
    settings = (F"Input file                 : {args.fasta}\n"
                F"Prefix                     : {args.prefix}\n"
                F"Minimum TRC total length   : {args.min_total_length}\n"
                F"Minimum array length       : {args.min_length}\n"
                F"Dust filter                : {'no' if args.no_dust else 'yes'}\n"
                F"TideHunter arguments       : {args.tidehunter_arguments}\n"
                F"CPU                        : {args.cpu}\n"
                F"Library                    : {args.library}\n"
                F"TideCluster version        : {version}\n")

    # create summar with number of clusters, number of omitted clusters, number of SSRs
    # total length all clusters, total length of reference
    summary = (
        F"Number of TRCs                  : {l_debug}\n"
        F"Number of TRCs above threshold  : {len(omitted_clusters)}\n"
        F"Number of SSRs in TRCs          : {len(ssr)}\n"
        F"Total length of TRCs            : {int(trc_total_length)} nt\n"
        F"Number of TRAs                  :"
        F" {sum([len(i) for i in fasta_dict.values()])}\n"
        F"Input sequence length           : {input_fasta_length} nt\n"
    )

    # replace all PREFIX_PLACEHOLDER with prefix value and save to new file
    # replace all SETTINGS_PLACEHOLDER with settings value and save to new file
    with open(html_src, "r") as f, open(html_dst, "w") as f2:
        for line in f:
            new_line = line.replace("PREFIX_PLACEHOLDER", prefix)
            new_line = new_line.replace("SETTINGS_PLACEHOLDER", settings)
            new_line = new_line.replace("SUMMARY_PLACEHOLDER", summary)
            f2.write(new_line)

    if len(cmd_list) == 0:
        print("No remaining sequences for TAREAN analysis; exiting")
        with open(F"{prefix}_tarean_report.html", "w") as f:
            f.write("No TRC passed the minimum total length threshold for TAREAN "
                    "analysis. ")
        return

    print("running TAREAN")
    with Pool(cpu) as p:
        total_jobs = len(cmd_list)
        completed_jobs = 0
        for _ in p.imap(tc.run_cmd, cmd_list):
            completed_jobs += 1
            print(F"completed {completed_jobs} of {total_jobs}")

    print("TAREAN finished")
    # get SSR info for tarean report from gff3 file

    # export SSR info to csv file
    with open(F"{tarean_dir}/SSRS_summary.csv", "w") as f:
        for k, v in ssr.items():
            f.write(F"{k}\t{v}\n")

    # final tarean report:
    cmd = (F"{script_path}/tarean/tarean_report.R -i {tarean_dir} -o"
           F" {prefix}_tarean_report -g {gff}")
    print("Making final tarean report.")
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
    gff_short = None
    gff_short_annot = None
    if gff is None:
        gff = prefix + "_clustering.gff3"
        gff_short = prefix + "_tidehunter_short.gff3"
        gff_short_annot = prefix + "_tidehunter_short_annotation.gff3"
    gff_out = prefix + "_annotation.gff3"
    gff3_dir_split_files = prefix + "_annotation_split_files"
    # get list consensus sequences from consensus directory
    # naming scheme is TRC_10_dimer.fasta
    # use glob to get all files in directory
    consensus_files = glob.glob(consensus_dir + "/TRC*dimers.fasta")
    # it is possible that consensus_files does not exist
    if len(consensus_files) > 0:
        print(F"Annotating based on consensus sequences in {consensus_dir}")
        # conncatenate all consensus sequences to one file
        consensus_files_concat = consensus_dir + "/consensus_sequences_all.fasta"
        if not os.path.exists(consensus_files_concat):
            with open(consensus_files_concat, "w") as f:
                for consensus_file in consensus_files:
                    with open(consensus_file, "r") as f_in:
                        for line in f_in:
                            f.write(line)
        seq_lengths = tc.read_fasta_sequence_size(consensus_files_concat)
        # run RepeatMasker on concatenated consensus sequences
        cmd = (F"RepeatMasker -pa {cpu} -lib {library} -e ncbi -s -no_is -norna "
               F"-nolow "
               F"-dir {consensus_dir} {consensus_files_concat}")
        subprocess.check_call(cmd, shell=True)

        rm_file = F"{consensus_files_concat}.out"
        # parse RepeatMasker output
        rm_annotation = tc.get_repeatmasker_annotation(rm_file, seq_lengths, prefix)
        # add annotation to gff3 file, only if gff3 file exists
        if os.path.exists(gff):
            tc.add_attribute_to_gff(gff, gff_out, "Name", "annotation", rm_annotation)
            tc.split_gff3_by_cluster_name(gff_out, gff3_dir_split_files)
        else:
            print(F"gff3 file {gff} does not exist, no annotation added to gff3 file")
    else:
        # when consensus is not available, it is expected that gff3 file contains
        # consensus sequences, which can be annotated
        print("No consensus sequences found in {consensus_dir}")
        print(F"Annotating based on consensus sequences in stored in {gff}")
        tc.annotate_gff(gff, gff_out, library, cpu=cpu)
    if gff_short is not None:
        if os.path.exists(gff_short):
            print('Running annotation of omitted short regions from TideHunter')
            tc.annotate_gff(gff_short, gff_short_annot, library, cpu=cpu)


def clustering(fasta, prefix, gff3=None, min_length=None, dust=True, cpu=4):
    """
    Run clustering on sequences defined in gff3 file and fasta file
    produce gff3 file with cluster information
    :param fasta: fasta file with sequences
    :param prefix: prefix - base naame for input and output files
    :param gff3: gff3 file with tidehunter results
    :param min_length: minimal length of repeat to be included in clustering
    :param dust: use dust filter in blast search
    :param cpu: number of cpu cores to use
    :return:

    """
    gff3_out = prefix + "_clustering.gff3"
    gff3_dir_split_files = prefix + "_clustering_split_files"
    fasta = fasta
    if gff3 is None:
        gff3 = prefix + "_tidehunter.gff3"
    if min_length is not None:
        print('running filtering on gff3 file')
        gff3 = tc.filter_gff_by_length(
                gff3,
                gff_short=prefix + "_tidehunter_short.gff3",
                min_length=min_length
                )

    # filtering on duplicates in gff3 file
    gff3 = tc.filter_gff_remove_duplicates(gff3)
    # check if gff3 has more than 1 sequence, it is enough to read
    # just beginning of the file
    with open(gff3, "r") as f:
        count = 0
        for i in f:
            if i.startswith("#"):
                continue
            count += 1
            if count > 1:
                break
    if count == 0:
        print("No tandem repeats found in gff3 file after filtering, exiting")
        exit(0)
    # get consensus sequences for clustering
    consensus_file = tempfile.NamedTemporaryFile(delete=False).name
    consensus_dimers_file = tempfile.NamedTemporaryFile(delete=False).name
    with open(consensus_file, "w") as f, open(consensus_dimers_file, "w") as f2:
        for seq_id, seq, cons in tc.gff3_to_fasta(gff3, fasta, "consensus_sequence"):
            mult = round(1 + 10000 / len(cons))
            consensus = cons * mult
            consensus_dimers = cons * 4
            if len(consensus) > 10000:
                consensus = consensus[0:10000]
            # write consensus sequence to file
            f.write(F">{seq_id}\n{consensus}\n")
            f2.write(F">{seq_id}\n{consensus_dimers}\n")
    # run dustmasker first, sequences which are completely masked
    # will not be used in clustering.
    mask_prop = tc.get_ssrs_proportions(consensus_dimers_file)
    # if count mask_prop above 0.9
    ssrs_id = [k for k, v in mask_prop.items() if v > 0.9]
    # remove ssrs from consensus sequences, they will be added back later
    # but not used for clustering
    ssrs_description = {}
    ssrs_seq = {}
    ssrs_dimers = {}
    # iterate over consensus_file and consensus_dimers_file
    consensus_file_filtered = tempfile.NamedTemporaryFile(delete=False).name
    consensus_dimers_file_filtered = tempfile.NamedTemporaryFile(delete=False).name
    with open(consensus_file, "r") as f, open(consensus_dimers_file, "r") as f2:
        for id, seq in tc.read_single_fasta_as_generator(f):
            if id not in ssrs_id:
                with open(consensus_file_filtered, "a") as f_out:
                    f_out.write(F">{id}\n{seq}\n")
    with open(consensus_dimers_file, "r") as f, open(consensus_dimers_file_filtered, "a") as f2:
        for id, seq in tc.read_single_fasta_as_generator(f):
            if id not in ssrs_id:
                f2.write(F">{id}\n{seq}\n")
            else:
                # NOTE if there are high proportion is simple
                # repeats, this could use a lot of memory!
                ssrs_dimers[id] = seq
                ssrs_description[id] = tc.get_ssrs_description(seq)
                ssrs_seq[id] = " ".join(
                        [i.split(" ")[0] for i in ssrs_description[id].split(
                                "\n"
                                )]
                        )

    # find unique ssrs seq
    ssrs_clusters = {}
    ssrs_representative = {}
    for k, ssrs in ssrs_seq.items():
        if ssrs not in ssrs_representative:
            ssrs_representative[ssrs] = k
    for k, ssrs in ssrs_seq.items():
        ssrs_clusters[k] = ssrs_representative[ssrs]
    # recalculate description for each ssrs_cluster
    dimers_ssrs_clusters = {}
    for n, repre_id in ssrs_clusters.items():
        if repre_id not in dimers_ssrs_clusters:
            dimers_ssrs_clusters[repre_id] = []
        dimers_ssrs_clusters[repre_id].append(ssrs_dimers[n])
    # recalculating ssrs description
    ssrs_cluster_description = {}
    for k, v in dimers_ssrs_clusters.items():
        ssrs_cluster_description[k] = tc.get_ssrs_description_multiple(v)
    if os.path.getsize(consensus_file_filtered) == 0:
        print("No tandem repeats left after dustmasking, skipping clustering")
    # first round of clustering by mmseqs2
    clusters1 = tc.find_cluster_by_mmseqs2(consensus_file_filtered, cpu=cpu)
    representative_id = set(clusters1.values())
    consensus_fasta_representative = tempfile.NamedTemporaryFile(delete=False).name
    tc.filter_fasta_file(consensus_dimers_file_filtered, consensus_fasta_representative,
                          representative_id
                          )
    # second round of clustering by blastn
    clusters2 = tc.find_clusters_by_blast_connected_component(
            consensus_fasta_representative, dust=dust, cpu=cpu
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
        clusters_final[k] = ssrs_clusters[k]
        clusters1[k] = k
    # get total size of each cluster, store in dict
    cluster_size = tc.get_cluster_size2(gff3, clusters_final)

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
    ssrs_description_final = {}  # with TRC names
    for k, v in ssrs_cluster_description.items():
        ssrs_description_final[cluster_names[k]] = v

    cons_cls, cons_cls_dimer = tc.add_cluster_info_to_gff3(gff3, gff3_out, clusters_final)

    tc.merge_overlapping_gff3_intervals(gff3_out, gff3_out)

    gff_tmp = gff3_out + "_tmp"

    tc.add_attribute_to_gff(gff3_out, gff_tmp, "Name", "repeat_type", ssrs_info)
    os.rename(gff_tmp, gff3_out)
    tc.add_attribute_to_gff(gff3_out, gff_tmp, "Name", "ssr", ssrs_description_final)
    os.rename(gff_tmp, gff3_out)

    # save also first round of clustering for debugging
    cons_cls1, cons_cls_dimer1_ = tc.add_cluster_info_to_gff3(
            gff3, gff3_out + "_1.gff3", clusters1
            )
    tc.merge_overlapping_gff3_intervals(gff3_out + "_1.gff3", gff3_out + "_1.gff3")
    # split gff3 file to parts by cluster
    tc.split_gff3_by_cluster_name(gff3_out, gff3_dir_split_files)

    # for debugging
    # write consensus sequences by clusters to directory
    #  used gff3_out as base name for directory
    consensus_dir = prefix + "_consensus"
    tc.save_consensus_files(consensus_dir, cons_cls, cons_cls_dimer, )
    consensus_dir = prefix + "_consensus_1"
    tc.save_consensus_files(consensus_dir, cons_cls1, cons_cls_dimer1_)
    # remove all temporary files
    os.remove(consensus_file)
    os.remove(consensus_dimers_file)
    os.remove(consensus_file_filtered)
    os.remove(consensus_dimers_file_filtered)
    os.remove(consensus_fasta_representative)


def tidehunter(fasta, tidehunter_arguments, prefix, cpu=4):
    """
    run tidehunter on fasta file
    :param fasta: file with sequences
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
    results = tc.run_tidehunter(
            fasta_file_chunked, tidehunter_arguments
            )
    with open(output, "w") as out:
        # write GFF3 header
        out.write("##gff-version 3\n")

        with open(results) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                feature = tc.TideHunterFeature(line)
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

    parser.add_argument("-v", "--version", action="version", version=__version__)

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
                  ', default value: %(default)s)'),
            )
    parser_tidehunter.add_argument(
            "-c", "--cpu", type=int, default=4,
            help="Number of CPUs to use"
            )
    # Clustering
    parser_clustering = subparsers.add_parser(
            'clustering', help='Run clustering on TideHunter output'
            )
    parser_clustering.add_argument(
            "-f", "--fasta", help="Reference fasta", required=True
            )

    parser_clustering.add_argument(
            "-m", "--min_length", help="Minimum length of tandem repeat array to be "
                                       "included in clustering step. Shorter arrays are "
                                       "discarded, default (%(default)s)",
            required=False,
            default=5000, type=int
            )

    parser_clustering.add_argument(
            "-pr", "--prefix", help=("Prefix is used as a base name for output files."
                                     "If --gff is not provided, prefix will be also used"
                                     "to identify GFF file from previous tidehunter "
                                     "step"),
            required=True
            )
    parser_clustering.add_argument(
            "-g", "--gff", help=("GFF3 output file from tidehunter step. If not provided "
                                 "the file named 'prefix_tidehunter.gff3' will be used"),
            required=False, default=None
            )
    parser_clustering.add_argument(
            '-nd', '--no_dust', required=False, default=False, action='store_true',
            help='Do not use dust filter in blastn when clustering'
            )
    parser_clustering.add_argument(
            "-c", "--cpu", type=int, default=4, help="Number of CPUs to use"
            )

    # Annotation
    parser_annotation = subparsers.add_parser(
            'annotation', help=('Run annotation on output from clustering step'
                                ' using reference library of tandem repeats')
            )
    parser_annotation.add_argument(
            "-pr", "--prefix", help=("Prefix is used as a base name for output files."
                                     "If --gff is not provided, prefix will be also used"
                                     "to identify GFF3 file from previous clustering "
                                     "step"),
            required=True
            )

    parser_annotation.add_argument(
            "-g", "--gff", help=("GFF3 output file from clustering step. If not provided "
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

    parser_annotation.add_argument(
            "-c", "--cpu", type=int, default=4, help="Number of CPUs to use"
            )

    # tarean
    parser_tarean = subparsers.add_parser(
            'tarean', help='Run TAREAN on clusters to extract representative sequences'
            )
    parser_tarean.add_argument(
            "-g", "--gff", help=("GFF3 output file from annotation or clustering step"
                                 "If not provided the file named "
                                 "'prefix_annotation.gff3' "
                                 "will be used instead. If 'prefix_annotation.gff3' is "
                                 "not "
                                 "found, 'prefix_clustering.gff3' will be used"
                                 ),
            required=False, default=None
            )
    parser_tarean.add_argument(
            "-f", "--fasta", help="Reference fasta", required=True
            )
    parser_tarean.add_argument(
            "-pr", "--prefix", help=("Prefix is used as a base name for output files."
                                     "If --gff is not provided, prefix will be also used"
                                     "to identify GFF3 files from previous clustering/"
                                     "annotation step"),
            required=True
            )
    parser_tarean.add_argument(
            "-c", "--cpu", type=int, default=4, help="Number of CPUs to use"
            )
    parser_tarean.add_argument(
            "-M", "--min_total_length", type=int, default=50000,
            help=("Minimum combined length of tandem repeat arrays within a single "
                  "cluster, required for inclusion in TAREAN analysis."
                  "Default (%(default)s)")
            )

    parser_run_all = subparsers.add_parser(
            'run_all', help='Run all steps of TideCluster'
            )

    parser_run_all.add_argument(
            "-f", "--fasta", help="Reference fasta", required=True, type=str
            )
    parser_run_all.add_argument(
            "-pr", "--prefix", help="Base name used for input and output files",
            required=True, type=str
            )
    parser_run_all.add_argument(
            "-l", "--library", help="Path to library of tandem repeats", required=False,
            type=str, default=None
            )
    parser_run_all.add_argument(
            "-m", "--min_length", help=("Minimum length of tandem repeat"
                                        " (%(default)s)"), required=False,
            default=5000, type=int
            )
    parser_run_all.add_argument(
            '-T', '--tidehunter_arguments', type=str, nargs="?", required=False,
            default="-p 40 -P 3000 -c 5 -e 0.25",
            help=('additional arguments for TideHunter in quotes'
                  ', default value: %(default)s)'),
            )
    parser_run_all.add_argument(
            "-nd", "--no_dust", help="Do not use dust filter in blastn when clustering",
            action="store_true", required=False, default=False
            )
    parser_run_all.add_argument(
            "-c", "--cpu", type=int, default=4, help="Number of CPUs to use"
            )

    parser_run_all.add_argument(
            "-M", "--min_total_length", type=int, default=50000,
            help=("Minimum combined length of tandem repeat arrays within a single "
                  "cluster, required for inclusion in TAREAN analysis."
                  "Default (%(default)s)")
            )

    parser.description = """Wrapper of TideHunter
    This script enable to run TideHunter on large fasta files in parallel. It splits
    fasta file into chunks and run TideHunter on each chunk. Identified tandem repeat 
    are then clustered, annotated and representative consensus sequences are extracted.
    
     
    """

    # make epilog, in epilog keep line breaks as preformatted text

    parser.epilog = ('''
    Example of usage:
    
    # first run tidehunter on fasta file to generate raw GFF3 output
    # TideCluster.py tidehunter -c 10 -f test.fasta -pr prefix 
    
    # then run clustering on the output from previous step to cluster similar tandem 
    repeats
    TideCluster.py clustering -c 10 -f test.fasta -pr prefix -m 5000
    
    # then run annotation on the clustered output to annotate clusters with reference
    # library of tandem repeats in RepeatMasker format
    TideCluster.py annotation -c 10 -pr prefix -l library.fasta
    
    # then run TAREAN on the annotated output to extract representative consensus
    # and generate html report
    TideCluster.py tarean -c 10 -f test.fasta -pr prefix
    
    Recommended parameters for TideHunter:
    short monomers: -T "-p 10 -P 39 -c 5 -e 0.25"
    long monomers: -T "-p 40 -P 3000 -c 5 -e 0.25" (default)
    
    For parallel processing include -c option before command name. 
    
    For more information about TideHunter parameters see TideHunter manual.
    
    Library of tandem repeats for annotation step are sequences in RepeatMasker format
    where header is in format:
    
    >id#clasification
    
    ''')

    cmd_args = parser.parse_args()
    if cmd_args.command == "tidehunter":
        tidehunter(
                cmd_args.fasta, cmd_args.tidehunter_arguments, cmd_args.prefix,
                cmd_args.cpu
                )
    elif cmd_args.command == "clustering":
        clustering(
                cmd_args.fasta, cmd_args.prefix, cmd_args.gff, cmd_args.min_length,
                not cmd_args.no_dust, cmd_args.cpu
                )
    elif cmd_args.command == "annotation":
        annotation(
                cmd_args.prefix, cmd_args.library, cmd_args.gff,
                cmd_args.consensus_directory,
                cmd_args.cpu
                )
    elif cmd_args.command == "tarean":
        tarean(
                prefix=cmd_args.prefix,
                gff=cmd_args.gff,
                fasta=cmd_args.fasta,
                cpu=cmd_args.cpu,
                min_total_length=cmd_args.min_total_length,
                args=cmd_args,
                version=__version__
                )
    elif cmd_args.command == "run_all":
        tidehunter(
                cmd_args.fasta, cmd_args.tidehunter_arguments, cmd_args.prefix,
                cmd_args.cpu
                )
        clustering(
                cmd_args.fasta, cmd_args.prefix,
                min_length=cmd_args.min_length,
                dust=not cmd_args.no_dust,
                cpu=cmd_args.cpu
                )
        if cmd_args.library:
            annotation(
                cmd_args.prefix, cmd_args.library,
                cpu=cmd_args.cpu
                )
        tarean(
                prefix=cmd_args.prefix,
                fasta=cmd_args.fasta,
                gff=None,
                cpu=cmd_args.cpu,
                min_total_length=cmd_args.min_total_length,
                args=cmd_args,
                version=__version__
                )

    else:
        parser.print_help()
        sys.exit(1)
