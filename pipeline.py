#!/usr/bin/env python2
__author__ = 'german'

import sys
import logging
from logging.config import fileConfig
import os
import argparse
import time

loginipath = ('logging_config.ini')
fileConfig(loginipath, defaults={'logfilename': 'pipeline.log'})
logger = logging.getLogger('CONVector_logger')


def get_quantiles(mode, CONVector_dirpath):
    """
    Reads files with hard-coded quantiles and gets quantiles, returns dictionaries, {degree of freedom : quantile}.
    """
    # Choose quantile based on mode.
    if mode == 'hard':
        quantile = '099'
    elif mode == 'normal':
        quantile = '098'
    elif mode == 'soft':
        quantile = '095'
    
    quantiles_dirpath = os.path.join(CONVector_dirpath, 'quantiles')
    quantiles_file = 'quantiles{}.txt'.format(quantile)
    quantiles_filepath = os.path.join(quantiles_dirpath, quantiles_file)

    quantiles_dict = {}
    with open(quantiles_filepath) as handle:
        quantiles_table = handle.readline()
        quantiles_array = quantiles_table.split()
        for i in xrange(2, 500):
            quantiles_dict[i] = float(quantiles_array[i - 2])

    return quantiles_dict


def get_tmp_del_files(amplicons_filepath, del_quantiles_dict, dup_quantiles_dict, out_prfx, out_dirpath, min_corr, control_is_clean):
    """
    Launch CONVector using predefined parameters, output - files with CNVs
    (outputId_of_task_before_LDA and outputId_of_task_after_LDA)
    """
    amplicon_coverage_file = out_prfx + '_AmplCov_filtered.tsv'
    amplicon_coverage_filepath = os.path.join(out_dirpath, amplicon_coverage_file)
    num_of_samples = 0
    with open(amplicon_coverage_filepath) as f:
        num_of_samples = len(f.readline().split()) - 2
    #
    lower_bound = '-' + str(del_quantiles_dict[num_of_samples])
    upper_bound = str(dup_quantiles_dict[num_of_samples])
    #
    compile_java_command = 'javac ./deletionsAnalysis/Main.java'
    os.system(compile_java_command)
    #
    call_variants_command = 'java  deletionsAnalysis/Main  -d {} -b {} -f {} -mc {} -lb {} -ub {} -mnm 5 -nne 4 -nod 4 -dist 1000000 -lcb 25000 -lca 50'.format(amplicon_coverage_filepath, amplicons_filepath, out_prfx, min_corr, lower_bound, upper_bound)
    if control_is_clean:
        call_variants_command += ' -c'
    # Options in deletionsAnalysis/OptionsParse.java:
    ## -b    Path to amplicon coordinates table.
    ## -d    Path to amplicon covarage table.
    ## -f    Output file prefix.
    ## -lb   Lower bound.
    ## -ub   Upper bound.
    logger.warn(call_variants_command)
    os.system(call_variants_command)


def get_results(out_prfx, out_dirpath):
    """Launch finalizer and creates a report about CNVs, .xls file with name
    result_before_LDA_id_of_task.xls and result_after_LDA_id_of_task.xls"""
    string_to_get_final_results_unsupervised = ("").join(["python2 finalizer.py -i output", out_prfx, ".xls -o result_before_LDA_", out_prfx, ".xls"])
    string_to_get_final_results_supervised = ("").join(["python2 finalizer.py -i output", out_prfx, "_after_LDA.xls -o result_after_LDA_", out_prfx, ".xls"])
    os.system(string_to_get_final_results_unsupervised)
    os.system(string_to_get_final_results_supervised)


def main():
    # Create arguments parser
    parser = argparse.ArgumentParser(description="""This is a pipeline for the detection of CNV in the data
    obtained with parallel target sequencing and AmpliSeq. For changing the parameters you can look through this Python script""")
    
    # Add arguments to the arguments parser
    ## Input
    parser.add_argument('--ampls', '-a', action='store', dest='amplicons_filepath', required=True,
                        help='Path to TSV table with amplicons coordinates.')
    parser.add_argument('--bam_dir', '-b', action='store', dest='bam_dirpath', required=True,
                        help='Path to directory with input BAM files.')
    #parser.add_argument('--control','-cs', action="store", dest = "control_dataset",
    #                    default = "", required=False,
    #                    help = 'Id of control dataset (to perform quality control and to merge them)')    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    parser.add_argument('--control_cov', action='store', dest="ControlCoverage_filepath", default='', required=False,
                        help='Path to TSV table with control data set amplicons coverage.')
    parser.add_argument('--clean_control', action='store_true', dest='control_is_clean', default = False, required=False,
                        help='Set this flag if a control data set is free of CNVs.')
    ## Output
    parser.add_argument('--out_prfx', action='store', dest='out_prfx', default='CONVector_result', required=False,
                        help='Output file name prefix; default "CONVector_result".')
    parser.add_argument('--out_dir', action='store', dest='out_dirpath', required=True,
                        help='Path to output directory.')
    ## Algorithm parameters
    parser.add_argument('--ampls_frac', action='store', dest='amplicons_frac', default='0.8', required=False,
                        help='Fraction of amplicons which will be taken into account when determining if a sample has normal coverage; default 0.8.')
    parser.add_argument('--min_corr', action='store', dest='min_corr', default='0.7', required=False,
                        help='Minimum correlation threshold for CNV calling; default 0.7.')
    parser.add_argument('--mode', '-m', action='store', dest='mode', default='normal', required=False,
                        help='CNV calling mode: "hard", "normal", "soft"; default "normal".')

    # Import arguments from the arguments parser
    args = parser.parse_args()
    ## Input
    amplicons_filepath = args.amplicons_filepath
    bam_dirpath = args.bam_dirpath
    ControlCoverage_filepath = args.ControlCoverage_filepath
    control_is_clean = args.control_is_clean
    ## Output
    out_prfx = args.out_prfx
    out_dirpath = args.out_dirpath
    ## Algorithm parameters
    amplicons_frac = args.amplicons_frac
    min_corr = args.min_corr
    mode = args.mode
    
    logger.info('Task with ID {} started'.format(out_prfx))

    CONVector_filepath = os.path.abspath(__file__)
    CONVector_dirpath = os.path.dirname(CONVector_filepath)

    del_quantiles_dict = get_quantiles(mode, CONVector_dirpath)
    dup_quantiles_dict = get_quantiles('soft', CONVector_dirpath)

    out_dirpath = os.path.abspath(out_dirpath)
    tmp_dirpath = os.path.join(out_dirpath, 'tmp')
    if not os.path.exists(tmp_dirpath):
        os.makedirs(tmp_dirpath)
    tmpChimeras_filepath = os.path.join(tmp_dirpath, 'tmp_chimeras.txt')
    tmpOutput_filepath = os.path.join(tmp_dirpath, 'tmp_output.txt')
    tmpReferences_filepath = os.path.join(tmp_dirpath, 'tmp_references.txt')
    if not os.path.exists(tmpChimeras_filepath):
        f = open(tmpChimeras_filepath, 'w')
        f.close()
    if not os.path.exists(tmpOutput_filepath):
        f = open(tmpOutput_filepath, 'w')
        f.close()
    if not os.path.exists(tmpReferences_filepath):
        f = open(tmpReferences_filepath, 'w')
        f.close()

    bam_dirpath = os.path.abspath(bam_dirpath)
    amplicons_filepath = os.path.abspath(amplicons_filepath)

    # Calculate amplicon coverages. (Run chimeric_solver.py).
    chimeric_solver_filepath = os.path.join(CONVector_dirpath, 'chimeric_solver.py')
    coverage_file = '{}_AmplCov.tsv'.format(out_prfx)
    chimeric_solver_command = 'python2 {} --convert --ampls {} --bam_dir {} --out_file {} --out_dir {}'.format(chimeric_solver_filepath,
                                                                                                               amplicons_filepath,
                                                                                                               bam_dirpath,
                                                                                                               coverage_file,
                                                                                                               out_dirpath)
    os.system(chimeric_solver_command)
    
    #output_file_name = out_prfx
    #if not out_prfx == ControlCoverage_filepath:
    #    output_file_name = out_prfx + "_" + ControlCoverage_filepath
    #    if float(amplicons_frac) < 0.9:
    #        amplicons_frac = "0.9"
    #        logger.info("Too low quality control for merging test and control datasets. Value at least 0.9 should be used.")

    qc_filepath = os.path.join(CONVector_dirpath, 'qc.py')
    coverage_filepath = os.path.join(out_dirpath, coverage_file)
    qc_command = 'python2 {} --ampls {} --cov {} --control_cov {} --out_prfx {} --out_dir {} --ampls_frac {}'.format(qc_filepath, amplicons_filepath, coverage_filepath,
                                                                                                                     ControlCoverage_file, out_prfx, out_dirpath, amplicons_frac)
    os.system(qc_command)

    get_tmp_del_files(amplicons_filepath, del_quantiles_dict, dup_quantiles_dict, out_prfx, out_dirpath, min_corr, control_is_clean)
    get_results(output_file_name, out_dirpath)
    string_to_cmd = (" ").join(["./visualisation.R", amplicons_filepath, "./visualisation", out_prfx])
    os.system(string_to_cmd)
    logger.info(" ".join(["Task with id", out_prfx, "finished! ;-)"]))


if __name__ == "__main__":
    main()



