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

def get_coverages(bam_dirpath, amplicons_filepath, out_prfx, out_dir):
    """Launch chimeric_solver and counts coverages. Save the output .xls file to the directory ./result by default"""
    string_to_cmd = ("").join(["python2 chimeric_solver.py --dir ", bam_dirpath, " --bed ", amplicons_filepath, " --resFile ", out_prfx, ".xls", " --conv --out ", out_dir])
    os.system(string_to_cmd)

def get_quantiles():
    """Reads files with hard coded quantiles and gets quantiles, returns dictionaries, {degree of freedom : quantile} """
    quant_0995 = {}
    quant_099 = {}
    quant_098 = {}
    quant_095 = {}
    with open(os.path.abspath("./quantiles/quantiles0995.txt")) as f:
        table = f.readline()
        array_of_quantiles = table.split()
        for i in xrange(2, 500):
            quant_0995[i] = float(array_of_quantiles[i - 2])
    with open(os.path.abspath("./quantiles/quantiles099.txt")) as f:
        table = f.readline()
        array_of_quantiles = table.split()
        for i in xrange(2, 500):
            quant_099[i] = float(array_of_quantiles[i - 2])
    with open(os.path.abspath("./quantiles/quantiles098.txt")) as f:
        table = f.readline()
        array_of_quantiles = table.split()
        for i in xrange(2, 500):
            quant_098[i] = float(array_of_quantiles[i - 2])
    with open(os.path.abspath("./quantiles/quantiles095.txt")) as f:
        table = f.readline()
        array_of_quantiles = table.split()
        for i in xrange(2, 500):
            quant_095[i] = float(array_of_quantiles[i - 2])

    return quant_0995, quant_099, quant_098, quant_095

def get_tmp_del_files(amplicons_filepath, quant_dict_del, quant_dict_dup, out_prfx, min_corr, mode):
    """Launch CONVector using predefined parameters, output - files with CNVs
    (outputId_of_task_before_LDA and outputId_of_task_after_LDA)"""
    num_of_samples = 0
    with open("./result/" + out_prfx + "_qc.xls") as f:
        num_of_samples = len(f.readline().split()) - 2
    x = str(quant_dict_del[num_of_samples])
    y = str(quant_dict_dup[num_of_samples])
    string_to_cmd = "javac ./deletionsAnalysis/Main.java"
    os.system(string_to_cmd)
    string_to_cmd = ("").join(["java  deletionsAnalysis/Main ", " -d ./result/", out_prfx,
                               "_qc.xls", " -b ", amplicons_filepath, " -f output", out_prfx, " -mc ", min_corr, " -mnm 5 -nne 4 -nod 4 -lb -",
                               x, " -ub ", y , " -dist 1000000 -lcb 25000 -lca 50"])
    if mode:
        string_to_cmd += " -c"
    logger.warn(string_to_cmd)
    os.system(string_to_cmd)

def get_results(out_prfx, out_dir):
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
    parser.add_argument('--out_dir', action='store', dest='out_dir', required=True,
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
    out_dir = args.out_dir
    ## Algorithm parameters
    amplicons_frac = args.amplicons_frac
    min_corr = args.min_corr
    mode = args.mode

    bam_dirpath = os.path.abspath(bam_dirpath)
    
    #if out_prfx == "CONVector_result":
    #    out_prfx += bam_dirpath
    #ControlCoverage_filepath = out_prfx
    #if args.ControlCoverage_filepath:
    #    ControlCoverage_filepath = args.ControlCoverage_filepath

    logger.info((" ").join(["Task with id", out_prfx, "started"]))
    quant_dict = {}
    quant_0995, quant_099, quant_098, quant_095 = get_quantiles()
    if mode == "hard":
        quant_dict_del = quant_099
    elif mode == "normal":
        quant_dict_del = quant_098
    elif mode == "soft":
        quant_dict_del = quant_095
    quant_dict_dup = quant_095

    if not os.path.exists("tmp_chimeras.txt"):
        f = open('tmp_chimeras.txt', 'w')
        f.close()
    if not os.path.exists("tmp_output.txt"):
        f = open('tmp_output.txt', 'w')
        f.close()
    if not os.path.exists("tmp_references.txt"):
        f = open('tmp_references.txt', 'w')
        f.close()

    amplicons_filepath = os.path.abspath(amplicons_filepath)
    if args.bam_dirpath:
        get_coverages(bam_dirpath, amplicons_filepath, out_prfx, out_dir)
    output_file_name = out_prfx
    if not out_prfx == ControlCoverage_filepath:
        output_file_name = out_prfx + "_" + ControlCoverage_filepath
        if float(amplicons_frac) < 0.9:
            amplicons_frac = "0.9"
            logger.info("Too low quality control for merging test and control datasets. Value at least 0.9 should be used.")

    string_to_cmd = ("").join(["python2 qc.py ", amplicons_filepath, " ./result/", out_prfx, ".xls ./result/", ControlCoverage_filepath, ".xls ",  output_file_name, " ", amplicons_frac])
    os.system(string_to_cmd)

    get_tmp_del_files(amplicons_filepath, quant_dict_del, quant_dict_dup, output_file_name, min_corr, control_is_clean)
    get_results(output_file_name, out_dir)
    string_to_cmd = (" ").join(["./visualisation.R", amplicons_filepath, "./visualisation", out_prfx])
    os.system(string_to_cmd)
    logger.info(" ".join(["Task with id", out_prfx, "finished! ;-)"]))


if __name__ == "__main__":
    main()



