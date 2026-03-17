#!/usr/bin/env python2
__author__ = 'german'

import argparse
import chimeric_solver
import sys
from collections import defaultdict
from math import log
import statistics
import os
import copy
import logging
from logging.config import fileConfig

loginipath = ('logging_config.ini')
fileConfig(loginipath, defaults={'logfilename': 'pipeline.log'})
logger = logging.getLogger('CONVector_logger')


# parse_file_with_coverages() ===============================================================================================================================================================

def parse_file_with_coverages(coverage_filepath):
    """
    param:        coverage_filepath: path to TSV file with amplicons' coverage; format described in manual
    
    return:               ampl_cov_dict: a dictionary with structure {sample: {amplicon: coverage}};
                 ampl_cov_NoLowCov_dict: a dictionary with structure {sample: {amplicon: coverage}} without low covered amplicons;
             log_ampl_cov_NoLowCov_dict: a dictionary with structure {sample: {amplicon: log(coverage)}};
                          low_cov_ampl_names: a list of low covered amplicons names.
    """
    # Read lines from TSV amplicon coverage table, including header.
    with open(coverage_filepath) as handle:
        coverage_table_lines = handle.readlines()

    # Get sample names from header.
    coverage_table_header = coverage_table_lines.pop(0)
    coverage_table_header = coverage_table_header.split()
    sample_names = coverage_table_header[2:]
    sample_num = len(sample_names)

    log_ampl_cov_dict = defaultdict(dict)
    ampl_cov_dict = defaultdict(dict)
    low_cov_ampl_names = []
    homozygous_deletion_threshold = 10
    
    # Iterate through table rows => amplicons
    for i in xrange(len(coverage_table_lines)):
        line = coverage_table_lines[i]
        splitted_line = line.split()
        ampl_name = splitted_line[1]
        coverages = splitted_line[2:]
        coverage_sum = 0

        # Iterate through columns in the row => coverage of an amplicon in each sample
        for j in xrange(len(sample_names)):
            sample_name = sample_names[j]
            try:
                coverage = (int(coverages[j]))
            except ValueError:
                break
            ampl_cov_dict[sample_name][ampl_name] = coverage
            if coverage > homozygous_deletion_threshold:
                # Append logarithm of coverage to the dictionary.
                log_ampl_cov_dict[sample_name][ampl_name] = log(float(coverage))
            else:
                log_ampl_cov_dict[sample_name][ampl_name] = 0.01
                # 0.01 is just a very low coverage; it indicates a homozygous deletion.
            # Total coverage of an amplicon in all samples
            coverage_sum += coverage

        if coverage_sum / sample_num < 50:
            low_cov_ampl_names.append(ampl_name)

    # If average amplicon coverage in a sample is less than 50, exclude the amplicon from all samples
    log_ampl_cov_NoLowCov_dict = log_ampl_cov_dict
    for sample_name in sample_names:
        for low_cov_ampl_name in low_cov_ampl_names:
            log_ampl_cov_NoLowCov_dict[sample_name].pop(low_cov_ampl_name, None)
    
    return ampl_cov_dict, log_ampl_cov_NoLowCov_dict, low_cov_ampl_names

# parse_file_with_coverages() ===============================================================================================================================================================


# return_quantiles_chisq() ==================================================================================================================================================================

def return_quantiles_chisq():
    """
    Reads files with quantiles of chi squared distribution, degrees of freedom: from 1 to 500 (depending on number
    of samples).
    return: tuple of dictionaries, {degree of freedom: corresponding quantile}
    """
    
    qc_filepath = os.path.abspath(__file__)
    CONVector_dirpath = os.path.dirname(qc_filepath)

    quantiles99_filepath = os.path.join(CONVector_dirpath, 'quantiles/quantile_chisq_99.txt')
    quantiles99 = {}
    with open(quantiles99_filepath) as f:
        counter = 1
        for line in f:
            quantiles99[counter] = float(line.split()[0])
            counter += 1

    quantiles95_filepath = os.path.join(CONVector_dirpath, 'quantiles/quantile_chisq_95.txt')
    quantiles95 = {}
    with open(quantiles95_filepath) as f:
        counter = 1
        for line in f:
            quantiles95[counter] = float(line.split()[0])
            counter += 1
    return quantiles99, quantiles95

# return_quantiles_chisq() ==================================================================================================================================================================


# normalize_by_chromosome_coverage() =======================================================================================================================================================

def normalize_by_chromosome_coverage(norm_cov_ampl_names, log_ampl_cov_NoLowCov_dict):
    """
     param:        norm_cov_ampl_names: dictionary of structure {chromosome_name: amplicon_name}
     param: log_ampl_cov_NoLowCov_dict: a dictionary with structure {sample: {amplicon: log(coverage)}}
    """
    for sample, ampl_cov_dict in log_ampl_cov_NoLowCov_dict.iteritems():
        for chrom, ampl_names_lst in norm_cov_ampl_names.iteritems():
            ampl_coverages = [ampl_cov_dict[ampl_name] for ampl_name in ampl_names_lst]
            for ampl_name in ampl_names_lst:
                ampl_cov_dict[ampl_name] -= statistics.mean(ampl_coverages)

# normalize_by_chromosome_coverage() =======================================================================================================================================================


# form_ellipsoid() ==========================================================================================================================================================================

def form_ellipsoid(log_ampl_cov_NoLowCov_dict, ampl_names_lst):
    """
     param: log_ampl_cov_NoLowCov_dict: a dictionary with structure {sample: {amplicon: log(coverage)}}
     param:             ampl_names_lst: list of names of amplicons from one chromosome
    return:                  ellipsoid: a dictionary of structure {amplicon: (stats_estimator_1, stats_estimator_2)}
    """
    ellipsoid = {}
    for ampl_name in ampl_names_lst:
        ampl_cov_lst = []
        for sample, ampl_cov_dict in log_ampl_cov_NoLowCov_dict.iteritems():
            ampl_cov = ampl_cov_dict[ampl_name]
            ampl_cov_lst.append(ampl_cov)
        # statistics.medianW() and statistics.sn_estimator() are imported from <CONVector directory>/statistics.py.
        ## statistics.medianW() calculates Hodges-Lehmann estimator (according to Google Gemini).
        ## statistics.sn_estimator calculates S_n estimator of scale developed by Rousseeuw and Croux (according to Google Gemini).
        ellipsoid[ampl_name] = (statistics.medianW(ampl_cov_lst), statistics.sn_estimator(ampl_cov_lst) ** 2)
    return ellipsoid

# form_ellipsoid() ==========================================================================================================================================================================


# diagnose_chromosome_ellipsoid() ===========================================================================================================================================================

def diagnose_chromosome_ellipsoid(log_ampl_cov_NoLowCov_dict, ampl_ellipsoid_dict, chrom_ampl_names, qChisq, accepted_ampls_num):
    """
     param: log_ampl_cov_NoLowCov_dict: a dictionary with structure {sample: {amplicon: log(coverage)}}
     param:        ampl_ellipsoid_dict: a dictionary of structure {amplicon: (stats_estimator_1, stats_estimator_2)}
     param:           chrom_ampl_names: list of amplicons from one chromosome without low covered amplicons
     param:                     qChisq: corresponding chi square quantile
     param:         accepted_ampls_num: number of amplicons that will be taken into account
    return:      sample_is_normal_dict: a dictionary of structure {sample: int}, where int is 1 if sample's coverage is normal and int is 0 if sample's coverage is irregular
    """
    
    sample_is_normal_dict = defaultdict(int)
    sample_ampl_dist_dict = defaultdict(list)
    
    for sample, ampl_cov_dict in log_ampl_cov_NoLowCov_dict.iteritems():
        for ampl_name in chrom_ampl_names:
            amplicon_coverage = ampl_cov_dict[ampl_name]
            amplicon_medianW = ampl_ellipsoid_dict[ampl_name][0]
            amplicon_Sn_estimator = ampl_ellipsoid_dict[ampl_name][1]
            ampl_dist_to_mean = (amplicon_coverage - amplicon_medianW)
            ampl_dist = (ampl_dist_to_mean ** 2) / (amplicon_Sn_estimator)
            sample_ampl_dist_dict[sample].append(ampl_dist)
    
    for sample, ampl_dist_lst in sample_ampl_dist_dict.iteritems():
        if sum(sorted(ampl_dist_lst)[:accepted_ampls_num]) < qChisq:
            sample_is_normal_dict[sample] = 1
        else:
            sample_is_normal_dict[sample] = 0
    return sample_is_normal_dict

# diagnose_chromosome_ellipsoid() ===========================================================================================================================================================


# form_list_of_qc_negative() ================================================================================================================================================================

def form_list_of_qc_negative(sample_is_normal_dict_lst, samples_to_test_qc):
    """
    :param sample_is_normal_dict_lst: dict {sample : list of 0 and 1, determining the irregularity of coverage inside
           chromomsome}
    :param samples_to_test_qc: list of sample in test dataset, names
    :return: list of samples that did not passed our QC algorithm
    """
    dict_of_sums = {}
    num_of_required_good_chromosomes = len(sample_is_normal_dict_lst) - 1
    list_of_negatives = []
    for sample in samples_to_test_qc:
        dict_of_sums[sample] = 0
        for dictionary in sample_is_normal_dict_lst:
            dict_of_sums[sample] += dictionary[sample]
    counter_of_negative = 0
    counter_of_positive = 0
    for sample in sorted(dict_of_sums.iterkeys()):
        if dict_of_sums[sample] < num_of_required_good_chromosomes:
            logger.info("QC Negative " + sample)
            counter_of_negative += 1
            list_of_negatives.append(sample)
        else:
            counter_of_positive += 1
            logger.info("QC Positive " + sample)
    logger.warn(" ".join(["Overall:", str(counter_of_negative), "samples was filtered out and", str(counter_of_positive), "were accepted"]))
    return list_of_negatives

# form_list_of_qc_negative() ================================================================================================================================================================


# output_result_file() ======================================================================================================================================================================

def output_result_file(test_ampl_cov_dict, train_ampl_cov_dict, qc_negative_list, norm_cov_ampl_names, out_prfx, mode):
    """
    :param test_ampl_cov_dict: coverages before normalization
    :param train_ampl_cov_dict: coverages before normalization
    :param qc_negative_list: list of samples that did not passed QC
    :param norm_cov_ampl_names: amplicons from clean samples
    :param out_prfx: name of task
    :param mode: merging control and test sample or not
    :return: output file with coverages and QC report
    """
    ordered_list_of_samples_test = list(sorted(test_ampl_cov_dict.iterkeys()))
    ordered_list_of_samples_control = list(sorted(train_ampl_cov_dict.iterkeys()))
    with open("./result/" + out_prfx + "_qc.xls", "wb") as f:
        top_string = "Gene\tTarget\t"
        for sample in ordered_list_of_samples_test:
            if sample not in qc_negative_list:
                top_string += "Case_" + sample + "\t"
        if mode == 1:
            for sample in ordered_list_of_samples_control:
                top_string += "Control_" + sample + "\t"
        top_string += '\n'
        f.write(top_string)
        for chr, ampls in norm_cov_ampl_names.iteritems():
            for ampl in ampls:
                new_result_string = "N/A\t" + ampl + "\t"
                for sample in ordered_list_of_samples_test:
                    if sample not in qc_negative_list:
                        new_result_string += str(test_ampl_cov_dict[sample][ampl]) + "\t"
                if mode == 1:
                    for sample in ordered_list_of_samples_control:
                        new_result_string += str(train_ampl_cov_dict[sample][ampl]) + "\t"

                new_result_string += "\n"
                f.write(new_result_string)

# output_result_file() ======================================================================================================================================================================


# main() ====================================================================================================================================================================================

def main():
    parser = argparse.ArgumentParser(description='Control quality of amplicon coverage tables.')
    
    # Add arguments to the arguments parser
    ## Input
    parser.add_argument('--ampls', '-a', action='store', dest='amplicons_filepath', required=True,
                        help='Path to TSV table with amplicons coordinates.')
    parser.add_argument('--cov', action='store', dest='coverage_filepath', required=True,
                        help='Path to TSV table with data set amplicons coverage.')
    parser.add_argument('--control_cov', action='store', dest="ControlCoverage_filepath", default='', required=True,
                        help='Path to TSV table with control data set amplicons coverage.')
    ## Output
    parser.add_argument('--out_prfx', action='store', dest='out_prfx', default='CONVector_result', required=False,
                        help='Output file name prefix; default "CONVector_result".')
    parser.add_argument('--out_dir', action='store', dest='out_dirpath', required=True,
                        help='Path to output directory.')
    ## Algorithm parameters
    parser.add_argument('--ampls_frac', action='store', dest='amplicons_frac', default='0.8', required=False,
                        help='Fraction of amplicons which will be taken into account when determining if a sample has normal coverage; default 0.8.')

    # Import arguments from the arguments parser
    args = parser.parse_args()
    ## Input
    amplicons_filepath = args.amplicons_filepath
    coverage_filepath = args.coverage_filepath
    ControlCoverage_filepath = args.ControlCoverage_filepath
    ## Output
    out_prfx = args.out_prfx
    out_dirpath = args.out_dirpath
    ## Algorithm parameters
    amplicons_frac = float(args.amplicons_frac)
    
    # Import amplicons' coordinates table
    directory = '/'.join(amplicons_filepath.split('/')[:-1])
    panel_of_amplicons = chimeric_solver.parse_bed_file(amplicons_filepath)

    # Import amplicons' coverage for "train" and "test" sets, gather statistics.
    ## samples_to_<test/train> is a dictionary with structure: {sample : {amplicon : coverages}}.
    ## low_covered_ampls_<test/train> is a list with low covered amplicons.
    ## clean_coverages_of_samples_to_<test/train> contains coverages of samples without low covered amplicons for statistical analysis.
    ## true_coverages_of_samples_to_<test/train> contains coverages of all samples (for output).
    test_ampl_cov_dict, test_log_ampl_cov_NoLowCov_dict, test_low_cov_ampl_names = parse_file_with_coverages(coverage_filepath)
    train_ampl_cov_dict, train_log_ampl_cov_NoLowCov_dict, train_low_cov_ampl_names = parse_file_with_coverages(ControlCoverage_filepath)

    # Create a list of amplicons with low coverage and a list of all amplicons.
    low_cov_ampl_names = set(test_low_cov_ampl_names + train_low_cov_ampl_names)
    norm_cov_ampl_names = defaultdict(list)
    ampl_names = defaultdict(list)
    for chrom in panel_of_amplicons:
        for amplicon in panel_of_amplicons[chrom]:
            if amplicon.ID not in low_cov_ampl_names:
                norm_cov_ampl_names[chrom].append(amplicon.ID)
            ampl_names[chrom].append(amplicon.ID)

    # Normalize amplicons' coverages by total chromosome coverage.
    normalize_by_chromosome_coverage(norm_cov_ampl_names, test_log_ampl_cov_NoLowCov_dict)
    normalize_by_chromosome_coverage(norm_cov_ampl_names, train_log_ampl_cov_NoLowCov_dict)

    # Calculate some statistics based on amplicon's coverages, write into the dictionaries.
    train_ellipsoids = {}
    test_ellipsoids = {}
    for chrom, ampl_names_lst in norm_cov_ampl_names.iteritems():
        train_ellipsoids[chrom] = form_ellipsoid(train_log_ampl_cov_NoLowCov_dict, ampl_names_lst)
        test_ellipsoids[chrom] = form_ellipsoid(test_log_ampl_cov_NoLowCov_dict, ampl_names_lst)

    # Create robust_variances_test_vs_test_lst variable.
    robust_variances_test_vs_test_lst = []
    for chrom, ampl_ellipsoid_dict in test_ellipsoids.iteritems():
        chrom_ampl_names = norm_cov_ampl_names[chrom]
        for ampl_name, ellipsoid in ampl_ellipsoid_dict.iteritems():
            if ampl_name in chrom_ampl_names:
                Sn_estimator = ellipsoid[1]
                robust_variances_test_vs_test_lst.append(Sn_estimator)

    # Create a sample_is_normal_dict_lst list.
    # sample_is_normal_dict_lst contains sample_is_normal_dict dictionaries.
    # sample_is_normal_dict is a dictionary of structure {sample: int}, where int=1 if sample's coverage is regular and int=0 if sample's coverage is irregular
    quantiles99, quantiles95 = return_quantiles_chisq()
    robust_variances_control_vs_control_lst = []
    sample_is_normal_dict_lst = []
    for chrom, ampl_ellipsoid_dict in train_ellipsoids.iteritems():
        chrom_ampl_names = norm_cov_ampl_names[chrom]
        accepted_ampls_num = int(len(chrom_ampl_names) * amplicons_frac)
        qChisq = quantiles99[(accepted_ampls_num)]
        sample_is_normal_dict = diagnose_chromosome_ellipsoid(test_log_ampl_cov_NoLowCov_dict, ampl_ellipsoid_dict, chrom_ampl_names, qChisq, accepted_ampls_num)
        sample_is_normal_dict_lst.append(sample_is_normal_dict)
        for ampl_name, ellipsoid in ampl_ellipsoid_dict.iteritems():
            if ampl_name in chrom_ampl_names:
                Sn_estimator = ellipsoid[1]
                robust_variances_control_vs_control_lst.append(Sn_estimator)

    # Get list of samples with irregular chromosomes' coverage.
    samples_to_test_qc = sorted(list(test_log_ampl_cov_NoLowCov_dict.iterkeys()))
    qc_negative_list = form_list_of_qc_negative(sample_is_normal_dict_lst, samples_to_test_qc)

    # Calculate ARV statistics for "train" (control) and test data sets.
    avrcc = statistics.mean(robust_variances_control_vs_control_lst)
    avrtt = statistics.mean(robust_variances_test_vs_test_lst)

    # Write ARV statistics and sample filtration results to qc_control_log.txt.
    with open("qc_control_log.txt","wb") as qc_report:
        arvc = (" ").join(["Average Robust Variance Of Control Dataset with filename", ControlCoverage_filepath, ":", str(avrcc), "\n"])
        qc_report.write(arvc + "\n")
        arvt = (" ").join(["Average Robust Variance Of Test Dataset with filename", coverage_filepath, ":", str(avrtt), "\n"])
        qc_report.write(arvt + "\n")
        logger.info(" ".join(["ARVc =", str(arvc)]))
        logger.info(" ".join(["ARVt =", str(arvt)]))
        qc_report.write((" ").join(["Total amount of samples filtered out using Quality Control algorithm:", str(len(qc_negative_list)), "\n\n"]) )
        for negative_sample in qc_negative_list:
            qc_report.write(("").join([negative_sample, " \n - did not passed QC Control", "\n\n"]))

    # Set mode.
    ## Mode is equal to 1 if the control (train) and test samples differ. Otherwise, mode = 0.
    mode = 0
    if coverage_filepath != ControlCoverage_filepath:
        mode = 1
    
    # Write amplicons' coverage without filtered out samples.
    # Add control samples to the table if merging of test and control samples is allowed by mode variable.
    output_result_file(test_ampl_cov_dict, train_ampl_cov_dict, qc_negative_list, ampl_names, out_prfx, mode)

# main() ====================================================================================================================================================================================


main()

