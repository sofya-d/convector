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

def parse_file_with_coverages(filename):
    """
    :param filename: name of file with coverages, .xls, format described in manual
    :return: dictionaries 'samples' {sample : {amplicon : coverages}};
             list with low covered amplicons;
             coverages of clean samples (without low covered amplicons, for statistical analysis);
             coverages of all samples (for output)
    """
    samples = defaultdict(dict)
    true_coverages_clean = defaultdict(dict)
    true_coverages = defaultdict(dict)
    total_number_of_samples = 0
    low_covered_ampls = []
    with open(filename) as f:
        array_of_lines = f.readlines()
        samples_names = []
    for i in xrange(len(array_of_lines)):
        line = array_of_lines[i]
        splitted_line = line.split()
        if not i:
            samples_names = splitted_line[2:]
            total_number_of_samples = len(samples_names)
        else:
            summa_of_coverages = 0
            ampl_name = splitted_line[1]
            for i in xrange(len(samples_names)):
                value = 0
                try:
                    value = (int(splitted_line[i + 2]))
                except ValueError:
                    break
                true_coverages_clean[samples_names[i]][ampl_name] = value
                true_coverages[samples_names[i]][ampl_name] = value
                # 10 - threshold for homozygous deletion
                if value > 10:
                    samples[samples_names[i]][ampl_name] = log(float(value))
                else:
                    samples[samples_names[i]][ampl_name] = 0.01
                    # 0.01 - just very low value, indicates a homo del
                summa_of_coverages += value
            if summa_of_coverages / total_number_of_samples < 50:
                for name in samples_names:
                    samples[name].pop(ampl_name, None)
                    true_coverages_clean[name].pop(ampl_name, None)
                # exclude low covered (in all samples) amplicons using threshold 50
                low_covered_ampls.append(ampl_name)
    return samples, low_covered_ampls, true_coverages_clean, true_coverages

# parse_file_with_coverages() ===============================================================================================================================================================


# return_quantiles_chisq() ==================================================================================================================================================================

def return_quantiles_chisq():
    """
    Reads files with quantiles of chi squared distribution, degrees of freedom: from 1 to 500 (depending on number
    of samples).
    :return: tuple of dictionaries, {degree of freedom: corresponding quantile}
    """
    quantiles99 = {}
    quantiles95 = {}
    with open("./quantiles/quantile_chisq_99.txt") as f:
        counter = 1
        for line in f:
            quantiles99[counter] = float(line.split()[0])
            counter += 1
    with open("./quantiles/quantile_chisq_95.txt") as f:
        counter = 1
        for line in f:
            quantiles95[counter] = float(line.split()[0])
            counter += 1
    return quantiles99, quantiles95

# return_quantiles_chisq() ==================================================================================================================================================================


# normalize_data_inside_chromosomes() =======================================================================================================================================================

def normalize_data_inside_chromosomes(clean_chromosomes_amplicons, samples):
    """
    :param clean_chromosomes_amplicons: coverages of amplicons inside one chromosome
    :param samples: list with names
    :return: normalized coverages
    """
    total_amount_sample_chromosome = defaultdict(dict)
    for sample, data in samples.iteritems():
        for chromosome, amplicons in clean_chromosomes_amplicons.iteritems():
            try:
                summa_of_coverages = [data[amplicon] for amplicon in amplicons]
            except KeyError:
                logger.warn("Smth went wrong...\nThe most probable thing - you are using the wrong file with coverages!")
                sys.exit(1)
            for amplicon in amplicons:
                data[amplicon] -= statistics.mean(summa_of_coverages)
            total_amount_sample_chromosome[sample][chromosome] = sum(summa_of_coverages)
    return total_amount_sample_chromosome

# normalize_data_inside_chromosomes() =======================================================================================================================================================


# form_ellipsoid() ==========================================================================================================================================================================

def form_ellipsoid(samples_to_train, amplicons_from_chromosome):
    """
    :param samples_to_train: dictionary {amplicon name : coverages in samples from training samples}
    :param amplicons_from_chromosome: coverages of samples from control dataset, list
    :return:
    """
    ellipsoid = {}
    for amplicon in amplicons_from_chromosome:
        amplicon_values = []
        for name, info in samples_to_train.iteritems():
            amplicon_values.append(info[amplicon])
        ellipsoid[amplicon] = (statistics.medianW(amplicon_values), statistics.sn_estimator(amplicon_values) ** 2)
    return ellipsoid

# form_ellipsoid() ==========================================================================================================================================================================


# diagnose_chromosome_ellipsoid() ===========================================================================================================================================================

def diagnose_chromosome_ellipsoid(samples_to_test, ellipsoid, list_of_amplicons_to_test, qChisq, num_of_accepted):
    """

    :param samples_to_test: coverages in samples and amplicons, dict {sample name : {ampl name : coverage} }
    :param ellipsoid: pairs for each amplicon, (estimation of mean, estimation of standard deviation)
    :param list_of_amplicons_to_test: list of amplicons from one chromosome without low covered amplicons
    :param qChisq: corresponding chi square quantile
    :param num_of_accepted: number of amplicons that will be taken into account
    :return: dict of lists with normal or not coverages inside one chromosome ( 1 = normal, 0 = irregular);
             average robust residuals for calculation of ARV
    """
    normal_or_not = defaultdict(int)
    statistic_for_sample_and_chromosome = defaultdict(list)
    avtc_residuals_for_amplicons = defaultdict(list)

    for sample, coverages_of_amplicons in samples_to_test.iteritems():
        for ampl in list_of_amplicons_to_test:
            distance_to_mean = (coverages_of_amplicons[ampl] - ellipsoid[ampl][0])
            dist = (distance_to_mean ** 2) / (ellipsoid[ampl][1])
            statistic_for_sample_and_chromosome[sample].append(dist)
            avtc_residuals_for_amplicons[ampl].append(distance_to_mean ** 2)
    for sample, statistic_values in statistic_for_sample_and_chromosome.iteritems():
        if sum(sorted(statistic_values)[:num_of_accepted]) < qChisq:
            normal_or_not[sample] = 1
        else:
            normal_or_not[sample] = 0
    return normal_or_not, avtc_residuals_for_amplicons

# diagnose_chromosome_ellipsoid() ===========================================================================================================================================================


# form_list_of_qc_negative() ================================================================================================================================================================

def form_list_of_qc_negative(list_of_normal_or_not_dicts, samples_to_test_qc):
    """
    :param list_of_normal_or_not_dicts: dict {sample : list of 0 and 1, determining the irregularity of coverage inside
           chromomsome}
    :param samples_to_test_qc: list of sample in test dataset, names
    :return: list of samples that did not passed our QC algorithm
    """
    dict_of_sums = {}
    num_of_required_good_chromosomes = len(list_of_normal_or_not_dicts) - 1
    list_of_negatives = []
    for sample in samples_to_test_qc:
        dict_of_sums[sample] = 0
        for dictionary in list_of_normal_or_not_dicts:
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

def output_result_file(true_coverages_of_samples_to_test, true_coverages_of_samples_to_train, qc_negative_list, clean_chromosomes_amplicons, out_prfx, mode):
    """
    :param true_coverages_of_samples_to_test: coverages before normalization
    :param true_coverages_of_samples_to_train: coverages before normalization
    :param qc_negative_list: list of samples that did not passed QC
    :param clean_chromosomes_amplicons: amplicons from clean samples
    :param out_prfx: name of task
    :param mode: merging control and test sample or not
    :return: output file with coverages and QC report
    """
    ordered_list_of_samples_test = list(sorted(true_coverages_of_samples_to_test.iterkeys()))
    ordered_list_of_samples_control = list(sorted(true_coverages_of_samples_to_train.iterkeys()))
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
        for chr, ampls in clean_chromosomes_amplicons.iteritems():
            for ampl in ampls:
                new_result_string = "N/A\t" + ampl + "\t"
                for sample in ordered_list_of_samples_test:
                    if sample not in qc_negative_list:
                        new_result_string += str(true_coverages_of_samples_to_test[sample][ampl]) + "\t"
                if mode == 1:
                    for sample in ordered_list_of_samples_control:
                        new_result_string += str(true_coverages_of_samples_to_train[sample][ampl]) + "\t"

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
    parser.add_argument('--cov', action='store', dest='coverage_filepath', default = False, required=True,
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
    amplicons_frac = args.amplicons_frac
    
    # Import amplicons' coordinates table
    directory = '/'.join(amplicons_filepath.split('/')[:-1])
    panel_of_amplicons = chimeric_solver.parse_bed_file(amplicons_filepath)

    # Import amplicons' coverage for "train" and "test" sets, gather statistics.
    ## samples_to_<test/train> is a dictionary with structure: {sample : {amplicon : coverages}}.
    ## low_covered_ampls_<test/train> is a list with low covered amplicons.
    ## clean_coverages_of_samples_to_<test/train> contains coverages of samples without low covered amplicons for statistical analysis.
    ## true_coverages_of_samples_to_<test/train> contains coverages of all samples (for output).
    samples_to_test, low_covered_ampls_test, clean_coverages_of_samples_to_test, true_coverages_of_samples_to_test = parse_file_with_coverages(coverage_filepath)
    samples_to_train, low_covered_ampls_train, clean_coverages_of_samples_to_train, true_coverages_of_samples_to_train = parse_file_with_coverages(ControlCoverage_filepath)

    # Compute a union of low-covered amplicons in "train" and "test" sets.
    set_of_low_covered_ampls = set(low_covered_ampls_test + low_covered_ampls_train)

    # Set mode.
    ## Mode is equal to 1 if the control (train) and test samples differ. Otherwise, mode = 0.
    mode = 0
    if coverage_filepath != ControlCoverage_filepath:
        mode = 1

    # Create a list of amplicons with low coverage and a list of all amplicons.
    clean_chromosomes_amplicons = defaultdict(list)
    all_amplicon_names = defaultdict(list)
    for key in panel_of_amplicons:
        for amplicon in panel_of_amplicons[key]:
            if amplicon.ID not in set_of_low_covered_ampls:
                clean_chromosomes_amplicons[key].append(amplicon.ID)
            all_amplicon_names[key].append(amplicon.ID)

    # Normalize amplicons' coverages by total chromosome coverage.
    total_amount_sample_chromosome_test = normalize_data_inside_chromosomes(clean_chromosomes_amplicons, samples_to_test)
    total_amount_sample_chromosome_train = normalize_data_inside_chromosomes(clean_chromosomes_amplicons, samples_to_train)

    # Calculate some statistics based on amplicon's coverages, write into the dictionaries.
    ellipsoids_for_chromosomes = {}
    ellipsoids_for_chromosomes_test = {}
    for chr, amplicons_from_chromosome in clean_chromosomes_amplicons.iteritems():
        if chr.startswith("chr"):
            ellipsoids_for_chromosomes[chr] = form_ellipsoid(samples_to_train, amplicons_from_chromosome)
            ellipsoids_for_chromosomes_test[chr] = form_ellipsoid(samples_to_test, amplicons_from_chromosome)

    # Create list_of_robust_variances_test_against_test variable.
    list_of_robust_variances_test_against_test = []
    for chr, ellipsoid_test in ellipsoids_for_chromosomes_test.iteritems():
        list_of_amplicons_to_test = clean_chromosomes_amplicons[chr]
        for amplicon, element in ellipsoid_test.iteritems():
            if amplicon in list_of_amplicons_to_test:
                list_of_robust_variances_test_against_test.append(element[1])

    # Create variables normal_or_not and avtc_residuals_for_amplicons.
    # normal_or_not contains values that signify whether chromosome coverage is regular or not.
    # avtc_residuals_for_amplicons contains statistics for ARV computing.
    quantiles99, quantiles95 = return_quantiles_chisq()
    list_of_robust_variances_control_against_control = []
    list_of_normal_or_not_dicts = []
    for chr, ellipsoid in ellipsoids_for_chromosomes.iteritems():
        list_of_amplicons_to_test = clean_chromosomes_amplicons[chr]
        for amplicon, element in ellipsoid.iteritems():
            if amplicon in list_of_amplicons_to_test:
                list_of_robust_variances_control_against_control.append(element[1])
            num_of_accepted = int(len(list_of_amplicons_to_test) * amplicons_frac)
            qChisq = quantiles99[(num_of_accepted)]
            normal_or_not, avtc_residuals_for_amplicons = diagnose_chromosome_ellipsoid(samples_to_test, ellipsoid, list_of_amplicons_to_test, qChisq, num_of_accepted)
            list_of_normal_or_not_dicts.append(normal_or_not)

    # Get list of samples with irregular chromosomes' coverage.
    samples_to_test_qc = sorted(list(samples_to_test.iterkeys()))
    qc_negative_list = form_list_of_qc_negative(list_of_normal_or_not_dicts, samples_to_test_qc)

    # Calculate ARV statistics for "train" (control) and test data sets.
    avrcc = statistics.mean(list_of_robust_variances_control_against_control)
    avrtt = statistics.mean(list_of_robust_variances_test_against_test)

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
    
    # Write amplicons' coverage without filtered out samples.
    # Add control samples to the table if merging of test and control samples is allowed by mode variable.
    output_result_file(true_coverages_of_samples_to_test, true_coverages_of_samples_to_train, qc_negative_list, all_amplicon_names, out_prfx, mode)

# main() ====================================================================================================================================================================================


main()

