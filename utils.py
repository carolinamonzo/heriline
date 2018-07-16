#! /usr/bin/env python3.6

# Utils for heriline.py - Pipeline for analysis of hereditary variants
# Date started - 2018-01-24

import argparse
import os
import logging
import coloredlogs
import json
import pprint
import datetime
import subprocess

coloredlogs.install()

def setLogger(mode, module):
    """
    Define the log that we will use
    Input: mode selected on the parsed arguments
    Output: logger mode and info
    """

    # Create logger with "spam_application"
    logger = logging.getLogger(module)

    log = logging.StreamHandler()

    # Select the mode

    if mode is "Debug":
        log.setLevel(level=logging.DEBUG)

    elif mode is "Production":
        log.setLevel(level=logging.INFO)

    # read database
    logger.info("Entering " + mode + " mode\n")

    return(logger)

def setMode(args):
    """
    Set mode on production or debug
    Input: arguments obtained from argparser
    Output: mode to run the program
    """

    if args.debug is True:
        mode = "Debug"
    else:
        mode = "Production"

    return (mode)

def parseArguments():
    """
    Function to parse arguments
    Input: mode of execution and files
    Output: parsed arguments from the command line
    """

    # Create argument parser class
    parser = argparse.ArgumentParser(description = "Pipeline for hereditary \
    variant analysis 'heriline'" + '/n' + '/n' + "python3.5 utils.py -fq fastq_path -b bed_file -p project_path")

    # Define arguments to parse
    parser.add_argument("--test", "-t", required = False, action = "store_true",
    help = "Run pipeline on test mode")

    parser.add_argument("--debug", "-d", required = False, action = "store_true",
    help = "Run pipeline on debug mode")

    parser.add_argument("--verbose", "-v", required = False, action = "store_true",
    help = "Run pipeline with verbosity")

    parser.add_argument("--ori_fastq", "-fq", required = False, type = str,
    help = "Set path to the original fastq directory for symlink to project")

    parser.add_argument("--ori_bed", "-b", required = False, type = str,
    help = "Set path to the original bed/manifest file for symlink to project")

    parser.add_argument("--project_path", "-p", required = False, type = str,
    help = "Set path to the project directory")

    parser.add_argument("--show_steps", "-s", required = False,
    action = "store_true", help = "Show steps of the pipeline")

    # Call for arguments
    args = parser.parse_args()

    return (args)

def load_config(args):
    """
    Load config dictionary
    """
    with open("/nfs/production2d/cmc_projects_tmp/heriline/config.json", "r") as jconfig:
        # Since we are loading from file, json.load
        config = json.load(jconfig)

    # Set specific variables for our analysis
    config["global_config"]["test"] = args.test
    config["global_config"]["verbose"] = args.verbose
    config["global_config"]["debug"] = args.debug

    # Chose directory with files of interest
    if args.test is True:
        config["global_config"]["ori_fastq"] = os.path.join(config["global_config"]["test_directory"], "ori_fastq_files/")
        config["global_config"]["ori_bed"] = os.path.join(config["global_config"]["test_directory"], "ori_bed_file/testing.bed")
        config["global_config"]["working_directory"] = config["global_config"]["test_directory"]
        config["global_config"]["project_directory"] = os.path.abspath(config["global_config"]["test_directory"])
    else:
        config["global_config"]["ori_fastq"] = os.path.abspath(args.ori_fastq)
        config["global_config"]["ori_bed"] = os.path.abspath(args.ori_bed)
        config["global_config"]["working_directory"] = args.project_path 
        config["global_config"]["project_directory"] = os.path.abspath(args.project_path)

    return(config)

def run_cmd(cmd_str, config, logger):
    """
    Run commands passed to the function by a string
    Input: command, config dictionary and logger
    Output: executed commands
    """
    # print command
    logger.info("[CMD]: {}\n".format(cmd_str))

    if config["global_config"]["verbose"] is True:
        subprocess.call(cmd_str, shell = True)
    else:
        # If not verbose, move standard output
        subprocess.call(cmd_str + " 2> /dev/null", shell = True)

    return None

def create_directories(config, logger):
    """
    Function to create structure of directories
    Input: config and logger
    Output: directories and paths set on the config dictionary
    """
    # Create directories
    logger.info("Creating directories\n")

    # Set paths on config dictionary
    if config['global_config']['test'] is False:
        config['paths']['analysis'] = os.path.join(config["global_config"]["working_directory"], 'analysis/')
        config['paths']['bed'] = os.path.join(config["global_config"]["working_directory"], 'bed/')
        config['paths']['fastq_files'] = os.path.join(config["global_config"]["working_directory"], 'fastq_files/')
    else:
        config['paths']['analysis'] = os.path.join(config['global_config']['test_directory'], 'analysis/')
        config['paths']['bed'] = os.path.join(config['global_config']['test_directory'], 'bed/')
        config['paths']['fastq_files'] = os.path.join(config['global_config']['test_directory'], 'fastq_files/')

    for dire in ['plots', 'coverage', 'fastq', 'heritage', 'mapping', 'mapping_QC',
    'fastq_QC', 'variant_calling', 'variant_annotation', 'cmd_files', 'fof_files']:
        config['paths'][dire] = os.path.join(config['paths']['analysis'], dire + '/')

    for dire in ['fastq_merged', 'fastq_trimmed']:
        config['paths'][dire] = os.path.join(config['paths']['fastq'], dire + '/')

    for dire in ['fastq_QC_raw', 'fastq_QC_merged', 'fastq_QC_trimmed']:
        config['paths'][dire] = os.path.join(config['paths']['fastq_QC'], dire + '/')

    config['paths']['variant_calling_decomposed'] = os.path.join(config['paths']['variant_calling'],
        'variant_calling_decomposed/')

    for dire in ['mapping_mkdup', 'mapping_perchr', 'mapping_merged']:
        config['paths'][dire] = os.path.join(config['paths']['mapping'], dire + '/')

    # Get list of directories
    directories = list(config['paths'].keys())

    # Check if directories are created, if not, create them
    for d in directories:
        if os.path.isdir(config['paths'][d]):
            continue
        else:
            logger.debug('Creating directory: {}'.format(config["paths"][d]))
            cmd_mkdir = ('mkdir -p {}'.format(config["paths"][d]))
            run_cmd(cmd_mkdir, config, logger)

    logger.info('Creating symlink to original bed file\n')

    # Symlink the original bed file to our directory of interest
    cmd_symlink = ('ln -s {} {}'.format(config["global_config"]["ori_bed"],
        config["paths"]["bed"]))
    run_cmd(cmd_symlink, config, logger)


    logger.info("Creating symlink to original fastq files using 'symlink_fastq.sh script from heriline'\n")

    # Symlink the original fastq files to the project
    cmd_fqsymlink = ("bash symlink_fastq.sh -i {} -p {}".format(config["global_config"]["ori_fastq"], config["global_config"]["working_directory"]))

    run_cmd(cmd_fqsymlink, config, logger)

    return(config)

def check_args(args, config, logger):
    '''
    Check that arguments have the format of interest
    Input: parsed arguments
    Output: check
    '''

    if os.path.isdir(config["global_config"]["ori_fastq"]) is False:
        logger.error('Input directory "{}" does not exist'.format(config["global_config"]["ori_fastq"]))
        exit(1)

    elif (config['paths']['bed']+ "*.bed").endswith('.bed') is False:
        logger.error('Input file "{}" is not a .bed file'.format(config['paths']['bed']))
        exit(1)

    elif os.path.isfile(config['global_config']['reference_genome']) is False:
        logger.error('Reference genome not found in {}'.format(config['global_config']['reference_genome']))
        exit(1)

    #elif not any(fastqs.endswith('.fastq.gz') for fastqs in os.listdir(args.ori_fastq)):
        #logger.error('The original fastq directory contains files with unexpected format, expecting ".fastq.gz"')
        #exit(1)

def get_software_versions(config, logger):
    '''
    Function to get current software versions and store them in a file to keep documentation
    Input: config dictionary and logger
    Output: software_versions file with information of interest
    '''

    versions = datetime.datetime.now().strftime("software_versions_%Y%m%d_%H_%M_%S.txt")

    # Check current software and store it in a file of interest
    cmd = "ls -l /software/ | grep current > {}{}".format(config["global_config"]["working_directory"], versions)

    subprocess.call(cmd, shell = True)

    logger.info('Software versions: {}{}'.format(config["global_config"]["working_directory"], versions))


def export_config(config, logger):
    '''
    Function to export the config dictionary with paths to json file
    Input: config dictionary and logger
    Output: config file with metadata from the study
    '''
    result = config["global_config"]["working_directory"] + datetime.datetime.now().strftime("config_%Y%m%d_%H-%M-%S.json")

    with open(result, 'w') as res:
        json.dump(config, res)

    logger.info('Exporting config log: ' + result)

def heriline_steps(config, logger):
    '''
    Function to show the steps of the pipeline and which ones have been ran
    '''

    print("Heriline steps: \n \
        step1_merge_fastq\n \
        step2_trimm_fastq\n \
        step3_map\n \
        step4_mark_duplicates\n \
        step5_split_chr\n \
        step6_merge_bam\n \
        step7_call_variants\n")


def main():
    # When running only utils, processes.

    args = parseArguments()

    logger = setLogger(setMode(args), "__main__")

    config = load_config(args)

    check_args(args, config, logger)

    if args.show_steps is True:
        heriline_steps(config, logger)

    config = create_directories(config, logger)

    get_software_versions(config, logger)

    # This must be the last thing to do
    export_config(config, logger)

    pprint.pprint(config)


if __name__ == "__main__":

    main()
