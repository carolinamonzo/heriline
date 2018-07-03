#! usr/bin/python3.6

# 2018-06-06, Carolina Monzo

import os
import pandas as pd
import re
import argparse
import json
import datetime
import glob
import subprocess

def parseArguments():
    '''
    Function to parse arguments
    Input: path to project
    Output: parsed arguments from command line
    '''

    # Create argument parser class
    parser = argparse.ArgumentParser(description = "Mark duplicates")

    # Define arguments to parse
    parser.add_argument("--project_path", "-p", required = False, type = str, help = "Argument to set the path to the project")

    # Call for arguments
    args = parser.parse_args()

    return(args)

def get_config(project_path):
    """
    Load config dictionary
    """

    config_list = []

    os.chdir(project_path)
    for file in glob.glob("config_*.json"):
        config_list.append(file)

    config_pwd = project_path + sorted(config_list)[-1]

    with open(config_pwd, "r") as jconfig:
        config = json.load(jconfig)

    return(config)

def read_input_fof(config):
    '''
    Function to get bam files with marked duplicates
    Input: path to the project of interest
    Output: fof file to do the mapping
    '''
    # Get fof file

    fof_list = []

    for file in glob.glob(config["paths"]["fof_files"] + "marked_duplicates_bam_*.fof"):
        fof_list.append(file)

    fof_pwd = os.path.normpath(sorted(fof_list)[-1])

    # Read fof file

    with open(fof_pwd, "r") as fi:
        fof = fi.read().splitlines()
    return(fof)

def cmd_mapping_statistics(config, fof):
    '''

    '''

    cmd_sh = datetime.datetime.now().strftime("cmd_mapping_statistics_%Y%m%d_%H-%M-%S.sh")

    with open(config["paths"]["cmd_files"] + cmd_sh, "a") as cmd_file:

        for fi in fof:
            sample_name = fi.split("/")[-1].split("_")[0]
            cmd_str = "samtools flagstat {} --threads 5 > {}/{}_flagstat.txt".format(fi, config["paths"]["mapping_QC"], sample_name)

            cmd_file.write(cmd_str + '\n')

        for fi in fof:
            sample_name = fi.split("/")[-1].split("_")[0]
            cmd_str = "java -jar /software/bin/picard.jar CollectAlignmentSummaryMetrics I={} O={}/{}_AlignmentSummaryMetrics.txt R={} ASSUME_SORTED=TRUE".format(fi, config["paths"]["mapping_QC"], sample_name, config["global_config"]["reference_genome"])

            cmd_file.write(cmd_str + '\n')



    print('[INFO]: CMD_FILE - {}{}'.format(config['paths']['cmd_files'], cmd_sh))

    return(cmd_sh)


def run_parallel(config, cmd_sh):
    '''

    '''

    log_str = datetime.datetime.now().strftime("mapping_statistics_%Y%m%d_%H-%M-%S.log")

    cmd = "parallel --joblog {}{} -j 4 :::: {}{}".format(config["paths"]["mapping_QC"], log_str, config["paths"]["cmd_files"], cmd_sh)

    print("[CMD]: " + cmd)

    subprocess.call(cmd + " 2> /dev/null", shell = True)

def write_output_fof(config):
    '''

    '''

    fof1 = datetime.datetime.now().strftime("mapping_statistics_flagstat_%Y%m%d_%H-%M-%S.fof")

    cmd_fof1 = 'find {} -name "*_flagstat.txt" > {}/{}'.format(config["paths"]["mapping_QC"], config["paths"]["fof_files"], fof1)

    print("[CMD]: " + cmd_fof1)
    print("[INFO]: FOF_FILE - {}{}".format(config["paths"]["fof_files"], fof1))

    # Execute command for fof 1
    os.system(cmd_fof1)



    #<TODO> We dont know yet the extensions of the other metrics files

def main():
    '''

    '''

    args = parseArguments()

    config = get_config(args.project_path)

    fof = read_input_fof(config)

    cmd_sh = cmd_mapping_statistics(config, fof)

    run_parallel(config, cmd_sh)

    write_output_fof(config)


if __name__ == '__main__':
    main()

