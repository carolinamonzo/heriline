#! usr/bin/python3.5

# 2018-06-20, Carolina Monzo

# Script: calculating fastqc stats and creating a small 
# report with general statistics

import os
import json
import argparse
import glob
import datetime
import subprocess

# GENERAL VARIABLES
merged_fastq = "merged_fastq_{}.fof"
trimmed_fastq = "merged_trimmed_fastq_files_{}.fof"
cmd_file = "cmd_fastq_QC_{}.sh"
parallel_log = "fastq_QC_{}.log"
output_fof = "merged_and_trimmed_fastqQC_{}.fof"

def parseArguments():
    '''
    Function to parse arguments
    Input: path to project
    Output: parsed arguments from command line
    '''

    # Create argument parser class
    parser = argparse.ArgumentParser(description = "Calculate fastq QC statistics")

    # Define arguments to parse
    parser.add_argument("--project_path", "-p", required = False, type = str, help = "Argument to set the path to the project")

    # Call for arguments
    args = parser.parse_args()

    return(args)

def get_config(project_path):
    '''
    Function to get the last config prepared for the analysis
    Input: path to project
    Output: config dictionary with paths for the project
    '''

    config_list = []

    for file in glob.glob(project_path + "config_*.json"):
        config_list.append(file)

    config_pwd = project_path + sorted(config_list)[-1]

    with open(config_pwd, 'r') as jconfig:
        config = json.load(jconfig)


    return(config)

def read_input_fof(config, input_fof):
    '''
    Function to create fastq files fof
    Input: path to the project of interest
    Output: fof file to do the merging
    '''
    # Get fof file

    fof_list = []

    for file in glob.glob(config["paths"]["fof_files"] + input_fof.format("*")):
        fof_list.append(file)

    fof_pwd = os.path.normpath(sorted(fof_list)[-1])

    with open(fof_pwd, "r") as fi:
        fof = fi.read().splitlines()

    return(fof)

def cmd_fastqQC(config, merged_fof, trimmed_fof):
    """
    Create cmd file to run fastqc analysis on all fastq files
    """
    # Concatenate raw_fof and trimmed_fof to run fastqc analysis on all files

    fof = merged_fof + trimmed_fof

    # Create cmd file with commands for running fastqc analysis

    cmd_time = datetime.datetime.now().strftime("%Y%m%d_%H-%M-%S")

    cmd_sh = cmd_file.format(cmd_time)

    for fastq in fof:
        samplename = fastq.split('/')[-1].split('.')[0]

        # If files are trimmed, send output to trimmed QC directory
        if "trimmed" in samplename:
            cmd_str = "fastqc {}{} -o {} --extract".format(config["paths"]["fastq_trimmed"], samplename, config["paths"]["fastq_QC_trimmed"])

        # If files are raw merged, send output to merged QC directory
        else:
            cmd_str = "fastqc {}{} -o {} --extract".format(config["paths"]["fastq_merged"], samplename, config["paths"]["fastq_QC_merged"])


        # Write the cmd.sh file
        with open("{}{}".format(config["paths"]["cmd_files"], cmd_sh), "a") as cmd_file:
            cmd_file.write(cmd_str + '\n')

    print("[INFO]: CMD_FILE - {}{}".format(config["paths"]["cmd_files"], cmd_sh))

    return(cmd_sh)


def run_parallel(config, cmd_sh):
    """

    """

    log_time = datetime.datetime.now().strftime("%Y%m%d_%H-%M-%S")

    log_str = parallel_log.format(log_time)

    cmd = "parallel --joblog {}{} -j 10 :::: {}{}".format(config["paths"]["fastq_QC"], log_str, config["paths"]["cmd_files"], cmd_sh)

    print("[CMD]: " + cmd)

    subprocess.call(cmd, shell = True)

def write_output_fof(config):
    """

    """

    # Create fof file for fastqc.html reports

    fof_time = datetime.datetime.now().strftime("%Y%m%d_%H-%M-%S")

    fof = output_fof.format(fof_time)

    cmd_str = 'find {} -name "*_fastqc.html" > {}{}'.format(config["paths"]["fastq_QC"], config["paths"]["fof_files"], fof)

    print("[CMD]: " + cmd_str)
    print("[INFO]: FOF_FILE - {}{}".format(config["paths"]["fof_files"], fof))
    subprocess.call(cmd_str, shell = True)

def clean_directory(config):
    """
    Fastqc creates many files we dont have interest in
        - Zip file with all files of the analysis
        - Inside each directory of files:fastqc_report.html, fastqc.fo, Icons/, Images
    """

    cmd_str = 'for file in $(find {} -name "fastqc_report.html" -o -name "fastqc.fo" -o -name "*.zip"); do rm $file; done'.format(config["paths"]["fastq_QC"])
    subprocess.call(cmd_str, shell = True)

    cmd_str = 'for dir in $(find {} -type d -name "Icons" -o -name "Images"); do rm -r $dir; done'.format(config["paths"]["fastq_QC"])

    subprocess.call(cmd_str, shell = True)


def main():


    args = parseArguments()

    config = get_config(args.project_path)

    merged_fof = read_input_fof(config, merged_fastq)

    trimmed_fof = read_input_fof(config, trimmed_fastq)

    cmd_sh = cmd_fastqQC(config, merged_fof, trimmed_fof)

    run_parallel(config, cmd_sh)

    write_output_fof(config)

    # Since fastqc creates extra files that we dont want we have to remove them
    clean_directory(config)



if __name__ == "__main__":
    main()

