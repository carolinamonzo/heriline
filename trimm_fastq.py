#! usr/bin/python3.6

# 2018-03-14, Carolina Monzo

# Script: trimming merged fastq files by quality

import os
import json
import argparse
import glob
import datetime
import subprocess

def parseArguments():
    '''
    Function to parse arguments
    Input: path to project
    Output: parsed arguments from command line
    '''

    # Create argument parser class
    parser = argparse.ArgumentParser(description = "Trimm fastq files")

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

    for file in glob.glob(project_path + "config_*.json"):
        config_list.append(file)

    config_pwd = os.path.normpath(project_path + sorted(config_list)[-1])

    with open(config_pwd, "r") as jconfig:
        config = json.load(jconfig)

    return(config)

def read_input_fof(config):
    """
    Get fof file for merged fastq files
    """
    # Get fof file

    fof_list = []

    for file in glob.glob(config["paths"]["fof_files"] + "merged_fastq_*.fof"):
        fof_list.append(file)

    fof_pwd = os.path.normpath(sorted(fof_list)[-1])

    with open(fof_pwd, "r") as fi:
        fof = fi.read().splitlines()
    return(fof)


def cmd_trimm_fastq(config, fof):
    """
    Create cmd file
    """
    # Create cmd file with commands for trimming fastq files

    cmd_sh = datetime.datetime.now().strftime("cmd_trimm_fastq_%Y%m%d_%H-%M-%S.sh")

    for fastq in fof:
        samplename = fastq.split('/')[-1].split('.')[0]
        cmd_str = "seqtk trimfq {} | gzip > {}/{}-trimmed.fastq.gz".format(fastq, config["paths"]["fastq_trimmed"], samplename)
        
        # Write the cmd.sh file
        with open("{}{}".format(config["paths"]["cmd_files"], cmd_sh), "a") as cmd_file:
            cmd_file.write(cmd_str + '\n')

    print("[INFO]: CMD_FILE - {}{}".format(config["paths"]["cmd_files"], cmd_sh))

    return(cmd_sh)


def run_parallel(config, cmd_sh):
    """
    Run the cmd_trimm_fastq.sh in parallel
    """
    log_str = datetime.datetime.now().strftime("trimming_fastq_%Y%m%d_%H-%M-%S.log")

    cmd_str = "parallel --joblog {}{} -j15 :::: {}{}".format(config["paths"]["fastq_trimmed"], log_str, config["paths"]["cmd_files"], cmd_sh)

    print ("[CMD]: " + cmd_str)
    subprocess.call(cmd_str + " 2> /dev/null", shell = True)

def write_output_fof(config):
    '''
    Function to create merged fastq files fof
    Input: path to the project of interest
    Output: fof file to do the mapping
    '''

    fof = datetime.datetime.now().strftime("merged_trimmed_fastq_files_%Y%m%d_%H-%M-%S.fof")

    # Create fof file
    cmd_fof = 'find {} -name "*.fastq.gz" > {}{}'.format(config["paths"]["fastq_trimmed"], config["paths"]["fof_files"], fof)

    # Write command on the command line
    print("[CMD]: " + cmd_fof)
    print("[INFO]: FOF_FILE - {}{}".format(config["paths"]["fof_files"], fof))

    # Execute command
    subprocess.call(cmd_fof, shell = True)


def main():

    args = parseArguments()

    config = get_config(args.project_path)

    fof = read_input_fof(config)

    cmd_sh = cmd_trimm_fastq(config, fof)

    run_parallel(config, cmd_sh)

    write_output_fof(config)

if __name__ == '__main__':
    main()
