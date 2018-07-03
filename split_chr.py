#! usr/bin/python3.6

# 2018-04-09, Carolina Monzo

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
    parser = argparse.ArgumentParser(description = "Split bam files per chromosome")

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
    Function to get bam files from individuals
    Input: path to the project of interest
    Output: fof file to do the splitting
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

def cmd_split_bam(config, fof):
    '''
    Create cmd file
    '''

    cmd_sh = datetime.datetime.now().strftime("cmd_split_bam_%Y%m%d_%H-%M-%S.sh")

    chr_list = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"]

    with open(config["paths"]["cmd_files"] + cmd_sh, "a") as cmd_file:

        for fi in fof:
            sample_name = fi.split("/")[-1].split("_")[0]
            
            for i in chr_list:
                bam_name = sample_name + "_chr" + i + ".bam"
                cmd_str = "samtools view -h {} {} | samtools sort -O bam -o {}{}".format(fi, i, config["paths"]["mapping_perchr"], bam_name)

                cmd_file.write(cmd_str + "\n")

    print("[INFO]: CMD_file - {}{}".format(config["paths"]["cmd_files"], cmd_sh))
    return(cmd_sh)

def run_parallel(config, cmd_sh):
    '''
    Function to run the cmd file in parallel
    Input: cmd file and config
    Output: splitted bam files
    '''

    log_str = datetime.datetime.now().strftime("split_bam_%Y%m%d_%H-%M-%S.log")

    cmd = "parallel --joblog {}{} -j 18 :::: {}{}".format(config["paths"]["mapping_perchr"], log_str, config["paths"]["cmd_files"], cmd_sh)

    print("[CMD] " + cmd)

    subprocess.call(cmd + " 2> /dev/null", shell = True)

def write_output_fof(config):
    '''

    '''

    fof = datetime.datetime.now().strftime("splitted_bams_%Y%m%d_%H-%M-%S.fof")

    # Create fof file
    cmd_fof = 'find {} -name "*chr*" > {}{}'.format(config["paths"]["mapping_perchr"], config["paths"]["fof_files"], fof)

    # Write command on the command line

    print("[CMD] " + cmd_fof)
    print("[INFO]: FOF_FILE - {}{}".format(config["paths"]["fof_files"], fof))

    # Execute command
    subprocess.call(cmd_fof, shell = True)


def main():

    args = parseArguments()

    config = get_config(args.project_path)

    fof = read_input_fof(config)

    cmd_sh = cmd_split_bam(config, fof)

    run_parallel(config, cmd_sh)

    write_output_fof(config)

if __name__ == "__main__":

    main()
