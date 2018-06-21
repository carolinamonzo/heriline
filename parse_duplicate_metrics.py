#! usr/bin/python3.5

# 2018-06-21; Carolina Monzo

# Script to parse mapping qc results and get general statistics

import os
import json
import argparse
import glob
import datetime
import subprocess
import re
import pandas as pd


def parseArguments():
    '''
    Function to parse arguments
    Input: path to project
    Output: parsed arguments from command line
    '''

    # Create argument parser class
    parser = argparse.ArgumentParser(description = "Parse mapping QC files")

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

def get_stat_files(config, str_file):
    """
    Function to get a list of the files to parse
    Input: config dictionary to obtain path to project
    Output: list of files to parse
    """
    stat_files = []

    # get a list of paths and files to parse
    for f in glob.glob(config["paths"]["mapping_QC"] + str_file):
        stat_files.append(f)

    return(stat_files)


def get_header(config, stat_files):

    header = ["STATISTIC"]
    for f in stat_files:
        header.append(f.split("/")[-1].split("_")[0])

    return("\t".join(header))

def paste_files_duplicate_metrics(config, stat_files):

 
    cmd_str = ["paste -d '\\t'", "<(cat {} | head -8 | tail -2 | sed 's/ /_/g' | bash /nfs/production2d/cmc_projects_tmp/transpose_first_two_lines.sh | cut -f 1)".format(stat_files[0])]
    for f in stat_files:

        string_files = "<(cat {} | head -8 | tail -2 | sed 's/ /_/g' | bash /nfs/production2d/cmc_projects_tmp/transpose_first_two_lines.sh | cut -f 2)".format(f)

        cmd_str.append(string_files)

    return(cmd_str)

def paste_files_AlignmentSummary(config, stat_files):

    cmd_str = ["paste -d '\\t'", "<(grep -v -E '(^#|^$)' {} | bash /nfs/production2d/cmc_projects_tmp/transpose_first_two_lines.sh | cut -f 1)".format(stat_files[0])]

    for f in stat_files:
        string_files = "<(grep -v -E '(^#|^$)' {} | bash /nfs/production2d/cmc_projects_tmp/transpose_first_two_lines.sh | cut -f 2)".format(f)

        cmd_str.append(string_files)

    return(cmd_str)



def write_output_file(config, header, cmd_fish, str_f):

    str_time = datetime.datetime.now().strftime("%Y%m%d_%H-%M-%S")

    str_file = str_f.format(str_time)

    cmd_fish.append(">> {}{}".format(config["paths"]["mapping_QC"], str_file))
    
    cmd_str = " ".join(cmd_fish)

    with open("{}{}".format(config["paths"]["mapping_QC"], str_file), "a") as fi:
        fi.write(header + "\n")
    
    subprocess.call(["bash", "-c", cmd_str])

    print("[INFO]: STATS_FILE - {}{}".format(config["paths"]["mapping_QC"], str_file))


def main():
    args = parseArguments()

    config = get_config(args.project_path)

    # Get stats for duplicate_metrics from picard
    str_f = "*_duplicate_metrics.txt"

    stat_files = get_stat_files(config, str_f)

    header = get_header(config, stat_files)    

    cmd_fish = paste_files_duplicate_metrics(config, stat_files)
    
    str_f = "duplicate_metrics_stats_{}.tsv"
    write_output_file(config, header, cmd_fish, str_f)

    # Get stats for AlignmentSummaryMetrics
    str_f = "*_AlignmentSummaryMetrics.txt"

    stat_files = get_stat_files(config, str_f)

    header = get_header(config, stat_files)

    cmd_fish = paste_files_AlignmentSummary(config, stat_files)

    str_f = "AlignmentSummaryMetrics_stats_{}.tsv"
    write_output_file(config, header, cmd_fish, str_f)

if __name__ == '__main__':
    main()
