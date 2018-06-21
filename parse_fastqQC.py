#! usr/bin/python3.5

# 2018-06-20; Carolina Monzo

# Script to parse fastqc results and get general statistics

import os
import json
import argparse
import glob
import datetime
import subprocess
import re
import pandas as pd

# Dictionary of modules of fastqc
_HEADER = ["Basic_statistics",
"Per_base_sequence_quality",
"Per_tile_sequence_quality",
"Per_sequence_quality_scores",
"Per_base_sequence_content",
"Per_sequence_GC_content",
"Per_base_N_content",
"Sequence_Length_Distribution",
"Sequence_Duplication_Levels",
"Overrepresented_sequences",
"Adapter_Content",
"Kmer_Content",
"File_name"
]


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

def get_stat_files(config, path):
    """
    Function to get a list of the files to parse
    Input: config dictionary to obtain path to project
    Output: list of files to parse
    """
    stat_files = []
    
    # Walk the directories and get a list of paths and files to parse
    for root, dirs, f in os.walk(path):
        for fi in f:
            if fi == "summary.txt":
                stat_files.append(os.path.join(root, fi))

    return(stat_files)

def get_stats(config, files):
    """
    Function to read summary.txt files and obtain metrics
    Input: list of files
    Output: pandas dataframe with values and general file with stats per sample
    """

    # Create dataframe

    stats = []

    stats_todf = []

    st_time = datetime.datetime.now().strftime("%Y%m%d_%H-%M-%S")

    stats_file = files[0].split("/")[-2].split("_")[-2] + "_general_stats_{}.tsv".format(st_time)


    # Get stats per file
    for f in files:
        with open(f, "r") as fi:
            for line in fi:
                stats.append(line.split("\t")[0])
        stats.append(f.split("/")[-2].split("_merged")[0])

        stats_todf.append(stats)


        # Clean list
        stats = []      


    df = pd.DataFrame(stats_todf, columns = _HEADER)

    # Write general statistic files to a csv file to have all stats per file
    df.to_csv("{}{}".format(config["paths"]["fastq_QC"], stats_file), sep="\t", index=False)


    print("[INFO]: STATS_FILE - {}{}".format(config["paths"]["fastq_QC"], stats_file))

    return(df)


def summary_stats(config, df, fi):
    """
    Function to summarize statistics files per statistic
    Input: statistics per sample dataframe and file name to write to
    Output: temporal csv file with statistics
    """

    melted_data = pd.melt(df, value_vars = _HEADER[:-1], var_name = "statistics", value_name = "summary")

    ab = melted_data.groupby(by=["statistics", "summary"])["summary"].count()    
    ab.groupby(by=["statistics", "summary"]).size().reset_index(name="count")

    ab.to_csv(fi, mode = "a", sep = '\t')
    
    
def reshape_summary_file(config, df_init):
    """
    Function to reshape data and format the output csv file
    Input: dataframe to reshape
    Output: formated dataframe to print to final file
    """
    str_file = "temp_{}.txt".format(datetime.datetime.now().strftime("%Y%m%d_%H-%M-%S"))

    header = ["statistics", "status", "count"]
    df_list = []

    # We dont have a dataframe, we have a series object that stops modifications
    with open("{}{}".format(config["paths"]["fastq_QC"], str_file), "a") as fi:
        summary_stats(config, df_init, fi)

    with open("{}{}".format(config["paths"]["fastq_QC"], str_file), "r") as ca:
        lines = ca.read().splitlines()

    os.system("rm {}{}".format(config["paths"]["fastq_QC"], str_file))

    for line in lines:
        df_list.append(line.split("\t"))

    # Create a dataframe object with the data so it can be manipulated
    df = pd.DataFrame(df_list, columns = header)

    # Shape the data adequately
    df = df.pivot(index = "statistics", columns = "status", values = "count").fillna(0).reset_index()

    df = df[["PASS", "WARN", "FAIL", "statistics"]]

    return(df)


def main():
    """
    
    """

    args = parseArguments()

    config = get_config(args.project_path)

    # Get lists of files to parse

    stats_merged = get_stat_files(config, config["paths"]["fastq_QC_merged"])

    stats_trimmed = get_stat_files(config, config["paths"]["fastq_QC_trimmed"])

    df_merged = get_stats(config, stats_merged)

    df_trimmed = get_stats(config, stats_trimmed)

    # Reparse to get a simpler file with both results
    str_time = datetime.datetime.now().strftime("%Y%m%d_%H-%M-%S")
    str_file = "summary_test_{}.tsv".format(str_time)

    df1 = reshape_summary_file(config, df_merged)
    df2 = reshape_summary_file(config, df_trimmed)

    with open("{}{}".format(config["paths"]["fastq_QC"], str_file), "a") as fi:
        fi.write("STATISTICS_FASTQ_MERGED\n")
        df1.to_csv(fi, mode = "a", sep = '\t', index = False)
        fi.write("\n\nSTATISTICS_FASTQ_MERGED_TRIMMED\n")
        df2.to_csv(fi, mode = "a", sep = "\t", index = False)

    print("[INFO]: SUMMARY_STATS_FILE - {}{}".format(config["paths"]["fastq_QC"], str_file))

    

if __name__ == '__main__':
    main()

