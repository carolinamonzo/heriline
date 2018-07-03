#! usr/bin/python3.6

# 2017-12-31, Carolina Monzo

# Modified for exomes-imprinting 2018-03-12

import os
import re   
import pandas as pd
import argparse
import glob
import json
import datetime
import coloredlogs
import logging
import subprocess

# GENERAL VARIABLES
_INPUT_FOF = "fastq_files_ori_{}.fof"
_CMD_FILE = "cmd_zcat_fastq_{}.sh"
_PARALLEL_LOG = "merge_fastq_{}.log"
_OUTPUT_FOF = "merged_fastq_{}.fof"

def parseArguments():
    '''
    Function to parse arguments
    Input: path to project
    Output: parsed arguments from command line
    '''

    # Create argument parser class
    parser = argparse.ArgumentParser(description = "Merge fastq files")

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


def fastq_dataframe(config, fastq_files):
    '''
    Function to take lines from fof file and create a sorted pandas dataframe with info from fastqs
    Input: list of fastq files read from the fof file
    Output: sorted dataframe with the fastq name metadata
    '''
    # Regular expression to create an annonymous dictionary
    dict_list = []
    rex = re.compile(r"{}(?P<sample>\d+)_S(?P<sh>\d+)_L(?P<lane>\d+)_R(?P<read>\d+).+".format(config["paths"]["fastq_files"]))
    for fastq in fastq_files:
        m = rex.match(fastq)
        dicc = m.groupdict()
        dicc['fastq_file'] = fastq.split('/')[-1]
        dicc['fastq_path'] = str(fastq)
        dict_list.append(dicc)

    # Read fastq info array into a pandas dataframe
    df_fastq = pd.DataFrame(dict_list, columns = ['sample', 'sh', 'lane', 'read', 'fastq_file', 'fastq_path'])
    # Sort by sample, lane and read
    df_fastq_sorted = df_fastq.sort_values(by=["fastq_path"], ascending=[True]).reset_index(drop=True)

    return(df_fastq_sorted)


def cmd_zcat_fastq(config, cmd_file, df_fastq_sorted):
    '''
    Function to generate the zcat commands, it selects R1 and R2 from sorted samples
    Input: dataframe of fastq sorted names
    Output: cmd file with the zcat commands of ordered fastqs to merge
    '''
    cmd_time = datetime.datetime.now().strftime("%Y%m%d_%H-%M-%S")

    cmd_sh = cmd_file.format(cmd_time)

    # Create list by sample
    samples = list(df_fastq_sorted['sample'].unique())
    
    # Set zcat template

    for i in range(len(samples)):
            
            sample_fastqs = list(df_fastq_sorted[df_fastq_sorted['sample']==samples[i]]['fastq_path'])
            R1_list = sample_fastqs[0::2]
            R2_list = sample_fastqs[1::2]
            merged_fastq_name = df_fastq_sorted['sample'].unique()[i]
            R1_str = "zcat {} | gzip > {}{}_R1_merged.fastq.gz".format(' '.join(R1_list), config["paths"]["fastq_merged"], merged_fastq_name)
            R2_str = "zcat {} | gzip > {}{}_R2_merged.fastq.gz".format(' '.join(R2_list), config["paths"]["fastq_merged"], merged_fastq_name)
                                    
            # Write the cmd.sh file
            with open('{}{}'.format(config["paths"]["cmd_files"], cmd_sh), 'a') as cmd_file:
                cmd_file.write(R1_str + '\n')
                cmd_file.write(R2_str + '\n')
                
    print("[INFO]: CMD_FILE - {}{}".format(config["paths"]["cmd_files"], cmd_sh))

    return(cmd_sh)

def run_parallel(config, parallel_log, cmd_sh):
    '''
    Function to run the cmd file in parallel
    Input: cmd file and project path to find it
    Output: merged fastq files
    '''
    log_time = datetime.datetime.now().strftime("%Y%m%d_%H-%M-%S")

    log_str = parallel_log.format(log_time)

    cmd = "parallel --joblog {}{} -j 15 :::: {}{}".format(config["paths"]["fastq_merged"], log_str, config["paths"]["cmd_files"], cmd_sh)

    print("[CMD]: " + cmd)

    subprocess.call(cmd + " 2> /dev/null", shell = True)


def write_output_fof(config, output_fof):
    """
    Create fof file for merged fastq files
    """
    # Create fof file for merged fastq files

    fof_time = datetime.datetime.now().strftime("%Y%m%d_%H-%M-%S")

    fof = output_fof.format(fof_time)

    cmd_str = 'find {} -name "*.fastq.gz" > {}{}'.format(config["paths"]["fastq_merged"],config["paths"]["fof_files"], fof)
    print("[CMD]: " + cmd_str)
    print("[INFO]: FOF_FILE - {}{}".format(config["paths"]["fof_files"], fof))
    subprocess.call(cmd_str, shell = True)

def main():
    '''
    Function to sort the orders
    '''

    args = parseArguments()

    config = get_config(args.project_path)

    fof = read_input_fof(config, _INPUT_FOF)

    df_fastq_sorted = fastq_dataframe(config, fof)
    
    cmd_sh = cmd_zcat_fastq(config, _CMD_FILE, df_fastq_sorted)

    run_parallel(config, _PARALLEL_LOG, cmd_sh)

    write_output_fof(config, _OUTPUT_FOF)

if __name__ == '__main__':
    main()
