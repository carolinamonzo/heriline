#! usr/bin/python3.6

# 2018-01-03, Carolina Monzo

# Modified for exomes-imprinting 2018-03-13

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
    parser = argparse.ArgumentParser(description = "Map fastq files")

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
    Function to get merged fastq files fof
    Input: path to the project of interest
    Output: fof file to do the mapping
    '''
    # Get fof file

    fof_list = []

    for file in glob.glob(config["paths"]["fof_files"] + "merged_trimmed_fastq_files_*.fof"):
        fof_list.append(file)

    fof_pwd = os.path.normpath(sorted(fof_list)[-1])

    # Read fof file

    with open(fof_pwd, "r") as fi:
        fof = fi.read().splitlines()
    return(fof)

def merged_fastq_dataframe(config, fastq_files):
    '''
    Function to take lines from fof file to create a sorted pandas dataframe with info from merged fastqs
    Input: list of fastq files read from the fof file
    Ouptut: sorted dataframe with the fastq name metadata
    '''
    #<TODO> join this function with the fastq_dataframe function from merge_fastq.py
        #Ideally we would only create the rex variable outside and pass it as input to the function
        #Then we create an if statement and chose the corresponding columns for the dataframe

    # Regular expression to create an annonymous dictionary
    dict_list = []

    rex = re.compile(r"{}(?P<sample>\d+)_R(?P<read>\d+)_merged-trimmed.fastq.gz".format(config["paths"]["fastq_trimmed"]))

    for fastq in fastq_files:
        m = rex.match(fastq)
        dicc = m.groupdict()
        dicc['fastq_file'] = fastq.split('/')[-1]
        dicc['fastq_path'] = str(fastq)
        dict_list.append(dicc)

    # Read fastq info array into a pandas dataframe
    
    df_fastq = pd.DataFrame(dict_list, columns = ['sample', 'read', 'fastq_file', 'fastq_path'])

    # Sort by sample and read

    df_fastq_sorted = df_fastq.sort_values(by=['fastq_path'], ascending=[True]).reset_index(drop=True)

    return(df_fastq_sorted)

def cmd_map_fastq(config, df_fastq_sorted):
    '''
    Function to create commands to map merged fastq files, we have
        a R1 and R2 fastq corresponding to each sample
    Input: dataframe with metadata from merged fastq files
    Output: files with commands for mapping with BWA
    '''

    cmd_sh = datetime.datetime.now().strftime("cmd_bwa_mem_%Y%m%d_%H-%M-%S.sh")

    # Create list by sample

    samples = list(df_fastq_sorted['sample'].unique())

    # Set templates

    for i in range(len(samples)):
        sample_fastqs = list(df_fastq_sorted[df_fastq_sorted['sample']==samples[i]]['fastq_path'])
        R1_list = sample_fastqs[0]
        R2_list = sample_fastqs[1]
        output_bam_name = df_fastq_sorted['sample'].unique()[i]

        # Create template for bwa
        
        cmd_bwa = """time bwa mem -M -L 5 -t 3 {}.gz -R "@RG\\tID:{}\\tPL:ILLUMINA\\tSM:{}\\tDS:ref=37d5\\tCN:UGDG\\tDT:{}\\tPU:{}" {} {} 2> {}{}_bwa_mem.err | samtools sort -O bam -o {}{}_sorted.bam; time samtools index {}{}_sorted.bam""".format(config["global_config"]["reference_genome"], output_bam_name, output_bam_name, datetime.datetime.now().strftime("%Y-%m-%d"), output_bam_name, R1_list, R2_list, config["paths"]["mapping"], output_bam_name, config["paths"]["mapping"], output_bam_name, config["paths"]["mapping"], output_bam_name)

        # Write the cmd.sh file

        with open('{}{}'.format(config["paths"]["cmd_files"], cmd_sh), 'a') as cmd_file:
            cmd_file.write(cmd_bwa + '\n')
    print("[INFO]: CMD_FILE - {}{}".format(config["paths"]["cmd_files"], cmd_sh))

    return(cmd_sh)


def run_parallel(config, cmd_sh):
    '''
    Function to run the cmd file in parallel
    Input: cmd file and project path to find it
    Output: merged fastq files
    '''
    log_str = datetime.datetime.now().strftime("bwa_mem_%Y%m%d_%H-%M-%S.log")

    cmd = "parallel --joblog {}{} -j 5 :::: {}{}".format(config["paths"]["mapping"], log_str, config["paths"]["cmd_files"], cmd_sh)

    print("[CMD]: " + cmd)

    subprocess.call(cmd, shell = True)

def write_output_fof(config):
    '''
    Function to generate a fof file for mapped files
    '''

    fof = datetime.datetime.now().strftime("mapped_sorted_bwa_mem_%Y%m%d_%H-%M-%S.fof")

    # Create fof file

    cmd_fof = 'find {} -name "*_sorted.bam" > {}{}'.format(config["paths"]["mapping"], config["paths"]["fof_files"], fof)

    # Write command on the command line

    print("[CMD]: " + cmd_fof)
    print("[INFO]: FOF_FILE - {}{}".format(config["paths"]["fof_files"], fof))

    # Execute command
    os.system(cmd_fof)


def main():
    '''
    Function to sort the order of functions to use
    '''

    args = parseArguments()

    config = get_config(args.project_path)

    fastq_files = read_input_fof(config)

    # Create fastq files dataframe

    df_fastq_sorted = merged_fastq_dataframe(config, fastq_files)

    # Create BWA MEM commands for mapping

    cmd_sh = cmd_map_fastq(config, df_fastq_sorted)

    run_parallel(config, cmd_sh)

    write_output_fof(config)

if __name__ == '__main__':
    main()
