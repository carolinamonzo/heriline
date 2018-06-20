#! /bin/bash

# Carolina Monzo - 2018-03-12

# Script: Symlink exome files to project

# Get path of original files

while getopts hp:i: option
do
    case "${option}"
    in
        h) help=${OPTARG}
        echo "Script to symlink files to project"
        echo " -i: path to directory of original files"
        echo " -p: path to project path"
        exit 1
        ;;
        i) ori_path=${OPTARG}
        ;;
        p) project_path=${OPTARG}
        ;;
    esac
done

# Set path variables

fastq_directory=${project_path}fastq_files/

analysis_directory=${project_path}analysis/

dt=`date '+%Y%m%d_%H-%M-%S'`

# Find files of interest (not empty) and create symlinks
# Write commands on cmd file as well

cmd=${analysis_directory}cmd_files/cmd_symlink_ori_fastq_${dt}.sh

for file in $(find $ori_path -name "*.fastq.gz"); do
    if [[ -s $file ]]; then 
    echo "ln -s $file $fastq_directory" >> $cmd;
    fi;
done

echo "[INFO]: CMD_FILE - $cmd"

bash $cmd

# Create fof with this symlinks

fof=${analysis_directory}fof_files/fastq_files_ori_${dt}.fof

# create fof file
cmd_fof="""find ${fastq_directory} -name "*.fastq.gz" > ${fof}"""

echo "[CMD]: $cmd_fof"
echo "[INFO]: FOF_FILE - $fof"
eval $cmd_fof
