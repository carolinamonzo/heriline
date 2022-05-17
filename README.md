# heriline
Scripts for hereditary variants analysis, can be run independently or:

### Arguments to run heriline:
   - Path to bed file
   - Path to original fastq files
   - Path to project

```
usage: utils.py [-h] [--test] [--debug] [--verbose] [--ori_fastq ORI_FASTQ]
                [--ori_bed ORI_BED] [--project_path PROJECT_PATH]
                [--show_steps]

Pipeline for hereditary variant analysis 'heriline'  
python3.5 utils.py -fq fastq_path -b bed_file -p project_path

optional arguments:
  -h, --help            show this help message and exit
  --test, -t            Run pipeline on test mode
  --debug, -d           Run pipeline on debug mode
  --verbose, -v         Run pipeline with verbosity
  --ori_fastq ORI_FASTQ, -fq ORI_FASTQ
                        Set path to the original fastq directory for symlink
                        to project
  --ori_bed ORI_BED, -b ORI_BED
                        Set path to the original bed/manifest file for symlink
                        to project
  --project_path PROJECT_PATH, -p PROJECT_PATH
                        Set path to the project directory
  --show_steps, -s      Show steps of the pipeline
  ```
  
  
