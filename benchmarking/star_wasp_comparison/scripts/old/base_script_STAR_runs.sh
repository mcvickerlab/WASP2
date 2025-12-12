# Running STAR with no variants-no WASP output (/STAR_Base). Sample: sample_id | xthreads   
# This bash file contains commands run for each sample. Base edits to STAR commands can be made to this file: "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR/BaseCode_STAR_Runs.sh". To apply these changes to \
# each sample and run, re-call "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR/STAR_Run_Dirs.py" to re-create all sample directories and copy BaseCode files and/edits to respective sample and thread directories. The "STAR_Run_Dirs.py" 
# file will also pass respective parameters such as the sample name and number of threads to sample-specific BaseCode files. 
# Specifically, parameter replacements (done through the python file (STAR_Run_Dirs.py)) per run are made to the SAMPLE_ID, XTHREADS and "RUNTHREADN X"
# Directory "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR/STAR_Runs" which receives all run output per sample needs to be created prior to running the python file: "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR/STAR_Run_Dirs.py" 


STAR=/usr/bin/STAR
genomeDirectory="/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/genome_directory/"
fastqFileDir="/scratch/asiimwe/STAR-WASP_FASTQs_VCFs/FASTQ"
ulimit -n 10000

STARdir=/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR/STAR_Runs/sample_id/xthreads

cd $STARdir

STARpar="--runThreadN x --genomeDir $genomeDirectory --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM MD jM jI rB MC  --alignEndsType EndToEnd --outSAMunmapped Within --outFilterMultimapNmax 1"

readFiles="--readFilesCommand gunzip -c --readFilesIn $fastqFileDir/sample_id/R1.fastq.gz  $fastqFileDir/sample_id/R2.fastq.gz"

$STAR $STARpar $readFiles
