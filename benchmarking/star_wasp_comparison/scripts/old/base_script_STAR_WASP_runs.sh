# Running STAR+WASP. Sample: sample_id | xthreads   
# This bash file contains commands run for each sample. Base edits to STAR+WASP commands can be made to this file: "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/BaseCode_STAR_WASP_Runs.sh". To apply these changes to \
# each sample and run, re-call "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/STAR_WASP_Run_Dirs.py" to re-create all sample directories and copy BaseCode files and/edits to respective sample and thread directories.
# The "STAR_WASP_Run_Dirs.py" file will also pass respective parameters such as the sample name and number of threads to sample-specific BaseCode files. 
# Specifically, parameter replacements (done through the python file (STAR_WASP_Run_Dirs.py)) per run are made to the SAMPLE_ID, XTHREADS, "RUNTHREADN X" and "SSD" which represents the unique snp directory per sample. (We have multiple\
# samples derived from sample NA12878. All NA12878-derived samples need to point to one snp directory to avoid mutiple copies of the same directory - this variable caters for the sample uniqueness we need).
# Directory "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/STAR_WASP_Runs" which receives all run output per sample needs to be created prior to running the python file:"/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/STAR_WASP_Run_Dirs.py"

STAR=/usr/bin/STAR
genomeDirectory="/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/genome_directory/"
fastqFileDir="/scratch/asiimwe/STAR-WASP_FASTQs_VCFs/FASTQ"
vcfFileDir="/scratch/asiimwe/STAR-WASP_FASTQs_VCFs/VCF"
ulimit -n 10000

STAR_WASPdir=/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/STAR_WASP_Runs/sample_id/xthreads

cd $STAR_WASPdir

STARpar="--runThreadN x --genomeDir $genomeDirectory  --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM MD jM jI rB MC vA vG vW --waspOutputMode SAMtag  --alignEndsType EndToEnd --outSAMunmapped Within --outFilterMultimapNmax 1"

readFiles="--readFilesCommand gunzip -c --readFilesIn $fastqFileDir/sample_id/R1.fastq.gz  $fastqFileDir/sample_id/R2.fastq.gz"

hetVcf="--varVCFfile $vcfFileDir/ssd_input_snp_dir/ssd.vcf.snv1het"


$STAR $STARpar $readFiles $hetVcf

