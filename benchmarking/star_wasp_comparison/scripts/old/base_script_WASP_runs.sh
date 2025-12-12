# Running WASP (mappability filtering for correcting allelic mapping biases) with STAR alignmnet. Sample: sample_id | xthreads   
# This bash file contains commands run for each sample. Base edits to WASP commands can be made to this file: "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/WASP/BaseCode_WASP_Runs.sh". To apply these changes to \
# each sample and run, re-call "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/WASP/WASP_Run_Dirs.py" to re-create all sample directories and copy BaseCode files and/edits to respective samples and thread directories.
# The "WASP_Run_Dirs.py" file will also pass respective parameters such as the sample name and number of threads to sample-specific BaseCode files. 
# Specifically, parameter replacements (done through the python file (WASP_Run_Dirs.py)) per run are made to the SAMPLE_ID, XTHREADS, "RUNTHREADN X" and "SSD" which represents the unique snp directory per sample. (We have multiple\
# samples derived from sample NA12878. All NA12878-derived samples need to point to one snp directory to avoid mutiple copies of the same directory - this variable caters for the sample uniqueness we need).
# Directory "/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/WASP/WASP_Runs" which receives all run output per sample needs to be created prior to running the python file:"/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/WASP/WASP_Run_Dirs.py"

WASP=/home/asiimwe/WASP
PYTHON=/home/asiimwe/miniconda3/bin/python3.9
export PATH=/usr/bin/samtools/:$PATH
STAR=/usr/bin/STAR
genomeDirectory="/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/genome_directory/"
vcfFileDir=/scratch/asiimwe/STAR-WASP_FASTQs_VCFs/VCF
fastqFileDir=/scratch/asiimwe/STAR-WASP_FASTQs_VCFs/FASTQ
ulimit -n 10000 

WASPdir=/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/WASP/WASP_Runs/sample_id/xthreads

cd $WASPdir

STARpar="--runThreadN x --genomeDir $genomeDirectory --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM MD jM jI  --alignEndsType EndToEnd --outSAMunmapped Within --outFilterMultimapNmax 1"

readFiles="--readFilesCommand gunzip -c --readFilesIn $fastqFileDir/sample_id/R1.fastq.gz  $fastqFileDir/sample_id/R2.fastq.gz"

hetVcf="--varVCFfile $vcfFileDir/ssd_input_snp_dir/ssd.vcf.snv1het" 

$STAR $STARpar $readFiles $hetVcf

mv $WASPdir/Aligned.sortedByCoord.out.bam $WASPdir/A_sorted.bam
samtools index $WASPdir/A_sorted.bam $WASPdir/A_sorted.bai
$PYTHON $WASP/mapping/find_intersecting_snps.py --is_paired_end --is_sorted --snp_dir $vcfFileDir/ssd_input_snp_dir/SNPdir --output_dir ./ A_sorted.bam 
$STAR $STARpar $hetVcf --readFilesCommand gunzip -c --readFilesIn  A_sorted.remap.fq1.gz A_sorted.remap.fq2.gz  
mv Aligned.sortedByCoord.out.bam Aligned.out.bam
samtools index Aligned.out.bam Aligned.out.bai
$PYTHON $WASP/mapping/filter_remapped_reads.py A_sorted.to.remap.bam  Aligned.out.bam A_sorted.keep.bam
