# Scripts below were modified by H.E. Rivera following the pipelines used by the Matz Lab at UT Austin
# See https://github.com/z0on/annotatingTranscriptomes for further documentation on the original scripts

# Directory layout
# Raw reads live in /DaviesData/Raw_reads/
# Contains 48 fq files (12 samples, 2 lanes/sample, with forward and reverse reads)

#######################################
# Preparing raw files, running Trinity and getting some long contigs

### Step 1
# Concatenate all the raw forward reads and all the raw reverse reads (in same order of files)

cd ./DaviesData/Raw_reads

cat *_1.fq >> all_1.fq
cat *_2.fq >> all_2.fq

# Move concatenated files into their own directory
mv all*.fq  ./DaviesData/concat_reads/

### Step 2 Run Trinity

# Read trimming and quality filtering was done within Trinity
# Using the trimmomatic option
# module loading will depend on cluster and environment specifications
# module load bio
# module load samtools
# module load jellyfish
# module load salmon
# module load bowtie2
# module load trinity

Trinity --seqType fq --max_memory 190G --left /DaviesData/concat_reads/all_1.fq --right /DaviesData/concat_reads/all_2.fq --CPU 19 --trimmomatic --quality_trimming_params "ILLUMINACLIP:/vortexfs1/omics/tarrant/hrivera/DaviesData/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25" >> log_run.txt

# the resulting assembled transcriptome file is Trinity.fasta

# Step 4 Remove short contigs (<500 bp) using noshort.pl (see additional file)
# module load bio

/Sarah_scripts/noshorts.pl Trinity.fasta 500

# retained:	 563,264
# discarded: 783,085 (lost 58%)

# The transcriptome file that will be used downstream is now called noshorts.fasta
