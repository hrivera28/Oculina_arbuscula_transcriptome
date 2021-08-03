# Scripts below were modified by H.E. Rivera following the pipelines used by the Matz Lab at UT Austin
# See https://github.com/z0on/annotatingTranscriptomes for further documentation on the original scripts

# Directory layout
# Raw reads live in /DaviesData/Raw_reads/
# Contains 48 fq files (12 samples, 2 lanes/sample, with forward and reverse reads)
# Concatenated reads live in /DaviesData/concat_reads/
# Holobiont transcriptome file is noshort.fasta and lives in /DaviesData/Trinity/trinity_out_dir/
# BLAST databases used downstream live in /DaviesData/Sarah_dbs
# Contains dirty host, dirty sym, clean host, clean sym, SILVA, and ITS dbs, and uniprot databases


########### Separate host and symbiont transcriptomes
# Using Davies's dirty host, dirty sym, clean host, clean sym fasta files (see Davies et al. 2016 Frontiers in Marine Science)


### Step 1 Make database fasta files into BLAST dbs
# module load bio
# module load blast

makeblastdb -in /DaviesData/Sarah_dbs/host_db/clean_host.fasta -dbtype nucl
makeblastdb -in /DaviesData/Sarah_dbs/sym_db/clean_sym.fasta -dbtype nucl
makeblastdb -in /DaviesData/Sarah_dbs/dirty_host_db/dirty_host.fasta -dbtype nucl
makeblastdb -in /DaviesData/Sarah_dbs/dirty_sym_db/dirty_sym.fasta -dbtype nucl

# copy the noshorts.fasta transcriptome file into directory

### Step 2 Split transcriptome into chucks to speed up blast
cp /DaviesData/Trinity/trinity_out_dir/noshorts.faasta /DaviesData/Sarah_dbs/dirty_host_db/

splitFasta.pl noshorts.fasta 120 # splits transcriptome into 120 files (see additional file)

### Step 3 Run blast for each db
# I wrote these as sbatch arrayed files in practice
# but the idea would be something like:

ls -1 subset[1-9]_*.fasta | while read line; do
	tblastx -query $line -db /Sarah_dbs/dirty_host_db/dirty_host.fasta -evalue 1e-3 -num_threads 36 -max_target_seqs 5 -outfmt "7 qseqid sseqid sgi evalue bitscore score length nident mismatch positive qstart qend sstart send qframe staxids stitle" -out $line.br
done
# This whole thing was repeated by directory (i.e. dirty host, dirty sym, etc).

# Each blast result was written to a subset[].br file
# These were then combined within each directory
cat subset*br > dirtyhost_blast.br
# etc for the four directories

### Step 4 Parse the blast results and produce clean host and sym transcript files
# Uses sarahs_blastparser.py script to get lists of which contigs belong to which set (see additional files)

# module load python2/2.7.14

# The 100 is the length of overlap and the 80 is the min percent identity to pull a match from the blast file

python /Sarah_scripts/sarahsblastparser.py 100 80 dirty_host_blast.br
python /Sarah_scripts/sarahsblastparser.py 100 80 dirty_sym_blast.br
python /Sarah_scripts/sarahsblastparser.py 100 80 clean_sym_blast.br
python /Sarah_scripts/sarahsblastparser.py 100 80 clean_host_blast.br

# ouput files will be names dirty_host_blast_goodmatches.txt, etc. 

### Step 5 Obtain lists of sequences to remove
# Final host contigs = dirty host contigs - clean sym contigs
# Final sym contigs = dirty sym contigs - clean host contigs

# module load anaconda
# source activate biopython

# This gets you a list of the genes that you want to keep for each set
comm -23 dirty_host_blast_goodmatches.txt clean_sym_blast_goodmatches.txt > clean_host_contigs_list.txt
comm -23 dirty_sym_blast_goodmatches.txt clean_host_blast_goodmatches.txt > clean_sym_contigs_list.txt

# use Sarah's get_seq2.py script to get clean fastas using BioPython
# args are input_fasta seq_ids_to_keep output_fasta
python /Sarah_scripts/get_seq2.py /DaviesData/Trinity/trinity_out_dir/noshorts.fasta clean_sym_contigs_list.txt sym_contigs_for16Sblast.fasta
python /Sarah_scripts/get_seq2.py /DaviesData/Trinity/trinity_out_dir/noshorts.fasta clean_host_contigs_list.txt host_contigs_for16Sblast.fasta


### Step 6 Remove bacterial contaminants/ribosomal seqs
# Using SILVA databases

# SILVA dbs: SILVA_132_LSUParc and SILVA_132_SSURef_Nr99 from SILVA website
# Combine the two fasta files to make one SILVA sequence file (SILVA_SSU_LSU_combined.fasta)

# module load bio
# module load blast

# Make blastdb for it
makeblastdb -in SILVA_SSU_LSU_combined.fasta -dbtype nucl

blastn -query host_contigs_for16Sblast.fasta -db SILVA_SSU_LSU_combined.fasta -outfmt 5 -evalue 0.05 -num_threads 36 -num_descriptions 20 -gapopen 5 -gapextend 2 -penalty -3 -num_alignments 20 -reward 2 -out SILVA_host_match.br
blastn -query sym_contigs_for16Sblast.fasta -db SILVA_SSU_LSU_combined.fasta -outfmt 5 -evalue 0.05 -num_threads 36 -num_descriptions 20 -gapopen 5 -gapextend 2 -penalty -3 -num_alignments 20 -reward 2 -out SILVA_sym_match.br

### Step 7 Remove SILVA contaminants

# parse the blast results
# module load bio
# module load blast
# module load python2/2.7.14

python2 ./Sarah_scripts/parse_blastn.py shared_hits.txt blastn SILVA_shared_seqs_match.br
python2 ./Sarah_scripts/parse_blastn.py host_hits.txt blastn SILVA_host_match.br
python2 ./Sarah_scripts/parse_blastn.py sym_hits.txt blastn SILVA_sym_match.br

## Then I used the same approach as with the dirty host dirty sym dbs to get a list of contigs without any of the 16S matches
## and then used the get_seq script to crate a new fasta file with just the final host, sym, shared seqs

cut -f 1 host_hits.txt | grep -E "TRINITY_\w+_\w+_\w+_\w+" -o | sort | uniq > host_16S_matches_toremove.txt
cut -f 1 sym_hits.txt | grep -E "TRINITY_\w+_\w+_\w+_\w+" -o | sort | uniq > sym_16S_matches_toremove.txt

# There were 78 sequences in the host files that had a good blast hit
# and 48 seqs in the sym file that had a good hit

# Generate full list of sequences from the host and sym "cleaned" transcriptomes produced in Step 5
grep 'TRINITY_\w+_\w+_\w+_\w+' -o host_contigs_for16Sblast.fasta > host_removed_symbiont_contigs.txt
grep 'TRINITY_\w+_\w+_\w+_\w+' -o sym_contigs_for16Sblast.fasta > sym_removed_symbiont_contigs.txt
# Gets list of genes to keep
comm -23 host_removed_symbiont_contigs.txt sym_16S_matches_toremove.txt > clean_sym_contigs_list.txt
comm -23 sym_removed_host_contigs.txt host_16S_matches_toremove.txt > clean_host_contigs_list.txt

# Produces clean fasta files
python /Sarah_scripts/get_seq2.py sym_contigs_for16Sblast.fasta clean_sym_contigs_list.txt sym_transcripts_for_anno.fasta
python /Sarah_scripts/get_seq2.py host_contigs_for16Sblast.fasta clean_host_contigs_list.txt host_transcripts_for_anno.fasta
