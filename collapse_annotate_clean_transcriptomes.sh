# Scripts below were modified by H.E. Rivera following the pipelines used by the Matz Lab at UT Austin
# See https://github.com/z0on/annotatingTranscriptomes for further documentation on the original scripts

### Step 1 Collapse really similar isoforms that Trinity didn't call the same thing

# module load anaconda/5.2.1
# source activate cd-hit

cd-hit-est -i host_transcripts_for_anno.fasta -o host_collapsed_transcripts.fasta -c 0.99 -G 0 -aL 0.3 -aS 0.3
cd-hit-est -i sym_transcripts_for_anno.fasta -o sym_collapsed_transcripts.fasta -c 0.99 -G 0 -aL 0.3 -aS 0.3

## Edit the headers so that it's what Misha's scripts want:
# >TRINITY_XXX gene=isogroupXXX
# I'm making mine >Oarb_XXX gene=isogroupXXX and similarly for symbiont >Sym..

cat host_collapsed_transcripts.fasta | sed -E 's/len=.*//'| sed -E 's/^>TRINITY_DN((.*)_i[0-9]+)/>Oarb_\1 gene=isogroupOarb_\2/' > host_collapsed_short_header_forblast.fasta
cat sym_collapsed_transcripts.fasta | sed -E 's/len=.*//' | sed -E 's/^>TRINITY_DN((.*)_i[0-9]+)/>Sym_\1 gene=isogroupSym_\2/' > sym_collapsed_short_header_forblast.fasta


# subset fastas to speed up blasts
/Sarah_scripts/splitFasta.pl sym_collapsed_short_header_forblast.fasta 8
/Sarah_scripts/splitFasta.pl host_collapsed_short_header_forblast.fasta 8

#### blast transcriptomes #as job scripts
# module load bio
# module load blast

# line below works if you're using slurm arrays for parallelization
blastx -query subset${SLURM_ARRAY_TASK_ID}_host_collapsed_short_header_forblast.fasta -db uniprot_sprot.fasta -evalue 1e-3 -num_threads 18 -num_descriptions 5 -num_alignments 5 -out subset${SLURM_ARRAY_TASK_ID}.fasta.br
blastx -query subset${SLURM_ARRAY_TASK_ID}_sym_collapsed_short_header_forblast.fasta -db uniprot_sprot.fasta -evalue 1e-3 -num_threads 18 -num_descriptions 5 -num_alignments 5 -out subset${SLURM_ARRAY_TASK_ID}.fasta.br

# to concatenate the blast results cut the first 14 lines of the subset files
ls -1 subset[1-8]*.fasta.br | while read line; do tail -n +15 $line >> host_blastx.br; done
ls -1 subset[1-8]*.fasta.br | while read line; do tail -n +15 $line >> sym_blastx.br; done

# Generate seq2iso files
grep ">" host_collapsed_short_header_forblast.fasta | sed -E 's/>((.*)_i[0-9]+)/\1\t\2/' |sed -E 's/\sgene=.*//'>host_seq2iso.tab
grep ">" sym_collapsed_short_header_forblast.fasta | sed -E 's/>((.*)_i[0-9]+)/\1\t\2/' |sed -E 's/\sgene=.*//'>sym_seq2iso.tab

# Get GO terms and gene names from the uniprot database
# Scripts from https://github.com/z0on/annotatingTranscriptomes
# For coral
getGOfromUniProtKB.pl blast=host_blastx.br prefix=host fastaQuery=host_collapsed_short_header_forblast.fasta db=uniprot_sprot.fasta >> getgo_host
getGeneNamefromUniProtKB.pl blast=host_blastx.br prefix=host fastaQuery=host_collapsed_short_header_forblast.fasta db=uniprot_sprot.fasta >> getgn_host
# For Sym
getGOfromUniProtKB.pl blast=sym_blastx.br prefix=sym fastaQuery=sym_collapsed_short_header_forblast.fasta db=../../uniprotdb/uniprot_sprot.fasta >> getgo_sym
getGeneNamefromUniProtKB.pl blast=sym_blastx.br prefix=sym fastaQuery=sym_collapsed_short_header_forblast.fasta db=../../uniprotdb/uniprot_sprot.fasta >> getgn_sym

###Removing some Sym contamination from the host transcriptome:
# Grep gene names for [Cc]hlorop and [pP]hotop and remove those contigs from the host transcriptome
grep '[cC]hlorop' host_iso2gene.tab | cut -f 1 | sed -E 's/isogroup//' >> genes_to_remove.txt
# wc -l
# 51 genes
grep '[pP]hotop' host_iso2gene.tab |  cut -f 1 >> genes_to_remove.txt
# wc -l
# No genes woot!
# add back in the gene=isogroupXXX that's part of the seq header in the fasta file
cat genes_to_remove.txt | while read line; do grep $line  host_seq2iso.tab |  sed -E 's/\t(Oarb.*)/ gene=isogroup\1/' >> seqIDs_toss.txt ; done
# There are 93 seq IDs here because there are still a few isogroups that didn't get collapsed but which got assigned the same gene

# use mothur to filter the fasta
# module load boost/gcc/1.67
# module load mothur
mothur #enter interactive mode
mothur > remove.seqs(accnos=seqIDs_toss.txt, fasta=host_collapsed_short_header_forblast.fasta)

mv host_collapsed_short_header_forblast.pick.fasta Brevolium_psygmophillum_from_Oculina_arbuscula_transcriptome.fasta


### For sym transcriptome
grep -E 'OS=\w+\s\w+' -o sym_iso2gene.tab | sort| uniq -c  >>sym_species_hits.txt
# there were two genes that were annotated to acropora and 10 to nematostella
# I'll remove those from the sym transcriptome

grep 'Acropora' sym_iso2gene.tab | cut -f 1 | sed -E 's/isogroup//' >> sym_genes_to_remove.txt
grep 'Nematostella' sym_iso2gene.tab | cut -f 1 | sed -E 's/isogroup//' >> sym_genes_to_remove.txt

cat sym_genes_to_remove.txt | while read line; do grep $line sym_seq2iso.tab |  sed -E 's/\t(Sym.*)/ gene=isogroup\1/' >> sym_seqIDs_toss.txt ; done

#module load boost/gcc/1.67
#module load mothur
mothur #enter interactive mode
mothur > remove.seqs(accnos=sym_seqIDs_toss.txt, fasta=sym_collapsed_short_header_forblast.fasta)

mv sym_collapsed_short_header_forblast.pick.fasta Brevolium_psygmophillum_from_Oculina_arbuscula_transcriptome.fasta


#Remove the bad seqs and isos from the seq2iso, iso2go and iso2gene files
cat genes_to_remove.txt | while read line; do sed -i "/isogroup$line.*/d" host_iso2gene.tab ; done
cat genes_to_remove.txt | while read line; do sed -i "/isogroup$line.*/d" host_iso2go.tab; done
sed -E 's/gene=.*//' seqIDs_toss.txt > seqs_forclean_seq2iso_list.txt
cat seqs_forclean_seq2iso_list.txt | while read line; do sed -i "/$line.*/d" host_seq2iso.tab; done

cat sym_genes_to_remove.txt | while read line; do sed -i "/isogroup$line.*/d" sym_iso2gene.tab ; done
cat sym_genes_to_remove.txt | while read line; do sed -i "/isogroup$line.*/d" sym_iso2go.tab; done
sed -E 's/gene=.*//' sym_seqIDs_toss.txt > sym_seqs_forclean_seq2iso_list.txt
cat sym_seqs_forclean_seq2iso_list.txt | while read line; do sed -i "/$line.*/d" sym_seq2iso2.tab; done
