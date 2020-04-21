# module load anaconda/5.1
# module load blast/2.7.1
# source activate busco4

#BUSCO host
busco -i host_collapsed_short_header_forblast.fasta  -o host_busco_results -l metazoa_odb10 -m tran -c 12 -f
busco -i host_collapsed_short_header_forblast.fasta  -o host_busco_results -l eukaryota_odb10 -m tran -c 12

#BUSCO sym
busco -i sym_collapsed_short_header_forblast.fasta -o sym_busco_results -l embryophyta_odb10 -m tran -c 12 -f
busco -i sym_collapsed_short_header_forblast.fasta -o sym_busco_results -l eukaryota_odb10 -m tran -c 12 -f

#conda deactivate
