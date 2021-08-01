## Oculina arbuscula transcriptome
This repo contains the data and bioinformatics scripts associated with generating the Oculina arbuscula transcriptome along with a less complete transcriptome of its endosymbiont Breviolum psygmophilum. 

### Coral colony collection and sampling
Six colonies of the facultatively symbiotic coral Oculina arbuscula (Agassiz, 1864) were collected from 15-20 ft depth at Radio Island Jetty, North Carolina (34˚ 42.520’ N, 76˚ 40.796’ W) using a hammer and chisel on May 25, 2018. All O. arbuscula colonies were collected under the NC Division of Marine Fisheries Permits #706481 and #1627488. Corals were over-night transported to the Davies Marine Population Genomics Lab at Boston University and maintained in aquaria at 25˚C (+/-0.17, SD), 34.6 PSU (+/-0.76, SD), and pH of 8.0 (+/-0.08, SD). Because symbiont density varied substantially among branches of the same colonies, two tissue samples were collected from each colony (A-F), one from a branch with high symbiont density (brown in color, assumed symbiotic) and another from a low symbiont density branch (pale or white in color, assumed aposymbiotic). Colonies were preserved in 200 proof EtOH and maintained at -80˚C until RNA isolation. This sampling yielded a total of 12 unique transcriptomic libraries (colonies A-F, with symbiotic and aposymbiotic branches). 

### RNA extraction and sequencing
RNA was extracted using the Invitrogen™ RNAqueous™ Total RNA Isolation Kit, with the addition of a DNase step to remove genomic DNA. Extracted RNA was normalized across samples using a NanoDrop™ spectrophometer. Library prep and sequencing was performed by NovoGene on two lanes of an Illumina HiSeq 4000 sequencer as 150bp paired-end reads. Sequencing depth ranged from 9.8 to 13.3 million reads per sample (11.5 million average), with average quality ranging from 40-42 across the full read length.

### Holobiont transcriptome assembly (transcriptome_assembly.sh)  
A total of 277 million paired-end reads were quality filtered and trimmed using trimmomatic through Trinity (v. 2.6.6) with the following settings: Illumina TruSeq paired end adapters were removed (“TruSeq-PE-2” file), with default Trimmomatic parameters. Reads <25 bases long after filtering were removed. 
A total of 268,009,326 (96.72%) of reads were used for de novo assembly using Trinity (v.2.6.6) with default parameters and in-silico normalization of reads. Contigs <500 bases after assembly were removed. The final holobiont assembly contained 563,264 contigs with an average length of 1,092 and an N50 of 1,202. 

### Separation of coral host and algal symbiont transcripts (transcriptome_host_sym_fraction.sh)  
Four custom databases previously described in Davies et al. 2016 were used to filter and assign contigs to host or symbiont origin. Contigs were aligned to database reads using BLAST, in tblastx mode, E-value of 0.001, and maximum of five target sequences. BLAST results were parsed (blast_parser.py script in Sarah's scripts folder) to pull matches that were at least 100 bp in length and minimum of 80% sequence identity. Matches by read were then compared. Coral-assigned transcripts consisted of reads mapping to the holobiont database with any reads that also mapped to the clean symbiont database removed. The symbiont fraction was generated from reads mapping to the holobiont database with reads mapping to the clean cnidarian database removed. Any still reads assigned to both host and symbiont fractions after these steps were also subsequently removed. 

Any potential microbial contamination in both symbiont and host transcriptome fractions were identified by blasting against the SILVA database (132 LSUParc and SSURef_Nr99 combined – accessed Dec 2017). BLAST parameters were E-value of 0.05, a gap open penalty of 5, and gap extend penalty of 2, and penalty of 3. Reads assigning to the SILVA database were removed. 

### Further collapsing of isoforms and transcription annotation (collapse_annotate_clean_transcriptomes.sh)
Highly similar isoforms generated from the Trinity assembly were further collapsed to keep only the longest isoform as a representative using CD-HIT (Li and Godzik 2006) (parameters: -c 0.99 -G 0 -aL 0.3 -aS 0.3). 

Transcriptomes were then annotated through blastx against the unitprot database (downloaded June 12, 2019), with an E-value of 0.001. Reads that were annotated as “chroloplastic” or “chlorophyll” were removed from the coral transcriptome file (N=93). Genes annotated to Acropora millepora or Nematostella vectensis were removed from the symbiont transcriptome files (N=12). 

### BUSCO scoring of transcriptomes (busco.sh) 
To quantify transcriptome completeness, the coral and symbiont transcriptomes compared against the BUSCO Eukaryote_ob10, Metazoa_ob10 (coral only), and Alveolata_ob10 (symbiont only) databases using BUSCO v.4.0.5. 

Final transcriptome files are avaialable from the Davies Marine Population Genomics Lab Website: http://sites.bu.edu/davieslab/data-code/

![Graph of BUSCO scores for coral (host) and symbiont trasncritomes against BUSCO databases](https://github.com/hrivera28/Oculina_arbuscula_transcriptome/blob/master/BUSCO_collapsed_alv.png?raw=true)
