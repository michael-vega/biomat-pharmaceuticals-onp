### This workflow describes the generic metagenomics pipeline used for the manuscript, "Pharmaceutical biotransformation is influenced by photosynthesis and microbial nitrogen cycling in a benthic wetland biomat"

### Additional details and steps are described in the manuscript and supporting information, as well as on The Wrighton Lab GitHub page: https://github.com/TheWrightonLab/

### The majority of commands were performed using a Slurm job management system and commands are only generically described here 

### Packages Used ###
		#Sickle [v1.33]
		#idba_ud [v1.1.0]
		#EukRep [v0.6.5]
		#metaBAT2 [v2.12.1]
		#checkM [v1.1.2]
		#GTDB-tk [v1.5.0]
		#dRep [v2.6.2]
		#DRAM [v1.2]
		#Bowtie2 [v2.3.5]
	
###Download raw filtered metagenomic reads: NCBI and JGI Gold Accessions are available in Table S5 of the manuscript

####Example commands in this workflow use "FASTA.fastq.gz" as file names; FASTA should be replaced throughout by the appropriate filename downloaded from NCBI or JGI

###Unzip fastq.gz files and deinterleave paired end reads

gunzip -c FASTA.fastq.gz > FASTA.fastq

###Deinterleave forward and reverse reads from the fastq file 

deinterleave_fastq.sh < FASTA.fastq FASTA.fastq_R1.fastq FASTA.fastq_R2.fastq

###Trim deinterleaved reads with Sickle, plot quality and trimming info, then assemble with idba-ud using default parameters

box_trim_box_assemble.sh FASTA.fastq_R1.fastq FASTA.fastq_R2.fastq 10

###Pull out scaffolds/contigs less than 2500 sequences and calculate stats

pullseq.py -i scaffold.fa -m 2500 -o FASTA_contigs_2500.fa

contig_stats.pl -i FASTA_contigs_2500.fa -o FASTA_contigs_2500_stats

### For the 'bottom' sample only, map all 'bottom' reads to all 'bottom' bins, then perform a subassembly on the unbinned reads ###  

#Concatenating bins from the 'bottom' sample

cat /Bottom_bins/*.fa > Bottom_concatenated_bins.fa

#Mapping 'bottom' sample reads to bins from the 'bottom' sample 

bbmap.sh -Xmx48G threads=6 minid=100 overwrite=t ref=Bottom_concatenated_bins.fa in1=Bottom_R1_All_trimmed.fastq in2=Bottom_R2_All_trimmed.fastq outu1=Bottom_NotMappedToBins_R1_ToAssembly.fastq outu2=Bottom_NotMappedToBins_R2_ToAssembly.fastq

#Same script as above to trim deinterleaved reads with Sickle, plot quality and trimming info, then assemble ('subassemble') with idba-ud using default parameters

box_trim_box_assemble.sh Bottom_NotMappedToBins_R1_ToAssembly.fastq Bottom_NotMappedToBins_R2_ToAssembly.fastq 10

### EukRep to parse prokaryotic scaffolds ###

#Note: performed on all assemblies from all samples, including the 'bottom' subassembly 

EukRep -i FASTA_contigs_2500.fa -o FASTA_contigs_2500_eukrep.fa --min 2500 --prokarya FASTA_contigs_2500_prokaryote.fa  

### Binning ### 

#Note: binning was performed with prokaryote contigs from EukRep (e.g., FASTA_contigs_2500_prokaryote.fa) but also unparsed contigs (e.g., FASTA_contigs_2500_eukrep.fa)

#Build index file

bowtie2-build FASTA_contigs_2500_prokaryote.fa FASTA_contigs_2500_prokaryote_index 

#Align index file against trimmed reads to create sequence alignment map (SAM) file

bowtie2 --fast -p 6 -x FASTA_contigs_2500_prokaryote_index  -S All_mappedtoall_paired.sam -1 R1_All_trimmed.fastq -2 R2_All_trimmed.fastq --un unmapped_paired.fq --al mapped_paired.fq > bowtie_log  

#Convert SAM to BAM with samtools

samtools view -@ 6 -bS All_mappedtoall_paired.sam > FASTA_contigs_2500_prokaryote.bam 

#Sort BAM file with samtools

samtools sort -@ 6 -T FASTA_contigs_2500_prokaryote.bam.sorted -o FASTA_contigs_2500_prokaryote.sorted.bam FASTA_contigs_2500_prokaryote.bam 

#Use sorted BAM file to estimate coverage and depth, then create bins with metaBAT2

runMetaBat.sh -t 6 --verysensitive FASTA_contigs_2500_prokaryote.fa FASTA_contigs_2500_prokaryote.sorted.bam

Executes: 'jgi_summarize_bam_contig_depths --outputDepth FASTA_contigs_2500_prokaryote.fa.depth.txt --pairedContigs FASTA_contigs_2500_prokaryote.fa.paired.txt --minContigLength 1000 --minContigDepth 1 FASTA_contigs_2500_prokaryote.sorted.bam' 

Executes: 'metabat2  -t 6 --verysensitive --inFile FASTA_contigs_2500_prokaryote.fa --outFile FASTA_contigs_2500_prokaryote.fa.metabat-bins--verysensitive/bin --abdFile FASTA_contigs_2500_prokaryote.fa.depth.txt' 

### checkM ### 

checkm lineage_wf -t 6 -x fa FASTA_contigs_2500_prokaryote.fa.metabat-bins--verysensitive /checkm

### Bin Dereplication ###

dRep dereplicate -p 6 -comp 50 -con 10 --genomeInfo dRep_checkm_input_final.csv dRep_out_final -g MAG_Database/*.fa

### Taxonomy Assignment ###

gtdbtk classify_wf --genome_dir MAG_Database/ --out_dir classify_wf_output_RS202 --extension fa --cpus 5

### Annotate MAG database with DRAM ###

DRAM.py annotate -i ${filelist} -o ${filelist}_DRAMOUT --threads 16 --min_contig_size 2499

### Concatenating all MAGs to a single fasta for later steps ###

cat /MAG_Database/*.fa > Concatenated_dRep_MAGs_MediumHighOnly.fa

### Read mapping with BowTie2 ###

#Base script below with subsequent lines run internally as part of map_reads_fasta_abundance.py

python map_reads_fasta_abundance.py -f Concatenated_dRep_MAGs_MediumHighOnly.fa -r reads.txt -n names.txt -t MapReads2Bins.txt -p 30 -m 3 --min_coverage 5 --min_percent_contig_coverage 80 -u F -d F

#Build index
bowtie2-build Concatenated_dRep_MAGs_MediumHighOnly.fa Concatenated_dRep_MAGs_MediumHighOnly.fa

#Then, for each FASTA in "names.txt" with corresponding filtered reads file names listed in "reads.txt", the below three lines are repeated: 

bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 30 -x Concatenated_dRep_MAGs_MediumHighOnly.fa -S FASTA_to_Concatenated_dRep_MAGs_MediumHighOnly.fa.sam -1 FASTA/JGI_nonfiltered_reads/R1_ALL_trimmed.fastq -2 FASTA/JGI_nonfiltered_reads/R2_ALL_trimmed.fastq

python sam_file.py -i FASTA_to_Concatenated_dRep_MAGs_MediumHighOnly.fa.sam -v 3 -o mismatches_3_FASTA_to_Concatenated_dRep_MAGs_MediumHighOnly.fa.sam

rm FASTA_to_Concatenated_dRep_MAGs_MediumHighOnly.fa.sam

### Parse unbinned contigs/scaffolds from MAG database ###

pullseq -i FASTA_contigs_2500_prokaryote.fa -n Concatenated_dRep_MAGs_MediumHighOnly.fa -e > FASTA_unbinned_scaffolds.fa

### Annotate unbinned scaffolds with DRAM ###

DRAM.py annotate -i FASTA_unbinned_scaffolds.fa -o FASTA_unbinned_scaffolds.fa_DRAMOUT --threads 16 --min_contig_size 2499

