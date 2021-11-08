# P._stutzeri_Complex_Manuscript
Serveral Perlscript are used in this Manuscript titled "Comparative Genomics of Pseudomonas stutzeri Complex: Taxonomic Assignments and Genetic Diversity".

#### 1. Pan-genome analysis

* extract_protein_dir.pl: Run this command to enble users to extract protein sequences for genomes using multiple threads.  

Usage: extract_protein_dir.pl -dir genbank_file_directory [options]

Firstly, put all genomes used for pan genome analysis under the same directory, then run using the following command: "perl stutzeri_123 -m 3". The "stutzeri_123" is the directory containing annotated genomes as Genbank format file, and the parameter "-m " is used to set thread number. The detailed usage instructions is included in extract_protein_dir.pl.  
    
After all predicted protein-coding sequences (CDSs) were extracted from each of 123 genomes separately using extract_protein_dir.pl, the output directory is used as input for Orthofinder to infer orthogroups. When orthofinder is done, Orthogroups.tsv and Orthogroups_UnassignedGenes.tsv are obtained. The former is a tab separated text file, in which each row contains the genes belonging to a single orthogroup. The latter is a tab separated text file that is identical in format to Orthogroups.csv but contains all of the genes that were not assigned to any orthogroup. By runing the command "cat" in Unix/Linux system, the users can merge these two files into a single file (refer to merged_Orthogroups.txt), which contains all gene families for all analyzed genomes. This merged file is modified using a Perlscript "format_modification.pl", then PGAP can use the modified file and "genomic_name.txt" as inputs to perform Pan-genome analysis. 

* format_modification.pl: This Perlscript is used to modify the format of merged file, and allows PGAP to use this file as input file.

Usage: format_modification.pl genomic_name.txt merged_Orthogroups.txt OUT.txt
"genomic_name.txt" is a tab separated text containing only one row, which contains all analyzed genomic names. 

"merged_Orthogroups.txt" is genarated by merging the Orthogroups.tsv and Orthogroups_UnassignedGenes.tsv.

"OUT.txt" is the output file, which is used as input file for PGAP to perfomr Pan-genome analysis. 

Please refer to PGAP manual (https://sourceforge.net/projects/pgap/) for Pan-genome analysis. 
   

#### 2. Concatenate alignment files
* Concat_Seq.pl: Run this command to enble users to concatenate alignment files into a pseudo-DNA fasta file.   

FOR EXAMPLE: perl Concat_Seq.pl -dir alignment_dir

"alignment_dir" is a directory containing alignment files as FASTA format file.

After run Concat_Seq.pl, a concatenation sequence file called "concatenation.fasta" is produced in directory where Concat_Seq.pl is located. The detailed usage instructions is included in extract_protein_dir.pl.

#### 3. COG annotation analysis
* COG_annotation.pl: A Perlscript is used to complete the whole process of cog annotation using COG database and protein sequences as FASTA format file as inputs. 

Usage: COG_annotation.pl cog-20.fa cog-20.cog.csv cog-20.def.tab fun-20.tab merged_cogs.fa Pan_repersnet_seq.fasta COG_out 16

The latest version of COG database can be download in website https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/, which contains cog-20.fa, cog-20.cog.csv, cog-20.def.tab, fun-20.tab and merged_cogs.fa files. 

"Pan_repersnet_seq.fasta" is the protein sequences as FASTA format file.

"16" refers to 16 thread number used when runing COG_annotation.pl

"COG_out" is the COG annotation result, which a tab separated text file, in which each row contains the protein ID, COG functional category (could include multiple letters in the order of importance), and  COG functional category as single letter. If a gene was assigned to more than one COG category, each COG category is shown as separate row. 

    For example, the content of "COG_out" file is as follows: 
    Shell|6(323_genes|109_taxa)	LX	L
    Shell|6(323_genes|109_taxa)	LX	X
    Softcore|7(320_genes|119_taxa)	G	G
