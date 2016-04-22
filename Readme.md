************************************************************************
PanTools is a disk-based java application for computational pan-genomics
developed by Siavash Sheikhizadeh et. al in Bioinformatics group of 
Wageningen university and research center, the Netherlands.  
************************************************************************

Requirements
------------
- KMC2: is a disk-based programm for counting k-mers from (possibly gzipped) FASTQ/FASTA files. 
       The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
       You need to unzipd kmc.zip and add the path of the corresponding version (linux, macos or windows) of kmc and kmc_tools executables to your OS path environment variable.
- Java Virtual Machine version 1.7 or higher: Add the path to the java executable to your OS path environment variable. 

How to run the program 
----------------------
java  [-server] [-XX:+UseConcMarkSweepGC]  [-Xmx(a number followed by g or m)] -jar ./pantools/dist/pantools.jar command arguments

List of commands
----------------
1. reconstruct:
   To reconstruct all or a set of genomes out of the pan-genome
   java -jar pantools.jar reconstruct all [or genome_names_file] database_path
2. build:
   To build a pan-genome out of a set of genomes
   java -jar pantools.jar build K database_path fasta_names_file
3. annotate:
   To add annotations to a pan-genome
   java -jar pantools.jar annotate database_path gff_names_file
4. add:
   To add new genomes to an available pan-genome
   java -jar pantools.jar add database_path fasta_names_file
5. retrieve genes:
   To extract sequence of some annotated genes
   java -jar pantools.jar retrieve_genes database_path annotation_records
6. retrieve regions:
   To extract region sequence
   java -jar pantools.jar retrieve_regions database_path regions_file
7. group:
   To group some genes by adding group nodes pointing to them
   java -jar pantools.jar group database_path groups_file
8. compare:
   To compare two pan-genomes
   java -jar pantools.jar compare database_1_path database_2_path
9. query
   To run Cypher queries.
   java -jar pantools.jar query database_path
   
Description of commands' arguments
----------------------------------
K :                 is size of K for construction of the de Bruijn graph.
genome_names_file : is a text file containing genome_number and genome_name seperated by a single space in each line.
fasta_names_file :  is a text file containing paths to FASTA files; each in one line.
gff_names_file :    is a text file containing paths to GFF files corresponding to the stored genomes,
                    in the same order. Missing annotations are shown by an empty line.
annotation_records :is a text file containing annotation titles as they appear in gff file.
regions_file:       is a text file containing genome_number, sequence_number, begin and
                    end of a region in each line seperated by one space. 
groups_file:        is a FASTA file with titles being names given to the groups followed by lines containing annotation titles as they appear in gff files.
