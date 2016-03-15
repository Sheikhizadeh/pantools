Requirements:

1- KMC: add KMC/bin to the path.
2- neo4j-community-2.3.1-unix
3- Java Virtual Machine jdk1.7 : add the corresponding path to java executable

To run the program:
java  [-server] [-XX:+UseConcMarkSweepGC]  [-Xmx??g] -jar ./pantools/dist/pantools.jar command arguments


List of commands:

build:    To build a pan-genome out of a set of genomes
          java -jar pantools.jar build K database_path fasta_names_file

annotate: To add annotations to a pan-genome
          java -jar pantools.jar annotate database_path gff_names_file

add:      To add new genomes to an available pan-genome
          java -jar pantools.jar add database_path fasta_names_file

retrieve_genes:    To extract sequence of some annotated genes
          java -jar pantools.jar retrieve_genes database_path annotation_records

retrieve_regions:  To extract region sequence
          java -jar pantools.jar retrieve_regions database_path regions_file
        
group:    To group some genes by adding group nodes pointing to them
          java -jar pantools.jar group database_path groups_file

compare:     To compare two pan-genomes
          java -jar pantools.jar compare database_1_path database_2_path

List of arguments:

K :  Size of k

fasta_names_file: A text file containing paths to FASTA files; each in one line

gff_names_file  : A text file containing paths to GFF files corresponding to the stored genomes,
 in the same order. Missing annotations are shown by an empty line.

annotation_records : A text file containing annotation titles as they appear in gff file.

regions_file    : A text file containing genome_number, sequence_number, begin and
 end of a region in each line seperated by one space 

groups_file : A FASTA file with titles being names given to the groups followed by lines containing annotation titles as they appear in gff files.
