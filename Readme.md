Requirements:

- KMC: add paths to kmc and kmc_tools executables to the OS path environment variable.
- neo4j-community-2.3.1-unix
- Java Virtual Machine jdk1.7 : add the corresponding path to java executable

To run the program:
java  [-server] [-XX:+UseConcMarkSweepGC]  [-Xmx(a number followed by g or m)] -jar ./pantools/dist/pantools.jar command arguments


List of commands:

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

Description of commands' arguments:

a) K :
Size of K for de Bruijn graph

b) genome_names_file:
A text file containing genome_number and genome_name seperated by a single space in each line

c) fasta_names_file:
A text file containing paths to FASTA files; each in one line

d) gff_names_file:
A text file containing paths to GFF files corresponding to the stored genomes,
in the same order. Missing annotations are shown by an empty line.

e) annotation_records:
A text file containing annotation titles as they appear in gff file.

f) regions_file:
A text file containing genome_number, sequence_number, begin and
end of a region in each line seperated by one space 

g) groups_file:
A FASTA file with titles being names given to the groups followed by lines containing annotation titles as they appear in gff files.
