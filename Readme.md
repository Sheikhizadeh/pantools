************************************************************************
PanTools is a disk-based java application for computational pan-genomics
developed by Siavash Sheikhizadeh et. al in Bioinformatics group of 
Wageningen university and research center, the Netherlands.  
************************************************************************

Requirements
------------
- KMC2: is a disk-based programm for counting k-mers from (possibly gzipped) FASTQ/FASTA files( http://sun.aei.polsl.pl/kmc ).
        You need to unzipd kmc.zip file provided in this release of PanTools and add the path of the corresponding version (linux, macos or windows) of kmc and kmc_tools executables to your OS path environment variable.

- Java Virtual Machine version 1.8 or higher: Add the path to the java executable to your OS path environment variable.

- MCL: The Markov Cluster Algorithm, is a fast and scalable unsupervised cluster algorithm for graphs which is needed for group functionality of PanTools.
       You need to download and compile it using the instruction provided by the publishers in http://micans.org/mcl, and also add the path to the mcl executable to your path environment variable.

How to run the program 
----------------------
java  [-server] [-XX:+UseConcMarkSweepGC]  [-Xmx(a number followed by g or m)] -jar ./pantools/dist/pantools.jar command arguments


List of commands and examples for the provided sample data :

1. build:
   To build a pan-genome out of a set of genomes.

   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  build  K_VALUE  PATH_TO_THE_PANGENOME_DATABASE  PATH_TO_THE_GENOMES_PATH_FILE

   K_VALUE :   size of K for construction of the de Bruijn graph which should be 6 <= K <= 255.
   PATH_TO_THE_GENOMES_PATH_FILE : a text file containing paths to FASTA files (genomes); each in a seperated line.
   PATH_TO_THE_PANGENOME_DATABASE : path where the resulting pangenome is stored.  

   Example: 
   
   java  -Xmx4g  -jar  /home/sheik005/pantools/dist/pantools.jar build  15  /home/sheik005/two_hiv_pangenome_database  /home/sheik005/pantools/example/sample_genomes_path.txt
             
2. annotate:
   To add annotations to a pan-genome. For each genome, this function also produce a FASTA file containing all the protein sequences in the same order as they have been annotated in the GFF file.

   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  annotate  PATH_TO_THE_PANGENOME_DATABASE  PATH_TO_THE_ANNOTATION_PATH_FILE

   PATH_TO_THE_ANNOTATION_PATH_FILE : a text file containing paths to the GFF files corresponding to the genomes in the same order as they apear in PATH_TO_THE_GENOMES_PATH_FILE.
                                      Missing annotations should be indicated by an empty line.

   Example: 

   java  -jar  /home/sheik005/pantools/dist/pantools.jar  annotate  /home/sheik005/two_hiv_pangenome_database  /home/sheik005/pantools/example/sample_annotations_path.txt

3. add:
   To add new genomes to an available pan-genome.

   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  add  PATH_TO_THE_PANGENOME_DATABASE  PATH_TO_THE_NEW_GENOMES_PATH_FILE
   
   PATH_TO_THE_NEW_GENOMES_PATH_FILE : a text file containing paths to FASTA files of the new genomes to be added to the pangeome; each in a seperated line.
                                       New genomes could also be annotated later in the same way; however, there should be empty lines in the annotations path file for pre-existing genomes.

4. retrieve genes:
   To extract sequence of some annotated genes. 

   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  retrieve  genes  PATH_TO_THE_PANGENOME_DATABASE  PATH_TO_THE_ANNOTATION_RECORDS_FILE

   PATH_TO_THE_ANNOTATION_RECORDS_FILE : a text file containing records of annotated genes, as they appear in GFF file, to be retrieved.
                                         The resulting FASTA file would have the same name as the PATH_TO_THE_ANNOTATION_RECORDS_FILE with an additional .fasta extention.

   Example: 

   java  -jar  /home/sheik005/pantools/dist/pantools.jar  retrieve  genes  /home/sheik005/two_hiv_pangenome_database  /home/sheik005/pantools/example/sample_annotaion_records.txt

5. retrieve regions:
   To extract sequence of some genomic regios.

   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  retrieve  regions  PATH_TO_THE_PANGENOME_DATABASE  PATH_TO_THE_GENOMIC_REGIONS_FILE

   PATH_TO_THE_GENOMIC_REGIONS_FILE : a text file containing records with genome_number, sequence_number, begin and end positions seperated by one space for each region.
                                      The resulting FASTA file would have the same name as the PATH_TO_THE_GENOMIC_REGIONS_FILE with an additional .fasta extention.

   Example: 

   java  -jar  /home/sheik005/pantools/dist/pantools.jar  retrieve  regions  /home/sheik005/two_hiv_pangenome_database  /home/sheik005/pantools/example/sample_genomic_regions.txt

6. reconstruct:
   To reconstruct all or a set of genomes out of the pan-genome.

   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  reconstruct all/PATH_TO_THE_GENOME_NAMES_FILE  PATH_TO_THE_PANGENOME_DATABASE

   PATH_TO_THE_GENOME_NAMES_FILE : a text file containing genome_number and a given name for that genome seperated by a single space in each line. 
                                   The resulting FASTA files are named as genome_X.fasta where X determines the number of the genome in the pangenome.

   Example: 

   java  -jar  /home/sheik005/pantools/dist/pantools.jar  reconstruct  all  /home/sheik005/two_hiv_pangenome_database

7. group:
   To add homology and orthology nodes which point to a groups of homologous or orthologous genes. This functionality needs MCL to be installed on your machine.

   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  group  PATH_TO_THE_PANGENOME_DATABASE 

   Example: 

   java  -jar  /home/sheik005/pantools/dist/pantools.jar  group  /home/sheik005/two_hiv_pangenome_database

Visualization in the Neo4j browser
----------------------------------
Neo4j browser allows you to run Cypher queries and receive the results in a Tabular or a graph-representation mode. To do that, you need to download the appropriate version of Neo4j from their website. 
There you will find an article demonstrating how to use the Neo4j browser for querying, visualization and data interaction. It is quite simple. 
For example, to visualize the pangenome of two HIV strains provided as a sample data, after downloading Neo4j I needed to take this actions, on my linux machine :

1. Add the path to the Neo4j /bin directory to the path environment variable.

2. Hard-code the path to your pangenome in the configuration file ( NEO4J-DIRECTORY/conf/neo4j.conf ) by : 
   dbms.directories.data = PATH_TO_THE_PANGENOME_DATABASE

3. Start the Neo4j database server by : 
   neo4j start

4. open an internet browser and Open the URL http://localhost:7474

5. To visualize the whole pangenome of two HIV strains type this simple Cypher command:
   MATCH (n) RETURN n

6. To stop the Neo4j server type :
   neo4j stop

- Windows users could also download a Neo4j desktop application for starting and stopping a server instead of doing it on the commandline.
