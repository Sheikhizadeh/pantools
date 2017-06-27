************************************************************************
PanTools is a disk-based java application for computational pan-genomics
developed by Siavash Sheikhizadeh et. al in Bioinformatics group of 
Wageningen university and research center, the Netherlands.  
************************************************************************

Requirements
------------
- KMC2: is a disk-based programm for counting k-mers from (possibly gzipped) FASTQ/FASTA files( http://sun.aei.polsl.pl/kmc ).
        You need to download it and add the path to the appropriate version (linux, macos or windows) of kmc and kmc_tools executables to your OS path environment variable.

- Java Virtual Machine version 1.8 or higher: Add the path to the java executable to your OS path environment variable.

- MCL: The Markov Cluster Algorithm, is a fast and scalable unsupervised cluster algorithm for graphs ( http://micans.org/mcl ) which is needed for group functionality of PanTools.
       You need to download, unzip and compile it (see README), and add the path to the mcl executable to your path environment variable.

How to run the program 
----------------------
java  [-server] [-XX:+UseConcMarkSweepGC]  [-Xmx(a number followed by g or m)] -jar ./pantools/dist/pantools.jar command arguments


List of commands and examples for the provided sample data :

1. build:
   To build a pan-genome out of a set of genomes or a pan-proteome out of a set of proteins.

   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  build  pangenome [or panproteome] PATH_TO_THE_DATABASE  PATH_TO_THE_GENOMES_PATH_FILE [or PATH_TO_THE_PROTEOMES_PATH_FILE] [K_SIZE]

   PATH_TO_THE_GENOMES_PATH_FILE : a text file containing paths to FASTA files of genomes; each in a seperated line.
   PATH_TO_THE_PROTEOMES_PATH_FILE : a text file containing paths to FASTA files of proteomes; each in a seperated line.
   PATH_TO_THE_PANGENOME_DATABASE : path where the resulting pangenome is stored. 
   K_SIZE : If it is not given or is out of range ( 6 <= K_SIZE <= 255 ), an optimal value would be calculated automatically.    

   Example: 
   
   java  -Xmx4g  -jar  /home/sheik005/pantools/dist/pantools.jar build  pangenome /home/sheik005/two_hiv_pangenome_database  /home/sheik005/pantools/example/sample_genomes_path.txt
             
2. add:
   To add new genomes and annotations to an available pan-genome. 

   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  add  genomes [or annotaions] PATH_TO_THE_PANGENOME_DATABASE  PATH_TO_THE_NEW_GENOMES_PATH_FILE [or PATH_TO_THE_ANNOTATION_PATH_FILE]
   
   PATH_TO_THE_NEW_GENOMES_PATH_FILE : a text file containing paths to FASTA files of the new genomes to be added to the pangeome; each in a seperated line.
                                       New genomes could also be annotated later in the same way; however, there should be empty lines in the annotations path file for pre-existing genomes.

   PATH_TO_THE_ANNOTATION_PATH_FILE : a text file each line of which contains genome number and path to the corresponding GFF file seperated by one space.

   Example: 

   java  -jar  /home/sheik005/pantools/dist/pantools.jar  add annotations /home/sheik005/two_hiv_pangenome_database  /home/sheik005/pantools/example/sample_annotations_path.txt

3. retrieve:
   To retrieve the sequence of annotated genes, genomic regios or constituent genomes. 

   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  retrieve  genes [or regions or genomes]  PATH_TO_THE_PANGENOME_DATABASE  PATH_TO_THE_ANNOTATION_RECORDS_FILE [or PATH_TO_THE_GENOMIC_REGIONS_FILE or PATH_TO_THE_GENOME_NUMBERS_FILE]

   PATH_TO_THE_ANNOTATION_RECORDS_FILE : a text file containing records of annotated genes, as they appear in GFF file, to be retrieved.
                                         The resulting FASTA file would have the same name as the PATH_TO_THE_ANNOTATION_RECORDS_FILE with an additional .fasta extention.

   PATH_TO_THE_GENOMIC_REGIONS_FILE : a text file containing records with genome_number, sequence_number, begin and end positions seperated by one space for each region.
                                      The resulting FASTA file would have the same name as the PATH_TO_THE_GENOMIC_REGIONS_FILE with an additional .fasta extention.

   PATH_TO_THE_GENOME_NUMBERS_FILE : a text file containing genome_numbers to be retrieved in each line. The resulting FASTA files are named as genome_X.fasta where X determines the number of the genome in the pangenome.

   Examples: 

   java  -jar  /home/sheik005/pantools/dist/pantools.jar  retrieve  genes  /home/sheik005/two_hiv_pangenome_database  /home/sheik005/pantools/example/sample_annotation_records.txt
   java  -jar  /home/sheik005/pantools/dist/pantools.jar  retrieve  regions  /home/sheik005/two_hiv_pangenome_database  /home/sheik005/pantools/example/sample_genomic_regions.txt
   java  -jar  /home/sheik005/pantools/dist/pantools.jar  retrieve  genomes  /home/sheik005/two_hiv_pangenome_database  /home/sheik005/pantools/example/sample_genome_numbers.txt

4. group:
   To add homology and orthology nodes which point to a groups of homologous or orthologous genes.

   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  group  PATH_TO_THE_PANGENOME_DATABASE THRESHOLD [K_SIZE] 

   THRESHOLD : The minimum similarity needed to call two proteins homologous. ( 0 <= THRESHOLD <= 100) 
   K_SIZE : If it is not given or is out of range ( 4 <= K_SIZE <= 6 ), it would be set to 5.    
   
Example: 

   java  -jar  /home/sheik005/pantools/dist/pantools.jar  group  /home/sheik005/two_hiv_pangenome_database 75 5

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
