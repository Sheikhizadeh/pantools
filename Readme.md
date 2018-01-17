****************************************************************
PanTools version 1.1,

is a java application based on Neo4j graph database community 
edition 3.3.1 for computational pan-genomics, developed by 
Siavash Sheikhizadeh in Wageningen university, the Netherlands.
If you use PanTools please do not forget to cite it :

doi: 10.1093/bioinformatics/btw455

https://github.com/Sheikhizadeh/pantools  

****************************************************************

Requirements
------------
- KMC: is a disk-based programm for counting k-mers from 
       (possibly gzipped) FASTQ/FASTA files
       (http://sun.aei.polsl.pl/kmc).
        You need to download it and add the path to the 
        appropriate version (linux, macos or windows) of kmc 
        and kmc_tools executables to your OS path environment 
        variable.

- Java Virtual Machine version 1.8 or higher: Add the path to 
       the java executable to your OS path environment variable.

- MCL: The Markov Cluster Algorithm, is a fast and scalable 
       unsupervised cluster algorithm for graphs 
       (http://micans.org/mcl ) which is needed for group 
       functionality of PanTools.
       You need to download, unzip and compile it (see README), 
       and add the path to the mcl executable to your path
       environment variable.

Running the program 
-------------------
java <JVM options> -jar pantools.jar <command> <arguments>

pantools.jar is available in folder pantools/dist/ 

JVM options
-----------
[-server] 
[-XX:+UseConcMarkSweepGC]  
[-Xmx(a number followed by g or m)]

PanTools commands
-----------------

<build pangenome>
   To build a pan-genome out of a set of genomes.
   <arguments>
   -d PATH_TO_THE_PANGENOME_DATABASE : 
      path to the pangenome database. 
   -g PATH_TO_THE_GENOMES_FILE : 
      a text file containing paths to FASTA files of genomes;
      each in a seperated line.
   -k K_SIZE : 
      If it is not given or is out of range ( 6 <= K_SIZE <= 255 ), 
      an optimal value would be calculated automatically.    

<build panproteome>
   To build a pan-proteome out of a set of proteins.
   <arguments>
   -d PATH_TO_THE_PANGENOME_DATABASE : 
      path to the pangenome database. 
   -p PATH_TO_THE_PROTEOMES_FILE : 
      a text file containing paths to FASTA files of proteomes; 
      each in a seperated line.
             
<add genomes>
   To add new genomes to an available pan-genome.    
   <arguments>
   -d PATH_TO_THE_PANGENOME_DATABASE : 
      path to the pangenome database. 
   -g PATH_TO_THE_GENOMES_FILE : 
      a text file containing paths to FASTA files of the new 
      genomes to be added to the pangeome; 
      each in a seperated line.

<add annotations>
   To add new annotations to an available pan-genome. 
   <arguments>
   -d PATH_TO_THE_PANGENOME_DATABASE : 
      path to the pangenome database. 
   -a PATH_TO_THE_ANNOTATIONS_FILE : 
      a text file each line of which contains genome number and 
      path to the corresponding GFF file seperated by one space.
      Genomes are numbered in the same order they have been added
      to the pangenome. The protein sequence of the annotated genes 
      will be also stored in the folder "proteins" in the same path 
      as the pangenome. 

<retrieve genes>
   To retrieve the sequence of annotated genes from the pangenome. 
   The results will be stored in the same folder as the pangenome.
   <arguments>
   -d PATH_TO_THE_PANGENOME_DATABASE : 
      path to the pangenome database. 
   -e PATH_TO_THE_GENE_RECORDS : 
      a text file containing records of annotated genes, 
      as they appear in GFF file, to be retrieved. The resulting 
      FASTA file would have the same name with an additional 
      .fasta extention.

<retrieve regions> 
   To retrieve the sequence of some genomic regios from the pangenome. 
   The results will be stored in the same folder as the pangenome.
   <arguments>
   -d PATH_TO_THE_PANGENOME_DATABASE : 
      path to the pangenome database. 
   -r PATH_TO_THE_REGIONS_FILE : 
      a text file containing records with genome_number, 
      sequence_number, begin and end positions seperated by one 
      space for each region. The resulting FASTA file would have 
      the same name with an additional .fasta extention.
<retrieve genomes>
   To retrieve the full sequence of some genomes. The results will be 
   stored in the same folder as the pangenome itself.
   <arguments>
   -d PATH_TO_THE_PANGENOME_DATABASE : 
      path to the pangenome database. 
   -n PATH_TO_THE_GENOME_NUMBERS_FILE : 
      a text file containing genome_numbers to be retrieved in each line. 
      The resulting FASTA files are named like Genome_x.fasta.

<group>
   To add homology nodes which point to a groups of homologous proteins.
   <arguments>
   -d PATH_TO_THE_PANGENOME_DATABASE : 
      path to the pangenome database. 
   -i INTERSECTION_RATE (default = 0.09): 
      determines the fraction of kmers needs to be shared by two 
      intersecting proteins. Should be in range [0.001, 0.1].
   -t THRESHOLD (default = 95): 
      the minimum similarity score. Should be in range [1-99]. 
   -m MCL_INFLATION (default = 9.6): 
      the MCL inflation. Should be in range ]1-19[.
   -c CONTRAST (default = 8): 
      the contrast factor. Should be in range ]0-10[.
   -r RELAXATION (default 1): 
      the relaxation about homology. Sould be in range [1, 8], 
      from strict to relaxed.

<version>
   To show the versions of PanTools and Neo4j.
   
Visualization in the Neo4j browser
----------------------------------
   Neo4j browser allows you to run Cypher queries and receive 
   the results in a tabular or a graph representation mode. 
   You need to download the appropriate version of Neo4j. 
   To visualize the pangenome of two HIV strains provided 
   as a sample data in pantools repositiory, take these actions 
   on a linux machine. Windows users could also download the
   Neo4j desktop application for starting and stopping a server 
   instead of usingn commandline.
1. Add the path to the Neo4j /bin directory to the path 
   environment variable.
2. Hard-code the path to your pangenome in the configuration file 
   ,NEO4J-DIRECTORY/conf/neo4j.conf, by: 
   dbms.directories.data = PATH_TO_THE_PANGENOME_DATABASE
3. Start the Neo4j database server by: 
   neo4j start
4. open an internet browser and Open the URL http://localhost:7474
5. To visualize the whole pangenome of two HIV strains, 
   type this simple Cypher command:
   MATCH (n) RETURN n
6. To stop the Neo4j server type:
   neo4j stop
