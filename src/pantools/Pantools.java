/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


package pantools;

import genome.SequenceDatabase;
import index.IndexDatabase;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.lang.management.MemoryUsage;
import org.neo4j.graphdb.DynamicLabel;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.RelationshipType;
import pangenome.AnnotationLayer;
import pangenome.SequenceLayer;

/**
 * Implements the main and shared functions. 
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class Pantools {
    public static String GRAPH_DATABASE_PATH = "/databases/graph.db/";
    public static String INDEX_DATABASE_PATH = "/databases/index.db/";
    public static String GENOME_DATABASE_PATH = "/databases/genome.db/";
    public static GraphDatabaseService graphDb;
    public static IndexDatabase indexDb;
    public static SequenceDatabase genomeDb;
    public static SequenceDatabase sequenceDb;
    public static int MAX_TRANSACTION_SIZE = 100;    //   The number of transactions to be committed in batch
    public static int cores = Runtime.getRuntime().availableProcessors() / 2 + 1;

    public static Label pangenome_label = DynamicLabel.label("pangenome");
    public static Label genome_label = DynamicLabel.label("genome");
    public static Label sequence_label = DynamicLabel.label("sequence");
    public static Label node_label = DynamicLabel.label("node");
    public static Label degenerate_label = DynamicLabel.label("degenerate");
//    public static Label panproteome_label = DynamicLabel.label("panproteome");
//    public static Label proteome_label = DynamicLabel.label("proteome");
    public static Label annotation_label = DynamicLabel.label("annotation");
    public static Label variation_label = DynamicLabel.label("variation");
    public static Label gene_label = DynamicLabel.label("gene");
    public static Label coding_gene_label = DynamicLabel.label("coding_gene");
    public static Label RNA_label = DynamicLabel.label("RNA");
    public static Label mRNA_label = DynamicLabel.label("mRNA");
    public static Label tRNA_label = DynamicLabel.label("tRNA");
    public static Label rRNA_label = DynamicLabel.label("rRNA");
    public static Label CDS_label = DynamicLabel.label("CDS");
    public static Label exon_label = DynamicLabel.label("exon");
    public static Label intron_label = DynamicLabel.label("intron");
    public static Label feature_label = DynamicLabel.label("feature");
    public static Label broken_protein_label = DynamicLabel.label("broken_protein");
    public static Label orthology_group_lable = DynamicLabel.label("orthology_group");
    public static Label homology_group_lable = DynamicLabel.label("homology_group");
    public static Label kmer_lable = DynamicLabel.label("kmer");

    
    public static enum RelTypes implements RelationshipType {
        FF, FR, RF, RR,
        has, // for pointing to genome and sequence nodes
        starts,
        stops,
        has_homolog, // for pointing to gene nodes from the homology group
        has_ortholog, // for pointing to gene nodes from the orthology group
        splits_into,
        codes_for,// for connecting genes to mRNAs
        is_parent_of,
        contributes_to,// for connecting CDSs to mRNA
        branches, //to connect tree nodes
        visits,
        precedes, 
        is_homolog_to
    }

    public static long startTime;
    public static long phaseTime;
    public static int num_nodes;
    public static int num_degenerates;
    public static int num_edges;
    public static int num_bases;

    public static SequenceLayer seqLayer;
    public static AnnotationLayer annLayer;

    /**
     * The starting point of the PanTools program
     * @param args Command line arguments
     */
    public static void main(String[] args) {
        int K;
        if (args.length < 2 || args[1].equals("--help") || args[1].equals("-h")) {
            print_help_comment();
            System.exit(1);
        }
        seqLayer = new SequenceLayer();
        annLayer = new AnnotationLayer();
        System.out.println("------------------------------- PanTools -------------------------------");
        switch (args[0]) {
            case "build":
                if (args.length < 4) {
                    print_help_comment();
                    System.exit(1);
                }
                if (args[1].equals("pangenome")){
                        if (args.length > 4){
                            K = Integer.parseInt(args[4]);
                            if (K < 6 || K > 255)
                                K = -1;
                        } else
                           K = -1;
                    seqLayer.initialize_pangenome(args[3],args[2], K);
                } else if (args[1].equals("panproteome")){
                        if (args.length > 4){
                            K = Integer.parseInt(args[4]);
                            if (K < 4 || K > 6)
                                K = 5;
                        } else
                           K = 5;
                        System.out.println("Kmer size set to " + K);
                        annLayer.initialize_panproteome(args[3], args[2], K);
                } else {
                    print_help_comment();
                    System.exit(1);
                }
                break;
            case "add":
                if (args.length < 4) {
                    print_help_comment();
                    System.exit(1);
                }
                if (args[1].equals("genomes"))
                    seqLayer.add_genomes(args[3],args[2]);
                else if (args[1].equals("annotations"))
                    annLayer.add_annotaions(args[3],args[2]);
                else if (args[1].equals("variations"))
                    seqLayer.add_variations(args[3],args[2]);
                else {
                    print_help_comment();
                    System.exit(1);
                }
                break;
            case "group":
                if (args.length < 2) {
                    print_help_comment();
                    System.exit(1);
                }
                annLayer.group(args);
                break;
            case "retrieve":
                if (args.length < 4) {
                    print_help_comment();
                    System.exit(1);
                }
                if (args[1].equals("genes"))
                    seqLayer.retrieve_genes(args[3],args[2]);
                else if (args[1].equals("regions"))
                    seqLayer.retrieve_regions(args[3],args[2]);
                else if (args[1].equals("genomes"))
                    seqLayer.retrieve_genomes(args[3],args[2]);
                        else {
                    print_help_comment();
                    System.exit(1);
                }
                break;
            default:
                print_help_comment();
                System.exit(1);
        }
        System.out.println("Total time : " + (System.currentTimeMillis() - startTime) / 1000 + "." + (System.currentTimeMillis() - startTime) % 1000 + " seconds");
        print_peak_memory();
        System.out.println("-----------------------------------------------------------------------");
    }

    /**
     * Print the manual of the software.
     */
    private static void print_help_comment() {
        System.out.println("************************************************************************\n" +
"PanTools is a disk-based java application for computational pan-genomics\n" +
"developed by Siavash Sheikhizadeh et. al in Bioinformatics group of \n" +
"Wageningen university and research center, the Netherlands.  \n" +
"************************************************************************\n" +
"\n" +
"Requirements\n" +
"------------\n" +
"- KMC2: is a disk-based programm for counting k-mers from (possibly gzipped) FASTQ/FASTA files( http://sun.aei.polsl.pl/kmc ).\n" +
"        You need to download it and add the path to the appropriate version (linux, macos or windows) of kmc and kmc_tools executables to your OS path environment variable.\n" +
"\n" +
"- Java Virtual Machine version 1.8 or higher: Add the path to the java executable to your OS path environment variable.\n" +
"\n" +
"- MCL: The Markov Cluster Algorithm, is a fast and scalable unsupervised cluster algorithm for graphs ( http://micans.org/mcl ) which is needed for group functionality of PanTools.\n" +
"       You need to download, unzip and compile it (see README), and add the path to the mcl executable to your path environment variable.\n" +
"\n" +
"How to run the program \n" +
"----------------------\n" +
"java  [-server] [-XX:+UseConcMarkSweepGC]  [-Xmx(a number followed by g or m)] -jar ./pantools/dist/pantools.jar command arguments\n" +
"\n" +
"\n" +
"List of commands and examples for the provided sample data :\n" +
"\n" +
"1. build:\n" +
"   To build a pan-genome out of a set of genomes or a pan-proteome out of a set of proteins.\n" +
"\n" +
"   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  build  pangenome [or panproteome] PATH_TO_THE_DATABASE  PATH_TO_THE_GENOMES_PATH_FILE [or PATH_TO_THE_PROTEOMES_PATH_FILE]\n" +
"\n" +
"   PATH_TO_THE_GENOMES_PATH_FILE : a text file containing paths to FASTA files of genomes; each in a seperated line.\n" +
"   PATH_TO_THE_PROTEOMES_PATH_FILE : a text file containing paths to FASTA files of proteomes; each in a seperated line.\n" +
"   PATH_TO_THE_PANGENOME_DATABASE : path where the resulting pangenome is stored.  \n" +
"\n" +
"   Example: \n" +
"   \n" +
"   java  -Xmx4g  -jar  /home/sheik005/pantools/dist/pantools.jar build  pangenome /home/sheik005/two_hiv_pangenome_database  /home/sheik005/pantools/example/sample_genomes_path.txt\n" +
"             \n" +
"2. add:\n" +
"   To add new genomes and annotations to an available pan-genome. \n" +
"\n" +
"   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  add  genomes [or annotaions] PATH_TO_THE_PANGENOME_DATABASE  PATH_TO_THE_NEW_GENOMES_PATH_FILE [or PATH_TO_THE_ANNOTATION_PATH_FILE]\n" +
"   \n" +
"   PATH_TO_THE_NEW_GENOMES_PATH_FILE : a text file containing paths to FASTA files of the new genomes to be added to the pangeome; each in a seperated line.\n" +
"                                       New genomes could also be annotated later in the same way; however, there should be empty lines in the annotations path file for pre-existing genomes.\n" +
"\n" +
"   PATH_TO_THE_ANNOTATION_PATH_FILE : a text file each line of which contains genome number and path to the corresponding GFF file seperated by one space.\n" +
"\n" +
"   Example: \n" +
"\n" +
"   java  -jar  /home/sheik005/pantools/dist/pantools.jar  add annotations /home/sheik005/two_hiv_pangenome_database  /home/sheik005/pantools/example/sample_annotations_path.txt\n" +
"\n" +
"3. retrieve:\n" +
"   To retrieve the sequence of annotated genes, genomic regios or constituent genomes. \n" +
"\n" +
"   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  retrieve  genes [or regions or genomes]  PATH_TO_THE_PANGENOME_DATABASE  PATH_TO_THE_ANNOTATION_RECORDS_FILE [or PATH_TO_THE_GENOMIC_REGIONS_FILE or PATH_TO_THE_GENOME_NUMBERS_FILE]\n" +
"\n" +
"   PATH_TO_THE_ANNOTATION_RECORDS_FILE : a text file containing records of annotated genes, as they appear in GFF file, to be retrieved.\n" +
"                                         The resulting FASTA file would have the same name as the PATH_TO_THE_ANNOTATION_RECORDS_FILE with an additional .fasta extention.\n" +
"\n" +
"   PATH_TO_THE_GENOMIC_REGIONS_FILE : a text file containing records with genome_number, sequence_number, begin and end positions seperated by one space for each region.\n" +
"                                      The resulting FASTA file would have the same name as the PATH_TO_THE_GENOMIC_REGIONS_FILE with an additional .fasta extention.\n" +
"\n" +
"   PATH_TO_THE_GENOME_NUMBERS_FILE : a text file containing genome_numbers to be retrieved in each line. The resulting FASTA files are named as genome_X.fasta where X determines the number of the genome in the pangenome.\n" +
"\n" +
"   Examples: \n" +
"\n" +
"   java  -jar  /home/sheik005/pantools/dist/pantools.jar  retrieve  genes  /home/sheik005/two_hiv_pangenome_database  /home/sheik005/pantools/example/sample_annotation_records.txt\n" +
"   java  -jar  /home/sheik005/pantools/dist/pantools.jar  retrieve  regions  /home/sheik005/two_hiv_pangenome_database  /home/sheik005/pantools/example/sample_genomic_regions.txt\n" +
"   java  -jar  /home/sheik005/pantools/dist/pantools.jar  retrieve  genomes  /home/sheik005/two_hiv_pangenome_database  /home/sheik005/pantools/example/sample_genome_numbers.txt\n" +
"\n" +
"4. group:\n" +
"   To add homology and orthology nodes which point to a groups of homologous or orthologous genes. This functionality needs MCL to be installed on your machine.\n" +
"\n" +
"   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  group  PATH_TO_THE_PANGENOME_DATABASE \n" +
"\n" +
"   Example: \n" +
"\n" +
"   java  -jar  /home/sheik005/pantools/dist/pantools.jar  group  /home/sheik005/two_hiv_pangenome_database\n" +
"\n" +
"Visualization in the Neo4j browser\n" +
"----------------------------------\n" +
"Neo4j browser allows you to run Cypher queries and receive the results in a Tabular or a graph-representation mode. To do that, you need to download the appropriate version of Neo4j from their website. \n" +
"There you will find an article demonstrating how to use the Neo4j browser for querying, visualization and data interaction. It is quite simple. \n" +
"For example, to visualize the pangenome of two HIV strains provided as a sample data, after downloading Neo4j I needed to take this actions, on my linux machine :\n" +
"\n" +
"1. Add the path to the Neo4j /bin directory to the path environment variable.\n" +
"\n" +
"2. Hard-code the path to your pangenome in the configuration file ( NEO4J-DIRECTORY/conf/neo4j.conf ) by : \n" +
"   dbms.directories.data = PATH_TO_THE_PANGENOME_DATABASE\n" +
"\n" +
"3. Start the Neo4j database server by : \n" +
"   neo4j start\n" +
"\n" +
"4. open an internet browser and Open the URL http://localhost:7474\n" +
"\n" +
"5. To visualize the whole pangenome of two HIV strains type this simple Cypher command:\n" +
"   MATCH (n) RETURN n\n" +
"\n" +
"6. To stop the Neo4j server type :\n" +
"   neo4j stop\n" +
"\n" +
"- Windows users could also download a Neo4j desktop application for starting and stopping a server instead of doing it on the commandline.\n" +
"");
    }

    /**
     * Estimates and prints the peak memory used during the execution of the program. 
     */
    public static void print_peak_memory() {
        long memoryUsage = 0;
        try {
            for (MemoryPoolMXBean pool : ManagementFactory.getMemoryPoolMXBeans()) {
                MemoryUsage peak = pool.getPeakUsage();
                memoryUsage += peak.getUsed();
            }
            System.out.println("Peak memory : " + memoryUsage / 1024 / 1024 + " MB");
        } catch (Throwable t) {
            System.err.println("Exception in agent: " + t);
        }
    }
    
    /**
     * Writes a sequence in a FASTA file with specified length for lines.
     * 
     * @param fasta_file The FASTA file object
     * @param seq The sequence to be written in the file.
     * @param length Length of the lines.
     */    
    public static void write_fasta(BufferedWriter fasta_file, StringBuilder seq, int length) {
        int i;
        try {
            for (i = 1; i <= seq.length(); ++i) {
                fasta_file.write(seq.charAt(i - 1));
                if (i % length == 0) {
                    fasta_file.write("\n");
                }
            }
            fasta_file.write("\n");
        } catch (IOException ioe) {

        }

    }    
    
    /**
     * Return reverse complement of the given string.
     * 
     * @param s The input string
     * @return The reverse complement of the input string
     */     
    public static void reverse_complement(StringBuilder s) {
        char ch;
        int i, j;
        for ( i=0, j = s.length() - 1; i < j; ++i, --j) {
            ch = s.charAt(i);
            s.setCharAt(i, complement(s.charAt(j)));
            s.setCharAt(j, complement(ch));
        }
        if (i == j)
            s.setCharAt(i, complement(s.charAt(i)));
    }

    public static char complement(char ch) {
        switch (ch) {
            case 'A':
                return 'T';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            case 'T':
                return 'A';
            case 'R':
                return 'Y';
            case 'Y':
                return 'R';
            case 'K':
                return 'M';
            case 'M':
                return 'K';
            case 'B':
                return 'V';
            case 'V':
                return 'B';
            case 'D':
                return 'H';
            case 'H':
                return 'D';
            default:
                return ch;
        }
    }
    
    
    /**
     * Executes a shell command. 
     * @param command The command
     * @return The output of the bash command
     */
    public static String executeCommand(String command) {
        StringBuilder exe_output = new StringBuilder();
        String line = "";
        Process p;
        try {
            p = Runtime.getRuntime().exec(command);
            p.waitFor();
            BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
            while ((line = reader.readLine()) != null) {
                exe_output.append(line + "\n");
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        return exe_output.toString();
    }    
}
