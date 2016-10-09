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
    public static String PATH;
    public static String GRAPH_DATABASE_PATH = "/databases/graph.db/";
    public static String INDEX_DATABASE_PATH = "/databases/index.db/";
    public static String GENOME_DATABASE_PATH = "/databases/genome.db/";
    public static GraphDatabaseService graphDb;
    public static IndexDatabase indexDb;
    public static SequenceDatabase genomeDb;
    public static SequenceDatabase sequenceDb;
    public static int MAX_TRANSACTION_SIZE = 1000;    //   The number of transactions to be committed in batch

    public static Label pangenome_label = DynamicLabel.label("pangenome");
    public static Label genome_label = DynamicLabel.label("genome");
    public static Label sequence_label = DynamicLabel.label("sequence");
    public static Label node_label = DynamicLabel.label("node");
    public static Label degenerate_label = DynamicLabel.label("degenerate");
    public static Label gene_label = DynamicLabel.label("gene");
    public static Label pseudogene_label = DynamicLabel.label("pseudogene");
    public static Label TEgene_label = DynamicLabel.label("TEgene");
    public static Label coding_gene_label = DynamicLabel.label("coding_gene");
    public static Label noncoding_gene_label = DynamicLabel.label("noncoding_gene");
    public static Label tRNA_gene_label = DynamicLabel.label("tRNA_gene");
    public static Label mRNA_label = DynamicLabel.label("mRNA");
    public static Label tRNA_label = DynamicLabel.label("tRNA");
    public static Label ncRNA_label = DynamicLabel.label("ncRNA");
    public static Label pseudogenic_transcript_label = DynamicLabel.label("pseudogenic_transcript_label");
    public static Label CDS_label = DynamicLabel.label("CDS");
    public static Label ortholog_lable = DynamicLabel.label("orthologs");
    public static Label homolog_lable = DynamicLabel.label("homologs");

    public static enum RelTypes implements RelationshipType {
        FF, FR, RF, RR,
        has, // for pointing to genome and sequence nodes
        visits, // for connecting genes to the nodes
        contains, // for pointing to gene nodes of the group
        codes_for,// for connecting genes to mRNAs
        contributes_to,// for connecting CDSs and LTRs to mRNA
        is_a, // for connecting genes to tRNAs, ncRNAs and pseusogenic_transcripts
        covers //to connect CDSs to the nodes
    }

    public static long startTime;
    public static long phaseTime;
    public static int K;
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
        if (args.length < 2 || args[1].equals("--help") || args[1].equals("-h")) {
            print_help_comment();
            System.exit(1);
        }
        seqLayer = new SequenceLayer();
        annLayer = new AnnotationLayer();
        System.out.println("------------------------------- PanTools -------------------------------");
        switch (args[0]) {
            case "reconstruct":
                PATH = args[2];
                seqLayer.reconstruct_genomes(args[1]);
                break;
            case "build":
                K = Integer.parseInt(args[1]);
                if (K < 6 || K > 256) {
                    System.out.println("Please enter a proper K value ( 6 <= K <= 256 ).");
                    System.exit(1);
                }
                PATH = args[2];
                seqLayer.build(args[3]);
                break;
            case "add":
                PATH = args[1];
                seqLayer.add(args[2]);
                break;
            case "annotate":
                PATH = args[1];
                annLayer.annotate(args[2]);
                break;
            case "group":
                PATH = args[2];
                if (args[1].equals("denovo"))
                    annLayer.denovo_homology_annotation();
                else
                    annLayer.group_ortholog_proteins(args[1]);
                break;
            case "compare":
                if (seqLayer.compare_pangenomes(args[1], args[2])) {
                    System.out.println("Databases are equal.");
                    System.out.println("....................");
                } else {
                    System.out.println("Databases are different");
                    System.out.println(".......................");
                }
                break;
            case "retrieve":
                PATH = args[2];
                if (args[1].equals("genes")) {
                    seqLayer.retrieve_genes(args[3]);
                } else if (args[1].equals("regions")) {
                    seqLayer.retrieve_regions(args[3]);
                } else {
                    print_help_comment();
                    System.exit(1);
                }
                break;
            case "query":
                PATH = args[1];
                seqLayer.run_query();
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
"        You need to unzipd kmc.zip file provided in this release of PanTools and add the path of the corresponding version (linux, macos or windows) of kmc and kmc_tools executables to your OS path environment variable.\n" +
"\n" +
"- Java Virtual Machine version 1.7 or higher: Add the path to the java executable to your OS path environment variable.\n" +
"\n" +
"\n" +
"How to run the program \n" +
"----------------------\n" +
"java  [-server] [-XX:+UseConcMarkSweepGC]  [-Xmx(a number followed by g or m)] -jar ./pantools/dist/pantools.jar command arguments\n" +
"\n" +
"\n" +
"List of commands and examples for the provided sample data :\n" +
"\n" +
"1. build:\n" +
"   To build a pan-genome out of a set of genomes.\n" +
"\n" +
"   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  build  K_VALUE  PATH_TO_THE_PANGENOME_DATABASE  PATH_TO_THE_GENOMES_PATH_FILE\n" +
"\n" +
"   K_VALUE :   size of K for construction of the de Bruijn graph.\n" +
"   PATH_TO_THE_GENOMES_PATH_FILE : a text file containing paths to FASTA files (genomes); each in a seperated line.\n" +
"   PATH_TO_THE_PANGENOME_DATABASE : path where the resulting pangenome is stored.  \n" +
"\n" +
"   Example: \n" +
"   \n" +
"   java  -Xmx4g  -jar  /home/pantools/dist/pantools.jar build  31  /home/two_hiv_pangenome_database  /home/pantools/example/sample_genomes_path.txt\n" +
"             \n" +
"2. annotate:\n" +
"   To add annotations to a pan-genome. This function also produce a FASTA file containing all the protein sequences in the same order as they have been annotated in the GFF file.\n" +
"\n" +
"   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  annotate  PATH_TO_THE_PANGENOME_DATABASE  PATH_TO_THE_ANNOTATION_PATH_FILE\n" +
"\n" +
"   PATH_TO_THE_ANNOTATION_PATH_FILE : a text file containing paths to the GFF files corresponding to the genomes in the same order as they apear in PATH_TO_THE_GENOMES_PATH_FILE.\n" +
"                                      Missing annotations are indicated by an empty line.\n" +
"\n" +
"   Example: \n" +
"\n" +
"   java  -jar  /home/pantools/dist/pantools.jar  annotate  /home/two_hiv_pangenome_database  /home/pantools/example/sample_annotations_path.txt\n" +
"\n" +
"3. add:\n" +
"   To add new genomes to an available pan-genome.\n" +
"\n" +
"   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  add  PATH_TO_THE_PANGENOME_DATABASE  PATH_TO_THE_NEW_GENOMES_PATH_FILE\n" +
"   \n" +
"   PATH_TO_THE_NEW_GENOMES_PATH_FILE : a text file containing paths to FASTA files of the new genomes to be added to the pangeome; each in a seperated line.\n" +
"                                       New genomes could also be annotated later in the same way; however, there should be an empty lines in the annotations path file for existing genomes.\n" +
"\n" +
"4. retrieve genes:\n" +
"   To extract sequence of some annotated genes. \n" +
"\n" +
"   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  retrieve  genes  PATH_TO_THE_PANGENOME_DATABASE  PATH_TO_THE_ANNOTATION_RECORDS_FILE\n" +
"\n" +
"   PATH_TO_THE_ANNOTATION_RECORDS_FILE : a text file containing records of annotated genes, as they appear in GFF file, to be retrieved.\n" +
"                                         The resulting FASTA file would have the same name as the PATH_TO_THE_ANNOTATION_RECORDS_FILE with an additional .fasta extention.\n" +
"\n" +
"   Example: \n" +
"\n" +
"   java  -jar  /home/pantools/dist/pantools.jar  retrieve  genes  /home/two_hiv_pangenome_database  /home/pantools/example/sample_annotaion_records.txt\n" +
"\n" +
"5. retrieve regions:\n" +
"   To extract sequence of some genomic regios.\n" +
"\n" +
"   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  retrieve  regions  PATH_TO_THE_PANGENOME_DATABASE  PATH_TO_THE_GENOMIC_REGIONS_FILE\n" +
"\n" +
"   PATH_TO_THE_GENOMIC_REGIONS_FILE : a text file containing records with genome_number, sequence_number, begin and end positions seperated by one space for each region.\n" +
"                                      The resulting FASTA file would have the same name as the PATH_TO_THE_GENOMIC_REGIONS_FILE with an additional .fasta extention.\n" +
"\n" +
"   Example: \n" +
"\n" +
"   java  -jar  /home/pantools/dist/pantools.jar  retrieve  regions  /home/two_hiv_pangenome_database  /home/pantools/example/sample_genomic_regions.txt\n" +
"\n" +
"6. reconstruct:\n" +
"   To reconstruct all or a set of genomes out of the pan-genome.\n" +
"\n" +
"   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  reconstruct all [or PATH_TO_THE_GENOME_NAMES_FILE]  PATH_TO_THE_PANGENOME_DATABASE\n" +
"\n" +
"   PATH_TO_THE_GENOME_NAMES_FILE : a text file containing genome_number and a given name for that genome seperated by a single space in each line. \n" +
"                                   The resulting FASTA files are named as genome_X.fasta where X determines the number of the genome in the pangenome.\n" +
"\n" +
"   Example: \n" +
"\n" +
"   java  -jar  /home/pantools/dist/pantools.jar  reconstruct  all  /home/two_hiv_pangenome_database\n" +
"\n" +
"7. group:\n" +
"   To group genes by adding group nodes pointing to them, either de novo (based on thier similarity) or using a group file.\n" +
"\n" +
"   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  group  denovo [or PATH_TO_THE_GROUP_FILE] PATH_TO_THE_PANGENOME_DATABASE \n" +
"\n" +
"   PATH_TO_THE_GROUP_FILE : a text file each line stars with name of the group follewed by space-seperated name of the group members (proteins) in each line.\n" +
"\n" +
"   Example: \n" +
"\n" +
"   java  -jar  /home/pantools/dist/pantools.jar  group  /home/two_hiv_pangenome_database  /home/pantools/example/sample_orthologous_groups.txt\n" +
"\n" +
"\n" +
"8. compare:\n" +
"   To compare topology of two pan-genomes.\n" +
"\n" +
"   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  compare  PATH_TO_THE_PANGENOME_DATABASE_1  PATH_TO_THE_PANGENOME_DATABASE_2\n" +
"\n" +
"9. query\n" +
"   To run Cypher queries.\n" +
"\n" +
"   java  -jar  PATH_TO_THE_JAR_FILE/pantools.jar  query  PATH_TO_THE_PANGENOME_DATABASE\n" +
"   \n" +
"\n" +
"Visualization in the Neo4j browser\n" +
"----------------------------------\n" +
"Neo4j browser allows you to run Cypher queries and receive the results in a Tabular or a graph-representation mode. To do that, you need to download the appropriate version of Neo4j from their website. \n" +
"There you will find an article demonstrating how to use the Neo4j browser for querying, visualization and data interaction. It is quite simple. \n" +
"For example, to visualize the pangenome of two HIV strains provided as a sample data, after download I needed to take this actions, on my linux machine :\n" +
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
"   \n" +
"   \n" +
"\n" +
" \n" +
""
                
);
    }

    /**
     * Estimates and prints the peak memory used during the execution of the program. 
     */
    private static void print_peak_memory() {
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
    public static void write_fasta(BufferedWriter fasta_file, String seq, int length) {
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
     * @param s    The input string
     * @return 
     */     
    public static String reverse_complement(String s) {
        StringBuilder rv = new StringBuilder();
        for (int i = s.length() - 1; i >= 0; --i) {
            switch (s.charAt(i)) {
                case 'A':
                    rv.append('T');
                    break;
                case 'C':
                    rv.append('G');
                    break;
                case 'G':
                    rv.append('C');
                    break;
                case 'T':
                    rv.append('A');
                    break;
                case 'R':
                    rv.append('Y');
                    break;
                case 'Y':
                    rv.append('R');
                    break;
                case 'K':
                    rv.append('M');
                    break;
                case 'M':
                    rv.append('K');
                    break;
                case 'B':
                    rv.append('V');
                    break;
                case 'V':
                    rv.append('B');
                    break;
                case 'D':
                    rv.append('H');
                    break;
                case 'H':
                    rv.append('D');
                    break;
                default:
                    rv.append(s.charAt(i));
            }
        }
        return rv.toString();
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
