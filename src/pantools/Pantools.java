/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


package pantools;

import genome.SequenceDatabase;
import index.IndexDatabase;
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
 *
 * @author sheik005
 */
public class Pantools {
    /*
     The main function
     */

    public static String PATH;
    public static String GRAPH_DATABASE_PATH = "/databases/graph.db/";
    public static String INDEX_DATABASE_PATH = "/databases/index.db/";
    public static String GENOME_DATABASE_PATH = "/databases/genome.db/";
    public static GraphDatabaseService graphDb;
    public static IndexDatabase indexDb;
    public static SequenceDatabase genomeDb;
    public static SequenceDatabase sequenceDb;
    public static int db_trsc_limit = 1000;    //   The number of transactions to be committed in batch

    /*
     There are following types of nodes:
     - pangenome   
     - genome
     - sequence
     - node
     - degenerate
     - gene
     - group
     */
    public static Label pangenome_label = DynamicLabel.label("pangenome");
    public static Label genome_label = DynamicLabel.label("genome");
    public static Label sequence_label = DynamicLabel.label("sequence");
    public static Label node_label = DynamicLabel.label("node");
    public static Label degenerate_label = DynamicLabel.label("degenerate");
    public static Label gene_label = DynamicLabel.label("gene");
    public static Label mRNA_label = DynamicLabel.label("mRNA");
    public static Label tRNA_label = DynamicLabel.label("tRNA");
    public static Label ncRNA_label = DynamicLabel.label("ncRNA");
    public static Label pgRNA_label = DynamicLabel.label("pgRNA");
    public static Label ortholog_lable = DynamicLabel.label("orthologs");
    public static Label homolog_lable = DynamicLabel.label("homologs");
    /*
     All possible relationship types between two nodes in the graph.
     */

    public static enum RelTypes implements RelationshipType {

        AA0, AA1, AA2, AA3,
        AC0, AC1, AC2, AC3,
        AG0, AG1, AG2, AG3,
        AT0, AT1, AT2, AT3,
        AN0, AN1, AN2, AN3,
        CA0, CA1, CA2, CA3,
        CC0, CC1, CC2, CC3,
        CG0, CG1, CG2, CG3,
        CT0, CT1, CT2, CT3,
        CN0, CN1, CN2, CN3,
        GA0, GA1, GA2, GA3,
        GC0, GC1, GC2, GC3,
        GG0, GG1, GG2, GG3,
        GT0, GT1, GT2, GT3,
        GN0, GN1, GN2, GN3,
        TA0, TA1, TA2, TA3,
        TC0, TC1, TC2, TC3,
        TG0, TG1, TG2, TG3,
        TT0, TT1, TT2, TT3,
        TN0, TN1, TN2, TN3,
        NA0, NA1, NA2, NA3,
        NC0, NC1, NC2, NC3,
        NG0, NG1, NG2, NG3,
        NT0, NT1, NT2, NT3,
        NN0, NN1, NN2, NN3,
        begin, // for pointing to start node of a gene
        end, // for pointing to end node of a gene
        has, // for pointing to genome and sequence nodes
        contains// for pointing to gene nodes of the group
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
                if (K < 1 || K > 256) {
                    System.out.println(" Please enter a proper K value ( 1 <= K <= 256 ).");
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
                } else {
                    System.out.println("Databases are different");
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
    /*
     To print a simple manual for the users
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
    /*
     To calculates and prints peak of memory usage of the program in mega bytes.
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

}
