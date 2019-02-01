/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


package pantools;

import index.IndexDatabase;
import index.IndexScanner;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.lang.management.MemoryUsage;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import org.apache.commons.compress.compressors.CompressorInputStream;
import org.apache.commons.compress.compressors.CompressorStreamFactory;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.RelationshipType;
import pangenome.AnnotationLayer;
import pangenome.ProteomeLayer;
import pangenome.GenomeLayer;
import sequence.SequenceDatabase;
import sequence.SequenceScanner;

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
    public static String OUTPUT_PATH = "";

    public static GraphDatabaseService graphDb;
    public static IndexDatabase indexDb;
    public static IndexScanner indexSc;    
    public static SequenceDatabase genomeDb;
    public static SequenceScanner genomeSc;    

    public static String PATH_TO_THE_PANGENOME_DATABASE;
    public static String PATH_TO_THE_GENOMES_FILE;
    public static String PATH_TO_THE_PROTEOMES_FILE;
    public static String PATH_TO_THE_ANNOTATIONS_FILE;
    public static String PATH_TO_THE_REGIONS_FILE;
    public static String PATH_TO_THE_GENOME_NUMBERS_FILE;
    public static String RAW_ABUNDANCE_FILE = "";
    public static String PATH_TO_THE_FIRST_SRA;
    public static String PATH_TO_THE_SECOND_SRA;
    public static String FEATURE = "gene";
    
    public static boolean CONNECT_ANNOTATIONS = false;
    
    public static double INTERSECTION_RATE = 0.08;
    public static double CONTRAST = 8;
    public static double MCL_INFLATION = 10.8;
    public static int MIN_PROTEIN_IDENTITY = 95;

    public static int K_SIZE = -1;
    public static int MAX_ALIGNMENT_LENGTH = 2000;
    public static int GAP_OPEN = -20;
    public static int GAP_EXT = -3;
    public static int ANCHORS_DISTANCE = 10000; // The distance between two anchor nodes
    public static int MAX_TRANSACTION_SIZE = 100;    //   The number of transactions to be committed in batch
    public static int cores = Runtime.getRuntime().availableProcessors();
    public static long heapSize = Runtime.getRuntime().maxMemory();
    public static boolean DEBUG;
    public static boolean SHOW_KMERS;
    public static int THREADS = 1;
    
    public static double MIN_IDENTITY = 0.5;
    public static int NUM_KMER_SAMPLES = 15;
    public static int MAX_NUM_LOCATIONS = 15;
    public static int MIN_HIT_LENGTH = 13;
    public static int MAX_FRAGMENT_LENGTH = 5000;
    public static int ALIGNMENT_BOUND = 5;    
    public static int SHOULDER = 100;    
    public static int ALIGNMENT_MODE = 2; // 0: all-hits    
                                          // -1: pan-genomic unique_best
                                          // -2: pan-genomic random_best
                                          // -3: pan-genomic all_bests
                                          // 1: genomic unique_best
                                          // 2: genomic random_best
                                          // 3: genomic all_bests
    public static int CLIPPING_STRINGENCY = 1; // 0: no-clipping
                                               // 1: low
                                               // 2: medium
                                               // 3: high    
    public static boolean BAMFORMAT = false;
    public static boolean INTERLEAVED = false;
    
    public static Label pangenome_label = Label.label("pangenome");
    public static Label genome_label = Label.label("genome");
    public static Label sequence_label = Label.label("sequence");
    public static Label nucleotide_label = Label.label("nucleotide");
    public static Label degenerate_label = Label.label("degenerate");
    public static Label annotation_label = Label.label("annotation");
    public static Label variation_label = Label.label("variation");
    public static Label gene_label = Label.label("gene");
    public static Label coding_gene_label = Label.label("coding_gene");
    public static Label mRNA_label = Label.label("mRNA");
    public static Label tRNA_label = Label.label("tRNA");
    public static Label rRNA_label = Label.label("rRNA");
    public static Label CDS_label = Label.label("CDS");
    public static Label exon_label = Label.label("exon");
    public static Label intron_label = Label.label("intron");
    public static Label feature_label = Label.label("feature");
    public static Label homology_group_label = Label.label("homology_group");
    
    public static enum RelTypes implements RelationshipType {
        FF, FR, RF, RR,
        has, // for pointing to genome and sequence nodes
        starts,
        stops,
        has_homolog, // for pointing to gene nodes from the homology group
        codes_for,// for connecting genes to mRNAs
        is_parent_of,
        contributes_to,// for connecting CDSs to mRNA
        is_similar_to,
        annotates,
        varies
    }

    public static char[] sym = new char[]{'A', 'C', 'G', 'T', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N'};
    public static int[] complement = new int[]{3, 2, 1, 0, 9, 8, 6, 7, 5, 4, 13, 12, 11, 10, 14};
    public static int[] binary = new int[256];
   
    public static long startTime;
    public static long phaseTime;
    public static long num_nodes;
    public static int num_degenerates;
    public static long num_edges;
    public static long num_bases;
    public static Node db_node;
    private static String[] label_strings;
    public static Map<String,Label> labels;

    public static GenomeLayer seqLayer;
    public static AnnotationLayer annLayer;
    public static ProteomeLayer proLayer;

    /**
     * The starting point of the PanTools program
     * @param args Command line arguments
     */
    public static void main(String[] args) {
        int x, i;
        double y;
        File theDir;
        if (args.length < 1) {
            print_help_message();
            System.exit(1);
        }
        binary['A'] = 0;
        binary['C'] = 1;
        binary['G'] = 2;
        binary['T'] = 3;
        binary['M'] = 4;
        binary['R'] = 5;
        binary['W'] = 6;
        binary['S'] = 7;
        binary['Y'] = 8;
        binary['K'] = 9;
        binary['V'] = 10;
        binary['H'] = 11;
        binary['D'] = 12;
        binary['B'] = 13;
        binary['N'] = 14; 
        seqLayer = new GenomeLayer();
        annLayer = new AnnotationLayer();
        proLayer = new ProteomeLayer();
        labels = new HashMap<String,Label>();
        label_strings = new String[]{
        "pangenome", "genome","sequence","nucleotide","degenerate",
        "annotation","variation","gene","coding_gene", "mRNA", 
        "tRNA", "rRNA", "CDS", "exon", "intron", "feature", 
        "broken_protein", "homology_group", "low_complexity"};        
        for (i = 0; i < label_strings.length; ++i)
            labels.put(label_strings[i], Label.label(label_strings[i]));
        System.out.println("\n------------------------------- PanTools ------------------------------");
        try{
            for (i = 1; i < args.length; i += 2){
                switch (args[i]){
                    case "--kmer-size": case "-ks":
                        x = Integer.parseInt(args[i + 1]);
                        if (x >= 6 && x <= 255)
                            K_SIZE = x;
                        else {
                            System.out.println("Choose K in the range [6..255] or do not specify it to be calculated automatically.");
                            System.exit(1);
                        }
                        System.out.println("K_SIZE = " + K_SIZE);
                        break;
                    case "--db-path": case "-dp":
                        PATH_TO_THE_PANGENOME_DATABASE = args[i + 1];
                        System.out.println("PATH_TO_THE_PANGENOME_DATABASE = " + PATH_TO_THE_PANGENOME_DATABASE);
                        break;
                    case "--out-path": case "-op":
                        OUTPUT_PATH = args[i + 1];
                        theDir = new File(OUTPUT_PATH);
                        if (!theDir.exists()) {
                                System.out.println(OUTPUT_PATH + " does not exist!");
                                System.exit(1);
                        } else
                            System.out.println("OUTPUT_PATH = " + OUTPUT_PATH);
                        break;
                    case "--genomes-file": case "-gf":
                        PATH_TO_THE_GENOMES_FILE = args[i + 1];
                        theDir = new File(PATH_TO_THE_GENOMES_FILE);
                        if (!theDir.exists()) {
                                System.out.println(PATH_TO_THE_GENOMES_FILE + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("PATH_TO_THE_GENOMES_FILE = " + PATH_TO_THE_GENOMES_FILE);
                        break;
                    case "--proteomes-file": case "-pf":
                        PATH_TO_THE_PROTEOMES_FILE = args[i + 1];
                        theDir = new File(PATH_TO_THE_PROTEOMES_FILE);
                        if (!theDir.exists()) {
                                System.out.println(PATH_TO_THE_PROTEOMES_FILE + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("PATH_TO_THE_PROTEOMES_FILE = " + PATH_TO_THE_PROTEOMES_FILE);
                        break;
                    case "--annotations-file": case "-af":
                        PATH_TO_THE_ANNOTATIONS_FILE = args[i + 1];
                        theDir = new File(PATH_TO_THE_ANNOTATIONS_FILE);
                        if (!theDir.exists()) {
                                System.out.println(PATH_TO_THE_ANNOTATIONS_FILE + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("PATH_TO_THE_ANNOTATIONS_FILE = " + PATH_TO_THE_ANNOTATIONS_FILE);
                        break;
                    case "--connect-annotations": case "-ca":
                        CONNECT_ANNOTATIONS = true;
                        --i;
                        System.out.println("CONNECT_ANNOTATIONS = true");
                        break;
                    case "--regions-file": case "-rf":
                        PATH_TO_THE_REGIONS_FILE = args[i + 1];
                        theDir = new File(PATH_TO_THE_REGIONS_FILE);
                        if (!theDir.exists()) {
                                System.out.println(PATH_TO_THE_REGIONS_FILE + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("PATH_TO_THE_REGIONS_FILE = " + PATH_TO_THE_REGIONS_FILE);
                        break;
                    case "--genome-numbers": case "-gn":
                        PATH_TO_THE_GENOME_NUMBERS_FILE = args[i + 1];
                        theDir = new File(PATH_TO_THE_GENOME_NUMBERS_FILE);
                        if (!theDir.exists()) {
                                System.out.println(PATH_TO_THE_GENOME_NUMBERS_FILE + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("PATH_TO_THE_GENOME_NUMBERS_FILE = " + PATH_TO_THE_GENOME_NUMBERS_FILE);
                        break;
                    case "--intersection-rate": case "-ir": 
                        y = Double.parseDouble(args[i + 1]);
                        if (y >= 0.001 && y <= 0.1)
                            INTERSECTION_RATE = y;
                        else {
                            System.out.println("Choose INTERSECTION_RATE in the range [0.001..0.1] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("INTERSECTION_RATE = " + INTERSECTION_RATE);
                        break;
                    case "--similarity-threshold": case "-st": 
                        x = Integer.parseInt(args[i + 1]);
                        if (x > 0 && x < 100)
                            MIN_PROTEIN_IDENTITY = x;
                        else {
                            System.out.println("Choose MIN_PROTEIN_IDENTITY in the range ]0..100[ or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("MIN_PROTEIN_IDENTITY = " + MIN_PROTEIN_IDENTITY);
                        break;
                    case "--mcl-inflation": case "-mi": 
                        y = Double.parseDouble(args[i + 1]);
                        if (y > 1 && y < 19)
                            MCL_INFLATION = y;
                        else {
                            System.out.println("Choose MCL_INFLATION in the range ]1.0..19.0[ or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("MCL_INFLATION = " + MCL_INFLATION);
                        break;
                    case "--contrast": case "-ct": 
                        y = Double.parseDouble(args[i + 1]);
                        if (y > 0 && y < 10)
                            CONTRAST = y;
                        else {
                            System.out.println("Choose CONTRAST in the range [1..9] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("CONTRAST = " + CONTRAST);
                        break;
                    case "--relaxation": case "-rn": 
                        x = Integer.parseInt(args[i + 1]);
                        if (x >= 1 && x <= 8){
                            INTERSECTION_RATE = new double[] {0, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01}[x];
                            MIN_PROTEIN_IDENTITY = new int[]   {0,95, 85, 75, 65, 55, 45, 35, 25}[x];
                            MCL_INFLATION = new double[]{0, 10.8, 9.6, 8.4, 7.2, 6.0, 4.8, 3.6, 2.4}[x];
                            CONTRAST = new double[] {0,8, 7, 6, 5, 4, 3, 2, 1 }[x];
                        }
                        else {
                            System.out.println("Choose RELAXATION in the range [1..8] or do not specify it to use the default values.");
                            System.exit(1);
                        }
                        System.out.println("INTERSECTION = " + INTERSECTION_RATE);
                        System.out.println("MIN_PROTEIN_IDENTITY = " + MIN_PROTEIN_IDENTITY);
                        System.out.println("MCL_INFLATION = " + MCL_INFLATION);
                        System.out.println("CONTRAST = " + CONTRAST);
                        break;
                    case "--threads_number": case "-tn":
                        x = Integer.parseInt(args[i + 1]);
                        if (x < cores)
                            THREADS = x;
                        else {
                            System.out.println("The maximum number of threads on this machine = " + cores + ".");
                            THREADS = cores;
                        }
                        System.out.println("THREADS = " + THREADS);
                        break;
                    case "--gap-open": case "-go":
                        x = Integer.parseInt(args[i + 1]);
                        if (x >= -50 && x <= -1)
                            GAP_OPEN = x;
                        else {
                            System.out.println("Choose GAP_OPEN in the range [-50..-1] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("GAP_OPEN = " + GAP_OPEN);
                        break;
                    case "--gap-extention": case "-ge":
                        x = Integer.parseInt(args[i + 1]);
                        if (x >= -5 && x <= -1)
                            GAP_EXT = x;
                        else {
                            System.out.println("Choose GAP_EXT in the range [-5..-1] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("GAP_EXT = " + GAP_EXT);
                        break;
                    case "--feature_type": case "-ft":
                        if (labels.containsKey(args[i + 1]))
                            FEATURE = args[i + 1];
                        else {
                            System.out.println(args[i + 1] + " is an unknown feature.");
                            System.exit(1);
                        }
                        System.out.println("FEATURE = " + FEATURE);
                        break;
                    case "--first_sra": case "-1":
                        PATH_TO_THE_FIRST_SRA = args[i + 1];
                        theDir = new File(PATH_TO_THE_FIRST_SRA);
                        if (!theDir.exists()) {
                                System.out.println(PATH_TO_THE_FIRST_SRA + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("PATH_TO_THE_FIRST_SRA = " + PATH_TO_THE_FIRST_SRA);
                        break;
                    case "--second_sra": case "-2":
                        PATH_TO_THE_SECOND_SRA = args[i + 1];
                        theDir = new File(PATH_TO_THE_SECOND_SRA);
                        if (!theDir.exists()) {
                                System.out.println(PATH_TO_THE_SECOND_SRA + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("PATH_TO_THE_SECOND_SRA = " + PATH_TO_THE_SECOND_SRA);
                        break;
                    case "--clip_strigency": case "-cs":
                        CLIPPING_STRINGENCY = Integer.parseInt(args[i + 1]);
                        switch (CLIPPING_STRINGENCY){
                            case 0:
                            System.out.println("CLIPPING_STRINGENCY = " + CLIPPING_STRINGENCY + " : no-clipping");
                                break;
                            case 1:
                            System.out.println("CLIPPING_STRINGENCY = " + CLIPPING_STRINGENCY + " : low");
                                break;
                            case 2:
                            System.out.println("CLIPPING_STRINGENCY = " + CLIPPING_STRINGENCY + " : medium");
                                break;
                            case 3:
                            System.out.println("CLIPPING_STRINGENCY = " + CLIPPING_STRINGENCY + " : high");
                                break;
                            default:
                                System.out.println("Choose CLIPPING_STRINGENCY 0, 1, 2, or 3, or do not specify it to use the default value of 0.");
                                System.exit(1);
                        }
                        break;
                    case "--interleaved": case "-il":
                        INTERLEAVED = true;
                        --i;
                        System.out.println("INTERLEAVED = true");
                        break;
                    case "--bam-format": case "-bf":
                        BAMFORMAT = true;
                        --i;
                        System.out.println("BAMFORMAT = true");
                        break;
                    case "--min_mapping-identity": case "-mmi":
                        x = Integer.parseInt(args[i + 1]);
                        if (x >= 0 && x < 100)
                           MIN_IDENTITY = x;
                        else {
                            System.out.println("Choose MIN_IDENTITY in the range [0..100[ or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("MIN_IDENTITY = " + MIN_IDENTITY);
                        break;
                    case "--num-kmer-samples": case "-nks":
                        x = Integer.parseInt(args[i + 1]);
                        if (x >= 1)
                           NUM_KMER_SAMPLES = x;
                        else {
                            System.out.println("Choose a non-zero NUM_KMER_SAMPLES or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("NUM_KMER_SAMPLES = " + NUM_KMER_SAMPLES);
                        break;
                    case "--max-alignment-length": case "-mal":
                        x = Integer.parseInt(args[i + 1]);
                        if (x >= 50 && x<=5000)
                           MAX_ALIGNMENT_LENGTH = x;
                        else {
                            System.out.println("Choose MAX_ALIGNMENT_LENGTH in the range [50..5000] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("MAX_ALIGNMENT_LENGTH = " + MAX_ALIGNMENT_LENGTH);
                        break;
                    case "--max-fragment-length": case "-mfl":
                        x = Integer.parseInt(args[i + 1]);
                        if (x >= 50 && x<=5000)
                           MAX_FRAGMENT_LENGTH = x;
                        else {
                            System.out.println("Choose MAX_FRAGMENT_LENGTH in the range [50..5000] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("MAX_FRAGMENT_LENGTH = " + MAX_FRAGMENT_LENGTH);
                        break;
                    case "--min-hit_length": case "-mhl":
                        x = Integer.parseInt(args[i + 1]);
                        if (x >= 10 && x <= 100)
                           MIN_HIT_LENGTH = x;
                        else {
                            System.out.println("Choose MIN_HIT_LENGTH in the range [10..100] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("MIN_HIT_LENGTH = " + MIN_HIT_LENGTH);
                        break;
                    case "--alignment_bound": case "-ab":
                        x = Integer.parseInt(args[i + 1]);
                        if (x >= 1 && x <= 100)
                           ALIGNMENT_BOUND = x;
                        else {
                            System.out.println("Choose ALIGNMENT_BOUND in the range [1..100] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("ALIGNMENT_BOUND = " + ALIGNMENT_BOUND);
                        break;
                    case "--max-num-locations": case "-mnl":
                        x = Integer.parseInt(args[i + 1]);
                        if (x >= 1 && x <= 100)
                           MAX_NUM_LOCATIONS = x;
                        else {
                            System.out.println("Choose MAX_NUM_LOCATIONS in the range [1..100] or do not specify it to use the default value.");
                            System.exit(1);
                        }
                        System.out.println("MAX_NUM_LOCATIONS = " + MAX_NUM_LOCATIONS);
                        break;
                    case "--alignment-mode": case "-am":
                        ALIGNMENT_MODE = Integer.parseInt(args[i + 1]);
                        switch (ALIGNMENT_MODE){
                            case 0:
                                System.out.println("ALIGNMENT_MODE = " + ALIGNMENT_MODE + " : all-hits");
                            break;
                            case -1:
                                System.out.println("ALIGNMENT_MODE = " + ALIGNMENT_MODE + " : unique pangenomic-best");
                            break;
                            case -2:
                                System.out.println("ALIGNMENT_MODE = " + ALIGNMENT_MODE + " : random pangenomic-best");
                            break;
                            case -3:
                                System.out.println("ALIGNMENT_MODE = " + ALIGNMENT_MODE + " : all pangenomic-bests");
                            break;
                            case 1:
                                System.out.println("ALIGNMENT_MODE = " + ALIGNMENT_MODE + " : unique genomic-best");
                            break;
                            case 2:
                                System.out.println("ALIGNMENT_MODE = " + ALIGNMENT_MODE + " : random genomic-best");
                            break;
                            case 3:
                                System.out.println("ALIGNMENT_MODE = " + ALIGNMENT_MODE + " : all genomic-bests");
                            break;
                            default:    
                                System.out.println("Choose ALIGNMENT_MODE in range [-3..3] or leave it to use the default value.");
                                System.exit(1);
                        }
                        break;
                    case "--raw-abundance-file": case "-raf":
                        RAW_ABUNDANCE_FILE = args[i + 1];
                        theDir = new File(RAW_ABUNDANCE_FILE);
                        if (!theDir.exists()) {
                                System.out.println(RAW_ABUNDANCE_FILE + " does not exist!");
                                System.exit(1);
                        }
                        System.out.println("RAW_ABUNDANCE_FILE = " + RAW_ABUNDANCE_FILE);
                        break;
                    case "--help": case "-h":
                        print_help_message();
                        System.exit(1);
                        break;
                }  
            }
        } catch (NumberFormatException ex){
            System.out.println("The given number is not in the correct format!");
            System.exit(1);
        }
        switch (args[0]) {
            case "bg": case "build_pangenome":
                seqLayer.initialize_pangenome();
                break;
            case "bp": case "build_panproteome":
                proLayer.initialize_panproteome();
                break;
            case "ag": case "add_genomes":
                seqLayer.add_genomes();
                break;
            case "aa": case "add_annotations":
                annLayer.add_annotaions();
                break;
            case "ra": case "remove_annotations":
                annLayer.remove_annotaions();
                break;
            case "g": case "group":
                proLayer.group();
                break;
            case "rf": case "retrieve_features":
                annLayer.retrieve_feature();
                break;
            case "rr": case "retrieve_regions":
                seqLayer.retrieve_regions();
                break;
            case "rg": case "retrieve_genomes":
                seqLayer.retrieve_genomes();
                break;
            case "rs": case "retrieve_synteny":
                seqLayer.retrieve_synteny(args[2]);
                break;
            case "m": case "map":
                seqLayer.map_reads();
                break;
            case "h": case "help":
                print_help_message();
                System.exit(1);
                break;
            case "v": case "version":
                System.out.println("PanTools version 1.1\nNeo4j community edition 3.3.1");
                System.exit(1);
            default:
                System.out.println(args[0] + " is not a valid PanTools command, type 'pantools.jar h [or help]' to see the mannual.");
                System.exit(1);
        }
        System.out.println("Total time : " + (System.currentTimeMillis() - startTime) / 1000 + "." + (System.currentTimeMillis() - startTime) % 1000 + " seconds");
        print_peak_memory();
        System.out.println("-----------------------------------------------------------------------");
    }

    /**
     * Print the manual of the software.
     */
    private static void print_help_message() {
        System.out.println("****************************************************************\n" +
"PanTools version 1.1,\n" +
"\n" +
"is a java application based on Neo4j graph database community \n" +
"edition 3.3.1 for computational pan-genomics, developed by \n" +
"Siavash Sheikhizadeh in Wageningen university, the Netherlands.\n" +
"If you use PanTools please do not forget to cite it :\n" +
"\n" +
"doi: 10.1093/bioinformatics/btw455\n" +
"\n" +
"https://github.com/Sheikhizadeh/pantools  \n" +
"\n" +
"****************************************************************\n" +
"\n" +
"Requirements\n" +
"------------\n" +
"- KMC: is a disk-based programm for counting k-mers from \n" +
"       (possibly gzipped) FASTQ/FASTA files\n" +
"       (http://sun.aei.polsl.pl/kmc).\n" +
"        You need to download it and add the path to the \n" +
"        appropriate version (linux, macos or windows) of kmc \n" +
"        and kmc_tools executables to your OS path environment \n" +
"        variable.\n" +
"\n" +
"- Java Virtual Machine version 1.8 or higher: Add the path to \n" +
"       the java executable to your OS path environment variable.\n" +
"\n" +
"- MCL: The Markov Cluster Algorithm, is a fast and scalable \n" +
"       unsupervised cluster algorithm for graphs \n" +
"       (http://micans.org/mcl ) which is needed for group \n" +
"       functionality of PanTools.\n" +
"       You need to download, unzip and compile it (see README), \n" +
"       and add the path to the mcl executable to your path\n" +
"       environment variable.\n" +
"\n" +
"Running the program \n" +
"-------------------\n" +
"java <JVM options> -jar pantools.jar <command> <arguments>\n" +
"\n" +
"pantools.jar is available in folder pantools/dist/ \n" +
"arguments is a list of key value pairs separated by whitespace.\n" +
"\n" +
"JVM options\n" +
"-----------\n" +
"[-server] \n" +
"[-XX:+UseConcMarkSweepGC]  \n" +
"[-Xmx(a number followed by g or m)]\n" +
"\n" +
"PanTools commands\n" +
"-----------------\n" +
"\n" +
"<build_pangenome or bg> \n" +
"   To build a pan-genome out of a set of genomes.\n" +
"\n" +
"   <argument keys>\n" +
"   --database_path or -dp\n" +
"      path to the pangenome database. \n" +
"   --genomes-file or -gf \n" +
"      a text file containing paths to FASTA files of genomes;\n" +
"      each in a seperated line.\n" +
"   --kmer-size or -ks\n" +
"      the size of k-mers, if not given or is out of range \n" +
"      (6 <= K_SIZE <= 255),an optimal value would be calculated automatically.    \n" +
"\n" +
"<build_panproteome or bp>\n" +
"   To build a pan-proteome out of a set of proteins.\n" +
"\n" +
"   <argument keys>\n" +
"   --database_path or -dp\n" +
"      path to the pangenome database. \n" +
"   --proteomes_file or -pf\n" +
"      a text file containing paths to FASTA files of proteomes; \n" +
"      each in a seperated line.\n" +
"             \n" +
"<add_genomes or ag>\n" +
"   To add new genomes to an available pan-genome.  \n" +
"  \n" +
"   <argument keys>\n" +
"   --database_path or -dp\n" +
"      path to the pangenome database. \n" +
"   --genomes-file or -gf\n" +
"      a text file containing paths to FASTA files of the new \n" +
"      genomes to be added to the pangeome; \n" +
"      each in a seperated line.\n" +
"\n" +
"<add_annotations or aa>\n" +
"   To add new annotations to an available pan-genome. \n" +
"\n" +
"   <argument keys>\n" +
"   --database_path or -dp \n" +
"      path to the pangenome database. \n" +
"   --annotations-file or -af\n" +
"      a text file each line of which contains genome number and \n" +
"      path to the corresponding GFF file seperated by one space.\n" +
"      Genomes are numbered in the same order they have been added\n" +
"      to the pangenome. The protein sequence of the annotated genes \n" +
"      will be also stored in the folder \"proteins\" in the same path \n" +
"      as the pangenome. \n" +
"   --connect_annotations or -ca\n" +
"      connect the annotated genomic features to the nodes of gDBG.\n" +
"\n" +
"<retrieve_features of rf>\n" +
"   To retrieve the sequence of annotated features from the pangenome. \n" +
"   For each genome the resulting FASTA file will be stored in the current \n" +
"   directory.\n" +
"\n" +
"   <argument keys>\n" +
"   --database_path or -dp\n" +
"      path to the pangenome database. \n" +
"   --genome-numbers or -gn\n" +
"      a text file containing genome_numbers for which the features will \n" +
"      be retrieved. The resulting FASTA files have two parts separated by a dot. \n" +
"      The first part determines the feature and the second determines the \n" +
"      genome number; for example, genes.1.fasta.\n" +
"   --feature-type or -ft (default = gene)\n" +
"      the feature name; for example gene, mRNA, exon, tRNA, ... \n" +
"\n" +
"<retrieve_regions or rr> \n" +
"   To retrieve the sequence of some genomic regios from the pangenome. \n" +
"   The results will be stored in the same folder as the pangenome.\n" +
"\n" +
"   <argument keys>\n" +
"   --database_path or -dp \n" +
"      path to the pangenome database. \n" +
"   --regions-file or -rf\n" +
"      a text file containing records with genome_number, \n" +
"      sequence_number, begin and end positions seperated by one \n" +
"      space for each region. The resulting FASTA file would have \n" +
"      the same name with an additional .fasta extention.\n" +
"\n" +
"<retrieve_genomes or rg>\n" +
"   To retrieve the full sequence of some genomes. The results will be \n" +
"   stored in the same folder as the pangenome itself.\n" +
"\n" +
"   <argument keys>\n" +
"   --database_path or -dp\n" +
"      path to the pangenome database. \n" +
"   --genome-numbers or -gn\n" +
"      a text file containing genome_numbers to be retrieved in each line. \n" +
"      The resulting FASTA files are named like Genome_x.fasta.\n" +
"\n" +
"<group or g>\n" +
"   To add homology nodes which point to a groups of homologous proteins.\n" +
"\n" +
"   <argument keys>\n" +
"   --database_path or -dp\n" +
"      path to the pangenome database. \n" +
"   --intersection-rate or -ir (default = 0.09)\n" +
"      the fraction of kmers needs to be shared by two \n" +
"      intersecting proteins. Should be in range [0.001, 0.1].\n" +
"   --min-protein-identity or -mpi (default = 95) \n" +
"      the minimum similarity score. Should be in range [1-99]. \n" +
"   --mcl-inflation or -mi (default = 9.6) \n" +
"      the MCL inflation. Should be in range ]1-19[.\n" +
"   --contrast or -ct (default = 8)\n" +
"      the contrast factor. Should be in range ]0-10[.\n" +
"   --relaxation or rn (default 1)\n" +
"      the relaxation in homology calls. Should be in range [1, 8], \n" +
"      from strict to relaxed.\n" +
"   --threads-number or -tn (default = 1) \n" +
"      the number of parallel working threads\n" +
"\n" +
"<map or m>\n" +
"   To map single or paired-end reads to all or a sebset of constituent genomes.\n" +
"\n" +
"   <argument keys>\n" +
"   --database_path or -dp\n" +
"      path to the pangenome database. \n" +
"   -1 \n" +
"      a text file containing path to the first short-read archive in FASTQ\n" +
"      or FASTA format. \n" +
"   -2 \n" +
"      optionally, a text file containing path to the second short-read \n" +
"      archive in FASTQ or FASTA format. \n" +
"   --genome-numbers or -gn\n" +
"      a text file containing genome_numbers to map reads against in \n" +
"      each line. \n" +
"   --output-path or -op (default: database path determined by -dp)\n" +
"      path to the output files.\n" +
"   --threads-number or -tn (default = 1) \n" +
"      the number of parallel working threads\n" +
"   --min-mapping-score or -mms (default = 20)\n" +
"      the minimum of read mapping score\n" +
"   --num-kmer-samples or -nks (default = 20)\n" +
"      the number of kmers sampled from read\n" +
"   --min-hit-length or -mhl (default = 17)\n" +
"      the minimum acceptable length of alignment after soft-clipping\n" +
"   --max-alignment-length or -mal (default = 1000)\n" +
"      the maximum acceptable length of alignment\n" +
"   --max-fragment-length or -mfl (default = 2000)\n" +
"      the maximum acceptable length of fragment\n" +
"   --max-num-locations or -mnl (default = 20)\n" +
"      the maximum number of location of candidate hits to examine\n" +
"   --alignment-bound or -ab (default = 7)\n" +
"      the length of bound of banded alignment\n" +
"   --clipping-stringency or -ci (default = 2)\n" +
"      the stringency of soft-clipping  \n" +
"      0 : no soft clipping\n" +
"      1 : low\n" +
"      2 : medium\n" +
"      3 : high\n" +
"   --bam-format or -bf (default = FALSE)\n" +
"      the alignment format (.sam or .bam)\n" +
"   --alignment-mode or -am (default = 0)\n" +
"      the alignment mode\n" +
"      0 : Competitive, only-best\n" +
"      1 : Competitive, all-bests\n" +
"      2 : Normal, only-best\n" +
"      3 : Normal, all-bests\n" +
"<version or v>\n" +
"   To show the versions of PanTools and Neo4j.\n" +
"   \n" +
"<help or h>\n" +
"   To see this document.\n" +
"\n" +
"Visualization in the Neo4j browser\n" +
"----------------------------------\n" +
"   Neo4j browser allows you to run Cypher queries and receive \n" +
"   the results in a tabular or a graph representation mode. \n" +
"   You need to download the appropriate version of Neo4j. \n" +
"   To visualize the pangenome of two HIV strains provided \n" +
"   as a sample data in pantools repositiory, take these actions \n" +
"   on a linux machine. Windows users could also download the\n" +
"   Neo4j desktop application for starting and stopping a server \n" +
"   instead of usingn commandline.\n" +
"1. Add the path to the Neo4j /bin directory to the path \n" +
"   environment variable.\n" +
"2. Hard-code the path to your pangenome in the configuration file \n" +
"   ,NEO4J-DIRECTORY/conf/neo4j.conf, by: \n" +
"   dbms.directories.data = PATH_TO_THE_PANGENOME_DATABASE\n" +
"3. Start the Neo4j database server by: \n" +
"   neo4j start\n" +
"4. open an internet browser and Open the URL http://localhost:7474\n" +
"5. To visualize the whole pangenome of two HIV strains, \n" +
"   type this simple Cypher command:\n" +
"   MATCH (n) RETURN n\n" +
"6. To stop the Neo4j server type:\n" +
"   neo4j stop\n" +
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
     * @param s The input string
     * @return The reverse complement of the input string
     */     
    public static void reverse_complement(StringBuilder s) {
        char ch;
        int i, j;
        for (i = 0, j = s.length() - 1; i < j; ++i, --j) {
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

    public static boolean executeCommand_for(String command, int seconds) {
        Process p;
        boolean success = false;
        try {
            p = Runtime.getRuntime().exec(command);
            success = p.waitFor(seconds, TimeUnit.SECONDS);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return success;
    }
    
    public static boolean is_fasta(String file_name){
        try{
            BufferedReader in = open_file(file_name);
            String line;
            while ((line = in.readLine()) != null){
                if (line.equals("")) 
                    continue;
                else {
                    in.close();
                    return line.charAt(0) == '>';
                }            
            }
        } catch (IOException ex){
            System.err.println(ex.getMessage());
        }
        return false;
    }

    public static boolean is_fastq(String file_name){
        try{
            BufferedReader in = open_file(file_name);
            String line;
            while ((line = in.readLine()) != null){
                if (line.equals("")) 
                    continue;
                else {
                    in.close();
                    return line.charAt(0) == '@';
                }
            }
        } catch (IOException ex){
            System.err.println(ex.getMessage());
        }
        return false;
    }
    
    public static BufferedReader open_file(String filename){
        try{        
            String[] fields = filename.split("\\.");
            String file_type = fields[fields.length - 1].toLowerCase();
            if (file_type.equals("gz") || file_type.equals("gzip") || file_type.equals("bz2") || file_type.equals("bzip2"))
                return getBufferedReaderForCompressedFile(filename);//BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filename)), "UTF-8"));                    
            else 
                return new BufferedReader(new BufferedReader(new FileReader(filename)));                    
        } catch (IOException ex){
            System.out.println(ex.getMessage());
            return null;
        }
    }
    
    public static int get_lines_count(String file_name, boolean skip_empty){
        int count = 0;
        try{
            BufferedReader in = open_file(file_name);
            String line;
            while ((line = in.readLine()) != null){
                if (skip_empty && line.equals("")) 
                    continue;
                ++count;
            }
            in.close();
        } catch (IOException ex){
            System.err.println(ex.getMessage());
        }
        return count;
    }
    
    public static BufferedReader getBufferedReaderForCompressedFile(String fileIn){
        try{
            FileInputStream fin = new FileInputStream(fileIn);
            BufferedInputStream bis = new BufferedInputStream(fin);
            CompressorInputStream input = new CompressorStreamFactory().createCompressorInputStream(bis);
            BufferedReader br2 = new BufferedReader(new InputStreamReader(input));
            return br2;
        } catch (Exception ex){
            System.err.println(ex.getMessage() + "\nFailed to open the compresse file!");
            return null;
        }
    }    


}
