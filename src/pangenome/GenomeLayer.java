/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import genome.SequenceDatabase;
import index.IndexPointer;
import index.IndexDatabase;
import index.kmer;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.graphdb.Direction;
import static org.neo4j.graphdb.factory.GraphDatabaseSettings.keep_logical_logs;
import org.neo4j.io.fs.FileUtils;
import static pantools.Pantools.GENOME_DATABASE_PATH;
import static pantools.Pantools.GRAPH_DATABASE_PATH;
import static pantools.Pantools.INDEX_DATABASE_PATH;
import pantools.Pantools.RelTypes;
import static pantools.Pantools.degenerate_label;
import static pantools.Pantools.gene_label;
import static pantools.Pantools.genomeDb;
import static pantools.Pantools.genome_label;
import static pantools.Pantools.graphDb;
import static pantools.Pantools.indexDb;
import static pantools.Pantools.nucleotide_label;
import static pantools.Pantools.num_bases;
import static pantools.Pantools.num_degenerates;
import static pantools.Pantools.num_edges;
import static pantools.Pantools.num_nodes;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.phaseTime;
import static pantools.Pantools.sequence_label;
import static pantools.Pantools.startTime;
import static pantools.Pantools.MAX_TRANSACTION_SIZE;
import static pantools.Pantools.complement;
import static pantools.Pantools.reverse_complement;
import static pantools.Pantools.write_fasta;
import static index.IndexDatabase.get_K;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.LinkedList;
import static pantools.Pantools.DEBUG;
import static pantools.Pantools.K;
import static pantools.Pantools.SHOW_KMERS;
import static pantools.Pantools.low_complexity_label;
import static pantools.Pantools.variation_label;

/**
 * Implements all the functionalities related to the sequence layer of the pangenome
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class GenomeLayer {

    private kmer fwd_kmer, rev_kmer, k_mer;
    private long seq_len;
    private static int genome, sequence, position;
    private long curr_index;
    private byte curr_side;
    private Node curr_node;
    private Node db_node;
    private Node new_node;
    private Node degenerate_node;
    private IndexPointer pointer;
    private int fwd_code, rev_code;
    private boolean finish;
    private static int ANCHORS = 10000; // The number of anchor nodes

    /**
     * The constructor of the class.
     */    
    public GenomeLayer() {
        pointer = new IndexPointer();
        finish = false;
    }

    /**
     * Constructs a pangenome database from given genomes.
     * 
     * @param genome_paths_file Path to the FASTA genome files. 
     * @param pangenome_path Path to the database folder
     */  
    public void initialize_pangenome(String genome_paths_file, String pangenome_path, int k) {
    // If a database folder is already exist in the specified path, removes all the content of it.    
        File theDir = new File(pangenome_path);
        if (theDir.exists()) {
            try {
                FileUtils.deleteRecursively(new File(pangenome_path));
            } catch (IOException ioe) {
                System.out.println("Failed to delete the database " + pangenome_path);
                System.exit(1);  
            }
        } else {
            try {
                theDir.mkdir();
            } catch (SecurityException se) {
                System.out.println("Failed to create " + pangenome_path);
                System.exit(1);
            }
        }
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(pangenome_path + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();  
        registerShutdownHook(graphDb);
        startTime = System.currentTimeMillis();
        num_nodes = 0;
        num_edges = 0;
        num_bases = 0;
        num_degenerates = 0;
        genomeDb = new SequenceDatabase(pangenome_path + GENOME_DATABASE_PATH, genome_paths_file);
        indexDb = new IndexDatabase(pangenome_path + INDEX_DATABASE_PATH, genome_paths_file, genomeDb, k);
        K = get_K();
        try (Transaction tx = graphDb.beginTx()) {
            db_node = graphDb.createNode(pangenome_label);
            db_node.setProperty("k_mer_size", K);
            db_node.setProperty("date", new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(new Date()));
            tx.success();
        }
        construct_pangenome(0);
        System.out.println("Number of kmers:   " + indexDb.length());
        System.out.println("Number of nodes:   " + num_nodes);
        System.out.println("Number of edges:   " + num_edges);
        System.out.println("Number of bases:   " + num_bases);
        System.out.println("Number of degenerate nodes:   " + num_degenerates);
        try (Transaction tx = graphDb.beginTx()) {
            db_node.setProperty("k_mer_size", K);
            db_node.setProperty("num_k_mers", indexDb.length());
            db_node.setProperty("num_nodes", num_nodes);
            db_node.setProperty("num_degenerate_nodes", num_degenerates);
            db_node.setProperty("num_edges", num_edges);
            db_node.setProperty("num_genomes", genomeDb.num_genomes);
            db_node.setProperty("num_bases", num_bases);
            tx.success();
        }
        graphDb.shutdown();
        genomeDb.close();
        indexDb.close();
        File directory = new File(pangenome_path + GRAPH_DATABASE_PATH);
        for (File f : directory.listFiles()) {
            if (f.getName().startsWith("neostore.transaction.db.")) {
                f.delete();
            }
        }
        System.out.println("graph.db size: " + getFolderSize(new File(pangenome_path + GRAPH_DATABASE_PATH)) + " MB");
        System.out.println("index.db size: " + getFolderSize(new File(pangenome_path + INDEX_DATABASE_PATH)) + " MB");
        System.out.println("genome.db size: " + getFolderSize(new File(pangenome_path + GENOME_DATABASE_PATH)) + " MB");
    }
    
    /**
     * Adds new genomes to an available pangenome.
     * 
     * @param genome_paths_file Path to the FASTA genome files. 
     * @param pangenome_path Path to the database folder
     */
    public void add_genomes(String genome_paths_file, String pangenome_path) {
        int j, len, previous_num_genomes;
        long byte_number = 0;
        int[] address = new int[4];
        Node start, seq_node;
        if (new File(pangenome_path + GRAPH_DATABASE_PATH).exists()) {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(pangenome_path + GRAPH_DATABASE_PATH))
                    .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
            registerShutdownHook(graphDb);
            startTime = System.currentTimeMillis();
            try (Transaction tx = graphDb.beginTx()) {
                db_node = graphDb.findNodes(pangenome_label).next();
                if (db_node == null) {
                    System.out.println("Can not locate database node!");
                    System.exit(1);
                }
            // Reads the properties of the pangenome    
                K = (int) db_node.getProperty("k_mer_size");
                num_nodes = (long) db_node.getProperty("num_nodes");
                num_edges = (long) db_node.getProperty("num_edges");
                num_degenerates = (int) db_node.getProperty("num_degenerate_nodes");
                num_bases = 0;
                previous_num_genomes = (int) db_node.getProperty("num_genomes");
            // if the genome database is not available, reconstruct it.    
                if (!Files.exists(Paths.get(pangenome_path + GENOME_DATABASE_PATH))) 
                {
                // read genomes information from the graph and rebuild the genomes database
                    genomeDb = new SequenceDatabase(pangenome_path + GENOME_DATABASE_PATH, graphDb);
                    StringBuilder seq = new StringBuilder();
                    for (address[0] = 1; address[0] <= genomeDb.num_genomes; ++address[0]) {
                        for (address[1] = 1; address[1] <= genomeDb.num_sequences[address[0]]; ++address[1]) {
                            seq_node = graphDb.findNode(sequence_label, "identifier", address[0] + "_" + address[1]);
                            start = seq_node.getRelationships(Direction.OUTGOING).iterator().next().getEndNode();
                            address[2] = 1;
                            address[3] = (int) genomeDb.sequence_length[address[0]][address[1]];
                            extract_sequence(seq, new IndexPointer(start.getId(), true, 0, -1l), address, K);
                            len = seq.length();
                            if (len % 2 == 1) {
                                --len;
                            }
                            for (j = 0; j < len; j += 2, ++byte_number) {
                                genomeDb.genomes_buff[(int) (byte_number / genomeDb.parts_size[0])].put((byte) ((genomeDb.binary[seq.charAt(j)] << 4) | genomeDb.binary[seq.charAt(j + 1)]));
                            }
                            if (len == seq.length() - 1) {
                                genomeDb.genomes_buff[(int) (byte_number / genomeDb.parts_size[0])].put((byte) (genomeDb.binary[seq.charAt(len)] << 4));
                                ++byte_number;
                            }
                        }
                    }
                    genomeDb = new SequenceDatabase(pangenome_path + GENOME_DATABASE_PATH); //Readable only
                } else {
                    genomeDb = new SequenceDatabase(pangenome_path + GENOME_DATABASE_PATH);
                }
                genomeDb.add_genomes(pangenome_path + GENOME_DATABASE_PATH, genome_paths_file);
                indexDb = new IndexDatabase(pangenome_path + INDEX_DATABASE_PATH, genome_paths_file, genomeDb, graphDb, previous_num_genomes);
                tx.success();
            }
        // the sequences should be dropped out as they will change and add_sequence_properties() function will rebuild them.    
            drop_nodes_property("sequence");
        // the edge colors should be dropped out as they will change and localize_nodes() function will rebuild them again.    
            drop_edges_colors();
            construct_pangenome(previous_num_genomes);
            System.out.println("Number of kmers:   " + indexDb.length());
            System.out.println("Number of nodes:   " + num_nodes);
            System.out.println("Number of edges:   " + num_edges);
            System.out.println("Number of bases:   " + num_bases);
            System.out.println("Number of degenerate nodes:   " + num_degenerates);
            try (Transaction tx = graphDb.beginTx()) {
                db_node.setProperty("k_mer_size", K);
                db_node.setProperty("num_k_mers", indexDb.length());
                db_node.setProperty("num_nodes", num_nodes);
                db_node.setProperty("num_degenerate_nodes", num_degenerates);
                db_node.setProperty("num_edges", num_edges);
                db_node.setProperty("num_genomes", genomeDb.num_genomes);
                db_node.setProperty("num_bases", num_bases);
                tx.success();
            }
            graphDb.shutdown();
            genomeDb.close();
            indexDb.close();
            File directory = new File(pangenome_path + GRAPH_DATABASE_PATH);
            for (File f : directory.listFiles()) {
                if (f.getName().startsWith("neostore.transaction.db.")) {
                    f.delete();
                }
            }
            System.out.println("graph.db size: " + getFolderSize(new File(pangenome_path + GRAPH_DATABASE_PATH)) + " MB");
            System.out.println("index.db size: " + getFolderSize(new File(pangenome_path + INDEX_DATABASE_PATH)) + " MB");
            System.out.println("genome.db size: " + getFolderSize(new File(pangenome_path + GENOME_DATABASE_PATH)) + " MB");
        } else {
            System.out.println("No database found in " + pangenome_path);
            System.exit(1);
        }
    }

    public void remove_genomes(String genome_paths_file, String pangenome_path) {
        
    }
    
    private void parse_vcf(int genome, BufferedWriter log_file, String vcf_file){
        int i, trsc, num_genes, num_mRNAs, num_tRNAs, num_rRNAs, feature_len, offset, node_len;
        String sequence_id, current_sequence_id=null, attribute;
        Node node, gene_node, rna_node, variation_node, feature_node, parent_node = null;
        Relationship rel;
        String[] fields,parent_ids;
        String strand, line;
        long variation_node_id;
        //LinkedList<Node> rna_nodes = new LinkedList();
        IndexPointer start_ptr, stop_ptr;
        variation_node = graphDb.createNode(variation_label);
        variation_node.setProperty("genome", genome);
        variation_node.setProperty("date", new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(new Date()));
        variation_node.createRelationshipTo(graphDb.findNodes(genome_label,"number",genome).next(), RelTypes.varies);
        variation_node_id = variation_node.getId();
        num_genes = num_mRNAs = num_tRNAs = num_rRNAs = 0;
        try (BufferedReader in = new BufferedReader(new FileReader(vcf_file))) {
            // for each record of gff file
            while (in.ready())
            {
                try (Transaction tx2 = graphDb.beginTx()) {
                    for (trsc = 0; trsc < MAX_TRANSACTION_SIZE && in.ready(); ++trsc) {
                        line = in.readLine().trim();
                        if (line.equals("") || line.charAt(0) == '#') // if line is empty or a comment skip it
                            continue;
                        fields = line.split("\\t");
                        if (fields.length < 8)
                            continue;
                        if (trsc % 500 == 1)
                            System.out.print("\r" + genome + "\t" + num_genes + "\t" + num_mRNAs + "\t" + num_tRNAs + "\t" + num_rRNAs + "\t");
                    }// for trsc
                    tx2.success();
                } // tx2
            } // while lines
            in.close();
            System.out.println("\r" + genome + "\t" + num_genes + "\t" + num_mRNAs + "\t" + num_tRNAs + "\t" + num_rRNAs + "\t");
            log_file.write("Genome "+genome + " : " + num_genes + " genes\t" + num_mRNAs + " mRNAs\t" + num_tRNAs + " tRNAs\t" + num_rRNAs + " rRNAs\n");
            log_file.write("----------------------------------------------------\n");
        }catch (IOException ioe) {
            System.out.println(ioe.getMessage());
            System.out.println("Could not open " + vcf_file + "!");
        }
    }
    
    /**
     * Retrieves the sequence of a number of genomic regions from the pangenome and stores them in a FASTA file.
     * 
     * @param region_records_file a text file with lines containing genome number, sequence number, start and stop positions
     *        of the genomic regions seperated by one space.
     * @param pangenome_path Path to the database folder
     */
    public void retrieve_regions(String region_records_file, String pangenome_path) {
        if (new File(pangenome_path + GRAPH_DATABASE_PATH).exists()) {
            String[] fields;
            String line, out_file_name;
            IndexPointer start_ptr;
            StringBuilder seq;
            int c, num_regions = 0, proper_regions = 0;
            int[] address = new int[4];
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(pangenome_path + GRAPH_DATABASE_PATH))
                    .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
            registerShutdownHook(graphDb);
            try (Transaction tx = graphDb.beginTx()) {
                K = (int) graphDb.findNodes(pangenome_label).next().getProperty("k_mer_size");
                tx.success();
            }
            seq = new StringBuilder();
            try {
                BufferedReader in = new BufferedReader(new FileReader(region_records_file));
                while (in.ready()) {
                    line = in.readLine().trim();
                    if (line.equals("")) {
                        continue;
                    }
                    ++num_regions;
                }
                in.close();
            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            }
            startTime = System.currentTimeMillis();
            genomeDb = new SequenceDatabase(pangenome_path + GENOME_DATABASE_PATH);
            try (Transaction tx = graphDb.beginTx()) {
                try (BufferedReader in = new BufferedReader(new FileReader(region_records_file))) {
                    fields = region_records_file.split("\\/");
                    out_file_name = pangenome_path + fields[fields.length - 1].split("\\.")[0] + ".fasta";
                    BufferedWriter out = new BufferedWriter(new FileWriter(out_file_name));
                    for (c = 0; in.ready();) {
                        line = in.readLine().trim();
                        if (line.equals("")) {
                            continue;
                        }
                        fields = line.trim().split("\\s");
                        address[0] = Integer.parseInt(fields[0]);
                        address[1] = Integer.parseInt(fields[1]);
                        address[2] = Integer.parseInt(fields[2]);
                        address[3] = Integer.parseInt(fields[3]);
                        if (address[0] <= genomeDb.num_genomes && address[1] <= genomeDb.num_sequences[address[0]] && address[2] >= 1 && address[3] <= genomeDb.sequence_length[address[0]][address[1]]){
                            start_ptr = locate(address, K);
                            proper_regions++;
                            extract_sequence(seq, start_ptr, address, K);
                            out.write(">genome:" + address[0] + " sequence:" + address[1] + " from:" + address[2] + " to:" + address[3] + " length:" + seq.length() + "\n");
                            write_fasta(out, seq.toString(), 70);
                            seq.setLength(0);
                            ++c;
                            if (c % (num_regions / 100 + 1) == 0) {
                                System.out.print((long) c * 100 / num_regions + 1 + "%\r");
                            }
                        } else
                            System.out.println(line + "is not a proper coordinate!");
                    }
                    in.close();
                    out.close();
                    System.out.println(proper_regions + " out of " + num_regions + " genomic regions found and retrieved successfully (See " + out_file_name + ")");
                } catch (IOException ioe) {
                    System.out.println("Failed to read file names!");
                    System.exit(1);
                }
                tx.success();
            }
            graphDb.shutdown();
        }
    }
    
    /**
     * Reconstructs all or some of the genomes in separated FASTA files.
     * 
     * @param genome_numbers_file A text file containing the genome numbers to be . 
     * @param pangenome_path Path to the database folder
     */
    public void retrieve_genomes(String genome_numbers_file, String pangenome_path) {
        if (new File(pangenome_path + GENOME_DATABASE_PATH).exists()) {
            BufferedReader in;
            BufferedWriter out;
            IndexPointer start;
            String genome_number;
            int[] address;
            StringBuilder seq;
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(pangenome_path + GRAPH_DATABASE_PATH))
                    .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
            registerShutdownHook(graphDb);
            startTime = System.currentTimeMillis();
            genomeDb = new SequenceDatabase(pangenome_path + GENOME_DATABASE_PATH);
            address = new int[4];
            seq = new StringBuilder();
            try (Transaction tx = graphDb.beginTx()) {
                db_node = graphDb.findNodes(pangenome_label).next();
                if (db_node == null) {
                    System.out.println("Can not locate database node!");
                    System.exit(1);
                }
                K = (int) db_node.getProperty("k_mer_size");
                try {
                    in = new BufferedReader(new FileReader(genome_numbers_file));
                    while (in.ready()) {
                        genome_number = in.readLine().trim();
                        if (genome_number.equals(""))
                            continue;
                        try{
                            address[0] = Integer.parseInt(genome_number);
                        }catch(NumberFormatException e){
                            System.out.println(genome_number + "is not a valid genome number.");
                            continue;
                        }
                        if (address[0] < 1 || address[0] > genomeDb.num_genomes){
                            System.out.println(genome_number + "is not a valid genome number.");
                            continue;
                        }
                        System.out.println("Reconstructing genome " + genome_number + "...");
                        try {
                            out = new BufferedWriter(new FileWriter(pangenome_path + "/Genome_" + genome_number + ".fasta"));
                            for (address[1] = 1; address[1] <= genomeDb.num_sequences[address[0]]; ++address[1]) {
                                System.out.println("Sequence " + address[1] + " length = " + genomeDb.sequence_length[address[0]][address[1]]);
                                address[2] = 1;
                                address[3] = (int)genomeDb.sequence_length[address[0]][address[1]];
                                start = locate(address, K);
                                out.write(">" + genomeDb.sequence_titles[address[0]][address[1]] + "\n");
                                extract_sequence(seq, start, address, K);
                                write_fasta(out, seq.toString(), 80);
                                seq.setLength(0);
                            }
                            out.close();
                        } catch (IOException e) {
                            System.out.println(e.getMessage());
                            System.exit(1);
                        }
                    }
                    in.close();
                } catch (IOException ioe) {
                    System.out.println("Failed to read file names!");
                    System.exit(1);
                }
                tx.success();
            }
            System.out.println("Genomes were stored in the database directory.");
            graphDb.shutdown();
            genomeDb.close();
        } else {
            System.out.println("No database found in " + pangenome_path);
            System.exit(1);
        }
    }

    public void retrieve_synteny(String genome, String pangenome_path) {
        if (new File(pangenome_path + GENOME_DATABASE_PATH).exists()) {
            int g;
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(pangenome_path + GRAPH_DATABASE_PATH))
                    .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
            registerShutdownHook(graphDb);
            startTime = System.currentTimeMillis();
            genomeDb = new SequenceDatabase(pangenome_path + GENOME_DATABASE_PATH);
            System.out.println( genomeDb.genome_names[2]);
            try (Transaction tx = graphDb.beginTx()) {
                db_node = graphDb.findNodes(pangenome_label).next();
                if (db_node == null) {
                    System.out.println("Can not locate database node!");
                    System.exit(1);
                }
                K = (int) db_node.getProperty("k_mer_size");
                tx.success();
            }
            try {
                try {
                    g = Integer.parseInt(genome.trim());
                    System.out.println("Reconstructing synteny map between genome " + g + " and the rest :");
                    extract_synteny(pangenome_path, g);
                }catch(NumberFormatException e){
                    System.out.println("Invalid genome number!");
                    System.exit(1);
                }
            } catch (IOException ioe) {
                System.out.println("Failed to read file names!");
                System.exit(1);
            }
            System.out.println("Synteny files are ready in " + pangenome_path);
            graphDb.shutdown();
            genomeDb.close();
        } else {
            System.out.println("No database found in " + pangenome_path);
            System.exit(1);
        }
    }
    
    public static void extract_synteny(String pangenome_path, int genome1) throws IOException {
        Relationship rel;
        Node node, neighbor;
        IndexPointer start;
        BufferedWriter[] out_fwd = new BufferedWriter[genomeDb.num_genomes + 1];
        BufferedWriter[] out_rev = new BufferedWriter[genomeDb.num_genomes + 1];
        String formatStr = "%8s%10s%10s\n";
        int[] loc2;
        int genome2, loc1, seq1, seq2, i, node_len;
        long offset1, offset2;
        char side1, side2;
        String origin1;
        for (seq1 = 1; seq1 <= genomeDb.num_sequences[genome1]; ++seq1) {
            try (Transaction tx = graphDb.beginTx()) {
                loc1 = 0;
                start = locate(new int[]{genome1, seq1, loc1 + 1}, K);
                origin1 = "a" + genome1 + "_" + seq1;
                node = graphDb.getNodeById(start.node_id);
                side1 = start.canonical ? 'F' : 'R';
                while (true) {
                    //System.out.println(loc1);
                    node_len = (int) node.getProperty("length");
                    for (Relationship r: node.getRelationships(Direction.INCOMING)){
                        side2 = r.getType().name().charAt(1);
                        for (String origin2: r.getPropertyKeys()){
                            loc2 = (int[])r.getProperty(origin2);
                            genome2 = Integer.parseInt(origin2.substring(1).split("_")[0]);
                            seq2 = Integer.parseInt(origin2.split("_")[1]);
                            offset2 = genomeDb.sequence_offset[genome2][seq2];
                            offset1 = genomeDb.sequence_offset[genome1][seq1];
                            if (side1 == side2){
                                if (out_fwd[genome2] == null){
                                    out_fwd[genome2] = new BufferedWriter(new FileWriter(pangenome_path + "/F_"+ genome1 + "_" + genome2 + ".smf"));
                                    out_fwd[genome2].write("> " + genomeDb.genome_names[genome2].split("\\s")[0] + "\n");
                                }
                                for (i = 0; i < loc2.length; ++i){
                                    out_fwd[genome2].write(String.format(formatStr, offset1 + loc1 + 1, offset2 + loc2[i] + 1, node_len));
                                    }
                            } else {
                                if (out_rev[genome2] == null){
                                    out_rev[genome2] = new BufferedWriter(new FileWriter(pangenome_path + "/R_"+ genome1 + "_" + genome2 + ".smf"));
                                    out_rev[genome2].write("> " + genomeDb.genome_names[genome2].split("\\s")[0] + " Reverse\n");
                                }
                                for (i = 0; i < loc2.length; ++i)
                                    out_rev[genome2].write(String.format(formatStr, offset1 + loc1 + (side1 == 'F' ? 1 : node_len), offset2 + loc2[i] + (side2 == 'F' ? 1 : node_len), node_len));
                            }
                        }
                    }
                    loc1 += node_len - K + 1;
                    rel = get_outgoing_edge(node, origin1, loc1);
                    if (rel == null)
                        break;
                    else
                        neighbor = rel.getEndNode();
                    node = neighbor;
                    side1 = rel.getType().name().charAt(1);
                } // while
                System.out.print("\rSequence " + seq1 + " / " + genomeDb.num_sequences[genome1] + " finished.");
                tx.success();
            }
        }
        System.out.println();
        for (genome2 = 1; genome2 <= genomeDb.num_genomes; ++genome2){
            if (out_fwd[genome2] != null)
                out_fwd[genome2].close();
            if (out_rev[genome2] != null)
                out_rev[genome2].close();  
        }
    }
    
    /**
     * Extracts the genomic region belonging to the specified sequence starting at th specified node.
     * 
     * @param seq Will contains the sequence after function ends.
     * @param start_ptr A pangenome pointer which points to the node where the sequence starts.
     * @param address An array determining {genome, sequence, begin, end}properties of the sequence.
     */
    public static void extract_sequence(StringBuilder seq, IndexPointer start_ptr, int[] address, int K) {
        Relationship rel;
        Node neighbor, node;
        int[] addr = new int[]{address[0],address[1],address[2],address[3]};
        int begin = addr[2] - 1, end = addr[3] - 1;
        int len = 0, node_len, neighbor_len, seq_len, position;
        String rel_name, origin = "G" + address[0] + "S" + address[1];
        seq_len = end - begin + 1;
        seq.setLength(0);
        position = start_ptr.offset;
        node = graphDb.getNodeById(start_ptr.node_id);
        node_len = (int) node.getProperty("length");
    // Takes the part of the region lies in the first node of the path that region takes in the graph    
        if (start_ptr.canonical) {
            if (position + seq_len - 1 <= node_len - 1) { // The whole sequence lies in this node
                len += append_fwd(seq, (String) node.getProperty("sequence"), position, position + seq_len - 1);
            } else {
                len += append_fwd(seq, (String) node.getProperty("sequence"), position, node_len - 1);
            }
        } else {
            if (position - (seq_len - 1) >= 0) { // The whole sequence lies in this node
                len += append_rev(seq, (String) node.getProperty("sequence"), position - (seq_len - 1), position);
            } else {
                len += append_rev(seq, (String) node.getProperty("sequence"), 0, position);
            }
        }
    //  traverse the path of the region   
        while (len < seq_len) {
            //System.out.println(node.getId()+" "+len + " " + seq_len);
            addr[2] = (begin + len) - K + 1;
            rel = get_outgoing_edge(node, origin, addr[2]);
            neighbor = rel.getEndNode();
            rel_name = rel.getType().name();
            neighbor_len = (int) neighbor.getProperty("length");
            if (rel_name.charAt(1) == 'F') {// Enterring forward side
                if (len + neighbor_len - K + 1 > seq_len) // neighbor is the last node of the path
                    len += append_fwd(seq, (String) neighbor.getProperty("sequence"), K - 1, seq_len - len + K - 2);
                else 
                    len += append_fwd(seq, (String) neighbor.getProperty("sequence"), K - 1, neighbor_len - 1);
            }else{ // Enterring reverse side
                if (len + neighbor_len - K + 1 > seq_len) // neighbor is the last node of the pat
                    len += append_rev(seq, (String) neighbor.getProperty("sequence"), neighbor_len - K - (seq_len - len) + 1, neighbor_len - K);
                else 
                    len += append_rev(seq, (String) neighbor.getProperty("sequence"), 0, neighbor_len - K);
            }
            node = neighbor;
        } // while
    }
  
    /**
     * Give the next node of the path to be traversed through.
     * 
     * @param current_node The current node of the path we are located at. 
     * @param origin 
     * @param address An array which determine the genome, sequence and position of the desirable outgoing edge.
     * @return The outgoing edge.
     */
    public static Relationship get_outgoing_edge(Node current_node, String origin, int pos) {
        int[] occurrence;
        for (Relationship r_out : current_node.getRelationships(Direction.OUTGOING)) {
            occurrence = (int[])r_out.getProperty(origin, null);
            if (occurrence != null) {
                if (Arrays.binarySearch(occurrence, pos) >= 0)
                    return r_out;
            }
        }
        return null;
    }
    
    /**
     * Appends substring s[from..to] to the string builder.
     * @param seq String builder.
     * @param s The string.
     * @param from Start position of the substring.
     * @param to Stop position of the substring.
     * @return The length of substring appended to the string builder.
     */
    public static int append_fwd(StringBuilder seq, String s, int from, int to) {
        for (int i = from; i <= to; ++i)
            seq.append(s.charAt(i));
        return to - from + 1;
    }
    
    /**
     * Appends the reverse complement of substring s[from..to] to the string builder.
     * @param seq String builder.
     * @param s The string.
     * @param from Start position of the substring.
     * @param to Stop position of the substring.
     * @return The length of substring appended to the string builder.
     */
    public static int append_rev(StringBuilder seq, String s, int from, int to) {
        for (int i = to; i >= from; --i)
            seq.append(complement(s.charAt(i)));
        return to - from + 1;
    }
    
    /**
     * Returns a pangenome pointer pointing to the specified genomic position.
     * 
     * @param address An integer array lile {genome_number, sequence_number, begin_position, end_position}
     * @return A pointer to the genomic position in the pangenome
     */
    public static IndexPointer locate(int[] addr, int K) {
        int node_start_pos, low, high, mid , node_len, genomic_pos;
        boolean forward;
        Node node, neighbor, seq_node;
        Relationship rel;
        String anchor_sides, origin = "G" + addr[0] + "S" + addr[1];
        long[] anchor_nodes;
        int[] anchor_positions;
        int[] address = Arrays.copyOf(addr,addr.length);
        genomic_pos = address[2] - 1;
        seq_node = graphDb.findNode(sequence_label, "identifier", address[0]+"_"+address[1]);
        anchor_nodes = (long[]) seq_node.getProperty("anchor_nodes");
        anchor_positions = (int[]) seq_node.getProperty("anchor_positions");
        anchor_sides = (String) seq_node.getProperty("anchor_sides");
    // Find the immediate preceding anchor_node, searching in the sorted array of anchor positions.      
        for (low = 0, high = anchor_sides.length() - 1, mid = (low + high) / 2; low <= high; mid = (low + high) / 2) {
            if (genomic_pos < anchor_positions[mid]) {
                high = mid - 1;
            } else if (genomic_pos > anchor_positions[mid]) {
                low = mid + 1;
            } else {
                break;
            }
        }
        if (genomic_pos < anchor_positions[mid]) {
            --mid;
        }
        forward = anchor_sides.charAt(mid) == 'F';
        try (Transaction tx = graphDb.beginTx()) {
            node = graphDb.getNodeById(anchor_nodes[mid]);
            node_start_pos = anchor_positions[mid];
            node_len = (int) node.getProperty("length");
        // Traverse the pangenome from the anchor node until reach to the target
            while (node_start_pos + node_len <= genomic_pos) 
            {
                address[2] = node_start_pos + node_len - K + 1;
                rel = get_outgoing_edge(node, origin, address[2]);
                if (rel == null){
                    System.out.println("Failed to locate address : " + address[0] + " " + address[1] + " "+ address[2]);
                    break;
                }
                neighbor = rel.getEndNode();
                forward = rel.getType().name().charAt(1) == 'F';
                node_start_pos += node_len - K + 1;
                node = neighbor;
                node_len = (int) node.getProperty("length");
            }
            tx.success();
        }
        return new IndexPointer(node.getId(), forward, forward ? genomic_pos - node_start_pos : node_len - 1 - (genomic_pos - node_start_pos), -1l);
    }
    
    /***
     * Creates an edge between source and destination nodes.
     * 
     * @param src Source node
     * @param des Destination node
     * @param edge_type One of the four possible edge types: FF, FR, RF, RR
     * @param address Specifies which genomic address the edge points to. 
     * @return The newly created edge
     */
    private void connect(Node src, Node des, RelationshipType edge_type) {
        if (DEBUG) System.out.println("connect "+src.getId()+" "+edge_type.name()+" "+des.getId());
        src.createRelationshipTo(des, edge_type);
    }
    
    /**
     * Splits a node at a specified position by creating a new node called split_node as a part separated from the node.
     * @param node The node which should be split.
     * @param pos The position of the split with respect to the start on the node.
     * @return The newly created split node.
     */
    private Node split(Node node, int pos) {
        int split_len, node_len;
        int i, s_id, gen, seq, loc,starts_at;
        long inx,split_first_kmer,node_last_kmer=0;
        int[] address;
        Node neighbor, split_node;
        Relationship rel;
        kmer fwd_k_mer, rev_k_mer;
        address = (int[]) node.getProperty("address");
        gen = address[0];
        seq = address[1];
        loc = address[2];
        node_len = (int) node.getProperty("length");
        ++num_nodes;
        split_node = graphDb.createNode(nucleotide_label);
        address[0] = gen;
        address[1] = seq;
        address[2] = loc + pos;
        split_node.setProperty("address", address);
        split_len = node_len - pos;
        split_node.setProperty("length", split_len);
    // Updates the edges comming from gene level to the node.    
        for (Relationship r : node.getRelationships(RelTypes.starts, Direction.INCOMING)) {
            starts_at = (int)r.getProperty("offset");
            if (starts_at >= pos) {
                rel = r.getStartNode().createRelationshipTo(split_node, RelTypes.starts);
                rel.setProperty("offset", starts_at - pos);
                rel.setProperty("forward", r.getProperty("forward"));
                rel.setProperty("genomic_position", r.getProperty("genomic_position"));
            r.delete();
            } 
        }        
        // Updating the Kmers chain in the index  
        node_last_kmer=indexDb.find(make_kmer(gen,seq,loc+pos-1));
        split_first_kmer=indexDb.get_next_index(node_last_kmer);
        indexDb.put_next_index(-1L, node_last_kmer); 
        split_node.setProperty("first_kmer",split_first_kmer);
        split_node.setProperty("last_kmer",node.getProperty("last_kmer"));
        s_id=(int)split_node.getId();
        for(i=0,inx=split_first_kmer;inx!=-1L;inx=indexDb.get_next_index(inx),++i) // update kmer coordinates
        {
            indexDb.put_node_id(s_id, inx);
            indexDb.put_position(i, inx);
        }  
        // Moving forward-outgoing and reverse-incoming edges from node to split node.    
        for (Relationship r : node.getRelationships(Direction.OUTGOING,RelTypes.FR,RelTypes.FF)) {
            neighbor = r.getEndNode();
            if (neighbor.equals(node)) 
                neighbor = r.isType(RelTypes.FF) ? node : split_node;
            connect(split_node, neighbor, r.getType());
            r.delete();
        }
        for (Relationship r : node.getRelationships(Direction.INCOMING,RelTypes.RR,RelTypes.FR)) {
            neighbor = r.getStartNode();
            if (neighbor.equals(node)) 
                neighbor = r.isType(RelTypes.RR) ? node : split_node;
            connect(neighbor, split_node, r.getType());
            r.delete();
        }
    //  Connecting node to split node
        if (node.hasRelationship(Direction.INCOMING, RelTypes.FF, RelTypes.RF)){
            connect(node, split_node, RelTypes.FF);
            ++num_edges;
        }
        if (split_node.hasRelationship(Direction.INCOMING, RelTypes.FR, RelTypes.RR)){
            connect(split_node ,node, RelTypes.RR);
            ++num_edges;
        }
        node.setProperty("last_kmer",node_last_kmer);
        node.setProperty("length", pos + K - 1);
        return split_node;
    }
    
    /**
     * Makes a kmer located at a specific genomic position.
     * 
     * @param genome The genome number
     * @param sequence The sequence number
     * @param position The position of the kmer in the sequence
     * @return kmer The canonical form of the kmer 
     */
    private kmer make_kmer(int genome, int sequence, int position) {
        int j,fwd_code,rev_code;
        kmer fwd_kmer=new kmer(K,indexDb.get_pre_len(),indexDb.get_suf_len());
        kmer rev_kmer=new kmer(K,indexDb.get_pre_len(),indexDb.get_suf_len());
        for(j=0;j<K;++j)
        {
            fwd_code=genomeDb.get_code(genome,sequence,position+j);
            rev_code=3-fwd_code;
            fwd_kmer.next_fwd_kmer(fwd_code);
            rev_kmer.next_rev_kmer(rev_code);
        }   
        fwd_kmer.set_canonical(fwd_kmer.compare(rev_kmer) == -1);
        return fwd_kmer.get_canonical() ? fwd_kmer : rev_kmer;             
    }
    
    /**
     * Extends a new node till reach to a previously visited K-mer or a degenerate region.
     * 
     * @param node The node to be extended.
     */
    private void extend(Node node) {
        if (DEBUG) System.out.println("extending node "+node.getId());
        int begin, len;
        long id,last_kmer=(long)node.getProperty("last_kmer");
        boolean broke, degenerate;
        
        len = (int) node.getProperty("length");
        id = node.getId();
        broke = false;
        degenerate = false;        
        while (position < seq_len - 1) { // Not reached to the end of the sequence
            //System.out.println("extend "+position);
            if (genomeDb.get_code(genome, sequence, position + 1) > 3) { // hit a degenerate region
                ++position;
                begin = position - K + 1;
                node.setProperty("length", len);
                node.setProperty("last_kmer",last_kmer);                
                jump();
                int[] add = new int[]{genome,sequence,begin};
                create_degenerate(add);
                connect(curr_node ,degenerate_node, RelTypes.FF); 
                ++num_edges;
                curr_node = degenerate_node;
                //curr_side = 0; // we have set it zero already in create()
                degenerate = true;
                break;
            }
            next_kmer(genomeDb);
            curr_index = indexDb.find(k_mer);
            indexDb.get_pointer(pointer, curr_index);
            if (pointer.node_id == -1L) {
                indexDb.put_next_index(curr_index,last_kmer);
                ++len;
                //node.setProperty("length",(int)node.getProperty("length")+1);
                pointer.node_id = id;
                pointer.canonical = fwd_kmer.get_canonical();
                pointer.offset = len - K;//(int)node.getProperty("length")-K;
                indexDb.put_pointer(pointer, curr_index);
                last_kmer=curr_index;
            } else {
                broke = true;
                break;
            }
        }
        if (!degenerate) {
            node.setProperty("length", len);
            node.setProperty("last_kmer",last_kmer);
        }
        if (!broke && position == seq_len - 1) {
            finish = true;
        }
    }
    
    /**
     * Initializes a new node.
     */
    private void create() {
        int[] address;
        address= new int[]{genome,sequence,position - K + 1};
        ++num_nodes;
        new_node = graphDb.createNode(nucleotide_label);
        if (DEBUG) System.out.println("create "+new_node.getId());
        new_node.setProperty("address", address);
        new_node.setProperty("length", K);
        new_node.setProperty("last_kmer",curr_index);
        new_node.setProperty("first_kmer",curr_index);
    // Set the pointer to the Kmer in the pointer database    
        pointer.node_id = new_node.getId();
        pointer.canonical = fwd_kmer.get_canonical();
        pointer.offset = 0;
        indexDb.put_pointer(pointer, curr_index);
        connect(curr_node ,new_node, RelTypes.values()[curr_side*2]);
        ++num_edges;
        curr_node = new_node;
        curr_side = 0;
    }

    /**
     * Enters the node in which the current Kmer is found in and performs zero, one or two splits.   
     */
    private void follow_forward() {
        int l, pos, begin, g, s, loc, side;
        Node node, split_node1, split_node2,des, src;
        RelationshipType rel_type;
        int[] address;
        boolean degenerated, loop, repeated_edge;
        pos = pointer.offset;
        node = graphDb.getNodeById(pointer.node_id);
        if (DEBUG) System.out.println("follow_forward "+pointer.node_id+" at "+pos);
    // The first split might be done to seperate the part we need to enter in.
        if (pos > 0) {  
            if (DEBUG) System.out.println("first_split "+node.getId()+" at "+pos);
            split_node1 = split(node, pos);
            if (loop = (curr_node.equals(node) && curr_side == 0)) 
                src = split_node1;
            else 
                src = curr_node;
            node = split_node1; // Note : assigning reference variables is dangerous! if split_node changes node will change as well.
        } else {
            split_node1 = node;
            if (loop = (curr_node.equals(node) && curr_side == 0)) 
                src = split_node1;
            else 
                src = curr_node;
        }
        des = split_node1;
        side=curr_side*2;
        curr_side = 0;
        l = (int) node.getProperty("length") - K;
        address = (int[]) node.getProperty("address");
        g = address[0];
        s = address[1];
        loc = address[2];
        degenerated = false;
    // Follow the shared part
        for (pos = 0; pos <= l && position <= seq_len - 1 && genomeDb.get_code(g, s, loc + pos + K - 1) == genomeDb.get_code(genome, sequence, position); ++pos) {
            ++position;
         // If hit a degenarate region aplit and branch to a degenerate node 
            if (position <= seq_len - 1 && genomeDb.get_code(genome, sequence, position) > 3) {
                begin = position - K + 1;
                jump();
                if (pos + 1 <= l) {
                    split_node2 = split(node, pos + 1);
                    if (loop)
                        src = split_node2;
                }
                degenerated = true;
                int[] add = new int[]{genome,sequence,begin};
                create_degenerate(add);
                connect(node ,degenerate_node, RelTypes.FF);
                ++num_edges;
                break;
            }
        }
        if (position == seq_len) {
            finish = true;
        } else if (!degenerated) { // build the Kmer of difference 
            position -= K;
            initial_kmers(genomeDb);
        }
    //  A second split might be needed   
        if (!degenerated && pos <= l) {
            if (DEBUG) System.out.println("second_split "+node.getId()+" at "+pos);
            split_node2 = split(node, pos);
            if (loop)
                src = split_node2;
        }
    // connect the current node before doing splits to the split_node1    
        rel_type = RelTypes.values()[side];
        repeated_edge = false;
        for (Relationship r: src.getRelationships(rel_type, Direction.OUTGOING))
            if (r.getEndNode().equals(des)){
                repeated_edge = true;
                break;
            }
        if (!repeated_edge){
            connect(src ,des, rel_type);
            ++num_edges;
        }
        if (degenerated) {
            curr_node = degenerate_node;
            curr_side = 0; // not really needed
        } else {
            curr_index = indexDb.find(k_mer);
            indexDb.get_pointer(pointer, curr_index);
            curr_node = node;
        }
    }
    
    /**
     * Enters the forward side of node in which the current Kmer found and performs zero, one or two splits.   
     */
    private void follow_reverse() {
        int pos, begin, g, s, loc, side;
        int[] address;
        Node node, split_node1, split_node2 ,des, src;
        boolean degenerated = false, loop, first_split = false, repeated_edge;
        pos = pointer.offset;
        node = graphDb.getNodeById(pointer.node_id);
        if (DEBUG) System.out.println("follow_reverse "+pointer.node_id+" at "+pos);
        RelationshipType rel_type;
        split_node2 = node; //if the second split does not happens remains unchanged
        if (pos < (int) node.getProperty("length") - K) {
            if (DEBUG) System.out.println("first_split "+node.getId()+" at "+(pos+1));
            first_split = true;
            split_node1 = split(node, pos+1);
            if (loop = curr_node.equals(node) && curr_side == 0) // might be in reverse side due to a follow reverse
                src = split_node1;
            else 
                src = curr_node;
        } else {
            split_node1 = node;
            if (loop = curr_node.equals(node) && curr_side == 0)
                src = split_node1;
            else 
                src = curr_node;
        }
        des = node;
        side=curr_side*2+1;
        curr_side = 1;
        address = (int[]) node.getProperty("address");
        g = address[0];
        s = address[1];
        loc = address[2];
        for (pos = (int) node.getProperty("length") - K; pos >= 0 && position <= seq_len - 1 && genomeDb.get_code(g, s, loc + pos) == genomeDb.get_complement_code(genome, sequence, position); --pos) {
            ++position;
            if (position <= seq_len - 1 && genomeDb.get_code(genome, sequence, position) > 3) {
                begin = position - K + 1;
                jump();
                if (pos > 0) {
                    split_node2 = split(node, pos);
                    des = split_node2;
                    if (!first_split && loop)
                        src = split_node2;
                }
                int[] add = new int[]{genome,sequence,begin};
                create_degenerate(add);
                connect(split_node2, degenerate_node, RelTypes.RF);
                ++num_edges;
                degenerated = true;
                break;
            }
        }
        if (position == seq_len) {
            finish = true;
        } else if (!degenerated) {
            position -= K;
            initial_kmers(genomeDb);
        }
        if (!degenerated && pos >= 0) {
            if (DEBUG) System.out.println("second_split "+node.getId()+" at "+ (pos + 1));
            split_node2 = split(node, pos+1);
            des = split_node2;
            if (!first_split && loop)
                src = split_node2;
        }
        rel_type = RelTypes.values()[side];
        repeated_edge = false;
        for (Relationship r: src.getRelationships(rel_type, Direction.OUTGOING))
            if (r.getEndNode().equals(des)){
                repeated_edge = true;
                break;
            }
        if (!repeated_edge){
            connect(src ,des, rel_type);
            ++num_edges;
        }
        if (degenerated) {
            curr_node = degenerate_node;
            curr_side = 0;
        } else {
            curr_index = indexDb.find(k_mer);
            indexDb.get_pointer(pointer, curr_index);
            curr_node = split_node2;
        }
    }

    /**
     * Jump over an ambiguous region; at first, position points to the first position which degenerate starts, 
     * after jumping it points to the last base of the first K-mer after the ambiguous region. 
     */
    private void jump() {
        int j;
        fwd_code = genomeDb.get_code(genome, sequence, position);
        rev_code = 3 - fwd_code;
        do {
            while (fwd_code > 3 && position < seq_len - 1) {
                ++position;
                fwd_code = genomeDb.get_code(genome, sequence, position);
                rev_code = 3 - fwd_code;
            }
            fwd_kmer.reset();
            rev_kmer.reset();
            fwd_kmer.next_fwd_kmer(fwd_code);
            rev_kmer.next_rev_kmer(rev_code);
            for (j = 0; j < K - 1 && position < seq_len - 1; ++j) {
                ++position;
                fwd_code = genomeDb.get_code(genome, sequence, position);
                rev_code = 3 - fwd_code;
                if (fwd_code > 3) {
                    break;
                }
                fwd_kmer.next_fwd_kmer(fwd_code);
                rev_kmer.next_rev_kmer(rev_code);
            }
            if (j == K - 1) {
                fwd_kmer.set_canonical ( fwd_kmer.compare(rev_kmer) == -1);
                k_mer = fwd_kmer.get_canonical() ? fwd_kmer : rev_kmer;
            }
        } while (fwd_code > 3 && position < seq_len - 1);
        if (position == seq_len - 1){
        // to acheive the right length for the degenerate node    
            finish = true;
            ++position;
        }
    }
    
    /**
     * creates a degenerate node starting at "begin" ending at position-1.
     * @param address The genomic position of the region
     */
    private void create_degenerate(int[] address) {
        ++num_degenerates;
        ++num_nodes;
        degenerate_node = graphDb.createNode(degenerate_label);
        degenerate_node.addLabel(nucleotide_label);
        if (DEBUG) System.out.println("create_degenerate:"+degenerate_node.getId()+" position:"+position+" begin:"+address[2]);
        degenerate_node.setProperty("address", address);
        degenerate_node.setProperty("length", position - address[2]);
        if (!finish) {
            curr_index = indexDb.find(k_mer);
            indexDb.get_pointer(pointer, curr_index);
        }
    }
    
    /**
     * Generates the next forward and reverse kmers; the canonical one would be the next kmer.
     * @param s_db The genome database we are reading nuleotides from.
     */
    private void next_kmer(SequenceDatabase s_db) {
        ++position;
        fwd_code = s_db.get_code(genome, sequence, position);
        rev_code = 3 - fwd_code;
        fwd_kmer.next_fwd_kmer(fwd_code);
        rev_kmer.next_rev_kmer(rev_code);
        fwd_kmer.set_canonical( fwd_kmer.compare(rev_kmer) == -1 );
        k_mer = fwd_kmer.get_canonical() ? fwd_kmer : rev_kmer;
        if (!DEBUG && position % (seq_len / 100 + 1) == 0) 
            System.out.print((long) position * 100 / seq_len + 1 + "%\r");
        if (SHOW_KMERS) System.out.println(k_mer.toString());
    }

    /**
     * initializes the first K-mer of the genome;  
     * It might jump over the degenerate regions creating a degenerate_node
     * @param s_db The genome database we are reading nuleotides from.
     */
    public void initial_kmers(SequenceDatabase s_db) {
        int i;
        fwd_kmer.reset();
        rev_kmer.reset();
        for (i = 0; i < K && position < seq_len - 1; ++i) {
            if (s_db.get_code(genome, sequence, position + 1) > 3) {
                ++position;
                jump();
                int[] add = new int[]{genome,sequence,0};
                try (Transaction tx = graphDb.beginTx()) {
                    create_degenerate(add);
                    connect(curr_node ,degenerate_node, RelTypes.values()[curr_side*2]);
                    tx.success();
                }
                ++num_edges;
                curr_node = degenerate_node;
                break;
            }
            next_kmer(s_db);
        }
        fwd_kmer.set_canonical( fwd_kmer.compare(rev_kmer) == -1 );
        k_mer = fwd_kmer.get_canonical() ? fwd_kmer : rev_kmer;
        //System.out.println(k_mer.toString());
    }

    /**
     * 
     Constructs the pangenome out of the provided input sequences.
     */
    void construct_pangenome(int previous_num_genomes) {
        int trsc = 0;
        Node genome_node;
        Node[][] sequence_nodes = new Node[genomeDb.num_genomes + 1][];
        phaseTime = System.currentTimeMillis();
        fwd_kmer = new kmer(K, indexDb.get_pre_len(), indexDb.get_suf_len());
        rev_kmer = new kmer(K, indexDb.get_pre_len(), indexDb.get_suf_len());
        Transaction tx = graphDb.beginTx();
        try {
            for (genome = previous_num_genomes + 1; genome <= genomeDb.num_genomes; ++genome) {
                System.out.println("Processing genome " + genome + " :             ");
                genome_node = graphDb.createNode(genome_label);
                genome_node.setProperty("path", genomeDb.genome_names[genome]);
                genome_node.setProperty("number", genome);
                genome_node.setProperty("num_sequences", genomeDb.num_sequences[genome]);
                genome_node.setProperty("date", new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(new Date()));
                db_node.createRelationshipTo(genome_node, RelTypes.has);
                sequence_nodes[genome] = new Node[genomeDb.num_sequences[genome] + 1];
                for (sequence = 1; sequence <= genomeDb.num_sequences[genome]; ++sequence) {
                    sequence_nodes[genome][sequence] = curr_node = graphDb.createNode(sequence_label);
                    sequence_nodes[genome][sequence].setProperty("genome", genome);
                    sequence_nodes[genome][sequence].setProperty("number", sequence);
                    sequence_nodes[genome][sequence].setProperty("identifier", genome + "_" + sequence);
                    sequence_nodes[genome][sequence].setProperty("title", genomeDb.sequence_titles[genome][sequence]);
                    sequence_nodes[genome][sequence].setProperty("length", genomeDb.sequence_length[genome][sequence]);
                    sequence_nodes[genome][sequence].setProperty("offset", genomeDb.sequence_offset[genome][sequence]);
                    genome_node.createRelationshipTo(sequence_nodes[genome][sequence], RelTypes.has);
                    finish = false;
                    System.out.println("sequence " + sequence + "/" + genomeDb.num_sequences[genome] + " of genome " + genome + "\tlength=" + genomeDb.sequence_length[genome][sequence]);
                    curr_node = sequence_nodes[genome][sequence];
                    curr_side = 0;
                    position = -1;
                    seq_len = genomeDb.sequence_length[genome][sequence];
                    initial_kmers(genomeDb);
                    if (!finish){
                        curr_index = indexDb.find(k_mer);
                        indexDb.get_pointer(pointer, curr_index);
                        while (!finish) {
                            if (pointer.node_id == -1L) // kmer is new
                            {
                                create();
                                extend(curr_node);
                            } else if (fwd_kmer.get_canonical() ^ pointer.canonical)// if sides don't agree
                            {
                                follow_reverse();
                            } else {
                                follow_forward();
                            }
                            ++trsc;
                            if (trsc >= MAX_TRANSACTION_SIZE){    
                                tx.success();
                                tx.close();
                                tx = graphDb.beginTx();
                                trsc = 0;
                            }
                        }
                    }
                    connect(curr_node, sequence_nodes[genome][sequence], RelTypes.values()[curr_side*2]);// to point to the last k-mer of the sequence located in the other strand
                    ++num_edges;
                }//sequences
                System.out.println((System.currentTimeMillis() - phaseTime) / 1000 + " seconds elapsed.");
            }//genomes
            tx.success();
        } finally {
            tx.close();
        }
        add_sequence_properties();
        localize_nodes();
    }
    
    /**
     * To add list of anchor nodes, anchor sides and anchor positions to each sequence_node.
     * These are used for locating genomic regions very quickly. 
     */
    void localize_nodes() {
        ResourceIterator<Node> sequence_iterator;
        LinkedList<Node> sequence_nodes;
        int anchors_distance, trsc = 0, i, len, m, neighbor_length = 0, count;
        char node_side, neighbor_side;
        long length;
        long[] anchor_nodes;
        int[] anchor_positions;
        int[] initial_coordinate = new int[1];
        Node node, neighbor = null, sequence_node;
        StringBuilder nds = new StringBuilder();
        StringBuilder pos = new StringBuilder();
        StringBuilder sds = new StringBuilder();
        String[] ids_list, posis_list;
        String rel_name,origin;
        int[] positions;
        int[] new_positions;
        int[] address = new int[3], addr = null;
        boolean is_node = false, is_degenerate = false, found = true;
        try (Transaction tx = graphDb.beginTx()){
            sequence_iterator = graphDb.findNodes(sequence_label);
            sequence_nodes = new LinkedList();
            while (sequence_iterator.hasNext())
                sequence_nodes.add(sequence_iterator.next());
            sequence_iterator.close();
            tx.success();
        }

        Transaction tx = graphDb.beginTx();
        try {
            while (!sequence_nodes.isEmpty()){
                sequence_node = sequence_nodes.remove();
                origin = "G" + ((String)sequence_node.getProperty("identifier")).replace('_', 'S');
                address[0] = (int)sequence_node.getProperty("genome");
                address[1] = (int)sequence_node.getProperty("number");
                System.out.println("\rLocalizing sequence "+address[1] + "/" + genomeDb.num_sequences[address[0]] + " of genome " + address[0] + "                        ");
                length = genomeDb.sequence_length[address[0]][address[1]] - 1;
                anchors_distance = (int)num_nodes / ANCHORS + 1;
                node = sequence_node;
                node_side = 'F';
                count = 0;
                for (address[2] = 0; address[2] + K - 1 <= length && found;){ // K-1 bases of the last node not added
                    //System.out.println((address[2] + K - 1)+" ? " + length);
                    found = false;
                    for (Relationship r : node.getRelationships(Direction.OUTGOING)) {
                        rel_name = r.getType().name();
                        if (rel_name.charAt(0) != node_side)
                            continue;
                        neighbor = r.getEndNode();
                        neighbor_side = rel_name.charAt(1);
                        is_node = neighbor.hasLabel(nucleotide_label) && !neighbor.hasLabel(degenerate_label);
                        is_degenerate = neighbor.hasLabel(degenerate_label);
                        if (is_node || is_degenerate){
                            addr = (int[]) neighbor.getProperty("address");
                            neighbor_length = (int) neighbor.getProperty("length");
                        }
                        //System.out.println(node.getId()+" "+address[2]+" "+node_side);
                        //System.out.println(neighbor.getId()+" "+addr[2]+" "+neighbor_side);
                        if ((is_node && genomeDb.compare(genomeDb, address, addr, K - 1,
                                neighbor_side == 'F' ? K - 1 : neighbor_length - K, 1, neighbor_side == 'F'))
                                || (is_degenerate && Arrays.equals(addr, address))) {
                            //System.out.println("found "+address[2]+" "+neighbor.getId());
                            found = true;
                            positions = (int[]) r.getProperty(origin, null);
                            if (positions != null) {
                                len = positions.length;
                                new_positions = new int[len + 1];
                                for (i = 0; i < len; ++i) {
                                    new_positions[i] = positions[i];
                                }
                                new_positions[i] = address[2];
                                r.setProperty(origin, new_positions);
                                if (len > 1000)
                                    neighbor.addLabel(low_complexity_label);
                            } else {
                                initial_coordinate[0] = address[2];
                                r.setProperty(origin, initial_coordinate);
                            }
                            if (count % anchors_distance == 0) {
                                nds.append(neighbor.getId()).append(" ");
                                sds.append(neighbor_side);
                                pos.append(address[2]).append(" ");
                            }
                            count++;
                            address[2] = address[2] + neighbor_length - K + 1;
                            node = neighbor;
                            node_side = neighbor_side;
                            break;
                        }
                    }
                    ++trsc;
                    if (trsc >= 1000 * MAX_TRANSACTION_SIZE){
                        tx.success();
                        tx.close();
                        tx = graphDb.beginTx();
                        trsc = 0;
                    }
                    System.out.print("%" + address[2] * 100 / length + "\t\r");
                }
                if (!found) {
                    System.out.println("Could not locate position " + address[2] + " from node ID=" + node.getId());
                    System.exit(1);
                }
                m = sds.length();
                ids_list = nds.toString().split("\\s");
                posis_list = pos.toString().split("\\s");
                anchor_nodes = new long[m];
                anchor_positions = new int[m];
                for (i = 0; i < m; ++i) {
                    anchor_nodes[i] = Long.valueOf(ids_list[i]);
                    anchor_positions[i] = Integer.valueOf(posis_list[i]);
                }
                sequence_node.setProperty("anchor_nodes", anchor_nodes);
                sequence_node.setProperty("anchor_positions", anchor_positions);
                sequence_node.setProperty("anchor_sides", sds.toString());
                nds.setLength(0);
                pos.setLength(0);
                sds.setLength(0);
            }//while
            System.out.println((System.currentTimeMillis() - phaseTime) / 1000 + " seconds elapsed.");
            tx.success();
        } finally {
           tx.close();
        }
        System.out.println();
    }

    /**
     * Extracts the sequence of the nodes from the genome database and store it as a property of the node.
     */
    void add_sequence_properties() {
        int trsc = 0, node_length;
        int[] addr;
        Node node;
        ResourceIterator<Node> nodes_iterator;
        LinkedList<Node> nodes = new LinkedList();
        System.out.println("Adding sequence to the nodes...");
        try(Transaction tx = graphDb.beginTx()){
            nodes_iterator = graphDb.findNodes(nucleotide_label);
            while (nodes_iterator.hasNext()){
                nodes.add(nodes_iterator.next());
            }
            tx.success();
            nodes_iterator.close();
        }
        Transaction tx = graphDb.beginTx();
        try {
        //num_bases = K - 1; // for the missed overlapped of the last node of each sequence which will not be stored 
            while (!nodes.isEmpty()){
                node = nodes.remove();
                addr = (int[]) node.getProperty("address");
                node_length = (int) node.getProperty("length");
                num_bases += node_length;
                //num_bases += node_length - K + 1;
                //node.setProperty("sequence", genomeDb.get_sequence(addr[0], addr[1], addr[2], node_length - K + 1, true).toString());
                node.setProperty("sequence", genomeDb.get_sequence(addr[0], addr[1], addr[2], node_length, true).toString());
                ++trsc;
                if (trsc >= MAX_TRANSACTION_SIZE){
                    tx.success();
                    tx.close();
                    tx = graphDb.beginTx();
                    trsc = 0;
                }
            }
            tx.success();
        } finally {
            tx.close();
        }
    }

    /**
     * Remove the given property from all the nodes and degenerates.
     * @param property 
     */    
    void drop_nodes_property(String property) {
        int i;
        ResourceIterator<Node> nodes_iterator;
        LinkedList<Node> nodes = new LinkedList();
        Node node;
        try(Transaction tx = graphDb.beginTx()){
            nodes_iterator = graphDb.findNodes(nucleotide_label);
            while (nodes_iterator.hasNext()){
                nodes.add(nodes_iterator.next());
            }
            tx.success();
            nodes_iterator.close();
        }
        while (!nodes.isEmpty()) {
            try (Transaction tx = graphDb.beginTx()) {
                for (i = 0; i < MAX_TRANSACTION_SIZE && !nodes.isEmpty(); ++i) {
                    node = nodes.remove();
                    node.removeProperty(property);
                }
                tx.success();
            }
        }
    }

    /**
     * Remove the occurrence arrays of the edges.
     */    
    void drop_edges_colors() {
        int i;
        ResourceIterator<Relationship> rels;
        Relationship r;
        try (Transaction tx = graphDb.beginTx()) {
            rels = graphDb.getAllRelationships().iterator();
            tx.success();
        }
        while (rels.hasNext()) {
            try (Transaction tx = graphDb.beginTx()) {
                for (i = 0; i < MAX_TRANSACTION_SIZE && rels.hasNext(); ++i) {
                    r = rels.next();
                    if (r.isType(RelTypes.FF) || r.isType(RelTypes.FR) || r.isType(RelTypes.RF) || r.isType(RelTypes.RR))
                        for(String p:r.getPropertyKeys())
                            r.removeProperty(p);
                }
                tx.success();
            }
        }
        rels.close();
    }

    /**
     * Shuts down the graph database if the program halts unexpectedly.
     * 
     * @param graphDb The graph database object 
     */
    private static void registerShutdownHook(final GraphDatabaseService graphDb) {
        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                graphDb.shutdown();
            }
        });
    }

    /**
     * Returns size of a given folder.
     * 
     * @param dir The folder File object.
     * @return Size of the folder in MegaBytes
     */
    public static long getFolderSize(File dir) {
        long size = 0;
        for (File file : dir.listFiles()) {
            if (file.isFile()) {
                // System.out.println(file.getName() + " " + file.length());
                size += file.length();
            } else {
                size += getFolderSize(file);
            }
        }
        return size / 1048576 + 1;
    }
}