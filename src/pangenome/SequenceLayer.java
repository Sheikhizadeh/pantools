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
import java.io.Console;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Map;
import java.util.Map.Entry;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.io.fs.FileUtils;
import org.neo4j.graphdb.Direction;
import org.neo4j.graphdb.QueryExecutionException;
import org.neo4j.graphdb.Result;
import static pantools.Pantools.GENOME_DATABASE_PATH;
import static pantools.Pantools.GRAPH_DATABASE_PATH;
import static pantools.Pantools.INDEX_DATABASE_PATH;
import static pantools.Pantools.K;
import static pantools.Pantools.PATH;
import pantools.Pantools.RelTypes;
import static pantools.Pantools.degenerate_label;
import static pantools.Pantools.gene_label;
import static pantools.Pantools.genomeDb;
import static pantools.Pantools.genome_label;
import static pantools.Pantools.graphDb;
import static pantools.Pantools.indexDb;
import static pantools.Pantools.node_label;
import static pantools.Pantools.num_bases;
import static pantools.Pantools.num_degenerates;
import static pantools.Pantools.num_edges;
import static pantools.Pantools.num_nodes;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.phaseTime;
import static pantools.Pantools.sequence_label;
import static pantools.Pantools.startTime;
import static pantools.Pantools.db_trsc_limit;

/**
 *
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class SequenceLayer {

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
    private byte[] byte_pointer;
    private int fwd_code, rev_code;
    private boolean finish;
    private static int anchor_distance = 100; // Distance between two anchor nodes

    public SequenceLayer() {
        pointer = new IndexPointer();
        byte_pointer = new byte[index.IndexDatabase.ptr_len];
        finish = false;
    }

    /**
     * @param args the command line arguments build : To build a pan-genome out
     * of a set of genomes add : To add new genomes to an available pan-genome
     * genes : To extract gene sequence according to their IDs group : To add
     * group nodes pointing to genes comp : To compare two pan-genomes
     * @throws java.io.IOException
     */
    /*
     To builds a pan-genome database out of a set of input genomes.
     file : a text file containing paths to FASTA files
     */
    /**
     *
     * @param genome_paths
     * @param args the command line arguments build : To build a pan-genome out
     * of a set of genomes add : To add new genomes to an available pan-genome
     * genes : To extract gene sequence according to their IDs group : To add
     * group nodes pointing to genes comp : To compare two pan-genomes
     * @throws java.io.IOException
     */
    public void build(String genome_paths) {
        int i, j;
        File theDir = new File(PATH);
        if (theDir.exists()) {
            try {
                FileUtils.deleteRecursively(new File(PATH)); // deletes old files
            } catch (IOException ioe) {
            }
        } else {
            try {
                theDir.mkdir();
            } catch (SecurityException se) {
                System.out.println("Failed to create " + PATH);
                System.exit(1);
            }
        }
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH + GRAPH_DATABASE_PATH))
                .setConfig("keep_logical_logs", "100M size").newGraphDatabase();
        registerShutdownHook(graphDb);
        startTime = System.currentTimeMillis();
        num_nodes = 0;
        num_edges = 0;
        try (Transaction tx = graphDb.beginTx()) {
            db_node = graphDb.createNode(pangenome_label);
            db_node.setProperty("k_mer_size", K);
            tx.success();
        }
        genomeDb = new SequenceDatabase(PATH + GENOME_DATABASE_PATH, genome_paths);
        indexDb = new IndexDatabase(PATH + INDEX_DATABASE_PATH, genome_paths, K, genomeDb);
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
        File directory = new File(PATH + GRAPH_DATABASE_PATH);
        for (File f : directory.listFiles()) {
            if (f.getName().startsWith("neostore.transaction.db.")) {
                f.delete();
            }
        }
        System.out.println("graph.db size: " + getFolderSize(new File(PATH + GRAPH_DATABASE_PATH)) + " MB");
        System.out.println("index.db size: " + getFolderSize(new File(PATH + INDEX_DATABASE_PATH)) + " MB");
        System.out.println("genome.db size: " + getFolderSize(new File(PATH + GENOME_DATABASE_PATH)) + " MB");
    }
    /*
     To adds new genomes to an available pan-genome database.
     genome_paths : a text file containing paths to FASTA files
     */

    public void add(String genome_paths) {
        int i, j, g, s, len, previous_num_genomes;
        long byte_number = 0;
        Node start, seq_node;
        if (new File(PATH + GRAPH_DATABASE_PATH).exists()) {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH + GRAPH_DATABASE_PATH))
                    .setConfig("keep_logical_logs", "100M size").newGraphDatabase();
            registerShutdownHook(graphDb);
            startTime = System.currentTimeMillis();
            try (Transaction tx = graphDb.beginTx()) {
                db_node = graphDb.findNodes(pangenome_label).next();
                if (db_node == null) {
                    System.out.println("Can not locate database node!");
                    System.exit(1);
                }
                K = (int) db_node.getProperty("k_mer_size");
                num_nodes = (int) db_node.getProperty("num_nodes");
                num_edges = (int) db_node.getProperty("num_edges");
                num_degenerates = (int) db_node.getProperty("num_degenerate_nodes");
                num_bases = (int) db_node.getProperty("num_bases");
                previous_num_genomes = (int) db_node.getProperty("num_genomes");
                if (!Files.exists(Paths.get(PATH + GENOME_DATABASE_PATH))) // if genome database has been lost
                {
                    // read genomes information from the graph and rebuild the genomes database
                    genomeDb = new SequenceDatabase(PATH + GENOME_DATABASE_PATH, graphDb);
                    StringBuilder seq = new StringBuilder();
                    for (g = 1; g <= genomeDb.num_genomes; ++g) {
                        for (s = 1; s <= genomeDb.num_sequences[g]; ++s) {
                            seq_node = graphDb.findNode(sequence_label, "number", g + "_" + s);
                            start = seq_node.getRelationships(Direction.OUTGOING).iterator().next().getEndNode();
                            extract_sequence(seq, new IndexPointer(start.getId(), true, 0, -1), new int[]{g,s,-1}, 0, (int) genomeDb.sequence_length[g][s] - 1);
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
                    genomeDb = new SequenceDatabase(PATH + GENOME_DATABASE_PATH); //Readable only
                } else {
                    genomeDb = new SequenceDatabase(PATH + GENOME_DATABASE_PATH);
                }
                genomeDb.add_genomes(PATH + GENOME_DATABASE_PATH, genome_paths);
                indexDb = new IndexDatabase(PATH + INDEX_DATABASE_PATH, genome_paths, genomeDb, graphDb, previous_num_genomes);
                tx.success();
            }
            drop_property("sequence");// these should be dropped out as include-genomes function recreate them
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
            File directory = new File(PATH + GRAPH_DATABASE_PATH);
            for (File f : directory.listFiles()) {
                if (f.getName().startsWith("neostore.transaction.db.")) {
                    f.delete();
                }
            }
            System.out.println("graph.db size: " + getFolderSize(new File(PATH + GRAPH_DATABASE_PATH)) + " MB");
            System.out.println("index.db size: " + getFolderSize(new File(PATH + INDEX_DATABASE_PATH)) + " MB");
            System.out.println("genome.db size: " + getFolderSize(new File(PATH + GENOME_DATABASE_PATH)) + " MB");
        } else {
            System.out.println("No database found in " + PATH);
            System.exit(1);
        }
    }

    /**
     * Retrieves sequence of the genes sequence from the pangenome and stores them in a FASTA file. 
     * 
     * @param annotation_records_file a text file containing annotation records of the genes to be retrieved.
     */
    public void retrieve_genes(String annotation_records_file) {
        if (new File(PATH + GRAPH_DATABASE_PATH).exists()) {
            BufferedReader in;
            ResourceIterator<Node> gene_nodes;
            Relationship rstart;
            Node start, gene;
            String origin, record;
            String line;
            boolean strand;
            int i, j, begin, end, num_genes;
            int[] address;
            String[] fields, origin_fields;
            StringBuilder gene_seq;

            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH + GRAPH_DATABASE_PATH))
                    .setConfig("keep_logical_logs", "100M size").newGraphDatabase();
            registerShutdownHook(graphDb);
            startTime = System.currentTimeMillis();
            try (Transaction tx = graphDb.beginTx()) {
                K = (int) graphDb.findNodes(pangenome_label).next().getProperty("k_mer_size");
                tx.success();
            }
            num_genes = 0;
            address = new int[3];
            gene_seq = new StringBuilder();
            try {
                in = new BufferedReader(new FileReader(annotation_records_file));
                while (in.ready()) {
                    line = in.readLine();
                    if (line.equals("")) {
                        continue;
                    }
                    num_genes++;
                }
                in.close();
            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            }
            String[] records = new String[num_genes];
            try (Transaction tx = graphDb.beginTx()) {
                try {
                    in = new BufferedReader(new FileReader(annotation_records_file));
                    // fill all the records in an array to be sorted    
                    for (i = 0; in.ready();) {
                        line = in.readLine();
                        if (line.equals("")) {
                            continue;
                        }
                        fields = line.split("\\t");
                        records[i] = fields[fields.length - 1]; // attribute field
                        ++i;
                    }
                    in.close();
                    Arrays.sort(records);
                } catch (IOException ioe) {
                    System.out.println("Failed to read file names!");
                    System.exit(1);
                }
                try {
                    BufferedWriter out = new BufferedWriter(new FileWriter(annotation_records_file + ".fasta"));
                    // for all the genes in the database    
                    for (i = j = 0, gene_nodes = graphDb.findNodes(gene_label); gene_nodes.hasNext();) {
                        gene = gene_nodes.next();
                        record = (String) gene.getProperty("attribute");
                        if (record != null && Arrays.binarySearch(records, record) >= 0) // gene is in the records
                        {
                            ++i;
                            //System.out.println(record);
                            start = graphDb.getNodeById((long)gene.getProperty("start_node_id"));
                            rstart = graphDb.getRelationshipById((long)gene.getProperty("start_edge_id"));
                            origin = (String) gene.getProperty("origin");
                            origin_fields = origin.split("_");
                            address[0] = Integer.parseInt(origin_fields[0]);
                            address[1] = Integer.parseInt(origin_fields[1]);
                            begin = (int) gene.getProperty("begin");
                            end = (int) gene.getProperty("end");
                            strand = gene.getProperty("strand").toString().equals("+");
                            extract_sequence(gene_seq, new IndexPointer(start.getId(), (boolean) rstart.getProperty("forward"), (int) rstart.getProperty("position"), -1l), address, begin - 1, end - 1);//
                            //genomeDb=new sequence_database(PATH+GENOME_DATABASE_PATH);
                            //if(gene_seq.toString().equals(genomeDb.get_sequence(genome, sequence, begin-1, end-begin+1, strand))
                            //|| gene_seq.toString().equals(genomeDb.get_sequence(genome, sequence, begin-1, end-begin+1, !strand)) )//gene_seq.length() == end-begin+1)//
                            if (gene_seq.length() == end - begin + 1) {
                                ++j;
                                out.write(">" + record + "\n");
                                if (strand) {
                                    write_fasta(out, gene_seq.toString(), 70);
                                } else {
                                    write_fasta(out, reverse_complement(gene_seq.toString()), 70);
                                }
                            } else {
                                System.out.println("Failed to assemble:\n" + record);
                            }
                            gene_seq.setLength(0);
                        }
                        if (i % (num_genes / 100 + 1) == 0) {
                            System.out.print((long) i * 100 / num_genes + 1 + "%\r");
                        }
                    }//for i
                    System.out.println(j + " out of " + i + " found genes retrieved successfully.");
                    out.close();
                } catch (IOException ioe) {
                    System.out.println("Failed to read file names!");
                    System.exit(1);
                }
                tx.success();
            }
            graphDb.shutdown();
        } else {
            System.out.println("No database found in " + PATH);
            System.exit(1);
        }
    }

    /**
     * Retrieves the sequence of a number of genomic regions from the pangenome and stores them in a FASTA file.
     * 
     * @param region_records a text file with lines containing genome number, sequence number, start and stop positions
     *        of the genomic regions seperated by one space.
     */
    public void retrieve_regions(String region_records) {
        if (new File(PATH + GRAPH_DATABASE_PATH).exists()) {
            String[] fields;
            String line;
            IndexPointer start_ptr;
            Node start_node;
            StringBuilder seq;
            int c, g, s, from, to, num_regions;
            
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH + GRAPH_DATABASE_PATH))
                    .setConfig("keep_logical_logs", "100M size").newGraphDatabase();
            registerShutdownHook(graphDb);
            try (Transaction tx = graphDb.beginTx()) {
                K = (int) graphDb.findNodes(pangenome_label).next().getProperty("k_mer_size");
                tx.success();
            }
            seq = new StringBuilder();
            num_regions = 0;
            try {
                BufferedReader in = new BufferedReader(new FileReader(region_records));
                while (in.ready()) {
                    line = in.readLine();
                    if (line.equals("")) {
                        continue;
                    }
                    num_regions++;
                }
                in.close();
            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            }
            startTime = System.currentTimeMillis();
            try (Transaction tx = graphDb.beginTx()) {
                try (BufferedReader in = new BufferedReader(new FileReader(region_records))) {
                    BufferedWriter out = new BufferedWriter(new FileWriter(region_records + ".fasta"));
                    for (c = 0; in.ready();) {
                        line = in.readLine();
                        if (line.equals("")) {
                            continue;
                        }
                        fields = line.split("\\s");
                        g = Integer.parseInt(fields[0]);
                        s = Integer.parseInt(fields[1]);
                        from = Integer.parseInt(fields[2]);
                        to = Integer.parseInt(fields[3]);
                        start_ptr = locate(g, s,from);
                        start_node = graphDb.getNodeById(start_ptr.node_id);
                        extract_sequence(seq, new IndexPointer(start_node.getId(), start_ptr.canonical, start_ptr.position, -1l), new int[]{g,s,-1}, from, to);
                        if (seq.length() != to - from + 1) {
                            System.out.println("Failed to assemble region!");
                        }
                        out.write(">genome:" + g + " sequence:" + s + " from:" + from + " to:" + to + " length:" + seq.length() + "\n");
                        write_fasta(out, seq.toString(), 70);
                        seq.setLength(0);
                        ++c;
                        if (c % (num_regions / 100 + 1) == 0) {
                            System.out.print((long) c * 100 / num_regions + 1 + "%\r");
                        }
                    }
                    in.close();
                    out.close();
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
     * @param genome_records a text file containing the number and a given name for each genome. 
     */
    public void reconstruct_genomes(String genome_records_file) {
        if (new File(PATH + GENOME_DATABASE_PATH).exists()) {
            BufferedReader in;
            BufferedWriter out;
            Node seq_node, start;
            String[] fields;
            String line;
            int[] address;
            StringBuilder seq;
            
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH + GRAPH_DATABASE_PATH))
                    .setConfig("keep_logical_logs", "100M size").newGraphDatabase();
            registerShutdownHook(graphDb);
            startTime = System.currentTimeMillis();
            genomeDb = new SequenceDatabase(PATH + GENOME_DATABASE_PATH, graphDb);
            address = new int[3];
            seq = new StringBuilder();
            try (Transaction tx = graphDb.beginTx()) {
                db_node = graphDb.findNodes(pangenome_label).next();
                if (db_node == null) {
                    System.out.println("Can not locate database node!");
                    System.exit(1);
                }
                K = (int) db_node.getProperty("k_mer_size");
                if (genome_records_file.equals("all")) {
                    for (address[0] = 1; address[0] <= genomeDb.num_genomes; ++address[0]) {
                        System.out.println("Reconstructing genome " + address[0] + "...");
                        try {
                            out = new BufferedWriter(new FileWriter(PATH + "/genome_" + address[0] + ".fasta"));
                            for (address[1] = 1; address[1] <= genomeDb.num_sequences[address[0]]; ++address[1]) {
                                out.write(">" + genomeDb.sequence_titles[address[0]][address[1]] + "\n");
                                seq_node = graphDb.findNode(sequence_label, "number", address[0] + "_" + address[1]);
                                start = seq_node.getRelationships(Direction.OUTGOING).iterator().next().getEndNode();
                                extract_sequence(seq, new IndexPointer(start.getId(), true, 0, -1l), address, 0, (int) genomeDb.sequence_length[address[0]][address[1]] - 1);
                                write_fasta(out, seq.toString(), 80);
                                seq.setLength(0);
                            }
                            out.close();
                        } catch (IOException e) {
                            System.out.println(e.getMessage());
                            System.exit(1);
                        }
                    }
                } else {
                    try {
                        in = new BufferedReader(new FileReader(genome_records_file));
                        while (in.ready()) {
                            line = in.readLine();
                            if (line.equals("")) {
                                continue;
                            }
                            fields = line.split(" ");
                            address[0] = Integer.parseInt(fields[0]);
                            System.out.println("Reconstructing genome " + fields[1] + "...");
                            try {
                                out = new BufferedWriter(new FileWriter(PATH + (fields.length > 1 ? "/" + fields[1] : "/genome_" + address[0]) + ".fasta"));
                                for (address[1] = 1; address[1] <= genomeDb.num_sequences[address[0]]; ++address[1]) {
                                    out.write(">" + genomeDb.sequence_titles[address[0]][address[1]] + "\n");
                                    seq_node = graphDb.findNode(sequence_label, "number", address[0] + "_" + address[1]);
                                    start = seq_node.getRelationships(Direction.OUTGOING).iterator().next().getEndNode();
                                    extract_sequence(seq, new IndexPointer(start.getId(), true, 0, -1l), address, 0, (int) genomeDb.sequence_length[address[0]][address[1]] - 1);
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
                }
                tx.success();
            }
            System.out.println("Genomes are ready in " + PATH);
            graphDb.shutdown();
            genomeDb.close();
        } else {
            System.out.println("No database found in " + PATH);
            System.exit(1);
        }
    }
    
    /**
     * Compares the topology of two pangenomes. 
     * @param path1 Path to the first pangenome
     * @param path2 Path to the second pangenome
     * @return TRUE if pangenomes are the same, FALSE otherwise.
     */
    public boolean compare_pangenomes(String path1, String path2) {
        int sq_num1, sq_num2, ed_num1, ed_num2, K1, K2, ng1, ng2, len1, len2;
        long k, knum1, knum2;
        IndexDatabase indexDb1;
        IndexDatabase indexDb2;
        SequenceDatabase genomeDb1;
        SequenceDatabase genomeDb2;
        Node db_node1, db_node2;
        
        if (new File(path1 + GRAPH_DATABASE_PATH).exists() && new File(path2 + GRAPH_DATABASE_PATH).exists()) {
            GraphDatabaseService g1 = new GraphDatabaseFactory().newEmbeddedDatabase(new File(path1 + GRAPH_DATABASE_PATH));
            GraphDatabaseService g2 = new GraphDatabaseFactory().newEmbeddedDatabase(new File(path2 + GRAPH_DATABASE_PATH));
            registerShutdownHook(g1);
            registerShutdownHook(g2);
            try (Transaction tx1 = g1.beginTx(); Transaction tx2 = g2.beginTx();) {
                db_node1 = g1.findNodes(pangenome_label).next();
                db_node2 = g2.findNodes(pangenome_label).next();
                if (db_node1 == null) {
                    System.out.println("Can not locate database node for " + path1);
                    System.exit(1);
                }
                if (db_node2 == null) {
                    System.out.println("Can not locate database node for " + path2);
                    System.exit(1);
                }
                K1 = (int) db_node1.getProperty("k_mer_size");
                K2 = (int) db_node2.getProperty("k_mer_size");
                sq_num1 = (int) db_node1.getProperty("num_nodes");
                sq_num2 = (int) db_node2.getProperty("num_nodes");
                ed_num1 = (int) db_node1.getProperty("num_edges");
                ed_num2 = (int) db_node2.getProperty("num_edges");
                ng1 = (int) db_node1.getProperty("num_genomes");
                ng2 = (int) db_node2.getProperty("num_genomes");
                knum1 = (long) db_node1.getProperty("num_k_mers");
                knum2 = (long) db_node2.getProperty("num_k_mers");
                if (K1 != K2 || sq_num1 != sq_num2 || sq_num1 != sq_num2 || ed_num1 != ed_num2 || ng1 != ng2 || knum1 != knum2) {
                    return false;
                }
                genomeDb1 = new SequenceDatabase(path1 + GENOME_DATABASE_PATH);
                genomeDb2 = new SequenceDatabase(path2 + GENOME_DATABASE_PATH);
                indexDb1 = new IndexDatabase(path1 + INDEX_DATABASE_PATH);
                indexDb2 = new IndexDatabase(path2 + INDEX_DATABASE_PATH);
                System.out.println("comparing pangenomes...");
                K = K1;
                Node n1, n2;
                byte[] s1 = new byte[21];
                byte[] s2 = new byte[21];
                kmer k_mer1 = new kmer(K, indexDb1.get_pre_len(), indexDb1.get_suf_len());
                kmer k_mer2 = new kmer(K, indexDb2.get_pre_len(), indexDb2.get_suf_len());
                IndexPointer p1 = new IndexPointer();
                IndexPointer p2 = new IndexPointer();
                for (k = 0; k < knum1; ++k) {
                    k_mer1 = indexDb1.get_kmer(k);
                    k_mer2 = indexDb2.get_kmer(k);
                    if (k_mer1.compare(k_mer2) != 0) {
                        return false;
                    }
                    indexDb1.get_pointer(s1, p1, k);
                    indexDb2.get_pointer(s2, p2, k);
                    if (p1.next_index == -1L && p2.next_index == -1L) {
                        n1 = g1.getNodeById(p1.node_id);
                        len1 = (int) n1.getProperty("length");
                        n2 = g2.getNodeById(p2.node_id);
                        len2 = (int) n2.getProperty("length");
                        if (len1 != len2 || (!genomeDb1.compare(genomeDb2, (int[]) n1.getProperty("address"), (int[]) n2.getProperty("address"), 0, 0, len1, true)
                                && !genomeDb1.compare(genomeDb2, (int[]) n1.getProperty("address"), (int[]) n2.getProperty("address"), 0, 0, len1, false))) {
                            return false;
                        }
                        if (n1.getDegree(Direction.BOTH) != n2.getDegree(Direction.BOTH)) {
                            return false;
                        }
                    }
                    if (k % (knum1 / 100 + 1) == 0) {
                        System.out.print(k * 100 / knum1 + 1 + "%\r");
                    }
                }
                System.out.println();
                tx1.success();
                tx2.success();
            }
            g1.shutdown();
            g2.shutdown();
            return true;
        }
        System.out.println("Error: a pair of databases are needed to be compared!");
        return false;
    }
    
    /**
     * To iteratively receive and run Cypher queries. Results will be written in output files.
     */
    public void run_query() {
        int i;
        String query;
        BufferedWriter out;
        
        if (new File(PATH + GRAPH_DATABASE_PATH).exists()) {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH + GRAPH_DATABASE_PATH))
                    .setConfig("keep_logical_logs", "100M size").newGraphDatabase();
            registerShutdownHook(graphDb);
            startTime = System.currentTimeMillis();
            for (i = 1;; ++i) {
                Console console = System.console();
                query = console.readLine("Enter a query or type exit : ");
                if (query.equals("exit")) {
                    break;
                }
                try (Transaction ignored = graphDb.beginTx(); Result result = graphDb.execute(query)) 
                {
                    try {
                        out = new BufferedWriter(new FileWriter(PATH + "/query" + i + ".result"));
                        while (result.hasNext()) {
                            Map<String, Object> row = result.next();
                            for (Entry<String, Object> column : row.entrySet()) {
                                out.write(column.getKey() + ": " + column.getValue() + "\n");
                                System.out.println(column.getKey() + ": " + column.getValue());
                            }
                        }
                        out.close();
                    } catch (IOException ioe) {
                        System.out.println("Failed to write the result of query!");
                        System.exit(1);
                    }
                    ignored.success();
                } catch (QueryExecutionException qee) {
                    System.out.println(qee.getMessage());
                }
            }
            graphDb.shutdown();
        }
    }
    
    /**
     * Extracts the genomic region belonging to the specified sequence starting at th specified node.
     * 
     * @param seq Will contains the sequence after function ends.
     * @param start_ptr A pangenome pointer which points to the node where the sequence starts.
     * @param address An array which its first two elements determine the genome and sequence number of the region, respectively. 
     * @param begin The start coordinate of the region.
     * @param end  The stop coordinate of the region.
     */
    public static void extract_sequence(StringBuilder seq, IndexPointer start_ptr, int[] address, int begin, int end) {
        Relationship rel;
        Node neighbor, node;
        int loc, node_len, neighbor_len, seq_len, position;
        String rel_name;
        
        seq_len = end - begin + 1;
        seq.setLength(0);
        loc = begin;
        position=start_ptr.position;
        node=graphDb.getNodeById(start_ptr.node_id);
        node_len = (int) node.getProperty("length");
    // Takes the part of the region lies in the first node of the path that region takes in the graph    
        if (start_ptr.canonical) {
            if (position + seq_len - 1 <= node_len - 1) {
                loc += append_fwd(seq, (String) node.getProperty("sequence"), position, position + seq_len - 1);
            } else {
                loc += append_fwd(seq, (String) node.getProperty("sequence"), position, node_len - 1);
            }
        } else {
            if (position - (seq_len - 1) >= 0) {
                loc += append_rev(seq, (String) node.getProperty("sequence"), position - (seq_len - 1), position);
            } else {
                loc += append_rev(seq, (String) node.getProperty("sequence"), 0, position);
            }
        }
    //  traverse the path of the region   
        while (seq.length() < seq_len ) {
            address[2] = loc - K + 1;
            rel = get_outgoing_edge(node, address);
            neighbor = rel.getEndNode();
            rel_name = rel.getType().name();
            neighbor_len = (int) neighbor.getProperty("length");
            if (rel_name.charAt(1) == 'F') // Enterring forward side
                if (seq.length() + neighbor_len - K + 1 > seq_len) // neighbor is the last node of the path
                    loc += append_fwd(seq, (String) neighbor.getProperty("sequence"), K - 1, seq_len - seq.length() + K - 2);
                else 
                    loc += append_fwd(seq, (String) neighbor.getProperty("sequence"), K - 1, neighbor_len - 1);
            else // Enterring reverse side
                if (seq.length() + neighbor_len - K + 1 > seq_len) // neighbor is the last node of the pat
                    loc += append_rev(seq, (String) neighbor.getProperty("sequence"), neighbor_len - K - (seq_len - seq.length()) + 1, neighbor_len - K);
                else 
                    loc += append_rev(seq, (String) neighbor.getProperty("sequence"), 0, neighbor_len - K);
            node = neighbor;
            node_len = (int) node.getProperty("length");
        } // while
    }
  
    /**
     * Give the next node of the path to be traversed through.
     * 
     * @param current_node The current node of the path we are located at. 
     * @param address An array which determine the genome, sequence and position of the desirable outgoing edge.
     * @return The outgoing edge.
     */
    public static Relationship get_outgoing_edge(Node current_node, int[] address) {
        int[] points_to;
        for (Relationship r_out : current_node.getRelationships(Direction.OUTGOING)) {
            if (r_out.hasProperty("points_to")) {
                points_to = (int[]) r_out.getProperty("points_to");
                if ( points_to[0]==address[0] && points_to[1]==address[1] && points_to[2]==address[2] )
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
        seq.append(s.substring(from, to + 1));
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
        seq.append(reverse_complement(s.substring(from, to + 1)));
        return to - from + 1;
    }
    
    /**
     * Returns a pangenome pointer pointing to the specified genomic position.
     * 
     * @param genome The genome number
     * @param sequence The sequence number 
     * @param position The position of interest
     * @return 
     */
    public static IndexPointer locate(int genome, int sequence, int position) {
        int loc, low, high, mid , node_len;
        boolean forward;
        Node node, neighbor, seq_node;
        Relationship rel;
        String origin, anchor_sides;
        long[] anchor_nodes;
        int[] anchor_positions, address;
        
        address = new int[]{genome,sequence,0};
        origin = genome + "_" + sequence;
        seq_node = graphDb.findNode(sequence_label, "number", origin);
        anchor_nodes = (long[]) seq_node.getProperty("anchor_nodes");
        anchor_positions = (int[]) seq_node.getProperty("anchor_positions");
        anchor_sides = (String) seq_node.getProperty("anchor_sides");
    // Find the immediate preceding anchor_node, searching in the sorted array of anchor positions.      
        for (low = 0, high = anchor_sides.length() - 1, mid = (low + high) / 2; low <= high; mid = (low + high) / 2) {
            if (position < anchor_positions[mid]) {
                high = mid - 1;
            } else if (position > anchor_positions[mid]) {
                low = mid + 1;
            } else {
                break;
            }
        }
        if (position < anchor_positions[mid]) {
            --mid;
        }
        try (Transaction tx = graphDb.beginTx()) {
            node = graphDb.getNodeById(anchor_nodes[mid]);
            loc = anchor_positions[mid];
            forward = anchor_sides.charAt(mid) == 'F';
            node_len = (int) node.getProperty("length");
        // Traverse the pangenome from the anchor node until reach to the target
            while (loc + node_len - K + 1 <= position) 
            {
                address[2] = loc + node_len - K + 1;
                rel = get_outgoing_edge(node, address);
                neighbor = rel.getEndNode();
                loc += node_len - K + 1;
                node = neighbor;
                node_len = (int) node.getProperty("length");
            }
            tx.success();
        }
        return new IndexPointer(node.getId(), forward, forward ? position - loc : node_len - 1 - (position - loc), -1L);
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
    private Relationship connect(Node src, Node des, RelationshipType edge_type, int[] address) {
        Relationship rel;
        //System.out.println("connect "+src.getId()+" "+edge_type.name()+" "+des.getId()+" "+address[2]);
        rel = src.createRelationshipTo(des, edge_type);
        rel.setProperty("points_to", address);
        ++num_edges;
        return rel;
    }
    

    /**
     * Splits a node at a specified position by creating a new node called split_node as a part separated from the node.
     * @param node The node which should be split.
     * @param pos The position of the split with respect to the start on the node.
     * @return The newly created split node.
     */
    private Node split(Node node, int pos) {
        //System.out.println("split "+node.getId()+" "+pos);
        int split_len, node_len;
        int i, s_id, gen, seq, loc;
        long inx, split_first_kmer, node_last_kmer;
        int[] address;
        Node neighbor, split_node;
        Relationship rel;
        
        address = (int[]) node.getProperty("address");
        gen = address[0];
        seq = address[1];
        loc = address[2];
        node_len = (int) node.getProperty("length");
        num_bases += (K - 1);
        ++num_nodes;
        split_node = graphDb.createNode(node_label);
        address[0] = gen;
        address[1] = seq;
        address[2] = loc + pos;
        split_node.setProperty("address", address);
        split_len = node_len - pos;
        split_node.setProperty("length", split_len);
    // Updates the edges comming from gene level to the node.    
        if (node.hasRelationship(Direction.INCOMING, RelTypes.visits)) {
            for (Relationship r : node.getRelationships(RelTypes.visits, Direction.INCOMING)) 
            {
                rel = r.getStartNode().createRelationshipTo(split_node, RelTypes.visits);
                if (r.hasProperty("position") && (int) r.getProperty("position") >= pos) {
                    rel.setProperty("position", (int) r.getProperty("position") - pos);
                    rel.setProperty("forward", r.getProperty("forward"));
                    r.delete();
                } 
            }
        }
    // Updating the Kmers chain in the index    
        node_last_kmer = indexDb.find(make_kmer(gen, seq, loc + pos - 1));
        split_first_kmer = indexDb.get_next_index(node_last_kmer);
        indexDb.put_next_index(-1L, node_last_kmer);
        split_node.setProperty("first_kmer", split_first_kmer);
        split_node.setProperty("last_kmer", node.getProperty("last_kmer"));
        s_id = (int) split_node.getId();
        for (i = 0, inx = split_first_kmer; inx != -1L; inx = indexDb.get_next_index(inx), ++i) 
        {
            indexDb.put_node_id(s_id, inx);
            indexDb.put_position(i, inx);
        }
    // Moving forward-outgoing and reverse-incoming edges from node to split node.    
        for (Relationship r : node.getRelationships(Direction.OUTGOING,RelTypes.FR,RelTypes.FF)) {
            neighbor = r.getEndNode();
            if (neighbor.equals(node)) 
                neighbor = r.isType(RelTypes.FF) ? node : split_node;
            connect(split_node,neighbor, r.getType(), (int[])r.getProperty("points_to"));
            --num_edges;
            r.delete();
        }
        for (Relationship r : node.getRelationships(Direction.INCOMING,RelTypes.RR,RelTypes.FR)) {
            neighbor = r.getStartNode();
            if (neighbor.equals(node)) 
                neighbor = r.isType(RelTypes.RR) ? node : split_node;
            connect(neighbor,split_node, r.getType(), (int[])r.getProperty("points_to"));
            --num_edges;
            r.delete();
        }
    //  Connecting node to split node
        for (Relationship r : node.getRelationships(Direction.INCOMING,RelTypes.RF,RelTypes.FF)) {
                address = (int[])r.getProperty("points_to");
                address[2] = address[2] + pos ;
                connect(node,split_node,RelTypes.FF,address);
        }
        for (Relationship r : split_node.getRelationships(Direction.INCOMING,RelTypes.FR,RelTypes.RR)) {
                address = (int[])r.getProperty("points_to");
                address[2] = address[2] + (split_len-K+1) ;
                connect(split_node,node,RelTypes.RR,address);
        }
        node.setProperty("last_kmer", node_last_kmer);
        node.setProperty("length", pos + K - 1);
        return split_node;
    }
    
    /**
     * Extends a new node till reach to a previously visited K-mer or a degenerate region.
     * 
     * @param node The node to be extended.
     */
    private void extend(Node node) {
        //System.out.println("extending node "+node.getId());
        int begin, len;
        long id, last_kmer;
        boolean broke, degenerate;
        
        len = (int) node.getProperty("length");
        id = node.getId();
        last_kmer = (long) node.getProperty("last_kmer");
        broke = false;
        degenerate = false;        
        while (position < seq_len - 1) { // Not reached to the end of the sequence
            //System.out.println("extend "+position);
            if (genomeDb.get_code(genome, sequence, position + 1) > 3) { // hit a degenerate region
                ++position;
                begin = position - K + 1;
                node.setProperty("length", len);
                node.setProperty("last_kmer", last_kmer);
                jump();
                int[] add = new int[]{genome,sequence,begin};
                create_degenerate(add);
                connect(curr_node, degenerate_node, RelTypes.FF, add); //.values()[curr_side*2]
                curr_node = degenerate_node;
                //curr_side = 0; // we have set it zero already in create()
                degenerate = true;
                break;
            }
            next_kmer(genomeDb);
            if (position % (seq_len / 100 + 1) == 0) 
                System.out.print((long) position * 100 / seq_len + 1 + "%\r");
            curr_index = indexDb.find(k_mer);
            indexDb.get_pointer(byte_pointer, pointer, curr_index);
            if (pointer.node_id == -1L) {
                indexDb.put_next_index(curr_index, last_kmer);
                ++len;
                //node.setProperty("length",(int)node.getProperty("length")+1);
                ++num_bases;
                pointer.node_id = id;
                pointer.canonical = fwd_kmer.canonical;
                pointer.position = len - K;//(int)node.getProperty("length")-K;
                pointer.next_index = -1L;
                indexDb.put_pointer(pointer, curr_index);
                last_kmer = curr_index;
                //node.setProperty("last_kmer",curr_index);
                if (position % (seq_len / 100 + 1) == 0) 
                    System.out.print((long) position * 100 / seq_len + "%\r");
            } else {
                broke = true;
                break;
            }
        }
        if (!degenerate) {
            node.setProperty("length", len);
            node.setProperty("last_kmer", last_kmer);
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
        num_bases += K;
        ++num_nodes;
        new_node = graphDb.createNode(node_label);
        //System.out.println("create "+new_node.getId());
        new_node.setProperty("address", address);
        new_node.setProperty("length", K);
        new_node.setProperty("last_kmer", curr_index); // Record number of the current Kmer in the Kmer table
        new_node.setProperty("first_kmer", curr_index);
    // Set the pointer to the Kmer in the pointer database    
        pointer.node_id = new_node.getId();
        pointer.canonical = fwd_kmer.canonical;
        pointer.position = 0;
        pointer.next_index = -1L;
        indexDb.put_pointer(pointer, curr_index);
        connect(curr_node, new_node, RelTypes.values()[curr_side*2], address);
        curr_node = new_node;
        curr_side = 0;
    }

    /**
     * Enters the node in which the current Kmer is found in and performs zero, one or two splits.   
     */
    private void follow_forward() {
        int l, pos, begin, g, s, loc, side, points_to;
        Node node, split_node1, split_node2,des, src;
        int[] address;
        boolean degenerated, loop;
        pos = pointer.position;
        node = graphDb.getNodeById(pointer.node_id);
        loop = false;
    // The first split might be done to seperate the part we need to enter in.
        if (pos > 0) {  
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
        points_to=position - K + 1;
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
                connect(node, degenerate_node, RelTypes.FF, add);//.values()[curr_side*2]
                break;
            }
            if (position % (seq_len / 100 + 1) == 0) {
                System.out.print((long) position * 100 / seq_len + 1 + "%\r");
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
            split_node2 = split(node, pos);
            if (loop)
                src = split_node2;
        }
    // connect the current node before doing splits to the split_node1    
        connect(src,des,RelTypes.values()[side],new int[]{genome,sequence,points_to});//graphDb.getNodeById(src),graphDb.getNodeById(des)
        if (degenerated) {
            curr_node = degenerate_node;
            curr_side = 0; // not really needed
        } else {
            curr_index = indexDb.find(k_mer);
            indexDb.get_pointer(byte_pointer, pointer, curr_index);
            curr_node = node;
        }
    }
    /*
     To enter a node and follow it in forward direction till reaching to a split position or end of the node
     */

    private void follow_reverse() {
        int pos, begin, g, s, loc, side, points_to;
        int[] address;
        Node node, split_node1, split_node2 ,des, src;
        boolean degenerated = false, loop = false, first_split = false;
        pos = pointer.position;
        node = graphDb.getNodeById(pointer.node_id);
        split_node2 = node; //if the second split does not happens remains unchanged
        if (pos < (int) node.getProperty("length") - K) {
            first_split = true;
            split_node1 = split(node, pos+1);
            if (loop = curr_node.equals(node))
                src = split_node1;
            else 
                src = curr_node;
        } else {
            split_node1 = node;
            if (loop = curr_node.equals(node))
                src = split_node1;
            else 
                src = curr_node;
        }
        des = node;
        side=curr_side*2+1;
        points_to=position - K + 1;
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
                connect(split_node2, degenerate_node, RelTypes.RF, add);
                degenerated = true;
                break;
            }
            if (position % (seq_len / 100 + 1) == 0) {
                System.out.print((long) position * 100 / seq_len + 1 + "%\r");
            }
        }
        if (position == seq_len) {
            finish = true;
        } else if (!degenerated) {
            position -= K;
            initial_kmers(genomeDb);
        }
        if (!degenerated && pos >= 0) {
            split_node2 = split(node, pos+1);
            des = split_node2;
            if (!first_split && loop)
                src = split_node2;
        }
        connect(src,des,RelTypes.values()[side],new int[]{genome,sequence,points_to});
        if (degenerated) {
            curr_node = degenerate_node;
            curr_side = 0;
        } else {
            curr_index = indexDb.find(k_mer);
            indexDb.get_pointer(byte_pointer, pointer, curr_index);
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
            fwd_kmer.next_up_kmer(fwd_code);
            rev_kmer.next_down_kmer(rev_code);
            for (j = 0; j < K - 1 && position < seq_len - 1; ++j) {
                ++position;
                fwd_code = genomeDb.get_code(genome, sequence, position);
                rev_code = 3 - fwd_code;
                if (fwd_code > 3) {
                    break;
                }
                fwd_kmer.next_up_kmer(fwd_code);
                rev_kmer.next_down_kmer(rev_code);
            }
            if (j == K - 1) {
                fwd_kmer.canonical = fwd_kmer.compare(rev_kmer) == -1;
                k_mer = fwd_kmer.canonical ? fwd_kmer : rev_kmer;
            } else if (position == seq_len - 1) {
                finish = true;
            }
        } while (fwd_code > 3 && position < seq_len - 1);
    }
    
    /**
     * creates a degenerate node starting at "begin" ending at position-1.
     * @param address The genomic position of the region
     */
    private void create_degenerate(int[] address) {
        ++num_degenerates;
        degenerate_node = graphDb.createNode(degenerate_label);
        //System.out.println("create_degenerate "+degenerate_node.getId()+" position "+position);
        degenerate_node.setProperty("address", address);
        degenerate_node.setProperty("length", position - address[2]);
        num_bases += (position - address[2]);
        if (!finish) {
            curr_index = indexDb.find(k_mer);
            indexDb.get_pointer(byte_pointer, pointer, curr_index);
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
        fwd_kmer.next_up_kmer(fwd_code);
        rev_kmer.next_down_kmer(rev_code);
        fwd_kmer.canonical = fwd_kmer.compare(rev_kmer) == -1;
        k_mer = fwd_kmer.canonical ? fwd_kmer : rev_kmer;
        //System.out.println(k_mer.toString());
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
                create_degenerate(add);
                connect(curr_node, degenerate_node, RelTypes.values()[curr_side*2], add);
                curr_node = degenerate_node;
                break;
            }
            next_kmer(s_db);
        }
        fwd_kmer.canonical = fwd_kmer.compare(rev_kmer) == -1;
        k_mer = fwd_kmer.canonical ? fwd_kmer : rev_kmer;
    }
    /*
     To give the K-mer in genomic position (genome,sequence,position) 
     fwd_kmer: Forward sequence of the K-mer
     rev_kmer: Reverse sequence of the K-mer
     kmer  : canonical k-mer
     canonical: a boolean determining if the K-mer is canonical  
     */

    private kmer make_kmer(int genome, int sequence, int position) {
        int j, fwd_code, rev_code;
        kmer fwd_kmer = new kmer(K, indexDb.get_pre_len(), indexDb.get_suf_len());
        kmer rev_kmer = new kmer(K, indexDb.get_pre_len(), indexDb.get_suf_len());
        for (j = 0; j < K; ++j) {
            fwd_code = genomeDb.get_code(genome, sequence, position + j);
            rev_code = 3 - fwd_code;
            fwd_kmer.next_up_kmer(fwd_code);
            rev_kmer.next_down_kmer(rev_code);
        }
        fwd_kmer.canonical = fwd_kmer.compare(rev_kmer) == -1;
        return fwd_kmer.canonical ? fwd_kmer : rev_kmer;
    }
    /*
     The main function constructs the pangenome of input sequences.
     */

    void construct_pangenome(int previous_num_genomes) {
        int i;
        Node genome_node, sequence_node;
        long[] sequence_ids;
        phaseTime = System.currentTimeMillis();
        fwd_kmer = new kmer(K, indexDb.get_pre_len(), indexDb.get_suf_len());
        rev_kmer = new kmer(K, indexDb.get_pre_len(), indexDb.get_suf_len());
        for (genome = previous_num_genomes + 1; genome <= genomeDb.num_genomes; ++genome) {
            try (Transaction tx = graphDb.beginTx()) {
                genome_node = graphDb.createNode(genome_label);
                genome_node.setProperty("number", genome);
                genome_node.setProperty("num_sequences", genomeDb.num_sequences[genome]);
                db_node.createRelationshipTo(genome_node, RelTypes.has);
                System.out.println("Processing genome " + genome + " :             ");
                tx.success();
            }
            sequence_ids = new long[genomeDb.num_sequences[genome] + 1];
            for (sequence = 1; sequence <= genomeDb.num_sequences[genome]; ++sequence) {
                System.out.println("sequence " + sequence + "/" + genomeDb.num_sequences[genome] + " of genome " + genome + "\tlength=" + genomeDb.sequence_length[genome][sequence]);
                try (Transaction tx = graphDb.beginTx()) {
                    sequence_node = curr_node = graphDb.createNode(sequence_label);
                    sequence_node.setProperty("number", genome + "_" + sequence);
                    sequence_node.setProperty("sequence_title", genomeDb.sequence_titles[genome][sequence]);
                    sequence_node.setProperty("sequence_length", genomeDb.sequence_length[genome][sequence]);
                    genome_node.createRelationshipTo(sequence_node, RelTypes.has);
                    sequence_ids[sequence] = curr_node.getId();
                    curr_side = 0;
                    position = -1;
                    seq_len = genomeDb.sequence_length[genome][sequence];
                    initial_kmers(genomeDb);
                    tx.success();
                }
                curr_index = indexDb.find(k_mer);
                indexDb.get_pointer(byte_pointer, pointer, curr_index);
                finish = false;
                while (!finish) {
                    try (Transaction tx = graphDb.beginTx()) {
                        for (i = 0; i < db_trsc_limit && !finish; ++i) {
                            if (pointer.node_id == -1L) // kmer is new
                            {
                                create();
                                extend(curr_node);
                            } else if (fwd_kmer.canonical ^ pointer.canonical)// if sides don't agree
                            {
                                follow_reverse();
                            } else {
                                follow_forward();
                            }
                        }
                        tx.success();
                    }
                }//while position<n
                try (Transaction tx = graphDb.beginTx()) {
                    connect(curr_node, sequence_node, RelTypes.values()[curr_side*2], new int[]{genome,sequence,(int)seq_len-K+1});// to point to the last k-mer of the sequence located in the other strand
                    tx.success();
                }
            }//sequences
            try (Transaction tx = graphDb.beginTx()) {
                genome_node.setProperty("sequence_ids", sequence_ids);
                tx.success();
            }
            System.out.println((System.currentTimeMillis() - phaseTime) / 1000 + " seconds elapsed.");
        }//genomes
        annotate_nodes();
    }
    /*
     To add list of "anchor_nodes", "anchor_sides" and "anchor_positions" to each sequence_node.
     These properties facilitate extracting genomic regions. 
     */

    void annotate_nodes() {
        int trsc, i, m, neighbor_length = 0, count;
        long seq_len;
        long[] anchor_nodes;
        int[] anchor_positions;
        Node node, neighbor, seq_node;
        Relationship rel;
        StringBuilder nds = new StringBuilder();
        StringBuilder pos = new StringBuilder();
        StringBuilder sds = new StringBuilder();
        String[] ids_list, posis_list;
        String rel_name;
        int[] addr = null, address = new int[3];
        for (address[0] = 1; address[0] <= genomeDb.num_genomes; ++address[0]) {
            for (address[1] = 1; address[1] <= genomeDb.num_sequences[address[0]]; ++address[1]) {
                System.out.println("Localizing sequence " + address[1] + "/" + genomeDb.num_sequences[address[0]] + " of genome " + address[0] + "...                        ");
                seq_len = genomeDb.sequence_length[address[0]][address[1]] - 1;
                try (Transaction tx = graphDb.beginTx()) {
                    node = seq_node = graphDb.findNode(sequence_label, "number", address[0] + "_" + address[1]);
                    tx.success();
                }
                count = 0;
                for (address[2] = 0; address[2] + K - 1 < seq_len ;) // K-1 bases of the last node not added
                {
                    try (Transaction tx = graphDb.beginTx()) {
                        for (trsc = 0; trsc < db_trsc_limit && address[2] + K - 1 < seq_len ; ++trsc) {
                            rel = get_outgoing_edge(node, address);
                            //System.out.println("node "+node.getId()+" next "+address[0]+"_"+address[1]+"_"+address[2]);
                            neighbor = rel.getEndNode();
                            rel_name = rel.getType().name();
                            neighbor_length = (int) neighbor.getProperty("length");
                            // To add "sequence" property to the nodes
                            if (!neighbor.hasProperty("sequence")) {
                                addr = (int[]) neighbor.getProperty("address");
                                neighbor.setProperty("sequence", genomeDb.get_sequence(addr[0], addr[1], addr[2], neighbor_length, true));
                            }
                            if (count % anchor_distance == 0) {
                                nds.append(neighbor.getId()).append(" ");
                                sds.append(rel_name.charAt(1));
                                pos.append(address[2]).append(" ");
                            }
                            count++;
                            address[2] = address[2] + neighbor_length - K + 1;
                            node = neighbor;
                            if (address[2] % (seq_len / 100 + 1) == 0) {
                                System.out.print((long) address[2] * 100 / seq_len + 1 + "%    \r");
                            }
                        }
                        tx.success();
                    }
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
                try (Transaction tx = graphDb.beginTx()) {
                    seq_node.setProperty("anchor_nodes", anchor_nodes);
                    seq_node.setProperty("anchor_positions", anchor_positions);
                    seq_node.setProperty("anchor_sides", sds.toString());
                    tx.success();
                }
                nds.setLength(0);
                pos.setLength(0);
                sds.setLength(0);
            }//sequences
        }//genomes
    }
    /*
     To add "sequence" property to the nodes
     */

    void drop_property(String side) {
        int i;
        byte[] sequence;
        ResourceIterator<Node> nodes;
        Node node;
        try (Transaction tx = graphDb.beginTx()) {
            nodes = graphDb.findNodes(node_label);
            tx.success();
        }
        while (nodes.hasNext()) {
            try (Transaction tx = graphDb.beginTx()) {
                for (i = 0; i < db_trsc_limit && nodes.hasNext(); ++i) {
                    node = nodes.next();
                    node.removeProperty(side);
                }
                tx.success();
            }
        }
        nodes.close();
        try (Transaction tx = graphDb.beginTx()) {
            nodes = graphDb.findNodes(degenerate_label);
            tx.success();
        }
        while (nodes.hasNext()) {
            try (Transaction tx = graphDb.beginTx()) {
                for (i = 0; i < db_trsc_limit && nodes.hasNext(); ++i) {
                    node = nodes.next();
                    node.removeProperty(side);
                }
                tx.success();
            }
        }
        nodes.close();
    }
    /*
     To provide a simple binary search
     */

    public static int b_search(int[] array, int key) {
        int low = 0, mid, high = array.length - 1;
        while (low <= high) {
            mid = (low + high) / 2;
            if (key < array[mid]) {
                high = mid - 1;
            } else if (key > array[mid]) {
                low = mid + 1;
            } else {
                return mid;
            }
        }
        return -1;
    }
    /*
     To return the reverse complement of "s".
     */

    private static String reverse_complement(String s) {
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

    private static StringBuilder reverse_complement(StringBuilder s) {
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
        return rv;
    }
    /*
     To register the action to be taken if the program halts unexpectedly
     */

    private static void registerShutdownHook(final GraphDatabaseService graphDb) {
        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                graphDb.shutdown();
            }
        });
    }
    /*
     To return size of a folder
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
    /*
     To write the "seq" in a FASTA file with lines justified to lines "length" long 
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

}
