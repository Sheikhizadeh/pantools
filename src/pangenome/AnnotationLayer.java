/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import alignment.ProteinAlignment;
import genome.SequenceDatabase;
import index.IndexPointer;
import protein.protein_builder;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import static java.lang.Integer.max;
import static java.lang.Integer.min;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;
import org.neo4j.io.fs.FileUtils;
import org.neo4j.graphdb.Direction;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import static org.neo4j.graphdb.factory.GraphDatabaseSettings.keep_logical_logs;
import static pangenome.SequenceLayer.extract_sequence;
import static pangenome.SequenceLayer.getFolderSize;

import static pantools.Pantools.GENOME_DATABASE_PATH;
import static pantools.Pantools.GRAPH_DATABASE_PATH;
import static pantools.Pantools.RelTypes;
import static pantools.Pantools.gene_label;
import static pantools.Pantools.genomeDb;
import static pantools.Pantools.graphDb;
import static pantools.Pantools.num_edges;
import static pantools.Pantools.num_nodes;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.startTime;
import static pangenome.SequenceLayer.locate;
import static pantools.Pantools.CDS_label;
import static pantools.Pantools.K;
import static pantools.Pantools.MAX_TRANSACTION_SIZE;
import static pantools.Pantools.annotation_label;
import static pantools.Pantools.broken_protein_label;
import static pantools.Pantools.coding_gene_label;
import static pantools.Pantools.cores;
import static pantools.Pantools.executeCommand_for;
import static pantools.Pantools.exon_label;
import static pantools.Pantools.feature_label;
import static pantools.Pantools.genome_label;
import static pantools.Pantools.homology_group_lable;
import static pantools.Pantools.intron_label;
import static pantools.Pantools.mRNA_label;
import static pantools.Pantools.rRNA_label;
import static pantools.Pantools.reverse_complement;
import static pantools.Pantools.tRNA_label;
import static pantools.Pantools.write_fasta;

/**
 * Implements all the functionalities related to the annotation layer of the pangenome
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class AnnotationLayer {
    private double FRACTION;
    private double CONTRAST;
    private double INFLATION;
    private int THRESHOLD;
    private int MAX_ALIGNMENT_LENGTH;
    private int MAX_INTERSECTIONS;
    private int MAX_KMER_FREQ;
    private int THREADS;
    private AtomicInteger num_intersections;
    private AtomicInteger num_similarities;
    private int num_hexamers;
    private double[][] phylogeny_distance;
    private int[][] count;
    private int num_genomes;
    private ConcurrentLinkedQueue[] kmers_proteins_list;
    private ConcurrentLinkedQueue<Node> proteins;
    private ConcurrentLinkedQueue<Node> kmerized_proteins;
    private BlockingQueue<intersection> intersections;
    private BlockingQueue<intersection> similarities;
    private String pangenome_path;
    
    public AnnotationLayer(){
        FRACTION = 0.06;
        CONTRAST = 11;
        INFLATION = 16;
        THRESHOLD = 90;
        MAX_ALIGNMENT_LENGTH  = 1000;
        MAX_INTERSECTIONS  = 10000000;
        THREADS = cores;
    }
    
    public class intersection{
        public Node protein1;
        public Node protein2;
        public double similarity;
        public intersection(Node p1, Node p2, double s){
            protein1 = p1;
            protein2 = p2;
            similarity = s;
        }
    }
        
    public class Kmerize_proteins implements Runnable {
        int total;
        public Kmerize_proteins(int t) {
            total = t;
        }

        @Override
        public void run() {
            int i = 0, mask = (1 << 24) - 1, c, chunk = total / 40;
            Node protein_node;
            int protein_length, kmer_index;
            String protein;
            int[] code = new int[256];
            char[] aminoacids = new char[]
            {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
            for (i = 0; i < 20; ++i)
                code[aminoacids[i]] = i;
            try (Transaction tx = graphDb.beginTx()) {
                protein_node = proteins.poll();
                for (c = 0; protein_node != null; ++c){
                    protein = (String)protein_node.getProperty("protein", "");
                    protein_length = protein.length();
                    if (protein_length > 6){
                        kmerized_proteins.add(protein_node);
                        kmer_index = 0;
                        for (i = 0; i < 5; ++i)
                            kmer_index = (kmer_index << 4) | code[protein.charAt(i)];
                        for (; i < protein_length; ++i){// for each kmer of the protein
                            kmers_proteins_list[kmer_index].add(protein_node.getId());
                            kmer_index = ((kmer_index << 4) & mask) | code[protein.charAt(i)];
                        }                            
                    }
                    if (c % chunk == 0)
                        System.out.print("|");
                    protein_node = proteins.poll();
                }
                tx.success();
            }
        }
    }    

    public class Find_intersections implements Runnable {
        int total;
        double frac = FRACTION;
        int max_intersection = MAX_INTERSECTIONS;
        public Find_intersections(int t) {
            total = t;
        }

        @Override
        public void run() {
            int i, mask = (1 << 24) - 1, chunk = total / 40;
            int p, counter, num_ids, kmer_index;
            long[] crossing_protein_ids= new long[max_intersection];
            int[] code = new int[256];
            char[] aminoacids = new char[]
            {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
            for (i = 0; i < 20; ++i)
                code[aminoacids[i]] = i;
            Node protein_node, crossing_protein_node;
            Iterator<Long> proteins_list;
            long protein_id;
            int protein_length, shorter_len;
            String protein, crossing_protein;
            long crossing_protein_id, p_id;
            try (Transaction tx = graphDb.beginTx()) {
                protein_node = kmerized_proteins.poll();
                for (p = 0; protein_node != null; ++p) {
                    protein = (String)protein_node.getProperty("protein");
                    protein_length = protein.length();
                    protein_id = protein_node.getId();
                    num_ids = 0;
                    kmer_index = 0;
                    for (i = 0; i < 5; ++i)
                        kmer_index = (kmer_index << 4) | code[protein.charAt(i)];
                    for (; i < protein_length - 1 && num_ids < max_intersection; ++i){// for each kmer of the protein
                        proteins_list = kmers_proteins_list[kmer_index].iterator();
                        while(proteins_list.hasNext() && num_ids < max_intersection){
                            crossing_protein_id = proteins_list.next();
                            if (crossing_protein_id > protein_id)
                                crossing_protein_ids[num_ids++] = crossing_protein_id;
                        }
                        kmer_index = ((kmer_index << 4) & mask) | code[protein.charAt(i)];
                    }
                    if (num_ids > 0){
                        --num_ids;
                        Arrays.sort(crossing_protein_ids, 0, num_ids);
                        for (counter = 0, crossing_protein_id = crossing_protein_ids[num_ids]; num_ids >= 0; --num_ids){
                            p_id = crossing_protein_ids[num_ids];
                            if (crossing_protein_id != p_id){
                                if(counter > 1){
                                    crossing_protein_node = graphDb.getNodeById(crossing_protein_id);
                                    crossing_protein = (String)crossing_protein_node.getProperty("protein");
                                    shorter_len = Math.min(protein_length, crossing_protein.length());
                                    if (counter >= frac * (shorter_len - 5)){
                                        intersections.add(new intersection(protein_node, crossing_protein_node, 0));
                                        num_intersections.getAndIncrement();
                                    }
                                }
                                crossing_protein_id = p_id;
                                counter = 1; 
                            } else
                                ++counter;
                        }
                        if(counter > 1){
                            crossing_protein_node = graphDb.getNodeById(crossing_protein_id);
                            crossing_protein = (String)crossing_protein_node.getProperty("protein");
                            shorter_len = Math.min(protein_length, crossing_protein.length());
                            if (counter >= frac * (shorter_len - 5)){
                                intersections.add(new intersection(protein_node, crossing_protein_node, 0));
                                num_intersections.getAndIncrement();
                            }
                        }
                    }
                    if (p % chunk == 0)
                        System.out.print("|");
                    protein_node = kmerized_proteins.poll();
                }// for protein
                for (i = 0; i < THREADS; ++i)
                    intersections.add(new intersection(null,null,0));// end of queue
                tx.success();
            }
        }
    }    

    public class Find_similarities implements Runnable {
        int threshold = THRESHOLD;
        int max_alg_len = MAX_ALIGNMENT_LENGTH;
        ProteinAlignment aligner;
        public Find_similarities() {
            aligner = new ProteinAlignment(-10,-1,max_alg_len);
        }

        @Override
        public void run() {
            int i;
            Node protein_node1, protein_node2;
            String protein1, protein2;
            intersection ints;
            double similarity;
            try{
            try (Transaction tx = graphDb.beginTx()) {
                ints = intersections.take();
                for (i = 0; ints.protein1 != null; ++i) {
                    protein_node1 = ints.protein1;
                    protein_node2 = ints.protein2;
                    protein1 = (String)protein_node1.getProperty("protein");
                    protein2 = (String)protein_node2.getProperty("protein");
                    similarity = (double)protein_similarity(aligner, protein1, protein2) / perfect_score(aligner, protein1, protein2) * 100;
                    if (similarity > threshold){
                        num_similarities.getAndIncrement();
                        ints.similarity = similarity;
                        similarities.add(ints);
                    }
                    ints = intersections.take();
                }// for intersection
                similarities.add(new intersection(null,null,0));// end of queue
                tx.success();
            }
            }catch(InterruptedException e){
                System.err.println(e.getMessage());
            }
        }
    }    
  
    public class Write_similarities implements Runnable {
        public Write_similarities() {
        }

        @Override
        public void run() {
            intersection ints;
            try{
                ints = similarities.take();
                while (ints.protein1 != null){
                    try(Transaction tx = graphDb.beginTx()){
                        for (int trs = 0; trs < 10 * MAX_TRANSACTION_SIZE && ints.protein1 != null; ++trs) {
                            //System.out.println(ints.similarity);
                            ints.protein1.createRelationshipTo(ints.protein2, RelTypes.is_similar_to).setProperty("similarity", ints.similarity);
                            ints = similarities.take();
                        }
                        tx.success();
                    }
                }
            } catch(InterruptedException e){
                System.err.println(e.getMessage());
            }
        }
    }      
    
    /**
     * Implements a comparator for integer arrays of size two
     */
    public class IntPairComparator implements Comparator<int[]> {
        @Override
        public int compare(int[] x, int[] y) {
            if (x[2] < y[2]) 
                return -1;
            else if (x[2] > y[2]) 
                return 1;
            else if (x[3] < y[3]) 
                return -1;
            else if (x[3] > y[3]) 
                return 1;
            else
                return 0;
        }
    }

    public void initialize_panproteome(String protein_paths_file, String pangenome_path){
        String file_path, line, protein_ID;
        StringBuilder protein = new StringBuilder();
        Node protein_node = null, panproteome;
        int trsc, num_proteins = 0, genome;
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
        try(BufferedReader protein_paths = new BufferedReader(new FileReader(protein_paths_file))) {
            for (genome = 1; protein_paths.ready(); ++genome){
                file_path = protein_paths.readLine().trim();
                if (file_path.equals("")) // if line is empty
                    continue;
                BufferedReader in = new BufferedReader(new FileReader(file_path));
                protein_ID = in.readLine().trim().substring(1);
                ++num_proteins;
                while (in.ready()) {
                    try (Transaction tx = graphDb.beginTx()) {
                        for (trsc = 0; in.ready() && trsc < MAX_TRANSACTION_SIZE; ++trsc){
                            line = in.readLine().trim();
                            if (line.equals("")) // if line is empty
                                continue;
                            else if (line.charAt(0) == '>'){
                                ++num_proteins;
                                protein_node = graphDb.createNode(mRNA_label);
                                protein_node.setProperty("protein_ID", protein_ID);
                                protein_node.setProperty("protein", protein.toString());
                                protein_node.setProperty("protein_length", protein.length());
                                protein_node.setProperty("genome",genome);
                                protein.setLength(0);
                                protein_ID = line.substring(1);
                            }
                            else
                                protein.append(line);
                            if (num_proteins % 11 == 1)
                                System.out.print("\r" + num_proteins + " proteins ");
                        }
                        tx.success();
                    }
                }
                in.close();
                try (Transaction tx = graphDb.beginTx()) {
                    protein_node = graphDb.createNode(mRNA_label);
                    protein_node.setProperty("protein_ID", protein_ID);
                    protein_node.setProperty("protein", protein.toString());
                    protein_node.setProperty("protein_length", protein.length());
                    protein_node.setProperty("genome",genome);
                    protein.setLength(0);
                    tx.success();
                }
            }
            System.out.println("\r" + num_proteins + " proteins ");
            try (Transaction tx1 = graphDb.beginTx()) {
                panproteome = graphDb.createNode(pangenome_label);
                panproteome.setProperty("date", new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(new Date()));
                panproteome.setProperty("num_genomes",genome - 1);
                panproteome.setProperty("num_proteins",num_proteins);
                tx1.success();
            }
            protein_paths.close();
        } catch (IOException ex){
            System.out.println("\n" + ex.getMessage());
        }        
    }
    
    /**
     * Adds nodes related to different genomic features to the pangenome.
     * 
     * @param gff_paths_file A text file listing the paths to the annotation files
     * @param pangenome_path Path to the database folder
     */
    public void add_annotaions(String gff_paths_file, String pangenome_path) {
        if (! new File(pangenome_path + GRAPH_DATABASE_PATH).exists()) {
            System.out.println("No database found in " + pangenome_path);
            System.exit(1);
        }
        Node db_node, annotation_node;
        String[] fields;
        String gff_file, line;
        int[] address = new int[4];
        BufferedWriter log_file;
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(pangenome_path + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        try (Transaction tx1 = graphDb.beginTx()) {
            db_node = graphDb.findNodes(pangenome_label).next();
            if (db_node == null) {
                System.out.println("Can not locate database node!");
                System.exit(1);
            }
            // Read the pangenome information    
            K = (int) db_node.getProperty("k_mer_size");
            num_nodes = (int) db_node.getProperty("num_nodes");
            num_edges = (int) db_node.getProperty("num_edges");
            tx1.success();
        }
        startTime = System.currentTimeMillis();
        genomeDb = new SequenceDatabase(pangenome_path + GENOME_DATABASE_PATH);
        try (BufferedReader gff_paths = new BufferedReader(new FileReader(gff_paths_file))) {
            log_file = new BufferedWriter(new FileWriter(pangenome_path + "/annotation.log"));
            if (! new File(pangenome_path + "/proteins").exists())
                Files.createDirectory(Paths.get(pangenome_path + "/proteins"));
            System.out.println("genome\tgenes\tmRNAs\ttRNAs\trRNAs");
            while (gff_paths.ready()) // for each gff file
            {
                line = gff_paths.readLine();
                if (line.trim().equals(""))
                    continue;
                fields = line.split("\\s");
                address[0] = Integer.parseInt(fields[0]);
                gff_file = fields[1];
                if (! new File(gff_file).exists()){
                    log_file.write("Genome "+address[0]+"'s GFF file not found.");
                    continue;
                }
                try (Transaction tx = graphDb.beginTx()) {
                    annotation_node = graphDb.createNode(annotation_label);
                    annotation_node.setProperty("genome", address[0]);
                    annotation_node.setProperty("date", new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(new Date()));
                    graphDb.findNodes(genome_label,"number",address[0]).next().createRelationshipTo(annotation_node, RelTypes.has);
                    tx.success();
                }
                parse_gff(address, annotation_node.getId(), log_file, gff_file, pangenome_path);
            } // for genomes
            gff_paths.close();
            log_file.close();
        } catch (IOException ioe) {
            System.out.println(ioe.getMessage());
            System.out.println("Could not open " + gff_paths_file);
        }
        graphDb.shutdown();
        genomeDb.close();
        // delete the database transaction files
        File directory = new File(pangenome_path + GRAPH_DATABASE_PATH);
        for (File f : directory.listFiles())
            if (f.getName().startsWith("neostore.transaction.db."))
                f.delete();
        System.out.println("Annotation finished. See "+pangenome_path + "/annotation.log");
        System.out.println("Annotated proteins available in "+pangenome_path + "/proteins");
    }

    private void parse_gff(int[] address, long annotation_node_id, BufferedWriter log_file, String gff_file, String pangenome_path){
        int i, trsc, num_genes, num_mRNAs, num_tRNAs, num_rRNAs, feature_len, offset;
        long seq_len;
        String sequence_id, current_sequence_id=null, attribute;
        Node node, gene_node, rna_node, feature_node, parent_node = null;
        Relationship rel;
        String[] fields,parent_ids;
        String strand, line;
        LinkedList<Node> gene_nodes = new LinkedList();
        LinkedList<Node> rna_nodes = new LinkedList();
        IndexPointer start_ptr, stop_ptr;
        num_genes = num_mRNAs = num_tRNAs = num_rRNAs = 0;
        try (BufferedReader in = new BufferedReader(new FileReader(gff_file))) {
            // for each record of gff file
            while (in.ready())
            {
                try (Transaction tx2 = graphDb.beginTx()) {
                    for (trsc = 0; trsc < MAX_TRANSACTION_SIZE && in.ready(); ++trsc) {
                        line = in.readLine();
                        if (line.equals("") || line.charAt(0) == '#') // if line is empty or a comment skip it
                            continue;
                        fields = line.split("\\t");
                        if (fields.length < 8)
                            continue;
                        address[2] = Integer.parseInt(fields[3]); // begin
                        address[3] = Integer.parseInt(fields[4]); // end
                        strand = fields[6];
                        feature_len = address[3] - address[2] + 1;
                        attribute = fields[fields.length - 1];
                        sequence_id = fields[0]; // sequence id
                        // if a feature belongs to a new sequence
                        // find the sequence number according to the sequence id in the gff file
                        if ( ! sequence_id.equals(current_sequence_id) ) {
                            address[1] = find_sequence(sequence_id, address[0]);
                            current_sequence_id = sequence_id;
                        }
                        if (address[1] == -1){// if sequence is not uniquely determined
                            log_file.write("Sequence ID = "+sequence_id+" missed in genome "+address[0]+"\n"); // usually organal genes
                            continue;
                        }
                        seq_len = genomeDb.sequence_length[address[0]][address[1]];
                        if(address[2] > seq_len) {
                            log_file.write("Position "+address[2] + " is out of range 1-"+seq_len+".\n");
                            continue;
                        }
                        if(address[3] > seq_len) {
                            log_file.write("Position "+address[3] + " is out of range 1-"+seq_len+".\n");
                            continue;
                        }
                        feature_node = graphDb.createNode(feature_label);
                        feature_node.setProperty("address", address);
                        feature_node.setProperty("strand", strand);
                        feature_node.setProperty("length", feature_len);
                        feature_node.setProperty("ID", get_property(attribute,"ID"));
                        feature_node.setProperty("attribute", attribute);
                        feature_node.setProperty("type", fields[2]);
                        feature_node.setProperty("name", get_property(attribute,"Name"));
                        feature_node.setProperty("annotation_node_id", annotation_node_id);
                        feature_node.setProperty("genome",address[0]);
                        parent_node = get_node_by_id(gene_nodes,get_property(attribute,"Parent"));
                        if (parent_node != null)
                            parent_node.createRelationshipTo(feature_node, RelTypes.is_parent_of);
                        start_ptr = locate(address, K);
                        offset = start_ptr.offset;
                        node = graphDb.getNodeById(start_ptr.node_id);
                        rel = feature_node.createRelationshipTo(node, RelTypes.starts);
                        rel.setProperty("offset", offset);
                        rel.setProperty("genomic_position", address[2]);
                        rel.setProperty("forward", start_ptr.canonical);
                        stop_ptr = locate(new int[]{address[0], address[1], address[3]}, K);
                        rel = feature_node.createRelationshipTo(graphDb.getNodeById(stop_ptr.node_id), RelTypes.stops);
                        rel.setProperty("offset", stop_ptr.offset);
                        rel.setProperty("genomic_position", address[2] + feature_len - 1);
                        rel.setProperty("forward", stop_ptr.canonical);
                        if (fields[2].endsWith("gene")) {
                            // create new gene node
                            feature_node.addLabel(gene_label);
                            gene_nodes.add(feature_node);
                            // adding gene_node id to the sequence node
                            //genes_list[address[0]][address[1]].add(gene_node.getId());
                            ++num_genes;
                        } else  switch(fields[2]) {
                                    case "mRNA":
                                        ++num_mRNAs;
                                        feature_node.addLabel(mRNA_label);
                                        if (parent_node != null)
                                            parent_node.createRelationshipTo(feature_node, RelTypes.codes_for);
                                        rna_nodes.addFirst(feature_node);
                                        break;
                                    case "tRNA":
                                        ++num_tRNAs;
                                        feature_node.addLabel(tRNA_label);
                                        break;
                                    case "rRNA":
                                        ++num_rRNAs;
                                        feature_node.addLabel(rRNA_label);
                                        break;
                                    case "CDS": 
                                        feature_node.addLabel(CDS_label);
                                        parent_ids = get_property(attribute,"Parent").split(",");
                                        // connect CDS to its parent RNAs    
                                        for (i=0;i<parent_ids.length;++i){
                                            rna_node = get_node_by_id(rna_nodes,parent_ids[i]);
                                            // for CDSs with a parent of type gene    
                                            if (rna_node == null){
                                                gene_node = get_node_by_id(gene_nodes,parent_ids[i]);
                                                if (gene_node != null){
                                                    rna_node = graphDb.createNode(feature_label);
                                                    rna_node.addLabel(mRNA_label);
                                                    rna_node.setProperty("address", gene_node.getProperty("address"));
                                                    rna_node.setProperty("strand", gene_node.getProperty("strand"));
                                                    rna_node.setProperty("length", gene_node.getProperty("length"));
                                                    rna_node.setProperty("ID", gene_node.getProperty("ID"));
                                                    rna_node.setProperty("attribute", gene_node.getProperty("attribute"));
                                                    rna_node.setProperty("type", gene_node.getProperty("type"));
                                                    rna_node.setProperty("name", gene_node.getProperty("name"));
                                                    rna_node.setProperty("annotation_node_id", annotation_node_id);
                                                    rna_node.setProperty("genome",address[0]);
                                                    ++num_mRNAs;
                                                    rna_nodes.addFirst(rna_node);                                           
                                                    gene_node.createRelationshipTo(rna_node, RelTypes.is_parent_of);
                                                    gene_node.createRelationshipTo(rna_node, RelTypes.codes_for);
                                                    feature_node.createRelationshipTo(rna_node, RelTypes.contributes_to);
                                                }
                                            } else
                                                feature_node.createRelationshipTo(rna_node, RelTypes.contributes_to);
                                        }
                                        break;
                                    case "exon":
                                        feature_node.addLabel(exon_label);
                                        break;
                                    case "intron":
                                        feature_node.addLabel(intron_label);
                                        break;
                        }
                        if (trsc % 500 == 1)
                            System.out.print("\r" + address[0] + "\t" + num_genes + "\t" + num_mRNAs + "\t" + num_tRNAs + "\t" + num_rRNAs + "\t");
                    }// for trsc
                    tx2.success();
                } // tx2
            } // while lines
            in.close();
            System.out.println("\r" + address[0] + "\t" + num_genes + "\t" + num_mRNAs + "\t" + num_tRNAs + "\t" + num_rRNAs + "\t");
            set_protein_sequences(gene_nodes, address[0], log_file, pangenome_path);
            log_file.write("Genome "+address[0] + " : " + num_genes + " genes\t" + num_mRNAs + " mRNAs\t" + num_tRNAs + " tRNAs\t" + num_rRNAs + " rRNAs\n");
            log_file.write("----------------------------------------------------\n");
        }catch (IOException ioe) {
            System.out.println(ioe.getMessage());
            System.out.println("Could not open " + gff_file + "!");
        }
    }
    
    /**
     * Gives the node of an annotated feature, given its ID
     * @param nodes A list of annotated features
     * @param id ID of the query feature
     * @return The annotation node (gene or RNA) of the query feature
     */
    private Node get_node_by_id(LinkedList<Node> nodes, String id){
        ListIterator<Node> itr = nodes.listIterator();
        Node node;
        while (itr.hasNext()){
            node = itr.next();
            if (((String)node.getProperty("ID")).equals(id)){
                return node;
            }
        }
        return null;
    }    

    /**
     * Translate the protein sequence of the mRNA nodes 
     * @param gene_nodes A list of all the annotated genes
     * @param gff_file Name of the GFF file
     * @param log_file 
     */
    private void set_protein_sequences(LinkedList<Node> gene_nodes, int genome, BufferedWriter log_file, String pangenome_path){
        IntPairComparator comp = new IntPairComparator();
        PriorityQueue<int[]> pq = new PriorityQueue(comp);
        protein_builder pb = new protein_builder();
        Node cds_node, mrna_node, gene_node;
        Relationship start_edge;
        IndexPointer start_ptr = new IndexPointer();
        int trsc, isoforms_num, gene_start_pos, protein_num = 0;
        long annotaion;
        int[] address, begin_end;
        StringBuilder rna_builder = new StringBuilder();
        StringBuilder gene_builder = new StringBuilder(), protein;
        try (BufferedWriter out = new BufferedWriter(new FileWriter(pangenome_path + "/proteins/proteins_" + genome + ".fasta"))) {
            System.out.print("Adding protein sequences...");
            while (!gene_nodes.isEmpty()){
                try (Transaction tx = graphDb.beginTx()) {
                    for (trsc = 0; !gene_nodes.isEmpty()&& trsc < 2 * MAX_TRANSACTION_SIZE; ++trsc){
                        gene_node = gene_nodes.remove();
                        address = (int[])gene_node.getProperty("address");
                        gene_start_pos = address[2];
                        // extract gene sequence as appears in the sequence and connects it to its nodes
                        gene_builder.setLength(0);
                        start_edge = gene_node.getSingleRelationship(RelTypes.starts, Direction.OUTGOING);
                        start_ptr.node_id = start_edge.getEndNode().getId();
                        start_ptr.offset = (int)start_edge.getProperty("offset");
                        start_ptr.canonical = (boolean)start_edge.getProperty("forward");
                        extract_sequence(gene_builder, start_ptr, address);
                        if (gene_builder.length() == 0)
                            continue;
                        isoforms_num =0;
                        for (Relationship r1: gene_node.getRelationships(RelTypes.codes_for, Direction.OUTGOING)) {
                            mrna_node = r1.getEndNode();
                            if (mrna_node.hasRelationship(RelTypes.contributes_to, Direction.INCOMING)){
                                ++isoforms_num;
                                gene_node.addLabel(coding_gene_label);
                                for (Relationship r2: mrna_node.getRelationships(RelTypes.contributes_to, Direction.INCOMING)) {
                                    cds_node = r2.getStartNode();
                                    pq.add((int[])cds_node.getProperty("address"));
                                }
                                for (rna_builder.setLength(0);!pq.isEmpty();) {
                                    begin_end = pq.remove();
                                // some RNAs are misannotated    
                                    try{
                                        rna_builder.append(gene_builder.substring(begin_end[2] - gene_start_pos, begin_end[3] - gene_start_pos + 1));
                                    } catch(StringIndexOutOfBoundsException ex){
                                        log_file.write("Gene ID = " + gene_node.getProperty("ID")+" miss-annotated!\n");
                                        break;
                                    }
                                }
                                if (gene_node.getProperty("strand").equals("-"))
                                    reverse_complement(rna_builder);
                                protein = pb.translate(rna_builder);
                                if (protein.length() > 0){
                                    ++protein_num;
                                    if (protein_num % 11 == 1)
                                        System.out.print("\rAdding protein sequences... " + protein_num);
                                    annotaion = (long)mrna_node.getProperty("annotation_node_id");
                                    out.write(">G" + genome + "A" + annotaion + "P" + protein_num + "\n");
                                    write_fasta(out, protein, 70);
                                    mrna_node.setProperty("protein", protein.toString());
                                    mrna_node.setProperty("protein_ID", "G" + genome + "A" + annotaion + "P" + protein_num);
                                    mrna_node.setProperty("protein_length", protein.length());
                                    if (protein.charAt(0) != 'M'  || protein.charAt(protein.length() - 1) != '*'){
                                        log_file.write("Protein ID = " + mrna_node.getProperty("ID")+" miss-annotated!\n");
                                        mrna_node.addLabel(broken_protein_label);
                                    }
                                }  
                            }
                        }
                        gene_node.setProperty("isoforms_num", isoforms_num);
                    }
                    tx.success();
                }
            }
            System.out.print("\r                                                        \r"); //clear the line
        }catch (IOException ioe) {
            System.out.println(ioe.getMessage());
        }        
    }     
    
    /**
     * Extracts the value of the given property from the feature attribute 
     * @param attribute Attribute of the feature 
     * @param property Property name
     * @return The value of the property
     */
    private String get_property(String attribute, String property){
        String[] fields = attribute.split(";");
        for (int i=0;i<fields.length;++i)
            if (fields[i].startsWith(property))
                return fields[i].split("=")[1];
        return "";
    }
    
    /**
     * Finds the number of a sequence in a genome given its name.
     * 
     * @param name Name of the sequence
     * @param genome Number of the genome
     * @return The number of the sequence in the genome
     */
    private int find_sequence(String name, int genome) {
        int i, sequence = -1, number_of_matches = 0;
        for (i = 1; i <= genomeDb.num_sequences[genome]; ++i)
            if (genomeDb.sequence_titles[genome][i].contains(name)){
                ++number_of_matches;
                sequence = i;
            }
        if (number_of_matches > 1) // is not unique
            sequence = -1;
        return sequence;
    }

    /**
     * Groups the similar genes into homology and homology groups
     * @param args The command line arguments:
     * args[1] Path to the database folder
     * args[2] THRESHOLD to be replaced by the default value, if given
     * args[3] LEN_FACTOR to be replaced by the default value, if given
     */
    public void group(String[] args) {
        pangenome_path = args[1];
        int i, n, d, p, proteins_num, kmer_table_size;
        double x;
        ResourceIterator<Node>  proteins_iterator;
        startTime = System.currentTimeMillis();
        for (i = 2; i < args.length; ++i){
            switch (args[i]){
                case "-f":
                    x = Double.parseDouble(args[i + 1]);
                    if (x >= 0.001 && x <= 0.1)
                        FRACTION = x;
                    ++i;
                    break;
                case "-t":
                    n = Integer.parseInt(args[i + 1]);
                    if (n > 0 && n < 100)
                        THRESHOLD = n;
                    ++i;
                    break;
                case "-i":
                    x = Double.parseDouble(args[i + 1]);
                    if (x > 1 && x < 29)
                        INFLATION = x;
                    ++i;
                    break;
                case "-c":
                    x = Double.parseDouble(args[i + 1]);
                    if (x > 0 && x < 12)
                        CONTRAST = x;
                    ++i;
                    break;
                case "-d":
                    d = Integer.parseInt(args[i + 1]);
                    d = d < 1 ? 1 : d;
                    d = d > 6 ? 6 : d;
                    FRACTION = new double[] {0, 0.055, 0.045, 0.035, 0.025, 0.015, 0.005}[d];
                    THRESHOLD = new int[]   {0, 90, 80, 70, 60, 50, 40 }[d];
                    INFLATION = new double[]{0, 16, 13, 10, 7,  4,  1.2}[d];
                    CONTRAST = new double[] {0, 11, 9,  7,  5,  3,  1  }[d];
                    ++i;
                    break;
                case "-p":
                    p = Integer.parseInt(args[i + 1]);
                    if (p >= 1 && p <= cores)
                        THREADS = p;
                    ++i;
                    break;
            }
        }
        System.out.println("Running on " + THREADS + " CPU cores ...");
        System.out.println("FRACTION = " + FRACTION);
        System.out.println("THRESHOLD = " + THRESHOLD);
        System.out.println("INFLATION = " + INFLATION);
        System.out.println("CONTRAST = " + CONTRAST);
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(pangenome_path + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        
        kmer_table_size = (int)Math.round(Math.pow(2, 24));
        kmers_proteins_list = new ConcurrentLinkedQueue[kmer_table_size];
        proteins = new ConcurrentLinkedQueue<>();
        kmerized_proteins = new ConcurrentLinkedQueue<>();
        intersections = new LinkedBlockingQueue<>();
        similarities = new LinkedBlockingQueue<>();
        try(Transaction tx = graphDb.beginTx()){
            proteins_iterator = graphDb.findNodes(mRNA_label);
            while (proteins_iterator.hasNext())
                proteins.add(proteins_iterator.next());
            tx.success();
        }

        proteins_num = proteins.size();
        MAX_KMER_FREQ = proteins.size() / 100;
        num_hexamers = 0;
        for (i = 0; i < kmer_table_size; ++i)
            kmers_proteins_list[i] = new ConcurrentLinkedQueue();
        System.out.println("\nKmerizing proteins :");
        System.out.print("0 .......................................................... 100\n  "); 
        try{
            ExecutorService es = Executors.newCachedThreadPool();
            for(i = 0; i < THREADS; i++)
                es.execute(new Kmerize_proteins(proteins_num));
            es.shutdown();
            es.awaitTermination(10, TimeUnit.DAYS);        
        } catch (InterruptedException e){
            
        }
        for (i = 0; i < kmer_table_size; ++i){
            if (!kmers_proteins_list[i].isEmpty()){
                ++num_hexamers;
                if (kmers_proteins_list[i].size() >= MAX_KMER_FREQ)
                    kmers_proteins_list[i].clear();
            }
            else 
                kmers_proteins_list[i] = null;
        }
        System.gc();
        
        /*try(Transaction tx = graphDb.beginTx()){
            Node node;
            ResourceIterator<Node> itr = graphDb.findNodes(mRNA_label);
            while(itr.hasNext()){
                node = itr.next();
                proteins.add(node);
                node.removeProperty("group_id");
                for (Relationship r: node.getRelationships(RelTypes.has_homolog, Direction.INCOMING))
                    r.delete();
            }
            tx.success();
        }*/
        proteins_num = kmerized_proteins.size();
        num_intersections = new AtomicInteger(0);
        num_similarities = new AtomicInteger(0);
        System.out.println("\n\nFinding intersections and similarities:");
        System.out.print("0 .......................................................... 100\n  "); 
        try{
            ExecutorService es = Executors.newCachedThreadPool();
            for(i = 0; i < Math.max(1, THREADS / 8); i++)
                es.execute(new Find_intersections(proteins_num));
            TimeUnit.SECONDS.sleep(1); // To avoid a first block
            for(; i < THREADS - 1 ; i++)
                es.execute(new Find_similarities());;
            es.execute(new Write_similarities());
            es.shutdown();
            es.awaitTermination(10, TimeUnit.DAYS);        
        } catch (InterruptedException e){
            
        }
        kmers_proteins_list = null;   
        System.gc();

        System.out.println("\n\nBuilding homology groups : ");
        build_homology_groups(pangenome_path);
        System.out.println("Hexamers = " + num_hexamers);
        System.out.println("Proteins = " + proteins_num);
        System.out.println("Intersections = " + num_intersections.intValue());
        System.out.println("Similarities = " + num_similarities.intValue());
        System.out.println("Database size = " + getFolderSize(new File(pangenome_path + GRAPH_DATABASE_PATH)) + " MB");        

        File directory = new File(pangenome_path + GRAPH_DATABASE_PATH);
        for (File f : directory.listFiles()) {
            if (f.getName().startsWith("neostore.transaction.db.")) {
                f.delete();
            }
        }
    }
        
    /**
     * Retrieves the first random relationship of a specific type which connects two nodes in any direction. 
     * @param node1 The first node
     * @param node2 The second node
     * @param rt The relationship type
     * @return 
     */

    private Relationship get_edge(Node node1, Node node2, RelationshipType rt){
        Relationship rel = null;
        for (Relationship r: node1.getRelationships(rt))
            if (r.getOtherNode(node1).equals(node2)){
                rel = r;
                break;
            }
        return rel;
    }

    private long perfect_score(ProteinAlignment aligner, String p1, String p2) {
        char match;
        int i, len1, len2;
        long score;
        len1 = p1.length();
        len2 = p2.length();
        if (len1 < len2){
            for (score = 0, i = 0; i < len1; ++i) {
                match = p1.charAt(i);
                score += aligner.match[match][match];
            }
        } else {
            for (score = 0, i = 0; i < len2; ++i) {
                match = p2.charAt(i);
                score += aligner.match[match][match];
            }  
        }
        return score;
    }
   
    /**
     * Creates homology nodes which connect the homologous coding genes
     */
    private void build_homology_groups(String pangenome_path){
        int i, num_groups = 0, num_grouped_proteins;
        Node protein_node, homology_group;
        int[] copy_number;
        BufferedWriter groups_file;
        ResourceIterator<Node>  proteins_iterator;
        LinkedList<Node> group = new LinkedList();
        LinkedList<Node> homology_group_nodes = new LinkedList();
        try (Transaction tx = graphDb.beginTx()) {
            proteins.clear();
            proteins_iterator = graphDb.findNodes(mRNA_label);
            while (proteins_iterator.hasNext())
                proteins.add(proteins_iterator.next());
            num_genomes = (int)graphDb.findNodes(pangenome_label).next().getProperty("num_genomes");
            copy_number = new int[num_genomes + 1];   
            phylogeny_distance = new double[num_genomes + 1][];
            count = new int[num_genomes + 1][];
            for (i = 1; i < phylogeny_distance.length; ++i){
                phylogeny_distance[i] = new double[num_genomes + 1];
                count[i] = new int[num_genomes + 1];
            } 
            tx.success();
        }
        try{
            groups_file = new BufferedWriter(new FileWriter(pangenome_path + "/pantools_homologs.txt"));
            while (!proteins.isEmpty()) {
                try(Transaction tx = graphDb.beginTx()){
                protein_node = proteins.remove();
                num_grouped_proteins = breadth_first_search(group, protein_node);
                if (num_grouped_proteins > 0){ // has not been grouped before
                    num_groups += break_group(homology_group_nodes, group, pangenome_path);
                    while (!homology_group_nodes.isEmpty()){
                        homology_group = homology_group_nodes.remove();
                        groups_file.write(Long.toString(homology_group.getId()) + ":");
                        set_and_write_homology_groups(homology_group, copy_number, groups_file);
                    }  
                    group.clear();
                    homology_group_nodes.clear();
                }
                System.out.print("\r" + num_groups + " groups");
                tx.success();
                }
            } // while 
            groups_file.close();
        }catch (IOException ex){
            System.out.print(ex.getMessage());
        }                
        System.out.println("\r" + num_groups + " groups\n");

    }
    
    private void calculate_phylogeny_distances(ListIterator<Node> proteins_itr){
        int i, j, g1, g2;
        Node p1, p2;
        double similarity;
        for (i = 1; i < phylogeny_distance.length; ++i)
            for (j = 1; j < phylogeny_distance.length; ++j){
                phylogeny_distance[i][j] = phylogeny_distance[i][j] = 0;
                count[i][j] = count[j][i] = 0;
            }
        while (proteins_itr.hasNext()){
            p1 = proteins_itr.next();
            g1 = (int)p1.getProperty("genome");
            for (Relationship r: p1.getRelationships(Direction.OUTGOING, RelTypes.is_similar_to)){
                p2 = r.getEndNode();
                g2 = (int)p2.getProperty("genome");
                similarity = (double)r.getProperty("similarity");
                phylogeny_distance[g1][g2] += similarity;
                phylogeny_distance[g2][g1] += similarity;
                ++count[g1][g2];
                ++count[g2][g1];
            }
        }
        for (i = 1; i < phylogeny_distance.length; ++i)
            for (j = 1; j < phylogeny_distance.length; ++j){
                if (count[i][j] > 0)
                    phylogeny_distance[i][j] = phylogeny_distance[i][j]/count[i][j];
                else
                    phylogeny_distance[i][j] = THRESHOLD;
            }
        for (i = 1; i < phylogeny_distance.length; ++i){
            for (j = 1; j < phylogeny_distance.length; ++j){
                if (i == j)
                    phylogeny_distance[i][j] = 0;
                else
                    phylogeny_distance[i][j] = 100 - phylogeny_distance[i][j];
            }
        } 
    }    

    /**
     * Puts all the genes homologous to a query gene in one homology group
     * @param similar_groups An empty list to be used for storing groups with a member similar the the query_gene
     * @param crossing_genes An empty priority queue to be filled with the genes crossing the query_gene 
     * @param query_gene The gene to which similar ones should be found
     * @return The value by which the number of groups have been increased. (Could be a negative value if some groups get merged) 
     */
    private int breadth_first_search(LinkedList<Node> group, Node start_protein){
        int num_members = 0;
        long start_id = start_protein.getId();
        Node crossing_protein;
        try (Transaction tx1 = graphDb.beginTx()) {
        // To avoid having one protein in different groups    
            if (!start_protein.hasRelationship(RelTypes.has_homolog, Direction.INCOMING) ) { 
                Queue<Node> homologs = new LinkedList();
                homologs.add(start_protein);
                // for all the candidates with some shared node with the protein    
                while (!homologs.isEmpty()) {
                    start_protein = homologs.remove();
                    group.add(start_protein);
                    ++num_members;
                    for (Relationship crossing_edge: start_protein.getRelationships(RelTypes.is_similar_to)) {
                        crossing_protein = crossing_edge.getOtherNode(start_protein);
                        if(!crossing_protein.hasProperty("group_id")){
                            crossing_protein.setProperty("group_id", start_id);
                            homologs.add(crossing_protein);
                        }
                    }
                }// while
            }
            tx1.success(); 
        }
        return num_members;
    }  
    
    /**
     * Given two proteins calculates the normalized similarity score between them which is less or equal to 1.
     * Proteins longer than MAX_LENGTH will be broken in smaller parts to be compared correspondingly.  
     * @param p1 The first protein
     * @param p2 The second protein
     * @return The normalized similarity score which is less or equal to 1
     */
    long protein_similarity(ProteinAlignment aligner, String p1, String p2){
        int m = p1.length(), n = p2.length(), max_len = max(m,n);
        int i, parts_num = 1, part_len1, part_len2;
        long score;
        if (max_len > MAX_ALIGNMENT_LENGTH){
            parts_num = (max_len / MAX_ALIGNMENT_LENGTH) + (max_len % MAX_ALIGNMENT_LENGTH == 0 ? 0 : 1);
            part_len1 = m / parts_num;
            part_len2 = n / parts_num;
            for (score =0, i = 0; i < parts_num; ++i)
                score += aligner.get_similarity(p1.substring(i * part_len1, min(m, (i + 1) * part_len1)),
                                                    p2.substring(i * part_len2, min(n, (i + 1) * part_len2)) );
            return score;
        } else
            return aligner.get_similarity(p1, p2);
    }
   
    /**
     * divides the pre-calculated homology groups into homology groups using MCL algorithm 
     */   
    int break_group(LinkedList<Node> homology_group_nodes, LinkedList<Node> group, String pangenome_path){
        int i, num_groups = 0, num_members;
        double infl;
        Node homology_group;
        String graph_path, clusters_path, line, command;
        String[] fields;
        BufferedReader clusters_file;
        File tmp_file;
        if (group.size() == 1){
            homology_group = graphDb.createNode(homology_group_lable); 
            homology_group.createRelationshipTo(group.remove(), RelTypes.has_homolog);
            homology_group.setProperty("num_members", 1);
            homology_group_nodes.add(homology_group);
            num_groups++;
        } else {       
            graph_path = pangenome_path + "/" + group.getFirst().getId() + ".graph";
            clusters_path = pangenome_path + "/" + group.getFirst().getId() + ".clusters";
            write_similaity_matrix(group, graph_path);
            command = "mcl " + graph_path + " --abc -I " + INFLATION + " -o " + clusters_path;
            if(!executeCommand_for(command, group.size())){
                System.out.println("Failed to run MCL on group represented by " + group.getFirst().getId());
                homology_group = graphDb.createNode(homology_group_lable); 
                num_members = 0;
                while (!group.isEmpty()){
                    homology_group.createRelationshipTo(group.remove(), RelTypes.has_homolog);
                    ++num_members;
                }
                homology_group.setProperty("num_members", num_members);
                homology_group_nodes.add(homology_group);
                num_groups++;
            } else {
                try{
                    clusters_file = new BufferedReader(new FileReader(clusters_path));
                    while (clusters_file.ready()){
                        line = clusters_file.readLine();
                        fields = line.split("\\s");
                        if (fields.length > 1){
                            homology_group = graphDb.createNode(homology_group_lable);
                            homology_group_nodes.add(homology_group);
                            ++num_groups;
                            for (i = 0; i < fields.length; ++i)
                                homology_group.createRelationshipTo(graphDb.getNodeById(Long.parseLong(fields[i])), RelTypes.has_homolog);
                            homology_group.setProperty("num_members", fields.length);
                        } else 
                            num_groups += group_singleton(homology_group_nodes, fields[0]);
                    }
                    clusters_file.close();
                }catch (IOException ex){
                    System.out.print(ex.getMessage());
                }         
                new File(clusters_path).delete();
            }
            new File(graph_path).delete();
        }
        return num_groups;
    }   
    
    /**
     * Given a homology group complete all the pairwise similarities between its members and 
     * writes the weighted edges of the graph in SIMILARITY_GRAPH_FILE_NAME to be used by MCL clustering algorithm.
     * @param homology_group_node The homology group
     */
    private void write_similaity_matrix(LinkedList<Node> group, String graph_path){
        ListIterator<Node> itr;
        Node protein1_node, protein2_node;
        int genome1, genome2;
        double similarity;
        calculate_phylogeny_distances(group.listIterator());
        try (PrintWriter graph = new PrintWriter(graph_path)){
            for (itr = group.listIterator(); itr.hasNext(); ){
                protein1_node = itr.next();
                genome1 = (int)protein1_node.getProperty("genome");
                try (Transaction tx = graphDb.beginTx()) {
                    for (Relationship homology_edge: protein1_node.getRelationships(RelTypes.is_similar_to, Direction.OUTGOING)){
                        protein2_node = homology_edge.getEndNode();
                        similarity = (double)homology_edge.getProperty("similarity");
                        similarity -= THRESHOLD;
                        genome2 = (int)protein2_node.getProperty("genome");
                        similarity += phylogeny_distance[genome1][genome2];
                        graph.write(protein1_node.getId()+" "+protein2_node.getId()+" "+ Math.pow(similarity, CONTRAST) + "\n");
                    }
                    tx.success();
                }
            }
            graph.close();
        } catch (IOException ex){
        }
    }
    
    private int group_singleton(LinkedList<Node> homology_group_nodes, String singleton_id){
        int num_groups = 0;
        double best_score, similarity;
        Node homology_group, single_protein, neighbor, best_hit;
        single_protein = graphDb.getNodeById(Long.parseLong(singleton_id));
                try(Transaction tx = graphDb.beginTx()){
        if (!single_protein.hasRelationship(RelTypes.has_homolog, Direction.INCOMING))
        {
            best_score = 0;
            best_hit = null;
            for (Relationship r: single_protein.getRelationships(RelTypes.is_similar_to)){
                neighbor = r.getOtherNode(single_protein);
                similarity = (double)r.getProperty("similarity");
                if (similarity > best_score){
                    best_score = similarity;
                    best_hit = neighbor;
                }
            }
            if (best_score > 0){
                if(best_hit.hasRelationship(RelTypes.has_homolog, Direction.INCOMING)){
                    homology_group = best_hit.getSingleRelationship(RelTypes.has_homolog, Direction.INCOMING).getStartNode();
                    homology_group.createRelationshipTo(single_protein, RelTypes.has_homolog);
                    homology_group.setProperty("num_members", (int)homology_group.getProperty("num_members") + 1);
                } else {
                    homology_group = graphDb.createNode(homology_group_lable);
                    homology_group_nodes.add(homology_group);
                    ++num_groups;
                    homology_group.createRelationshipTo(single_protein, RelTypes.has_homolog);
                    homology_group.createRelationshipTo(best_hit, RelTypes.has_homolog);
                    homology_group.setProperty("num_members", 2);
                }
            } else {
                    homology_group = graphDb.createNode(homology_group_lable);
                    homology_group_nodes.add(homology_group);
                    ++num_groups;
                    homology_group.createRelationshipTo(single_protein, RelTypes.has_homolog);
                    homology_group.setProperty("num_members", 1);
            }
        }
        tx.success();}
        return num_groups;
    }    
    
    /**
     * Calculates the copy number of the genes in different geneomes to be stored in the homology nodes.
     */
    void set_and_write_homology_groups(Node homology_group_node, int[] copy_number, BufferedWriter groups_file) throws IOException {
        int i;
        Node protein_node;
        for (Relationship rel: homology_group_node.getRelationships(Direction.OUTGOING, RelTypes.has_homolog)){
            protein_node = rel.getEndNode();
            ++copy_number[(int)protein_node.getProperty("genome")];
            groups_file.write(" " + ((String)protein_node.getProperty("protein_ID")).replace(' ', '_'));
        }
        groups_file.write("\n");
        homology_group_node.setProperty("copy_number_variation", copy_number);
        for (i=0; i < copy_number.length; ++i)
            copy_number[i] = 0;
    }    
    
    /**
     * Shuts down the graph database if the program halts unexpectedly.
     * 
     * @param graphDb The graph database object 
     */
    private void registerShutdownHook(final GraphDatabaseService graphDb) {
        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                graphDb.shutdown();
            }
        });
    }
}
