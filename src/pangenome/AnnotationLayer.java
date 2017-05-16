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
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.PriorityQueue;
import java.util.Queue;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import org.neo4j.graphdb.Direction;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import static org.neo4j.graphdb.factory.GraphDatabaseSettings.keep_logical_logs;
import org.neo4j.graphdb.schema.IndexDefinition;
import org.neo4j.graphdb.schema.Schema;
import static pangenome.SequenceLayer.append_fwd;
import static pangenome.SequenceLayer.append_rev;

import static pantools.Pantools.GENOME_DATABASE_PATH;
import static pantools.Pantools.GRAPH_DATABASE_PATH;
import static pantools.Pantools.RelTypes;
import static pantools.Pantools.gene_label;
import static pantools.Pantools.genomeDb;
import static pantools.Pantools.graphDb;
import static pantools.Pantools.RNA_label;
import static pantools.Pantools.num_edges;
import static pantools.Pantools.num_nodes;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.sequence_label;
import static pantools.Pantools.startTime;
import static pangenome.SequenceLayer.locate;
import static pangenome.SequenceLayer.get_outgoing_edge;
import static pantools.Pantools.CDS_label;
import static pantools.Pantools.MAX_TRANSACTION_SIZE;
import static pantools.Pantools.annotation_label;
import static pantools.Pantools.broken_protein_label;
import static pantools.Pantools.coding_gene_label;
import static pantools.Pantools.executeCommand;
import static pantools.Pantools.exon_label;
import static pantools.Pantools.feature_label;
import static pantools.Pantools.genome_label;
import static pantools.Pantools.orthology_group_lable;
import static pantools.Pantools.homology_group_lable;
import static pantools.Pantools.intron_label;
import static pantools.Pantools.mRNA_label;
import static pantools.Pantools.kmer_lable;
import static pantools.Pantools.print_peak_memory;
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
    private static int K;
    private int L;
    private double THRESHOLD = 75;
    private int MAX_LENGTH  = 1000;
    private ProteinAlignment pro_aligner;
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

    public void initialize_panproteome(String protein_paths_file, String pangenome_path, int k_size){
        String file_path, line, protein_ID = "";
        StringBuilder protein = new StringBuilder();
        ResourceIterator<Node> kmers;
        Node kmer, panproteome, protein_node = null;
        int i,trsc, num_proteins = 1, genome;
        L = k_size;
        long[] kmer_node_ids = new long[(int)Math.pow(20,L)];
        int[] code = new int[256];
        char[] aminoacids = new char[]
        {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
        for (i = 0; i < 20; ++i)
            code[aminoacids[i]] = i;
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
                                if (protein.length() >= L){
                                    ++num_proteins;
                                    protein_node = graphDb.createNode(mRNA_label);
                                    protein_node.setProperty("protein_ID", protein_ID);
                                    protein_node.setProperty("protein", protein.toString());
                                    protein_node.setProperty("protein_length", protein.length());
                                    protein_node.setProperty("genome",genome);
                                    kmerize_protein(protein_node, protein, kmer_node_ids, code);
                                    protein.setLength(0);
                                }
                                protein_ID = line.substring(1);
                            }
                            else
                                protein.append(line);
                            if (num_proteins % 11 == 1)
                                System.out.print("\rnum_proteins : " + num_proteins);
                        }
                        tx.success();
                    }
                }
                in.close();
                if (protein.length() >= L){
                    try (Transaction tx = graphDb.beginTx()) {
                        protein_node = graphDb.createNode(mRNA_label);
                        protein_node.setProperty("protein_ID", protein_ID);
                        protein_node.setProperty("protein", protein.toString());
                        protein_node.setProperty("protein_length", protein.length());
                        protein_node.setProperty("genome",genome - 1);
                        kmerize_protein(protein_node, protein, kmer_node_ids, code);
                        protein.setLength(0);
                        tx.success();
                    }
                }
            }
            try (Transaction tx = graphDb.beginTx()) {
                panproteome = graphDb.createNode(pangenome_label);
                panproteome.setProperty("kmer_size", L);
                panproteome.setProperty("num_genomes", genome - 1);
                panproteome.setProperty("date", new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(new Date()));
                System.out.println("\rnum_proteins : " + num_proteins);
                kmers = graphDb.findNodes(kmer_lable);
                while (kmers.hasNext()){
                    kmer = kmers.next(); 
                    kmer.setProperty("degree", kmer.getDegree());
                }
                tx.success();
            }
            protein_paths.close();
        } catch (IOException ex){
            System.out.println("\n" + ex.getMessage());
        }        
    }
    
    private void kmerize_protein(Node protein_node, StringBuilder protein, long[] kmer_node_ids, int[] code){
        Node kmer_node;
        int i,start, protein_length, kmer_index, caa;
        protein_length = protein.length();
        //char[] kmer = new char[L];
        int[] seed = new int[L];
        for (i = 0; i < L - 1; ++i)
            seed[i] = protein.charAt(i);
        for (start = 0; start < protein_length - L; ++start){// for each penta-mer of the protein
            seed[(start + L - 1) % L] = protein.charAt(start + L - 1);
            kmer_index = 0;
            for (i = 0; i < L; ++i){
                caa = seed[(start + i) % L];
                kmer_index = kmer_index * 20 + code[caa];
               //kmer[i] = (char)caa;
            }
            if (kmer_node_ids[kmer_index] == 0){
                kmer_node = graphDb.createNode(kmer_lable);
                //kmer_node.setProperty("sequence", kmer);
                kmer_node_ids[kmer_index] = kmer_node.getId();
            } else 
                kmer_node = graphDb.getNodeById(kmer_node_ids[kmer_index]);
            protein_node.createRelationshipTo(kmer_node, RelTypes.visits);//.setProperty("position",String.valueOf(start+1));
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
                fields = line.split(" ");
                address[0] = Integer.parseInt(fields[0]);
                gff_file = fields[1];
                if (! new File(gff_file).exists()){
                    log_file.write("Annotaion file for genome " + address[0]+" not found.");
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
                        // if a feature belongs to a new sequence
                        // find the sequence number according to the sequence id in the gff file
                        if ( ! sequence_id.equals(current_sequence_id) ) {
                            address[1] = find_sequence(sequence_id, address[0]);
                            current_sequence_id = sequence_id;
                        }
                        // if sequence is uniquely determined
                        if (address[1] > 0) {
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
                        } else // if sequence not found
                            log_file.write("Sequence ID = "+sequence_id+" missed in genome "+address[0]+"\n"); // usually organal genes
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
                        get_sequence_of_gene(gene_builder, address, gene_node);
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
     * Connect the annotated gene to all the nodes its sequence passes through in the pangenome.
     * At the same time assembles the sequence of the gene
     * 
     * @param gene_builder String builder which will contain the sequence of the gene  
     * @param address An array containing {genome, sequence, begin, end}
     * @param gene_node The gene node itself
     */
    public void get_sequence_of_gene(StringBuilder gene_builder, int[] address, Node gene_node){
        Relationship start_edge, rel;
        Node neighbor, node;
        int genomic_pos, node_len, neighbor_len, seq_len = address[3] - address[2] + 1, offset;
        String rel_name;
        boolean fwd;
        gene_builder.setLength(0);
        genomic_pos = address[2]-1;
        start_edge = gene_node.getSingleRelationship(RelTypes.starts, Direction.OUTGOING);
        offset = (int)start_edge.getProperty("offset");
        node = start_edge.getEndNode();
        node_len = (int) node.getProperty("length");
        fwd = (boolean)start_edge.getProperty("forward");
        if (fwd) {
            if (offset + seq_len - 1 <= node_len - 1) // node covers the gene
                genomic_pos += append_fwd(gene_builder, (String) node.getProperty("sequence"), offset, offset + seq_len - 1);
            else
                genomic_pos += append_fwd(gene_builder, (String) node.getProperty("sequence"), offset, node_len - 1);
        } else {
            if (offset - (seq_len - 1) >= 0) // node covers the gene
                genomic_pos += append_rev(gene_builder, (String) node.getProperty("sequence"), offset - (seq_len - 1), offset);
            else
                genomic_pos += append_rev(gene_builder, (String) node.getProperty("sequence"), 0, offset);
        }
        while (gene_builder.length() < seq_len ) {
            try (Transaction tx = graphDb.beginTx()) {
                for (int trsc = 0; gene_builder.length() < seq_len && trsc < MAX_TRANSACTION_SIZE; ++trsc){
                    address[2] = genomic_pos - K + 1;
                    rel = get_outgoing_edge(node, address);
                    neighbor = rel.getEndNode();
                    rel_name = rel.getType().name();
                    neighbor_len = (int) neighbor.getProperty("length");
                    if (rel_name.charAt(1) == 'F') {
                        if (gene_builder.length() + neighbor_len - K + 1 > seq_len) 
                            genomic_pos += append_fwd(gene_builder, (String) neighbor.getProperty("sequence"), K - 1, seq_len - gene_builder.length() + K - 2);
                        else 
                            genomic_pos += append_fwd(gene_builder, (String) neighbor.getProperty("sequence"), K - 1, neighbor_len - 1);
                    } else {
                        if (gene_builder.length() + neighbor_len - K + 1 > seq_len) 
                        genomic_pos += append_rev(gene_builder, (String) neighbor.getProperty("sequence"), neighbor_len - K - (seq_len - gene_builder.length()) + 1, neighbor_len - K);
                    else 
                        genomic_pos += append_rev(gene_builder, (String) neighbor.getProperty("sequence"), 0, neighbor_len - K);
                    }
                    node = neighbor;
                }
                tx.success();
            }
        } // while forward, forward ? genomic_pos - node_start_pos : node_len - 1 - (genomic_pos - node_start_pos)
    }

    /**
     * Groups the similar genes into homology and orthology groups
     * @param args The command line arguments:
     * args[1] Path to the database folder
     * args[2] THRESHOLD to be replaced by the default value, if given
     * args[3] LEN_FACTOR to be replaced by the default value, if given
     */
    public void group(String[] args) {
        String pangenome_path = args[1];
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(pangenome_path + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        startTime = System.currentTimeMillis();
        try (Transaction tx = graphDb.beginTx()) {
            L = (int)graphDb.findNodes(pangenome_label).next().getProperty("kmer_size");
            tx.success();
        }
        pro_aligner = new ProteinAlignment(-10,-1,MAX_LENGTH);
        if (args.length > 2)
            THRESHOLD = Double.parseDouble(args[2]);
        System.out.println("Threshold = " + THRESHOLD);
        find_crossing_proteins();
        build_homology_groups();
        build_orthology_groups(pangenome_path);
    }
    
    public void find_crossing_proteins(){
        int i, degree, p, count, trsc, num_ids;
        double similarity, max_kmer_freq;
        LinkedList<Node> proteins = new LinkedList();
        ResourceIterator<Node> proteins_iterator;
        long[] crossing_protein_ids= new long[10000000];
        Node protein_node, kmer_node, crossing_protein_node;
        Relationship r;
        long other_node_id, protein_node_id;
        double fraction = (10.0 - L) / 100.0;
        int protein_length, crossing_protein_length, shorter_len, total_proteins, longer_len;
        char match;
        String protein, crossing_protein, shorter_protein;
        long perfect_score, crossing_protein_id, p_id;
        System.out.println("Finding crossing proteins...");
        try (Transaction tx = graphDb.beginTx()) {
            proteins_iterator = graphDb.findNodes(mRNA_label);
            while (proteins_iterator.hasNext())
                proteins.add(proteins_iterator.next());
            tx.success();
            proteins_iterator.close();
        }
        total_proteins = proteins.size();
        max_kmer_freq = 0.01 * total_proteins;
        for (p = 0; !proteins.isEmpty(); ) {
            try (Transaction tx = graphDb.beginTx()) {
                for (trsc = 0; !proteins.isEmpty()&& trsc < 2 * MAX_TRANSACTION_SIZE; ++trsc){
                    protein_node = proteins.remove();
                    protein = (String)protein_node.getProperty("protein", null);
                    if (protein != null){
                        ++p;
                        protein_length = protein.length();
                        if (protein_length > L){
                            protein_node_id = protein_node.getId();
                            num_ids = 0;
                            for (Relationship rel1: protein_node.getRelationships(Direction.OUTGOING)){
                                kmer_node = rel1.getEndNode();
                                degree = (int)kmer_node.getProperty("degree");
                                if (degree < max_kmer_freq){ // halfs the run-time
                                    for (Relationship rel2: kmer_node.getRelationships(Direction.INCOMING)){
                                        other_node_id = rel2.getStartNode().getId();
                                        if (other_node_id != protein_node_id)
                                            crossing_protein_ids[num_ids++] = other_node_id;
                                    }
                                }
                            }
                            if (num_ids > 0){
                                --num_ids;
                                Arrays.sort(crossing_protein_ids, 0, num_ids);
                                for (count = 0, crossing_protein_id = crossing_protein_ids[0]; num_ids >= 0; --num_ids){
                                    p_id = crossing_protein_ids[num_ids];
                                    if (crossing_protein_id != p_id){
                                        if (count > fraction * protein_length && get_edge((crossing_protein_node = graphDb.getNodeById(crossing_protein_id)), protein_node, RelTypes.is_homolog_to) == null){
                                            crossing_protein = (String)crossing_protein_node.getProperty("protein");
                                            crossing_protein_length = crossing_protein.length();
                                            if (protein_length < crossing_protein_length){
                                                shorter_protein = protein;
                                                shorter_len = protein_length;
                                                longer_len = crossing_protein_length;
                                            } else {
                                                shorter_protein = crossing_protein;
                                                shorter_len = crossing_protein_length;
                                                longer_len = protein_length;
                                            }
                                            if (Math.abs(protein_length - crossing_protein_length) <= (100 - THRESHOLD) / 100 * longer_len){
                                                for (perfect_score = 0, i = 0; i < shorter_len; ++i){
                                                    match = shorter_protein.charAt(i);
                                                    perfect_score += pro_aligner.match[match][match];
                                                }
                                                similarity = (double)protein_similarity(protein, crossing_protein)/ perfect_score * 100;
                                                if (similarity > THRESHOLD){
                                                    r = protein_node.createRelationshipTo(crossing_protein_node, RelTypes.is_homolog_to); 
                                                    r.setProperty("similarity",similarity);
                                                }
                                            }
                                        }
                                        crossing_protein_id = p_id;
                                        count = 1;
                                    }
                                    else
                                        ++count;
                                }
                                if (count > fraction * protein_length && get_edge((crossing_protein_node = graphDb.getNodeById(crossing_protein_id)), protein_node, RelTypes.is_homolog_to) == null){
                                    crossing_protein = (String)crossing_protein_node.getProperty("protein");
                                    crossing_protein_length = crossing_protein.length();
                                    if (protein_length < crossing_protein_length){
                                        shorter_protein = protein;
                                        shorter_len = protein_length;
                                        longer_len = crossing_protein_length;
                                    } else {
                                        shorter_protein = crossing_protein;
                                        shorter_len = crossing_protein_length;
                                        longer_len = protein_length;
                                    }
                                    if (Math.abs(protein_length - crossing_protein_length) <= (100 - THRESHOLD) / 100 * longer_len){
                                        for (perfect_score = 0, i = 0; i < shorter_len; ++i){
                                            match = shorter_protein.charAt(i);
                                            perfect_score += pro_aligner.match[match][match];
                                        }
                                        similarity = (double)protein_similarity(protein, crossing_protein)/ perfect_score * 100;
                                        if (similarity > THRESHOLD){
                                            r = protein_node.createRelationshipTo(crossing_protein_node, RelTypes.is_homolog_to); 
                                            r.setProperty("similarity",similarity);
                                        }
                                    }
                                }
                            }
                        }  
                        if (p % 11 == 1)
                            System.out.print("\r" + p + "/" + total_proteins);
                    }
                }
                tx.success();
            }
        } // for protein
    }
   
    /**
     * Creates homology nodes which connect the homologous coding genes
     */
    private void build_homology_groups(){
        int num_groups = 0, num_grouped_proteins;
        Node protein_node;
        LinkedList<Node> proteins = new LinkedList();
        ResourceIterator<Node> proteins_iterator;
        try (Transaction tx = graphDb.beginTx()) {
            proteins_iterator = graphDb.findNodes(mRNA_label);
            while (proteins_iterator.hasNext())
                proteins.add(proteins_iterator.next());
            tx.success();
            proteins_iterator.close();
        } 
        while (!proteins.isEmpty()) {
            protein_node = proteins.remove();
            num_grouped_proteins = breadth_first_search(protein_node);
            if (num_grouped_proteins > 0)
                num_groups += 1;
            System.out.print("\rHomology groups : " + num_groups);
        } // while 
        System.out.println("\rHomology groups : " + num_groups);
    }
    
    /**
     * Puts all the genes homologous to a query gene in one homology group
     * @param similar_groups An empty list to be used for storing groups with a member similar the the query_gene
     * @param crossing_genes An empty priority queue to be filled with the genes crossing the query_gene 
     * @param query_gene The gene to which similar ones should be found
     * @return The value by which the number of groups have been increased. (Could be a negative value if some groups get merged) 
     */
    private int breadth_first_search(Node start_protein){
        int num_members = 0;
        double similarity;
        Node homology_group_node, crossing_protein;
        Relationship crossing_edge;
        Iterator<Relationship> itr;
        try (Transaction tx1 = graphDb.beginTx()) {
        // To avoid having one protein in different groups    
            if (start_protein.hasProperty("protein") && !start_protein.hasRelationship(RelTypes.has_homolog, Direction.INCOMING)) { 
                Queue<Node> homologs = new LinkedList();
                homology_group_node = graphDb.createNode(homology_group_lable);
                homology_group_node.createRelationshipTo(start_protein, RelTypes.has_homolog);
                ++num_members;
                homologs.add(start_protein);
                // for all the candidates with some shared node with the protein    
                while (!homologs.isEmpty()) {
                    start_protein = homologs.remove();
                    itr = start_protein.getRelationships(RelTypes.is_homolog_to).iterator();
                    while (itr.hasNext()) {
                        try (Transaction tx2 = graphDb.beginTx()) {
                            for (int trsc = 0; itr.hasNext() && trsc < MAX_TRANSACTION_SIZE; ++trsc) {
                                crossing_edge = itr.next();
                                crossing_protein = crossing_edge.getOtherNode(start_protein);
                                if(!crossing_protein.hasRelationship(RelTypes.has_homolog, Direction.INCOMING) && crossing_protein.hasProperty("protein")){
                                    similarity = (double)crossing_edge.getProperty("similarity");
                                    if ( similarity > THRESHOLD ) {
                                        crossing_edge.setProperty("similarity", similarity);
                                        homology_group_node.createRelationshipTo(crossing_protein, RelTypes.has_homolog);
                                        num_members++;
                                        homologs.add(crossing_protein);
                                    }
                                }
                            }
                            tx2.success(); 
                        }
                    }
                }// while
                homology_group_node.setProperty("num_members", homology_group_node.getDegree(Direction.OUTGOING));
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
    long protein_similarity(String p1, String p2){
        int m = p1.length(), n = p2.length(), max_len = max(m,n);
        int i, parts_num = 1, part_len1, part_len2;
        long score;
        if (max_len > MAX_LENGTH){
            parts_num = (max_len / MAX_LENGTH) + (max_len % MAX_LENGTH == 0 ? 0 : 1);
            part_len1 = m / parts_num;
            part_len2 = n / parts_num;
            for (score =0, i = 0; i < parts_num; ++i)
                score += pro_aligner.get_similarity(p1.substring(i * part_len1, min(m, (i + 1) * part_len1)),
                                                    p2.substring(i * part_len2, min(n, (i + 1) * part_len2)) );
            return score;
        } else
            return pro_aligner.get_similarity(p1, p2);
    }
   
    /**
     * divides the pre-calculated homology groups into orthology groups using MCL algorithm 
     */   
    void build_orthology_groups(String pangenome_path) {
        int num_groups = 0, trsc;
        int[] copy_number;
        ResourceIterator<Node> nodes;
        LinkedList<Node> homology_group_nodes = new LinkedList();
        LinkedList<Node> orthology_group_nodes = new LinkedList();
        LinkedList<Node> members = new LinkedList();
        Node homology_group_node, orthology_group_node, single_protein;
        String graph_path, clusters_path;
        FileReader clusters_file;
        BufferedWriter groups_file;
        try (Transaction tx = graphDb.beginTx()) {
            copy_number = new int[(int)graphDb.findNodes(pangenome_label).next().getProperty("num_genomes")+1];                
            nodes = graphDb.findNodes(homology_group_lable);
            while (nodes.hasNext())
                homology_group_nodes.add(nodes.next());
            tx.success();
            nodes.close();
        }
        try{
            groups_file = new BufferedWriter(new FileWriter(pangenome_path + "/pantools_orthologs.txt"));
            while (!homology_group_nodes.isEmpty()){
                try (Transaction tx = graphDb.beginTx()) {
                    for (trsc = 0; !homology_group_nodes.isEmpty() && trsc < MAX_TRANSACTION_SIZE; ++trsc) {
                        homology_group_node = homology_group_nodes.remove();
                        if( homology_group_node.getDegree() > 1){
                            graph_path = pangenome_path + "/" + homology_group_node.getId() + ".graph";
                            clusters_path = pangenome_path + "/" + homology_group_node.getId() + ".clusters";
                            for (Relationship rel: homology_group_node.getRelationships(Direction.OUTGOING))
                                members.add(rel.getEndNode());
                            write_similaity_matrix(members, graph_path);
                            executeCommand("mcl " + graph_path + " --abc -I 5 -o " + clusters_path);
                            new File(graph_path).delete();
                            clusters_file = new FileReader(clusters_path);
                            while (clusters_file.ready()){
                                orthology_group_node = graphDb.createNode(orthology_group_lable);
                                orthology_group_nodes.add(orthology_group_node);
                                ++num_groups;
                                make_group(clusters_file, orthology_group_node);
                                homology_group_node.createRelationshipTo(orthology_group_node, RelTypes.splits_into);
                            }
                            clusters_file.close();
                            new File(clusters_path).delete();
                        } else {
                            single_protein = homology_group_node.getSingleRelationship(RelTypes.has_homolog, Direction.OUTGOING).getEndNode();
                            orthology_group_node = graphDb.createNode(orthology_group_lable); 
                            orthology_group_nodes.add(orthology_group_node);
                            ++num_groups;
                            orthology_group_node.createRelationshipTo(single_protein, RelTypes.has_ortholog);
                        }
                        System.out.print("\rOrthology groups : " + num_groups);
                    } // for
                    tx.success();
                }
            }
            System.out.println("\rOrthology groups : " + num_groups);
        // compute copy_number array in each group and delete singleton groups merged into the other ones and now empty 
            while (!orthology_group_nodes.isEmpty()){
                try (Transaction tx = graphDb.beginTx()) {
                    for (trsc = 0; !orthology_group_nodes.isEmpty() && trsc < MAX_TRANSACTION_SIZE; ++trsc) {
                        orthology_group_node = orthology_group_nodes.remove();
                        groups_file.write(Long.toString(orthology_group_node.getId()) + ":");
                        set_and_write_orthology_groups(orthology_group_node, copy_number, groups_file);
                    }
                    tx.success();
                }
            }       
            groups_file.close();
        }catch (IOException ex){
            System.out.print(ex.getMessage());
        }
    }   
    
    private void make_group(FileReader clusters_file, Node orthology_group_node){
        char ch;
        long id = 0;
        Node protein_node;
        try{
            while ((ch = (char)clusters_file.read()) != '\n'){
                if (ch == '\t'){
                    protein_node = graphDb.getNodeById(id);
                    orthology_group_node.createRelationshipTo(protein_node, RelTypes.has_ortholog);
                    id = 0;
                }
                else
                    id = id *10 + Character.getNumericValue(ch);
            }
            protein_node = graphDb.getNodeById(id);
            orthology_group_node.createRelationshipTo(protein_node, RelTypes.has_ortholog);
        } catch (IOException ex){

        }
    }
    
    /**
     * Given a homology group complete all the pairwise similarities between its members and 
     * writes the weighted edges of the graph in SIMILARITY_GRAPH_FILE_NAME to be used by MCL clustering algorithm.
     * @param homology_group_node The homology group
     */
    private void write_similaity_matrix(LinkedList<Node> members, String graph_path){
        ListIterator<Node> itr1, itr2;
        Node protein1_node, protein2_node;
        Relationship homology_edge;
        String protein1, protein2, shorter_protein;
        int i, protein_len1, protein_len2, shorter_len, longer_len, genome1, genome2;
        double similarity;
        long perfect_score;
        try (PrintWriter graph = new PrintWriter(graph_path)){
            for (itr1 = members.listIterator(); itr1.hasNext(); ){
                protein1_node = itr1.next();
                for (itr2 = members.listIterator(itr1.nextIndex()); itr2.hasNext(); ){
                    try (Transaction tx = graphDb.beginTx()) {
                        for (int trsc = 0; itr2.hasNext() && trsc < MAX_TRANSACTION_SIZE; ++trsc) {
                            protein2_node = itr2.next();
                            homology_edge = get_edge(protein1_node, protein2_node, RelTypes.is_homolog_to);
                            if (homology_edge != null){
                                similarity = (double)homology_edge.getProperty("similarity");
                                genome1 = (int)protein1_node.getProperty("genome");
                                genome2 = (int)protein2_node.getProperty("genome");
                            // to penalize paralogs
                                graph.write(protein1_node.getId()+" "+protein2_node.getId()+" "+ Math.pow(similarity - THRESHOLD, genome1 != genome2 ? 1.2 : 1.0) +"\n");
                            } else { // for detecting non_crossing homologs which are rare
                                protein1 = (String)protein1_node.getProperty("protein");
                                protein2 = (String)protein2_node.getProperty("protein");
                                protein_len1 = protein1.length();
                                protein_len2 = protein2.length();
                                if (protein_len1 < protein_len2){
                                    shorter_protein =  protein1;
                                    shorter_len = protein_len1;
                                    longer_len = protein_len2;
                                } else {
                                    shorter_protein =  protein2;
                                    shorter_len = protein_len2;
                                    longer_len = protein_len1;
                                }
                                if (Math.abs(protein_len1 - protein_len2) <= (100 - THRESHOLD) / 100 * longer_len ){
                                    for (perfect_score = 0, i = 0; i < shorter_len; ++i)
                                        perfect_score += pro_aligner.match[shorter_protein.charAt(i)][shorter_protein.charAt(i)];
                                    similarity = (double)protein_similarity(protein1, protein2) / perfect_score * 100;
                                    if (similarity > THRESHOLD){
                                        genome1 = (int)protein1_node.getProperty("genome");
                                        genome2 = (int)protein2_node.getProperty("genome");
                                    // to penalize paralogs
                                        graph.write(protein1_node.getId()+" "+protein2_node.getId()+" "+ Math.pow(similarity - THRESHOLD, genome1 != genome2 ? 1.2 : 1.0) +"\n");
                                        protein1_node.createRelationshipTo(protein2_node, RelTypes.is_homolog_to).setProperty("similarity", similarity);
                                    }
                                }
                            }
                        } 
                        tx.success();
                    }
                }
            }
            graph.close();
            members.clear();
        } catch (IOException ex){
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
        for (Relationship rel: node1.getRelationships(rt,Direction.OUTGOING))
            if (rel.getEndNode().equals(node2))
                return rel;
        for (Relationship rel: node1.getRelationships(rt,Direction.INCOMING))
            if (rel.getStartNode().equals(node2))
                return rel;
        return null;
    }

    /**
     * Calculates the copy number of the genes in different geneomes to be stored in the orthology nodes.
     */
    void set_and_write_orthology_groups(Node orthology_group_node, int[] copy_number, BufferedWriter groups_file) throws IOException {
        int i, num_members;
        Node protein_node;
        num_members = 0;
        for (Relationship rel: orthology_group_node.getRelationships(Direction.OUTGOING)){
            ++num_members;
            protein_node = rel.getEndNode();
            ++copy_number[(int)protein_node.getProperty("genome")];
            groups_file.write(" " + protein_node.getProperty("protein_ID"));
        }
        groups_file.write("\n");
        orthology_group_node.setProperty("copy_number_variation", copy_number);
        orthology_group_node.setProperty("num_members", num_members);
        for (i=0; i<copy_number.length; ++i)
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
