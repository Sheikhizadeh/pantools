/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import alignment.ProteinAlignment;
import alignment.SequenceAlignment;
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
import static pantools.Pantools.K;
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
import static pantools.Pantools.broken_protein_label;
import static pantools.Pantools.coding_gene_label;
import static pantools.Pantools.executeCommand;
import static pantools.Pantools.orthology_group_lable;
import static pantools.Pantools.homology_group_lable;
import static pantools.Pantools.mRNA_label;
import static pantools.Pantools.pentamer_lable;
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
    public double THRESHOLD = 75;
    public int MAX_LENGTH  = 1000;
    ProteinAlignment pro_aligner;
    SequenceAlignment seq_aligner;
    
    /**
     * Implements a comparator for integer arrays of size two
     */
    public class PairComparator implements Comparator<int[]> {
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

    public void build_protein_graph(String protein_paths_file, String pangenome_path){
        String file_path, line, protein_number = "";
        StringBuilder protein = new StringBuilder();
        Node mrna, pangenome;
        int g, trsc, num_proteins = 0;
        int[] address = new int[3];
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(pangenome_path + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        startTime = System.currentTimeMillis();
        try(BufferedReader protein_paths = new BufferedReader(new FileReader(protein_paths_file))) {
            for (g = 1; protein_paths.ready(); ++g){
                file_path = protein_paths.readLine().trim();
                BufferedReader in = new BufferedReader(new FileReader(file_path));
                address[0] = g;
                while (in.ready()) {
                    try (Transaction tx = graphDb.beginTx()) {
                        for (trsc = 0; in.ready() && trsc < MAX_TRANSACTION_SIZE; ++trsc){
                            line = in.readLine().trim();
                            if (line.equals("")) // if line is empty
                                continue;
                            else if (line.charAt(0) == '>'){
                                if (protein.length() > 0){
                                    ++num_proteins;
                                    mrna = graphDb.createNode(mRNA_label);
                                    mrna.setProperty("protein_number", protein_number);
                                    mrna.setProperty("protein", protein.toString());
                                    mrna.setProperty("protein_length", protein.length());
                                    mrna.setProperty("address",address);
                                    protein.setLength(0);
                                }
                                protein_number = line.substring(1);
                            } else
                                protein.append(line);
                            if (num_proteins % 11 == 1)
                                System.out.print("\rnum_proteins : " + num_proteins);
                        }
                        tx.success();
                    }
                }
                ++num_proteins;
                try (Transaction tx = graphDb.beginTx()) {
                    mrna = graphDb.createNode(mRNA_label);
                    mrna.setProperty("protein_number", protein_number);
                    mrna.setProperty("protein", protein.toString());
                    mrna.setProperty("protein_length", protein.length());
                    mrna.setProperty("address",address);
                    protein.setLength(0);
                    tx.success();
                }
            }
            System.out.println("\rnum_proteins : " + num_proteins);
            try (Transaction tx = graphDb.beginTx()) {
                pangenome = graphDb.createNode(pangenome_label);
                pangenome.setProperty("num_genomes", g - 1);
                tx.success();
            }
        } catch (IOException ex){
            System.out.print(ex.getMessage());
        }        
    }
    
    /**
     * Adds nodes related to different genomic features to the pangenome.
     * 
     * @param gff_paths_file A text file listing the paths to the annotation files
     * @param pangenome_path Path to the database folder
     */

    public void annotate(String gff_paths_file, String pangenome_path) {
        if (! new File(pangenome_path + GRAPH_DATABASE_PATH).exists()) {
            System.out.println("No database found in " + pangenome_path);
            System.exit(1);
        }
        int i, trsc, num_genes, num_mRNAs, num_tRNAs, num_rRNAs, feature_len;
        String sequence_id, current_sequence_id=null, attribute;
        Node db_node, gene_node = null, rna_node, cds_node;
        String[] fields,parent_ids;
        String strand, gff_file, line, ID;
        //long[] genes_array;
        int[] address = new int[4];
        LinkedList<Node> gene_nodes = new LinkedList();
        LinkedList<Node> rna_nodes = new LinkedList();
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
        /*genes_list = new ArrayList[genomeDb.num_genomes + 1][];
        for (i = 1; i <= genomeDb.num_genomes; ++i) {
            genes_list[i] = new ArrayList[genomeDb.num_sequences[i] + 1];
            for (j = 1; j <= genomeDb.num_sequences[i]; ++j) {
                genes_list[i][j] = new ArrayList();
            }
        }*/
        try (BufferedReader gff_paths = new BufferedReader(new FileReader(gff_paths_file))) {
            BufferedWriter log_file = new BufferedWriter(new FileWriter(pangenome_path + "/annotation.log"));
            Files.createDirectory(Paths.get(pangenome_path + "/proteins"));
            System.out.println("genome\tgenes\tmRNAs\ttRNAs\trRNAs");
            for (address[0] = 1; gff_paths.ready() && address[0] <= genomeDb.num_genomes; ++address[0]) // for each gff file
            {
                rna_nodes.clear();
                num_genes = num_mRNAs = num_tRNAs = num_rRNAs = 0;
                gff_file = gff_paths.readLine();
                if (gff_file.equals("")){
                    log_file.write("Genome " + address[0]+" Skipped.");
                    continue;
                }
                BufferedReader in = new BufferedReader(new FileReader(gff_file));
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
                        // if sequence uniquely determined    
                            if (address[1] > 0) 
                            {
                                if (fields[2].endsWith("gene")) {
                                // create new gene node    
                                    gene_node = graphDb.createNode(gene_label);
                                    gene_node.setProperty("address", address);
                                    gene_node.setProperty("strand", strand);
                                    gene_node.setProperty("length", feature_len);
                                    gene_node.setProperty("ID", get_property(attribute,"ID"));
                                    gene_node.setProperty("fields", fields);
                                    gene_node.setProperty("type", fields[2]);
                                    gene_node.setProperty("name", get_property(attribute,"Name"));
                                    
                                    gene_nodes.add(gene_node);
                                // adding gene_node id to the sequence node
                                    //genes_list[address[0]][address[1]].add(gene_node.getId());
                                    ++num_genes;
                                } else if(fields[2].endsWith("RNA") || fields[2].equals("pseudogenic_transcript")) {
                                    gene_node = get_node_by_id(gene_nodes,get_property(attribute,"Parent"));
                                    if (gene_node != null){
                                        rna_node = graphDb.createNode(RNA_label);
                                        rna_node.setProperty("ID", get_property(attribute,"ID"));
                                        rna_node.setProperty("address", address);
                                        rna_node.setProperty("length", feature_len);
                                        rna_nodes.addFirst(rna_node);
                                        gene_node.createRelationshipTo(rna_node, RelTypes.codes_for);
                                        switch (fields[2]){
                                            case "mRNA":
                                                ++num_mRNAs;
                                                rna_node.addLabel(mRNA_label);
                                                break;
                                            case "tRNA":
                                                ++num_tRNAs;
                                                rna_node.addLabel(tRNA_label);
                                                break;
                                            case "rRNA":
                                                ++num_rRNAs;
                                                rna_node.addLabel(rRNA_label);
                                                break;
                                        }
                                    }
                                } else if (fields[2].equals("CDS")) {
                                        cds_node = graphDb.createNode(CDS_label);
                                        ID = get_property(attribute,"ID");
                                        cds_node.setProperty("ID", ID);
                                        cds_node.setProperty("address", address);
                                        parent_ids = get_property(attribute,"Parent").split(",");
                                    // connect CDS to its parent RNA    
                                        for (i=0;i<parent_ids.length;++i){
                                            rna_node = get_node_by_id(rna_nodes,parent_ids[i]);
                                        // for CDSs with a parent of type gene    
                                            if (rna_node == null){
                                                gene_node = get_node_by_id(gene_nodes,parent_ids[i]);
                                                if (gene_node != null){
                                                    rna_node = graphDb.createNode(RNA_label);
                                                    rna_node.addLabel(mRNA_label);
                                                    rna_node.setProperty("ID", gene_node.getProperty("ID"));
                                                    rna_node.setProperty("address", gene_node.getProperty("address"));
                                                    rna_node.setProperty("length", gene_node.getProperty("length"));
                                                    ++num_mRNAs;
                                                    rna_nodes.addFirst(rna_node);                                           
                                                    gene_node.createRelationshipTo(rna_node, RelTypes.codes_for);
                                                }
                                            }
                                            if (rna_node != null)
                                                cds_node.createRelationshipTo(rna_node, RelTypes.contributes_to);
                                        }
                                }
                            } else // if sequence not found
                                log_file.write("Sequence ID = "+sequence_id+" missed in genome "+address[0]+"\n"); // usually organal genes
                            if (trsc % 500 == 1)
                                System.out.print("\r" + address[0] + "\t" + num_genes + "\t" + num_mRNAs + "\t" + num_tRNAs + "\t" + num_rRNAs + "\t");
                        }// for trsc
                        tx2.success();
                    } // tx2
                } // while lines
            // for the last gene in the file. Translate the proteins of the gene.
                in.close();
                System.out.println("\r" + address[0] + "\t" + num_genes + "\t" + num_mRNAs + "\t" + num_tRNAs + "\t" + num_rRNAs + "\t");
                log_file.write("Genome "+address[0] + " : " + num_genes + " genes\t" + num_mRNAs + " mRNAs\t" + num_tRNAs + " tRNAs\t" + num_rRNAs + " rRNAs\n");
                log_file.write("----------------------------------------------------\n");
                set_protein_sequences(gene_nodes, address[0], log_file, pangenome_path);
            } // for genomes
            gff_paths.close();
            log_file.close();
        } catch (IOException ioe) {
            System.out.println(ioe.getMessage());
            System.exit(1);
        }
    // setting the genes property of sequence nodes    
        /*try (Transaction tx4 = graphDb.beginTx()) {
            for (i = 1; i <= genomeDb.num_genomes; ++i) {
                for (j = 1; j <= genomeDb.num_sequences[i]; ++j) {
                    genes_array = new long[genes_list[i][j].size()];
                    Iterator<Long> genesIterator = genes_list[i][j].iterator();
                    for (int k = 0; genesIterator.hasNext(); ++k) {
                        genes_array[k] = genesIterator.next();
                    }
                    graphDb.findNode(sequence_label, "number", i + "_" + j).setProperty("genes", genes_array);
                }
            }
            db_node.setProperty("num_genes",total_genes);
            tx4.success();
        }*/
        graphDb.shutdown();
        genomeDb.close();
        // delete the database transaction files
        File directory = new File(pangenome_path + GRAPH_DATABASE_PATH);
        for (File f : directory.listFiles())
            if (f.getName().startsWith("neostore.transaction.db."))
                f.delete();
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
        PairComparator comp = new PairComparator();
        PriorityQueue<int[]> pq = new PriorityQueue(comp);
        protein_builder pb = new protein_builder();
        Node cds_node, rna_node, gene_node;
        int trsc, isoforms_num, gene_start_pos, protein_num = 0;
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
                            ++isoforms_num;
                            rna_node = r1.getEndNode();
                            if (rna_node.hasLabel(mRNA_label) && rna_node.hasRelationship(RelTypes.contributes_to, Direction.INCOMING)){
                                gene_node.addLabel(coding_gene_label);
                                for (Relationship r2: rna_node.getRelationships(RelTypes.contributes_to, Direction.INCOMING)) {
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
                                if (protein.length() > 0 && protein.charAt(0) == 'M'  && protein.charAt(protein.length() - 1) == '*'){
                                    ++protein_num;
                                    if (protein_num % 11 == 1)
                                        System.out.print("\rAdding protein sequences... " + protein_num);
                                    out.write(">" + genome + "_" + protein_num + ":" + gene_node.getProperty("name") + "\n");
                                    write_fasta(out, protein, 70);
                                    rna_node.setProperty("protein", protein.toString());
                                    rna_node.setProperty("protein_number", genome + "_" + protein_num);
                                    rna_node.setProperty("protein_length", protein.length());
                                    
                                } else {
                                    log_file.write("Protein ID = " + rna_node.getProperty("ID")+" miss-annotated!\n");
                                    rna_node.addLabel(broken_protein_label);
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
        Relationship rel;
        Node neighbor, node;
        int loc, node_len, neighbor_len, seq_len = address[3] - address[2] + 1, starts_at;
        String rel_name;
        gene_builder.setLength(0);
        loc = address[2]-1;
        IndexPointer start_ptr = locate(address);
        starts_at = start_ptr.position;
        node = graphDb.getNodeById(start_ptr.node_id);
        node_len = (int) node.getProperty("length");
        rel = gene_node.createRelationshipTo(node, RelTypes.starts);
        rel.setProperty("starts_at", starts_at);
        rel.setProperty("position", loc);
        rel.setProperty("forward", start_ptr.canonical);
        if (start_ptr.canonical) {
            if (starts_at + seq_len - 1 <= node_len - 1)
                loc += append_fwd(gene_builder, (String) node.getProperty("sequence"), starts_at, starts_at + seq_len - 1);
            else
                loc += append_fwd(gene_builder, (String) node.getProperty("sequence"), starts_at, node_len - 1);
        } else {
            if (starts_at - (seq_len - 1) >= 0)
                loc += append_rev(gene_builder, (String) node.getProperty("sequence"), starts_at - (seq_len - 1), starts_at);
            else
                loc += append_rev(gene_builder, (String) node.getProperty("sequence"), 0, starts_at);
        }
        while (gene_builder.length() < seq_len ) {
            try (Transaction tx = graphDb.beginTx()) {
                for (int trsc = 0; gene_builder.length() < seq_len && trsc < MAX_TRANSACTION_SIZE; ++trsc){
                    address[2] = loc - K + 1;
                    rel = get_outgoing_edge(node, address);
                    neighbor = rel.getEndNode();
                    rel_name = rel.getType().name();
                    neighbor_len = (int) neighbor.getProperty("length");
                    if (rel_name.charAt(1) == 'F') {
                        if (gene_builder.length() + neighbor_len - K + 1 > seq_len) 
                            loc += append_fwd(gene_builder, (String) neighbor.getProperty("sequence"), K - 1, seq_len - gene_builder.length() + K - 2);
                        else 
                            loc += append_fwd(gene_builder, (String) neighbor.getProperty("sequence"), K - 1, neighbor_len - 1);
                    } else {
                        if (gene_builder.length() + neighbor_len - K + 1 > seq_len) 
                        loc += append_rev(gene_builder, (String) neighbor.getProperty("sequence"), neighbor_len - K - (seq_len - gene_builder.length()) + 1, neighbor_len - K);
                    else 
                        loc += append_rev(gene_builder, (String) neighbor.getProperty("sequence"), 0, neighbor_len - K);
                    }
                    node = neighbor;
                }
                tx.success();
            }
        } // while
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
        seq_aligner = new SequenceAlignment(-5,-0.1,4.0,-2.0,MAX_LENGTH);
        pro_aligner = new ProteinAlignment(-10,-0.2,MAX_LENGTH);
        if (args.length > 2)
            THRESHOLD = Double.parseDouble(args[2]);
        System.out.println("THRESHOLD = " + THRESHOLD);
        find_crossing_proteins();
        build_homology_groups();
        build_orthology_groups(pangenome_path);
    }
    
    public void find_crossing_proteins(){
        long[] penta_mer_node_ids;
        int[] code = new int[256];
        char[] aminoacids = new char[]
        {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
        for (int i = 0; i < 20; ++i)
            code[aminoacids[i]] = i;
        LinkedList<Node> proteins = new LinkedList();
        ResourceIterator<Node> proteins_iterator;
        PriorityQueue<Long> crossing_protein_ids = new PriorityQueue();
        Node protein_node, pentamer_node, db_node;
        long other_node_id, protein_node_id;
        int protein_length, total_proteins, pentamer_index;
        String protein;
        //StringBuilder pentamer = new StringBuilder();
        int[] seed = new int[5];
        int i, p, count, start, trsc;
        long crossing_protein_id, p_id;
        System.out.println("Finding crossing proteins...");
        try (Transaction tx = graphDb.beginTx()) {
            db_node = graphDb.findNodes(pangenome_label).next();
            if (db_node.hasProperty("penta_mer_node_ids"))
                penta_mer_node_ids = (long[])db_node.getProperty("penta_mer_node_ids");
            else
                penta_mer_node_ids = new long[3200000]; //20 ^ 5 
            proteins_iterator = graphDb.findNodes(mRNA_label);
            while (proteins_iterator.hasNext())
                proteins.add(proteins_iterator.next());
            tx.success();
            proteins_iterator.close();
        }
        total_proteins = proteins.size();
        for (p = 0; !proteins.isEmpty(); ) {
            try (Transaction tx = graphDb.beginTx()) {
                for (trsc = 0; !proteins.isEmpty()&& trsc < 2 * MAX_TRANSACTION_SIZE; ++trsc){
                    protein_node = proteins.remove();
                    if (protein_node.hasProperty("protein")){
                        protein_node_id = protein_node.getId();
                        ++p;
                        protein = (String)protein_node.getProperty("protein");
                        protein_length = protein.length();
                        for (i = 0; i < 4; ++i)
                            seed[i] = protein.charAt(i);
                        for (start = 0; start < protein_length - 5; ++start){
                            seed[(start + 4) % 5] = protein.charAt(start + 4);
                            pentamer_index = 0;
                            //pentamer.setLength(0);
                            for (i = 0; i < 5; ++i){
                                pentamer_index = pentamer_index * 20 + code[seed[(start + i) % 5]];
                                //pentamer.append((char)seed[(start + i) % 5]);
                            }
                            if (penta_mer_node_ids[pentamer_index] == 0){
                                pentamer_node = graphDb.createNode(pentamer_lable);
                                //pentamer_node.setProperty("index", pentamer_index);
                                //pentamer_node.setProperty("peptide", pentamer.toString());
                                penta_mer_node_ids[pentamer_index] = pentamer_node.getId();
                            } else {
                                pentamer_node = graphDb.getNodeById(penta_mer_node_ids[pentamer_index]);
                                for (Relationship rel: pentamer_node.getRelationships()){
                                    other_node_id = rel.getStartNode().getId();
                                    if (other_node_id != protein_node_id)
                                        crossing_protein_ids.add(other_node_id);
                                }
                            }
                            protein_node.createRelationshipTo(pentamer_node, RelTypes.visits);
                        }
                        if (!crossing_protein_ids.isEmpty()){
                            for (count = 0, crossing_protein_id = crossing_protein_ids.peek(); !crossing_protein_ids.isEmpty();){
                                p_id = crossing_protein_ids.remove();
                                if (crossing_protein_id != p_id){
                                    if (count > protein_length/20 + 1)
                                        protein_node.createRelationshipTo(graphDb.getNodeById(crossing_protein_id), RelTypes.crosses);
                                    crossing_protein_id = p_id;
                                    count = 1;
                                }
                                else
                                    ++count;
                            }
                            if (count > protein_length/20 + 1) // for the last range of IDs in the queue
                                protein_node.createRelationshipTo(graphDb.getNodeById(crossing_protein_id), RelTypes.crosses);
                        }
                    }  
                    if (p % 11 == 1)
                        System.out.print("\r" + p + "/" + total_proteins);
                }
                tx.success();
            }
        } // for protein
        try (Transaction tx = graphDb.beginTx()) {
            db_node.setProperty("penta_mer_node_ids", penta_mer_node_ids);
            tx.success();
        }
        System.out.println("\r" + p + "/" + total_proteins);
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
                    itr = start_protein.getRelationships(RelTypes.crosses).iterator();
                    while (itr.hasNext()) {
                        try (Transaction tx2 = graphDb.beginTx()) {
                            for (int trsc = 0; itr.hasNext() && trsc < MAX_TRANSACTION_SIZE; ++trsc) {
                                crossing_edge = itr.next();
                                crossing_protein = crossing_edge.getOtherNode(start_protein);
                                if(!crossing_protein.hasRelationship(RelTypes.has_homolog, Direction.INCOMING) && crossing_protein.hasProperty("protein")){
                                    similarity = protein_similarity((String)start_protein.getProperty("protein"), (String)crossing_protein.getProperty("protein"));
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
    double protein_similarity(String p1, String p2){
        int m = p1.length(), n = p2.length(), max_len = max(m,n);
        int i, parts_num = 1, part_len1, part_len2;
        double score;
        if (max_len > MAX_LENGTH){
            parts_num = (max_len / MAX_LENGTH) + (max_len % MAX_LENGTH == 0 ? 0 : 1);
            part_len1 = m / parts_num;
            part_len2 = n / parts_num;
            for (score =0, i = 0; i < parts_num; ++i)
                score += pro_aligner.get_similarity(p1.substring(i * part_len1, min(m, (i + 1) * part_len1)),
                                                    p2.substring(i * part_len2, min(n, (i + 1) * part_len2)) );
            return score / parts_num;
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
                            compute_pairwise_similaities(members, graph_path);
                            executeCommand("mcl " + graph_path + " --abc -I 5 -o " + clusters_path);
                            new File(graph_path).delete();
                            clusters_file = new FileReader(clusters_path);
                            while (clusters_file.ready()){
                                orthology_group_node = graphDb.createNode(orthology_group_lable);
                                orthology_group_nodes.add(orthology_group_node);
                                ++num_groups;
                                make_group(clusters_file, orthology_group_node);
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
                        groups_file.write(Long.toString(orthology_group_node.getId()));
                        //degree = orthology_group_node.getDegree();
                        //if (degree > 0)
                            set_and_write_orthology_groups(orthology_group_node, copy_number, groups_file);
                        //else
                        //    orthology_group_node.delete();   
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
    private void compute_pairwise_similaities(LinkedList<Node> members, String graph_path){
        ListIterator<Node> itr1, itr2;
        Node protein1, protein2;
        Relationship crossing_edge;
        double similarity = 0;
        try (PrintWriter graph = new PrintWriter(graph_path)){
            for (itr1 = members.listIterator(); itr1.hasNext(); ){
                protein1 = itr1.next();
                for (itr2 = members.listIterator(itr1.nextIndex()); itr2.hasNext(); ){
                    try (Transaction tx = graphDb.beginTx()) {
                        for (int trsc = 0; itr2.hasNext() && trsc < MAX_TRANSACTION_SIZE; ++trsc) {
                            protein2 = itr2.next();
                            crossing_edge = get_edge(protein1, protein2, RelTypes.crosses);
                            if (crossing_edge != null){
                                similarity = (double)crossing_edge.getProperty("similarity", -1.0);
                                if(similarity == -1.0) // do not have similarity property
                                    similarity = protein_similarity((String)protein1.getProperty("protein"), (String)protein2.getProperty("protein"));
                                if (similarity > THRESHOLD){
                                    crossing_edge.setProperty("similarity", similarity); // to be usedful for group_singles()
                                // converts the similarity scores to a number between 0 and ( 100 - THRESHOLD) ^ 2.
                                // This increases the sensitivity of MCL clustering algorithm.
                                    similarity = similarity - THRESHOLD;
                                    graph.write(protein1.getId()+" "+protein2.getId()+" "+ (similarity * similarity) +"\n");
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
        Node protein_node, gene_node;
        num_members = 0;
        for (Relationship rel: orthology_group_node.getRelationships(Direction.OUTGOING)){
            ++num_members;
            protein_node = rel.getEndNode();
            ++copy_number[((int[])protein_node.getProperty("address"))[0]];
            gene_node = protein_node.getSingleRelationship(RelTypes.codes_for, Direction.INCOMING).getStartNode();
            groups_file.write("\t" + protein_node.getProperty("protein_number") + ":" + gene_node.getProperty("name"));
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
