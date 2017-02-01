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
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import org.neo4j.graphdb.Direction;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import static org.neo4j.graphdb.factory.GraphDatabaseSettings.keep_logical_logs;
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
import static pantools.Pantools.print_peak_memory;
import static pantools.Pantools.reverse_complement;
import static pantools.Pantools.tree_node_lable;
import static pantools.Pantools.tree_root_lable;
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
    public String PATH;
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

    /**
     * Adds nodes related to different genomic features to the pangenome.
     * 
     * @param gff_paths_file A text file listing the paths to the annotation files
     * @param PATH Path to the database folder
     */
    public void annotate(String gff_paths_file, String PATH) {
        if (! new File(PATH + GRAPH_DATABASE_PATH).exists()) {
            System.out.println("No database found in " + PATH);
            System.exit(1);
        }
        int i, j, num_genes, num_mRNAs, num_tRNAs, num_rRNAs, total_genes=0, feature_len;
        String sequence_id, current_sequence_id=null, attribute;
        Node db_node, gene_node = null, rna_node, cds_node;
        String[] fields,parent_ids;
        String strand, gff_file, line, ID;
        List<Long>[][] genes_list;
        int[] address = new int[4];
        long[] genes_array;
        LinkedList<Node> gene_nodes = new LinkedList();
        LinkedList<Node> rna_nodes = new LinkedList();
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH + GRAPH_DATABASE_PATH))
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
        genomeDb = new SequenceDatabase(PATH + GENOME_DATABASE_PATH);
        genes_list = new ArrayList[genomeDb.num_genomes + 1][];
        for (i = 1; i <= genomeDb.num_genomes; ++i) {
            genes_list[i] = new ArrayList[genomeDb.num_sequences[i] + 1];
            for (j = 1; j <= genomeDb.num_sequences[i]; ++j) {
                genes_list[i][j] = new ArrayList();
            }
        }
        try (BufferedReader gff_paths = new BufferedReader(new FileReader(gff_paths_file))) {
            BufferedWriter log_file = new BufferedWriter(new FileWriter(PATH + "/annotation.log"));
            System.out.println("genome\tgenes\tmRNAs\ttRNAs\trRNAs");
            for (address[0] = 1; gff_paths.ready() && address[0] <= genomeDb.num_genomes; ++address[0]) // for each gff file
            {
                gene_nodes.clear();
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
                        for (i = 0; i < MAX_TRANSACTION_SIZE && in.ready(); ++i) {
                            line = in.readLine();
                            //System.out.print(line);
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
                        // if sequence found    
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
                                    
                                    gene_nodes.addFirst(gene_node);
                                // adding gene_node id to the sequence node
                                    genes_list[address[0]][address[1]].add(gene_node.getId());
                                    ++total_genes;
                                    ++num_genes;
                                } else if(fields[2].endsWith("RNA") || fields[2].equals("pseudogenic_transcript")) {
                                    gene_node = get_node_by_id(gene_nodes,get_property(attribute,"Parent"));
                                    if (gene_node != null){
                                        rna_node = graphDb.createNode(RNA_label);
                                        rna_node.setProperty("ID", get_property(attribute,"ID"));
                                        rna_node.setProperty("address", address);
                                        rna_node.setProperty("length", feature_len);
                                        rna_node.setProperty("type", fields[2]);
                                        rna_nodes.addFirst(rna_node);
                                        gene_node.createRelationshipTo(rna_node, RelTypes.codes_for);
                                        switch (fields[2]){
                                            case "mRNA":
                                                ++num_mRNAs;
                                                break;
                                            case "tRNA":
                                                ++num_tRNAs;
                                                break;
                                            case "rRNA":
                                                ++num_rRNAs;
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
                                                    rna_node.setProperty("ID", gene_node.getProperty("ID"));
                                                    rna_node.setProperty("address", gene_node.getProperty("address"));
                                                    rna_node.setProperty("length", gene_node.getProperty("length"));
                                                    rna_node.setProperty("type", "mRNA");
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
                            if (i % 500 == 1)
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
                set_protein_sequences(gene_nodes, address[0], log_file, PATH);
            } // for genomes
            gff_paths.close();
            log_file.close();
        } catch (IOException ioe) {
            System.out.println(ioe.getMessage());
            System.exit(1);
        }
    // setting the genes property of sequence nodes    
        try (Transaction tx4 = graphDb.beginTx()) {
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
        }
        graphDb.shutdown();
        genomeDb.close();
        // delete the database transaction files
        File directory = new File(PATH + GRAPH_DATABASE_PATH);
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
    private void set_protein_sequences(LinkedList<Node> gene_nodes, int genome, BufferedWriter log_file, String PATH){
        LinkedList<Node> crossing_genes = new LinkedList();
        ListIterator<Node> itr = gene_nodes.listIterator();
        PairComparator comp = new PairComparator();
        PriorityQueue<int[]> pq = new PriorityQueue(comp);
        protein_builder pb = new protein_builder();
        Node cds_node, rna_node, gene_node;
        int count, isoforms_num, gene_start_pos, protein_num = 0;
        int[] address, begin_end;
        StringBuilder rna_builder = new StringBuilder();
        StringBuilder gene_builder = new StringBuilder();
        String protein, rna_type, strand;
        try (BufferedWriter out = new BufferedWriter(new FileWriter(PATH + "/Genome_" + genome + "_proteins.fasta"))) {
            System.out.print("Adding protein sequences...");
            for (count = 0; itr.hasNext(); count = 0){
                try (Transaction tx = graphDb.beginTx()) {
                    for (;itr.hasNext() && count < MAX_TRANSACTION_SIZE; ++count){
                        //System.out.print("count = " + count + "\r");
                        gene_node = itr.next();
                        address = (int[])gene_node.getProperty("address");
                        strand = (String)gene_node.getProperty("strand");
                        gene_start_pos = address[2];
                        // extract gene sequence as appears in the sequence and connects it to its nodes
                        gene_builder.setLength(0);
                        connect_gene_to_nodes(gene_builder, address, gene_node, crossing_genes);
                        isoforms_num =0;
                        for (Relationship r1: gene_node.getRelationships(RelTypes.codes_for, Direction.OUTGOING)) {
                            ++isoforms_num;
                            rna_node = r1.getEndNode();
                            rna_type = (String)rna_node.getProperty("type");
                            switch (rna_type){
                                case "mRNA":
                                    if (rna_node.hasRelationship(RelTypes.contributes_to, Direction.INCOMING)){
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
                                        protein = pb.translate(gene_node.getProperty("strand").equals("+")?rna_builder.toString():reverse_complement(rna_builder.toString()));
                                        if (protein.startsWith("M") && protein.endsWith("*")){
                                            ++protein_num;
                                            if (protein_num % 11 == 1){
                                                System.out.print("\rAdding protein sequences... " + protein_num);
                                            }
                                            out.write(">" + genome + "_" + protein_num + "\n");
                                            write_fasta(out, protein, 70);
                                            rna_node.setProperty("protein", protein);
                                            rna_node.setProperty("protein_number", genome + "_" + protein_num);
                                            rna_node.setProperty("protein_length", protein.length());
                                        } else {
                                            log_file.write("Protein ID = " + rna_node.getProperty("ID")+" miss-annotated!\n");
                                            rna_node.addLabel(broken_protein_label);
                                        }
                                    }
                                    break;
                                case "tRNA": case"rRNA":
                                    rna_node.setProperty("sequence", strand.equals("+")?gene_builder.toString():reverse_complement(gene_builder.toString()));
                                    break;
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
        boolean found = false;
        int sequence;
        for (found = false, sequence = 1; !found && sequence <= genomeDb.num_sequences[genome]; ++sequence) {
            if (genomeDb.sequence_titles[genome][sequence].contains(name)) {
                found = true;
            }
        }
        --sequence;
        if (found) {
            return sequence;
        } else {
            return -1;
        }
    }

    /**
     * Connect the annotated gene to all the nodes its sequence passes through in the pangenome.
     * At the same time assembles the sequence of the gene
     * 
     * @param gene_builder String builder which will contain the sequence of the gene  
     * @param address An array containing {genome, sequence, begin, end}
     * @param gene_node The gene node itself
     */
    public void connect_gene_to_nodes(StringBuilder gene_builder, int[] address, Node gene_node, LinkedList<Node> crossing_genes) {
        Relationship rel;
        Node neighbor, crossing_gene, node;
        int loc, node_len, neighbor_len, seq_len = address[3] - address[2] + 1, starts_at;
        String rel_name;
        gene_builder.setLength(0);
        loc = address[2]-1;
        IndexPointer start_ptr = locate(address);
        starts_at = start_ptr.position;
        node = graphDb.getNodeById(start_ptr.node_id);
        node_len = (int) node.getProperty("length");
        gene_node.createRelationshipTo(node, RelTypes.visits).setProperty("position", loc);
        rel = gene_node.createRelationshipTo(node, RelTypes.starts);
        rel.setProperty("starts_at", starts_at);
        rel.setProperty("position", loc);
        rel.setProperty("forward", start_ptr.canonical);
    //  finds new genes crossing gene_node at node   
        for (Relationship r: node.getRelationships(Direction.INCOMING,RelTypes.visits)) {
            crossing_gene = r.getStartNode();
            if (!crossing_genes.contains(crossing_gene))
                crossing_genes.add(crossing_gene);
        }
        
        if (start_ptr.canonical) {
            if (starts_at + seq_len - 1 <= node_len - 1) {
                loc += append_fwd(gene_builder, (String) node.getProperty("sequence"), starts_at, starts_at + seq_len - 1);
            } else {
                loc += append_fwd(gene_builder, (String) node.getProperty("sequence"), starts_at, node_len - 1);
            }
        } else {
            if (starts_at - (seq_len - 1) >= 0) {
                loc += append_rev(gene_builder, (String) node.getProperty("sequence"), starts_at - (seq_len - 1), starts_at);
            } else {
                loc += append_rev(gene_builder, (String) node.getProperty("sequence"), 0, starts_at);
            }
        }
        while (gene_builder.length() < seq_len ) {
            address[2] = loc - K + 1;
            rel = get_outgoing_edge(node, address);
            if (rel == null) // happened for miss-annotaion
                break;
            neighbor = rel.getEndNode();
            rel_name = rel.getType().name();
            neighbor_len = (int) neighbor.getProperty("length");
            if (rel_name.charAt(1) == 'F')
                if (gene_builder.length() + neighbor_len - K + 1 > seq_len) 
                    loc += append_fwd(gene_builder, (String) neighbor.getProperty("sequence"), K - 1, seq_len - gene_builder.length() + K - 2);
                else 
                    loc += append_fwd(gene_builder, (String) neighbor.getProperty("sequence"), K - 1, neighbor_len - 1);
            else if (gene_builder.length() + neighbor_len - K + 1 > seq_len) 
                loc += append_rev(gene_builder, (String) neighbor.getProperty("sequence"), neighbor_len - K - (seq_len - gene_builder.length()) + 1, neighbor_len - K);
            else 
                loc += append_rev(gene_builder, (String) neighbor.getProperty("sequence"), 0, neighbor_len - K);
            node = neighbor;
        //  finds new genes crossing gene_node at node   
            for (Relationship r: node.getRelationships(Direction.INCOMING,RelTypes.visits)) {
                crossing_gene = r.getStartNode();
                if (!crossing_genes.contains(crossing_gene))
                    crossing_genes.add(crossing_gene);
            }
            gene_node.createRelationshipTo(node, RelTypes.visits).setProperty("position", loc);
            
        } // while
        while (!crossing_genes.isEmpty()){
            crossing_gene = crossing_genes.remove();
            if (!crossing_gene.equals(gene_node))
                gene_node.createRelationshipTo(crossing_gene, RelTypes.corsses);
        }
        gene_node.setProperty("stop_node", node.getId());
    }

    /**
     * Groups the similar genes into homology and orthology groups
     * @param args The command line arguments:
     * args[1] Path to the database folder
     * args[2] THRESHOLD to be replaced by the default value, if given
     * args[3] LEN_FACTOR to be replaced by the default value, if given
     */
    public void group(String[] args) {
        PATH = args[1];
        if (args.length > 2)
            THRESHOLD = Double.parseDouble(args[2]);
        if (new File(PATH + GRAPH_DATABASE_PATH).exists()) {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH + GRAPH_DATABASE_PATH))
                    .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
            registerShutdownHook(graphDb);
            startTime = System.currentTimeMillis();
            seq_aligner = new SequenceAlignment(-5,-0.1,4.0,-2.0,MAX_LENGTH);
            pro_aligner = new ProteinAlignment(-10,-0.2,MAX_LENGTH);
            build_homology_groups();
            build_orthology_groups();
            drop_edges_of_type(RelTypes.resembles);
        } else {
            System.out.println("pangenome database not found!");
            System.exit(0);
        }
    }

    /**
     * Creates homology nodes which connect the homologous coding genes
     */
    private void build_homology_groups(){
        int num_groups = 0, total_genes = 0;
        Node gene_node;
        ResourceIterator<Node> genes_iterator;
        System.out.println("Grouping coding genes...");
        System.out.println("THRESHOLD = "+THRESHOLD);
        try (Transaction tx = graphDb.beginTx()) {
            genes_iterator = graphDb.findNodes(coding_gene_label);
            System.out.println("Genes\tHomology groups");
            while (genes_iterator.hasNext()) {
                gene_node=genes_iterator.next();
            // To avoid having one gene in different groups    
                if (!gene_node.hasRelationship(RelTypes.has_homolog, Direction.INCOMING)) { 
                    num_groups += 1;
                    total_genes += group_homologs(gene_node);
                    System.out.print("\r" + total_genes + "\t" + num_groups);
                    if (num_groups % 1000 == 0){
                        System.gc ();
                        System.runFinalization ();
                    }
                }
            } // while 
            System.gc ();
            System.runFinalization ();
            tx.success();
            genes_iterator.close();
        } // transaction 1
        System.out.println("\nTotal time : " + (System.currentTimeMillis() - startTime) / 1000 + "." + (System.currentTimeMillis() - startTime) % 1000 + " seconds");
        print_peak_memory();
        System.out.println();
    }
    
    /**
     * Puts all the genes homologous to a query gene in one homology group
     * @param similar_groups An empty list to be used for storing groups with a member similar the the query_gene
     * @param crossing_genes An empty priority queue to be filled with the genes crossing the query_gene 
     * @param query_gene The gene to which similar ones should be found
     * @return The value by which the number of groups have been increased. (Could be a negative value if some groups get merged) 
     */
    private int group_homologs(Node start_gene){
        int num_members = 1;
        double similarity;
        Node homology_group_node, crossing_gene;
        Iterator<Relationship> itr;
        Queue<Node> homologs = new LinkedList();
        try (Transaction tx1 = graphDb.beginTx()) {
            homology_group_node = graphDb.createNode(homology_group_lable);
            homology_group_node.createRelationshipTo(start_gene, RelTypes.has_homolog);
            homologs.add(start_gene);
            // for all the candidates with some shared node with the gene    
            while (!homologs.isEmpty()) {
                start_gene = homologs.remove();
                itr = start_gene.getRelationships(RelTypes.corsses).iterator();
                while (itr.hasNext()) {
                    try (Transaction tx2 = graphDb.beginTx()) {
                        for (int trsc = 0; itr.hasNext() && trsc < MAX_TRANSACTION_SIZE; ++trsc) {
                            crossing_gene = itr.next().getOtherNode(start_gene);
                            if(crossing_gene.hasLabel(coding_gene_label) && !crossing_gene.hasRelationship(RelTypes.has_homolog, Direction.INCOMING)){
                                similarity = compute_similarity(start_gene, crossing_gene);
                                start_gene.createRelationshipTo(crossing_gene, RelTypes.resembles).setProperty("similarity", similarity);
                                if ( similarity > THRESHOLD ) {
                                    homology_group_node.createRelationshipTo(crossing_gene, RelTypes.has_homolog);
                                    num_members++;
                                    homologs.add(crossing_gene);
                                }
                            }
                        }
                        tx2.success(); 
                    }
                }
            }// while
            homology_group_node.setProperty("num_members", homology_group_node.getDegree(Direction.OUTGOING));
            tx1.success(); 
        }
        return num_members;
    }  
    
    
    /**
     * Given two genes calculates the normalized similarity score between them. (score <= 100)
     * @param gene1 The first coding gene
     * @param gene2 The second coding gene
     * @return The normalized similarity score which is at most 1
     */
    private double compute_similarity(Node gene1 , Node gene2){
        Node mRNA_node1, mRNA_node2;
        String p1 = null, p2 = null;
        int mRNA_len1, mRNA_len2, len;
        mRNA_len1 = 0;
        for (Relationship r1: gene1.getRelationships(Direction.OUTGOING,RelTypes.codes_for)) {
            mRNA_node1 = r1.getEndNode();
            if ( mRNA_node1.hasRelationship(RelTypes.contributes_to, Direction.INCOMING) && !mRNA_node1.hasLabel(broken_protein_label) ){
                len = (int)mRNA_node1.getProperty("protein_length");
                if (len > mRNA_len1){
                    mRNA_len1 = len;
                    p1 = (String)mRNA_node1.getProperty("protein");
                }
            }
        }
        mRNA_len2 = 0;
        for (Relationship r2: gene2.getRelationships(Direction.OUTGOING,RelTypes.codes_for)) {
            mRNA_node2 = r2.getEndNode();
            if ( mRNA_node2.hasRelationship(RelTypes.contributes_to, Direction.INCOMING) && !mRNA_node2.hasLabel(broken_protein_label) ){
                len = (int)mRNA_node2.getProperty("protein_length");
                if (len > mRNA_len2){
                    mRNA_len2 = len;
                    p2 = (String)mRNA_node2.getProperty("protein");
                }
            }
        }
        if ( mRNA_len1 > 0 && mRNA_len2 > 0 && Math.abs(mRNA_len1 - mRNA_len2) <= Math.min(mRNA_len1, mRNA_len2))
            return protein_similarity(p1, p2);
        else
            return 0;
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
     * Given two nucleotide sequences, calculates the normalized similarity score between them which is less or equal to 1.
     * @param s1 The first sequence
     * @param s2 The second sequence
     * @return The normalized similarity score which is less or equal to 1
     */
    double sequence_similarity(String s1, String s2){
        int m = s1.length(), n = s2.length(), max_len = max(m,n);
        int i, parts_num = 1, part_len1, part_len2;
        double score = 0;
        if (max_len > MAX_LENGTH){
            parts_num = (max_len / MAX_LENGTH) + (max_len % MAX_LENGTH == 0 ? 0 : 1);
        }
        part_len1 = m / parts_num;
        part_len2 = n / parts_num;
        for (i = 0; i < parts_num; ++i)
            score += seq_aligner.get_similarity(s1.substring(i * part_len1, min(m, (i + 1) * part_len1)),
                                                s2.substring(i * part_len2, min(n, (i + 1) * part_len2)) );
        return score / parts_num;
    }
    
    /**
     * divides the pre-calculated homology groups into orthology groups using MCL algorithm 
     */   
    void build_orthology_groups() {
        int num_groups = 0, degree, trsc;
        int[] copy_number;
        ResourceIterator<Node> homology_group_nodes;
        LinkedList<Node> orthology_group_nodes = new LinkedList();
        LinkedList<Node> members = new LinkedList();
        Node homology_group_node, orthology_group_node, node, single_gene;
        String graph_path, clusters_path;
        System.out.println("Building orthology groups...");
        try (Transaction tx1 = graphDb.beginTx()) {
            copy_number = new int[(int)graphDb.findNodes(pangenome_label).next().getProperty("num_genomes")+1];
            homology_group_nodes = graphDb.findNodes(homology_group_lable);
            while (homology_group_nodes.hasNext()){
                try (Transaction tx2 = graphDb.beginTx()) {
                    for (trsc = 0; homology_group_nodes.hasNext() && trsc < MAX_TRANSACTION_SIZE; ++trsc) {
                        homology_group_node = homology_group_nodes.next();
                        if( homology_group_node.getDegree() > 1){
                            graph_path = PATH + "/" + homology_group_node.getId() + ".graph";
                            clusters_path = PATH + "/" + homology_group_node.getId() + ".clusters";
                            for (Relationship rel: homology_group_node.getRelationships(Direction.OUTGOING))
                                members.add(rel.getEndNode());
                            compute_pairwise_similaities(members, graph_path);
                            executeCommand("mcl " + graph_path + " --abc -I 5 -o " + clusters_path);
                            new File(graph_path).delete();
                            try (FileReader clusters_file = new FileReader(clusters_path)){
                                while (clusters_file.ready()){
                                    orthology_group_node = graphDb.createNode(orthology_group_lable);
                                    orthology_group_nodes.add(orthology_group_node);
                                    ++num_groups;
                                    make_group(clusters_file, orthology_group_node);
                                }
                                clusters_file.close();
                                new File(clusters_path).delete();
                            }catch (IOException ex){
                            }
                        } else {
                            single_gene = homology_group_node.getSingleRelationship(RelTypes.has_homolog, Direction.OUTGOING).getEndNode();
                            orthology_group_node = graphDb.createNode(orthology_group_lable); 
                            ++num_groups;
                            orthology_group_node.createRelationshipTo(single_gene, RelTypes.has_ortholog);
                            orthology_group_nodes.add(orthology_group_node);
                        }
                        System.out.print("\rOrthology groups : " + num_groups);
                        if (num_groups % 1000 == 0){
                            System.gc ();
                            System.runFinalization ();
                        }
                    } // for
                    tx2.success();
                    homology_group_nodes.close();
                }
            }
            tx1.success();
        }
        System.out.println("\rOrthology groups : " + num_groups);
        System.out.println("After housing singletons : " + (num_groups - group_singletons(orthology_group_nodes)));
    // compute copy_number array in each group and delete singleton groups merged into the other ones and now empty 
        while (!orthology_group_nodes.isEmpty()){
            try (Transaction tx = graphDb.beginTx()) {
                for (trsc = 0; !orthology_group_nodes.isEmpty() && trsc < MAX_TRANSACTION_SIZE; ++trsc) {
                    orthology_group_node = orthology_group_nodes.remove();
                    degree = orthology_group_node.getDegree();
                    if (degree > 0)
                        set_group_properties(orthology_group_node, copy_number);
                    else
                        orthology_group_node.delete();   
                }
                tx.success();
            }
        }       
    }   

    int group_singletons(LinkedList<Node> orthology_group_nodes) {
        int num = 0;
        double best_similarity, similarity;
        Node orthology_group_node, gene_node, mate_gene_node, mate_homology_group, gene_homology_group;
        long best_mate_id;
        Relationship r;
        ListIterator<Node> itr = orthology_group_nodes.listIterator();
        while (itr.hasNext()){
            try (Transaction tx2 = graphDb.beginTx()) {
                for (int trsc = 0; itr.hasNext() && trsc < MAX_TRANSACTION_SIZE; ) {
                    orthology_group_node = itr.next();
                    if (orthology_group_node.getDegree() == 1){
                        ++trsc;
                        r = orthology_group_node.getSingleRelationship(RelTypes.has_ortholog, Direction.OUTGOING);
                        gene_node = r.getEndNode();
                        best_similarity = THRESHOLD/2;
                        best_mate_id = -1l;
                        for (Relationship rel: gene_node.getRelationships(RelTypes.resembles)){
                            mate_gene_node = rel.getOtherNode(gene_node);
                            similarity = (double)rel.getProperty("similarity");
                            if (similarity > best_similarity){
                                best_similarity = similarity;
                                best_mate_id = mate_gene_node.getId();

                            }
                        }
                        if (best_mate_id != -1l){
                            ++num;
                            mate_gene_node = graphDb.getNodeById(best_mate_id);
                        // delete single orthology edge, orthology node will be deleted in the main function    
                            r.delete();
                        // join the best orthology group    
                            mate_gene_node.getSingleRelationship(RelTypes.has_ortholog, Direction.INCOMING)
                                    .getStartNode().createRelationshipTo(gene_node, RelTypes.has_ortholog);

                            mate_homology_group = mate_gene_node.getSingleRelationship(RelTypes.has_homolog, Direction.INCOMING).getStartNode();
                            r = gene_node.getSingleRelationship(RelTypes.has_homolog, Direction.INCOMING);
                            gene_homology_group = r.getStartNode();
                            if (!mate_homology_group.equals(gene_homology_group)){
                            // delete the single homology node    
                                r.delete();
                                gene_homology_group.delete();
                            // join the best homology group    
                                mate_homology_group.createRelationshipTo(gene_node, RelTypes.has_homolog);
                                mate_homology_group.setProperty("num_members", 1 + (int) mate_homology_group.getProperty("num_members"));
                            }
                        }
                    }
                } 
                tx2.success();
            }
        }
        System.gc ();
        System.runFinalization ();
        return num;
    }   

    
    private void make_group(FileReader clusters_file, Node orthology_group_node){
        char ch;
        long id = 0;
        try{
            while ((ch = (char)clusters_file.read()) != '\n'){
                if (ch == '\t'){
                    orthology_group_node.createRelationshipTo(graphDb.getNodeById(id), RelTypes.has_ortholog);
                    id = 0;
                }
                else
                    id = id *10 + Character.getNumericValue(ch);
            }
            orthology_group_node.createRelationshipTo(graphDb.getNodeById(id), RelTypes.has_ortholog);
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
        Node gene1, gene2;
        double similarity;
        try (PrintWriter graph = new PrintWriter(graph_path)){
            for (itr1 = members.listIterator(); itr1.hasNext(); ){
                gene1 = itr1.next();
                for (itr2 = members.listIterator(itr1.nextIndex()); itr2.hasNext(); ){
                    try (Transaction tx = graphDb.beginTx()) {
                        for (int trsc = 0; itr2.hasNext() && trsc < MAX_TRANSACTION_SIZE; ++trsc) {
                            gene2 = itr2.next();
                            if (get_edge(gene1, gene2, RelTypes.resembles) == null)
                                similarity = compute_similarity(gene1, gene2);
                            else
                                similarity = read_similarity(gene1, gene2);
                            if (similarity > THRESHOLD){
                            // converts the similarity scores to a number between 0 and ( 100 - THRESHOLD) ^ 2.
                            // This increases the sensitivity of MCL clustering algorithm.
                                similarity = similarity - THRESHOLD;
                                graph.write(gene1.getId()+" "+gene2.getId()+" "+ (similarity * similarity) +"\n");
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
     * Retrieves the similarity score of two nodes. Nodes represent tips (genes) or internal nodes of the phylogeny tree. 
     * @param node1 The first node
     * @param node2 The second node
     * @return The similarity score between two nodes of the phylogeny tree.
     */
    private double read_similarity(Node node1, Node node2){
        Relationship rel = get_edge(node1, node2, RelTypes.resembles);
        return (double)rel.getProperty("similarity");
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
    void set_group_properties(Node orthology_group_node, int[] copy_number) {
        int i, num_members;
        Node gene_node;
        num_members = 0;
        for (Relationship rel: orthology_group_node.getRelationships(Direction.OUTGOING)){
            ++num_members;
            gene_node = rel.getEndNode();
            ++copy_number[((int[])gene_node.getProperty("address"))[0]];
        }
        orthology_group_node.setProperty("copy_number_variation", copy_number);
        orthology_group_node.setProperty("num_members", num_members);
        for (i=0; i<copy_number.length; ++i)
            copy_number[i] = 0;
    }    
    
    /**
     * Removes the relationships of a specific type from the pan-genome
     * @param type The relationship type of interest.
     */   
    void drop_edges_of_type(RelationshipType type) {
        int trsc;
        ResourceIterator<Relationship> rels;
        Relationship r;
        try (Transaction tx = graphDb.beginTx()) {
            rels = graphDb.getAllRelationships().iterator();
            tx.success();
                    System.gc ();
                    System.runFinalization ();
        }
        while (rels.hasNext()) {
            try (Transaction tx = graphDb.beginTx()) {
                for (trsc = 0; trsc < MAX_TRANSACTION_SIZE && rels.hasNext(); ) {
                    r = rels.next();
                    if (r.isType(type)){
                        ++trsc;
                        r.delete();
                    }
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
    private void registerShutdownHook(final GraphDatabaseService graphDb) {
        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                graphDb.shutdown();
            }
        });
    }
}
