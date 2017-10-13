/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import genome.SequenceDatabase;
import index.IndexPointer;
import protein.protein_builder;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.PriorityQueue;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.Date;
import org.neo4j.graphdb.Direction;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import static org.neo4j.graphdb.factory.GraphDatabaseSettings.keep_logical_logs;
import static pangenome.GenomeLayer.extract_sequence;

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
import static pangenome.GenomeLayer.locate;
import static pantools.Pantools.CDS_label;
import static pantools.Pantools.MAX_TRANSACTION_SIZE;
import static pantools.Pantools.annotation_label;
import static pantools.Pantools.broken_protein_label;
import static pantools.Pantools.coding_gene_label;
import static pantools.Pantools.exon_label;
import static pantools.Pantools.feature_label;
import static pantools.Pantools.assembly_label;
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
        int k;
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
            k = (int) db_node.getProperty("k_mer_size");
            num_nodes = (long) db_node.getProperty("num_nodes");
            num_edges = (long) db_node.getProperty("num_edges");
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
                    graphDb.findNodes(assembly_label,"number",address[0]).next().createRelationshipTo(annotation_node, RelTypes.has);
                    tx.success();
                }
                parse_gff(address, annotation_node.getId(), log_file, gff_file, pangenome_path, k);
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

    private void parse_gff(int[] address, long annotation_node_id, BufferedWriter log_file, String gff_file, String pangenome_path, int k){
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
                        start_ptr = locate(address, k);
                        offset = start_ptr.offset;
                        node = graphDb.getNodeById(start_ptr.node_id);
                        rel = feature_node.createRelationshipTo(node, RelTypes.starts);
                        rel.setProperty("offset", offset);
                        rel.setProperty("genomic_position", address[2]);
                        rel.setProperty("forward", start_ptr.canonical);
                        stop_ptr = locate(new int[]{address[0], address[1], address[3]}, k);
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
