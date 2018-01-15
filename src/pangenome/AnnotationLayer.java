/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import genome.SequenceDatabase;
import index.IndexPointer;

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
import java.util.Arrays;
import java.util.Date;
import org.neo4j.graphdb.Direction;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.ResourceIterator;
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
import static pantools.Pantools.GRAPH_DATABASE_PATH;
import static pantools.Pantools.K;
import static pantools.Pantools.MAX_TRANSACTION_SIZE;
import static pantools.Pantools.annotation_label;
import static pantools.Pantools.broken_protein_label;
import static pantools.Pantools.coding_gene_label;
import static pantools.Pantools.exon_label;
import static pantools.Pantools.feature_label;
import static pantools.Pantools.genome_label;
import static pantools.Pantools.db_node;
import static pantools.Pantools.gene_label;
import static pantools.Pantools.graphDb;
import static pantools.Pantools.intron_label;
import static pantools.Pantools.mRNA_label;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.rRNA_label;
import static pantools.Pantools.reverse_complement;
import static pantools.Pantools.startTime;
import static pantools.Pantools.tRNA_label;
import static pantools.Pantools.write_fasta;

/**
 * Implements all the functionalities related to the annotation layer of the pangenome
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class AnnotationLayer {
    
    private int num_proteins;
    private static char[] aminoacid_table;
    private static int[] binary;
    private static StringBuilder protein_builder;
    
    public AnnotationLayer(){
        aminoacid_table = new char[]
          {'K','N','K','N',
           'T','T','T','T',
           'R','S','R','S',
           'I','I','M','I',
           'Q','H','Q','H',
           'P','P','P','P',
           'R','R','R','R',
           'L','L','L','L',
           'E','D','E','D',
           'A','A','A','A',
           'G','G','G','G',
           'V','V','V','V',
           '*','Y','*','Y',
           'S','S','S','S',
           '*','C','W','C',
           'L','F','L','F'
          };
        binary=new int[256];
        binary['A'] = 0; 
        binary['C'] = 1; 
        binary['G'] = 2; 
        binary['T'] = 3; 
    // All degenerate bases will be replaces by 'A' in the translation 
        protein_builder = new StringBuilder();
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

    public class feature{
        Node node;
        String ID;
        public feature(Node n, String id){
            node = n;
            ID = id;
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
        Node annotation_node, genome_node;
        String[] fields;
        String gff_path, line;
        int[] address = new int[4];
        int k, degree;
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
        num_proteins = 0;
        try (BufferedReader gff_paths = new BufferedReader(new FileReader(gff_paths_file))) {
            log_file = new BufferedWriter(new FileWriter(pangenome_path + "/annotation.log"));
            if (! new File(pangenome_path + "/proteins").exists())
                Files.createDirectory(Paths.get(pangenome_path + "/proteins"));
            System.out.println("genome\tgenes\tmRNAs\ttRNAs\trRNAs");
            while (gff_paths.ready()) // for each gff file
            {
                line = gff_paths.readLine().trim();
                if (line.equals(""))
                    continue;
                fields = line.split("\\s");
                address[0] = Integer.parseInt(fields[0]);
                gff_path = fields[1];
                if (! new File(gff_path).exists()){
                    log_file.write("Genome "+address[0]+"'s GFF file not found.");
                    System.out.println("Genome "+address[0]+"'s GFF file not found.");
                    continue;
                }
                try (Transaction tx = graphDb.beginTx()) {
                    annotation_node = graphDb.createNode(annotation_label);
                    genome_node = graphDb.findNodes(genome_label,"number",address[0]).next();
                    annotation_node.createRelationshipTo(genome_node, RelTypes.annotates);
                    degree = genome_node.getDegree(RelTypes.annotates, Direction.INCOMING);
                    annotation_node.setProperty("date", new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(new Date()));
                    annotation_node.setProperty("path", gff_path);
                    annotation_node.setProperty("genome", address[0]);
                    annotation_node.setProperty("number", degree);
                    annotation_node.setProperty("identifier", address[0] + "_" + degree);
                    tx.success();
                }
                parse_gff(address, address[0] + "_" + degree, log_file, gff_path, pangenome_path, k);
            } // for genomes
            gff_paths.close();
            log_file.close();
        } catch (IOException ioe) {
            System.out.println(ioe.getMessage());
            System.out.println("Could not open " + gff_paths_file);
        }
        try (Transaction tx = graphDb.beginTx()) {
            db_node.setProperty("num_proteins", num_proteins);
            tx.success();
        }
        graphDb.shutdown();
        genomeDb.close();
        // delete the database transaction files
        File directory = new File(pangenome_path + GRAPH_DATABASE_PATH);
        for (File f : directory.listFiles())
            if (f.getName().startsWith("neostore.transaction.db."))
                f.delete();
        System.out.println("Annotation finished (see annotation.log).");
        System.out.println("Annotated proteins available in directory proteins.");
    }

    private void parse_gff(int[] address, String annotation_id, BufferedWriter log_file, String gff_path, String pangenome_path, int k){
        int i, trsc, num_genes, num_mRNAs, num_tRNAs, num_rRNAs, feature_len, offset;
        long seq_len;
        String sequence_id, current_sequence_id=null, attribute;
        Node node, gene_node, rna_node, feature_node, parent_node = null;
        Relationship rel;
        String[] fields,parent_ids;
        String strand, line, ID;
        LinkedList<feature> gene_nodes = new LinkedList();
        LinkedList<feature> rna_nodes = new LinkedList();
        StringBuilder attr = new StringBuilder();
        IndexPointer start_ptr, stop_ptr;
        num_genes = num_mRNAs = num_tRNAs = num_rRNAs = 0;
        try (BufferedReader in = new BufferedReader(new FileReader(gff_path))) {
            // for each record of gff file
            while (in.ready())
            {
                try (Transaction tx2 = graphDb.beginTx()) {
                    for (trsc = 0; trsc < MAX_TRANSACTION_SIZE && in.ready(); ++trsc) {
                        line = in.readLine().trim();
                        if (line.equals("") || line.charAt(0) == '#') // if line is empty or a comment skip it
                            continue;
                        if (line.contains("\t")){
                            fields = line.split("\\t");
                            attribute = fields[fields.length - 1];
                        } else {
                            log_file.write("This line is not tab-delimited and is split by whitespaces: \n" + line + "\n"); 
                            fields = line.split("\\s");
                            attr.setLength(0);
                            for (i = 7; i < fields.length; ++i)
                                attr.append(fields[i]).append(" ");
                            attribute = attr.toString().trim();
                        }
                        address[2] = Integer.parseInt(fields[3]); // begin
                        address[3] = Integer.parseInt(fields[4]); // end
                        strand = fields[6];
                        feature_len = address[3] - address[2] + 1;
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
                        ID = get_property(attribute,"ID");
                        feature_node.setProperty("ID", ID);
                        feature_node.setProperty("attribute", attribute);
                        feature_node.setProperty("type", fields[2]);
                        feature_node.setProperty("name", get_property(attribute,"Name"));
                        feature_node.setProperty("annotation_id", annotation_id);
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
                            gene_nodes.addFirst(new feature(feature_node, ID));
                            // adding gene_node id to the sequence node
                            //genes_list[address[0]][address[1]].add(gene_node.getId());
                            ++num_genes;
                        } else  switch(fields[2]) {
                                    case "mRNA":
                                        ++num_mRNAs;
                                        feature_node.addLabel(mRNA_label);
                                        if (parent_node != null)
                                            parent_node.createRelationshipTo(feature_node, RelTypes.codes_for);
                                        rna_nodes.addFirst(new feature(feature_node, ID));
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
                                                    rna_node.setProperty("annotation_id", annotation_id);
                                                    rna_node.setProperty("genome",address[0]);
                                                    ++num_mRNAs;
                                                    rna_nodes.addFirst(new feature(rna_node, ID));                                           
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
            set_protein_sequences(gene_nodes, address[0], log_file, pangenome_path, k);
            log_file.write("Genome "+address[0] + " : " + num_genes + " genes\t" + num_mRNAs + " mRNAs\t" + num_tRNAs + " tRNAs\t" + num_rRNAs + " rRNAs\n");
            log_file.write("----------------------------------------------------\n");
        }catch (IOException ioe) {
            System.out.println(ioe.getMessage());
            System.out.println("Could not open " + gff_path + "!");
        }
    }
    
    /**
     * Gives the node of an annotated feature, given its ID
     * @param nodes A list of annotated features
     * @param id ID of the query feature
     * @return The annotation node (gene or RNA) of the query feature
     */
    private Node get_node_by_id(LinkedList<feature> features, String id){
        ListIterator<feature> itr = features.listIterator();
        feature f;
        while (itr.hasNext()){
            f = itr.next();
            if (f.ID.equals(id.toLowerCase())){
                return f.node;
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
    private void set_protein_sequences(LinkedList<feature> gene_nodes, int genome, BufferedWriter log_file, String pangenome_path, int K){
        IntPairComparator comp = new IntPairComparator();
        PriorityQueue<int[]> pq = new PriorityQueue(comp);
        Node cds_node, mrna_node, gene_node;
        Relationship start_edge;
        IndexPointer start_ptr = new IndexPointer();
        int trsc, isoforms_num, gene_start_pos;
        String annotaion_id, protein;
        int[] address, begin_end;
        StringBuilder gene_builder = new StringBuilder();
        StringBuilder rna_builder = new StringBuilder();
        try (BufferedWriter out = new BufferedWriter(new FileWriter(pangenome_path + "/proteins/proteins_" + genome + ".fasta"))) {
            System.out.print("Adding protein sequences...");
            while (!gene_nodes.isEmpty()){
                try (Transaction tx = graphDb.beginTx()) {
                    for (trsc = 0; !gene_nodes.isEmpty()&& trsc < 2 * MAX_TRANSACTION_SIZE; ++trsc){
                        gene_node = gene_nodes.remove().node;
                        address = (int[])gene_node.getProperty("address");
                        gene_start_pos = address[2];
                        // extract gene sequence as appears in the sequence and connects it to its nodes
                        //gene_builder.setLength(0);
                        start_edge = gene_node.getSingleRelationship(RelTypes.starts, Direction.OUTGOING);
                        //start_ptr.node_id = start_edge.getEndNode().getId();
                        //start_ptr.offset = (int)start_edge.getProperty("offset");
                        //start_ptr.canonical = (boolean)start_edge.getProperty("forward");
                        //extract_sequence(gene_builder, start_ptr, address, K);
                        address[2] -= 1;
                        address[3] -= 1;
                        genomeDb.get_sequence(gene_builder, address, (boolean)start_edge.getProperty("forward"));
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
                                protein = translate(rna_builder);
                                if (protein.length() > 0){
                                    ++num_proteins;
                                    if (num_proteins % 11 == 1)
                                        System.out.print("\rAdding protein sequences... " + num_proteins);
                                    annotaion_id = (String)mrna_node.getProperty("annotation_id");
                                    out.write(">G" + annotaion_id.replace('_', 'A') + "P" + num_proteins + "\n");
                                    write_fasta(out, protein, 70);
                                    mrna_node.setProperty("protein", protein.toString());
                                    mrna_node.setProperty("protein_ID", "G" + annotaion_id.replace('_', 'A') + "P" + num_proteins);
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
            if (fields[i].toLowerCase().startsWith(property.toLowerCase())){
                if (fields[i].contains("="))
                    return fields[i].trim().split("=")[1].toLowerCase();
                else
                    return fields[i].trim().split("\\s")[1].toLowerCase();
            }                    
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
     * Retrieves sequence of the genes sequence from the pangenome and stores them in a FASTA file. 
     * 
     * @param annotation_records_file a text file containing annotation records of the genes to be retrieved.
     * @param pangenome_path Path to the database folder
     */
    public void retrieve_genes(String annotation_records_file, String pangenome_path) {
        if (new File(pangenome_path + GRAPH_DATABASE_PATH).exists()) {
            BufferedReader in;
            ResourceIterator<Node> gene_nodes;
            Relationship rstart;
            Node start, gene;
            String record;
            String line, out_file_name;
            boolean strand;
            int i, j, num_genes, begin,end;
            int[] address;
            String[] fields;
            StringBuilder gene_seq;
            StringBuilder attr = new StringBuilder();

            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(pangenome_path + GRAPH_DATABASE_PATH))
                    .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
            registerShutdownHook(graphDb);
            startTime = System.currentTimeMillis();
            genomeDb = new SequenceDatabase(pangenome_path + GENOME_DATABASE_PATH);
            try (Transaction tx = graphDb.beginTx()) {
                K = (int) graphDb.findNodes(pangenome_label).next().getProperty("k_mer_size");
                tx.success();
            }
            num_genes = 0;
            gene_seq = new StringBuilder();
            try {
                in = new BufferedReader(new FileReader(annotation_records_file));
                while (in.ready()) {
                    line = in.readLine().trim();
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
                        line = in.readLine().trim();
                        if (line.equals("")) {
                            continue;
                        }
                        if (line.contains("\t")){
                            fields = line.split("\\t");
                            records[i] = fields[fields.length - 1];
                        } else {
                            fields = line.split("\\s");
                            attr.setLength(0);
                            for (int x = 7; x < fields.length; ++x)
                                attr.append(fields[x]).append(" ");
                            records[i] = attr.toString().trim();
                        }
                        ++i;
                    }
                    in.close();
                    Arrays.sort(records);
                } catch (IOException ioe) {
                    System.out.println("Failed to read file names!");
                    System.exit(1);
                }
                try {
                    fields = annotation_records_file.split("\\/");
                    out_file_name = pangenome_path + "/" + fields[fields.length - 1] + ".fasta";
                    BufferedWriter out = new BufferedWriter(new FileWriter(out_file_name));
                    // for all the genes in the database    
                    for (i = j = 0, gene_nodes = graphDb.findNodes(gene_label); gene_nodes.hasNext();) {
                        gene = gene_nodes.next();
                        record = (String) gene.getProperty("attribute");
                        if (record != null && Arrays.binarySearch(records, record) >= 0) // gene is in the records
                        {
                            ++i;
                            //System.out.println(record);
                            rstart = gene.getSingleRelationship(RelTypes.starts, Direction.OUTGOING);
                            start = rstart.getEndNode();
                            address = (int[]) gene.getProperty("address");
                            begin = address[2];
                            end = address[3];
                            strand = gene.getProperty("strand").toString().equals("+");
                            //extract_sequence(gene_seq, new IndexPointer(start.getId(), (boolean) rstart.getProperty("forward"), (int) rstart.getProperty("offset"),-1l), address, K);//
                            address[2] -= 1;
                            address[3] -= 1;
                            genomeDb.get_sequence(gene_seq, address, strand);
                            //genomeDb=new sequence_database(pangenome_path+GENOME_DATABASE_PATH);
                            //if(gene_seq.toString().equals(genomeDb.get_sequence(genome, sequence, begin-1, end-begin+1, strand))
                            //|| gene_seq.toString().equals(genomeDb.get_sequence(genome, sequence, begin-1, end-begin+1, !strand)) )//gene_seq.length() == end-begin+1)//
                            if (gene_seq.length() == end - begin + 1) {
                                ++j;
                                out.write(">" + record + "\n");
                                if (strand) {
                                    write_fasta(out, gene_seq.toString(), 70);
                                } else {
                                    reverse_complement(gene_seq);
                                    write_fasta(out, gene_seq.toString(), 70);
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
                    System.out.println(j + " out of " + i + " genes found and retrieved successfully (See " + out_file_name + ").");
                    out.close();
                } catch (IOException ioe) {
                    System.out.println("Failed to read file names!");
                    System.exit(1);
                }
                tx.success();
            }
            graphDb.shutdown();
            genomeDb.close();
        } else {
            System.out.println("No database found in " + pangenome_path);
            System.exit(1);
        }
    }

    /**
     * Translates a coding nucleotide sequence to a protein sequence.
     * 
     * @param mRNA The sequence of the codinf RNA
     * @return The protein sequence
     */
    public static String translate(StringBuilder mRNA){
        protein_builder.setLength(0);
        int i;   
        for(i=0;i<=mRNA.length()-3;i+=3)
            protein_builder.append(aminoacid_table[binary[mRNA.charAt(i)]*16+binary[mRNA.charAt(i+1)]*4+binary[mRNA.charAt(i+2)]]);
        return protein_builder.toString();
    }    

    public void remove_annotaions(String gff_paths_file, String pangenome_path) {
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
