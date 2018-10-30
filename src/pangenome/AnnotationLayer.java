/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import index.IndexDatabase;
import sequence.SequenceDatabase;
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
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import static org.neo4j.graphdb.factory.GraphDatabaseSettings.keep_logical_logs;
import static pangenome.GenomeLayer.locate;

import static pantools.Pantools.GENOME_DATABASE_PATH;
import static pantools.Pantools.FEATURE;
import static pantools.Pantools.RelTypes;
import static pantools.Pantools.labels;
import static pantools.Pantools.CDS_label;
import static pantools.Pantools.CONNECT_ANNOTATIONS;
import static pantools.Pantools.GRAPH_DATABASE_PATH;
import static pantools.Pantools.INDEX_DATABASE_PATH;
import static pantools.Pantools.K_SIZE;
import static pantools.Pantools.MAX_TRANSACTION_SIZE;
import static pantools.Pantools.PATH_TO_THE_ANNOTATIONS_FILE;
import static pantools.Pantools.PATH_TO_THE_GENOME_NUMBERS_FILE;
import static pantools.Pantools.PATH_TO_THE_PANGENOME_DATABASE;
import static pantools.Pantools.annotation_label;
import static pantools.Pantools.broken_protein_label;
import static pantools.Pantools.coding_gene_label;
import static pantools.Pantools.db_node;
import static pantools.Pantools.exon_label;
import static pantools.Pantools.feature_label;
import static pantools.Pantools.gene_label;
import static pantools.Pantools.genome_label;
import static pantools.Pantools.intron_label;
import static pantools.Pantools.mRNA_label;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.rRNA_label;
import static pantools.Pantools.reverse_complement;
import static pantools.Pantools.sequence_label;
import static pantools.Pantools.startTime;
import static pantools.Pantools.tRNA_label;
import static pantools.Pantools.write_fasta;
import sequence.SequenceScanner;
import static pantools.Pantools.genomeDb;
import static pantools.Pantools.genome_scanner;
import static pantools.Pantools.graphDb;
import static pantools.Pantools.indexDb;

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
    public void add_annotaions() {
        if (PATH_TO_THE_ANNOTATIONS_FILE == null){
            System.out.println("PATH_TO_THE_ANNOTATIONS_FILE is empty.");
            System.exit(1);
        }  
        if (! new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH).exists()) {
            System.out.println("No database found in " + PATH_TO_THE_PANGENOME_DATABASE);
            System.exit(1);
        }
        Node annotation_node, genome_node;
        String[] fields;
        String annotation_file, line;
        int[] address = new int[4];
        int degree;
        BufferedWriter log_file;
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        try (Transaction tx1 = graphDb.beginTx()) {
            db_node = graphDb.findNodes(pangenome_label).next();
            if (db_node == null) {
                System.out.println("Can not locate database node!");
                System.exit(1);
            }
            K_SIZE = (int) db_node.getProperty("k_mer_size");
            tx1.success();
        }
        startTime = System.currentTimeMillis();
        genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
        indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, "sorted");
        genome_scanner = new SequenceScanner(genomeDb, 1, 1, K_SIZE, indexDb.get_pre_len());
        num_proteins = 0;
        try{
            BufferedReader paths = new BufferedReader(new FileReader(PATH_TO_THE_ANNOTATIONS_FILE));
            log_file = new BufferedWriter(new FileWriter(PATH_TO_THE_PANGENOME_DATABASE + "/annotation.log"));
            if (! new File(PATH_TO_THE_PANGENOME_DATABASE + "/proteins").exists())
                Files.createDirectory(Paths.get(PATH_TO_THE_PANGENOME_DATABASE + "/proteins"));
            System.out.println("genome\tgenes\tmRNAs\ttRNAs\trRNAs");
            while (paths.ready()) // for each gff file
            {
                line = paths.readLine();
                if (line.equals(""))
                    continue;
                fields = line.split("\\s+");
                address[0] = Integer.parseInt(fields[0]);
                annotation_file = fields[1];
                if (! new File(annotation_file).exists()){
                    log_file.write("Genome "+address[0]+"'s annotation file not found.");
                    System.out.println("Genome "+address[0]+"'s annotation file not found.");
                    continue;
                }
                try (Transaction tx = graphDb.beginTx()) {
                    annotation_node = graphDb.createNode(annotation_label);
                    genome_node = graphDb.findNode(genome_label,"number",address[0]);
                    //genome_node = graphDb.findNodes(genome_label,"number",address[0]).next();
                    annotation_node.createRelationshipTo(genome_node, RelTypes.annotates);
                    degree = genome_node.getDegree(RelTypes.annotates, Direction.INCOMING);
                    annotation_node.setProperty("date", new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(new Date()));
                    annotation_node.setProperty("path", annotation_file);
                    annotation_node.setProperty("genome", address[0]);
                    annotation_node.setProperty("number", degree);
                    annotation_node.setProperty("identifier", address[0] + "_" + degree);
                    if (annotation_file.endsWith(".gff") || annotation_file.endsWith(".gff3")){
                        parse_gff(address, address[0] + "_" + degree, log_file, annotation_file);
                        annotation_node.setProperty("type", "GFF");
                    }
                    else if (annotation_file.endsWith(".gbk") || annotation_file.endsWith(".gbff")){
                        parse_gbk(address, address[0] + "_" + degree, log_file, annotation_file, PATH_TO_THE_PANGENOME_DATABASE + "/proteins");
                        annotation_node.setProperty("type", "GenBank");
                    } else {
                        System.err.println("Invalid extension for the annotaion file. Valid extensions are .gff .gff3 .gbk .gbff");
                    }
                    tx.success();
                }
            } // for genomes
            paths.close();
            log_file.close();
        } catch (IOException ioe) {
            System.out.println("Failed to open " + PATH_TO_THE_ANNOTATIONS_FILE);
        } catch (NumberFormatException ioe) {
            System.out.println("Wrong number format in " + PATH_TO_THE_ANNOTATIONS_FILE);
        }
        try (Transaction tx = graphDb.beginTx()) {
            db_node.setProperty("num_proteins", num_proteins);
            tx.success();
        }
        graphDb.shutdown();
        genomeDb.close();
        // delete the database transaction files
        File directory = new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH);
        for (File f : directory.listFiles())
            if (f.getName().startsWith("neostore.transaction.db."))
                f.delete();
        System.out.println("Annotation finished (see annotation.log).");
        System.out.println("Annotated proteins available in directory proteins.");
    }

    private void parse_gff(int[] address, String annotation_id, BufferedWriter log_file, String gff_path){
        int i, trsc, num_genes, num_mRNAs, num_tRNAs, num_rRNAs, feature_len, offset;
        long seq_len = -1;
        String sequence_id, current_sequence_id=null, attribute;
        Node seq_node = null, gene_node, rna_node, feature_node, parent_node = null, node;
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
                        line = in.readLine();
                        if (line.equals("") || line.charAt(0) == '#') // if line is empty or a comment skip it
                            continue;
                        if (line.contains("\t")){
                            fields = line.split("\\t");
                            attribute = fields[fields.length - 1];
                        } else {
                            //log_file.write("This line is not tab-delimited and is split by whitespaces: \n" + line + "\n"); 
                            fields = line.split("\\s+");
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
                            if (address[1] == -1){// if sequence is not uniquely determined
                                log_file.write("Sequence ID = "+sequence_id+" missed in genome "+address[0]+"\n"); // usually organal genes
                                continue;
                            } else {
                                seq_node = graphDb.findNode(sequence_label, "identifier", address[0] + "_" + address[1]);
                                seq_len = genomeDb.sequence_length[address[0]][address[1]];
                            }
                            if(address[2] > seq_len) {
                                log_file.write("Position "+address[2] + " is out of range 1-"+seq_len+".\n");
                                continue;
                            }
                            if(address[3] > seq_len) {
                                log_file.write("Position "+address[3] + " is out of range 1-"+seq_len+".\n");
                                continue;
                            }
                            current_sequence_id = sequence_id;
                        }
                        feature_node = graphDb.createNode(feature_label);
                        feature_node.setProperty("address", address);
                        feature_node.setProperty("strand", strand);
                        feature_node.setProperty("length", feature_len);
                        ID = get_property(attribute,"ID");
                        feature_node.setProperty("id", ID);
                        feature_node.setProperty("attribute", attribute);
                        feature_node.setProperty("type", fields[2]);
                        feature_node.setProperty("name", get_property(attribute,"Name"));
                        feature_node.setProperty("annotation_id", annotation_id);
                        feature_node.setProperty("genome",address[0]);
                        parent_node = get_node_by_id(gene_nodes,get_property(attribute,"Parent"));
                        if (parent_node != null)
                            parent_node.createRelationshipTo(feature_node, RelTypes.is_parent_of);
                        if (CONNECT_ANNOTATIONS){
                            start_ptr = locate(address);
                            offset = start_ptr.offset;
                            node = graphDb.getNodeById(start_ptr.node_id);
                            rel = feature_node.createRelationshipTo(node, RelTypes.starts);
                            rel.setProperty("offset", offset);
                            rel.setProperty("genomic_position", address[2]);
                            rel.setProperty("forward", start_ptr.canonical);
                            stop_ptr = locate(new int[]{address[0], address[1], address[3]});
                            rel = feature_node.createRelationshipTo(graphDb.getNodeById(stop_ptr.node_id), RelTypes.stops);
                            rel.setProperty("offset", stop_ptr.offset);
                            rel.setProperty("genomic_position", address[2] + feature_len - 1);
                            rel.setProperty("forward", stop_ptr.canonical);
                        }
                        if (fields[2].endsWith("gene")) {
                            // create new gene node
                            feature_node.addLabel(gene_label);
                            seq_node.createRelationshipTo(feature_node, RelTypes.has);
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
                                                    rna_node.setProperty("id", gene_node.getProperty("id"));
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
                        if (trsc % 987 == 1)
                            System.out.print("\r" + address[0] + "\t" + num_genes + "\t" + num_mRNAs + "\t" + num_tRNAs + "\t" + num_rRNAs + "\t");
                    }// for trsc
                    tx2.success();
                } // tx2
            } // while lines
            in.close();
            System.out.println("\r" + address[0] + "\t" + num_genes + "\t" + num_mRNAs + "\t" + num_tRNAs + "\t" + num_rRNAs + "\t");
            set_protein_sequences(gene_nodes, address[0], log_file, PATH_TO_THE_PANGENOME_DATABASE);
            log_file.write("Genome "+address[0] + " : " + num_genes + " genes\t" + num_mRNAs + " mRNAs\t" + num_tRNAs + " tRNAs\t" + num_rRNAs + " rRNAs\n");
            log_file.write("----------------------------------------------------\n");
        }catch (IOException ioe) {
            System.out.println(ioe.getMessage());
            System.out.println("Could not open " + gff_path + "!");
        }
    }

    private void parse_gbk(int[] address, String annotation_id, BufferedWriter log_file, String gbk_path, String proteins_path){
        int i, trsc, num_genes, num_mRNAs, num_tRNAs, num_rRNAs, feature_len, offset;
        String sequence_id;
        Node seq_node = null, gene_node = null, feature_node = null;
        String[] fields, coordinates;
        String strand, line, protein_id = "";
        StringBuilder protein = new StringBuilder();
        IndexPointer start_ptr, stop_ptr;
        boolean complement;
        num_genes = num_mRNAs = num_tRNAs = num_rRNAs = 0;
        try (BufferedReader in = new BufferedReader(new FileReader(gbk_path))) {
             BufferedWriter out = new BufferedWriter(new FileWriter(proteins_path + "/proteins_" + address[0] + ".fasta"));
            // for each record of gff file
            while (in.ready()) {
                try (Transaction tx2 = graphDb.beginTx()) {
                    for (trsc = 0; trsc < MAX_TRANSACTION_SIZE && in.ready(); ++trsc) {
                        line = in.readLine();
                        if (line.equals("") || line.charAt(0) == '#') // if line is empty or a comment skip it
                            continue;
                        if (line.startsWith("LOCUS")){
                            fields = line.split("\\s+");
                            sequence_id = fields[1];
                            address[1] = find_sequence(sequence_id, address[0]);
                            if (address[1] == -1){
                                System.err.println("Failed to find Sequence ID = "+sequence_id+" in genome "+address[0]+"\n"); // usually organal genes
                                 System.exit(1);
                            } else
                                seq_node = graphDb.findNode(sequence_label, "number", address[1]);
                        } else {
                            if(!line.startsWith("                     ")){
                                if(line.startsWith("     gene") || line.startsWith("     mRNA") || line.startsWith("     tRNA") || line.startsWith("     rRNA")){
                                    complement = line.contains("complement");
                                    fields = line.replaceAll("[^0-9]+", " ").trim().split("\\s");
                                    address[2] = Integer.parseInt(fields[0]);
                                    address[3] = Integer.parseInt(fields[fields.length - 1]);
                                    strand = complement ?  "-" : "+";
                                    feature_len = address[3] - address[2] + 1;
                                    feature_node = graphDb.createNode(feature_label);
                                    seq_node.createRelationshipTo(feature_node, RelTypes.has);
                                    feature_node.setProperty("address", address);
                                    feature_node.setProperty("strand", strand);
                                    feature_node.setProperty("length", feature_len);
                                    feature_node.setProperty("annotation_id", annotation_id);
                                    feature_node.setProperty("genome",address[0]);
                                    if(line.startsWith("     gene")){
                                        ++num_genes;
                                        gene_node = feature_node;
                                        gene_node.addLabel(gene_label);
                                        feature_node.setProperty("type", "gene");
                                    } else if(line.startsWith("     mRNA")){
                                        ++num_mRNAs;
                                        feature_node.addLabel(mRNA_label);
                                        gene_node.createRelationshipTo(feature_node, RelTypes.codes_for);
                                        feature_node.setProperty("type", "mRNA");
                                    } else if(line.startsWith("     tRNA")){
                                        ++num_tRNAs;
                                        feature_node.addLabel(tRNA_label);
                                        feature_node.setProperty("type", "tRNA");
                                    } else if(line.startsWith("     rRNA")){
                                        ++num_rRNAs;
                                        feature_node.addLabel(rRNA_label);
                                        feature_node.setProperty("type", "rRNA");
                                    } 
                                } 
                            } else {
                                    if(line.startsWith("                     /gene"))
                                        feature_node.setProperty("id", line.substring(27));
                                    else if(line.startsWith("                     /locus_tag"))
                                        feature_node.setProperty("name", line.substring(32));
                                    else if(line.startsWith("                     /protein_id"))
                                        feature_node.setProperty("protein_ID", protein_id = line.substring(33));
                                    else if(line.startsWith("                     /product"))
                                        feature_node.setProperty("product", line.substring(30));
                                    else if(line.startsWith("                     /note"))
                                        feature_node.setProperty("note", line.substring(27));
                                    else if(line.startsWith("                     /translation")){
                                            protein.setLength(0);
                                            protein.append(line.split("/translation=")[1].replaceAll("\"", ""));
                                            while (in.ready() && !line.endsWith("\"")){
                                                line = in.readLine().trim();
                                                if (line.endsWith("\"")){
                                                    protein.append(line.replaceAll("\"", ""));
                                                    break;
                                                } else
                                                    protein.append(line);
                                            }
                                            //System.out.println(protein);
                                            feature_node.setProperty("protein", protein.toString());
                                            feature_node.setProperty("protein_length", protein.length());
                                            out.write(">" + protein_id + "\n" + protein.toString() + "\n");
                                    }
                            }
                        }
                        if (trsc % 987 == 1)
                            System.out.print("\r" + address[0] + "\t" + num_genes + "\t" + num_mRNAs + "\t" + num_tRNAs + "\t" + num_rRNAs + "\t");
                    }// for trsc
                    tx2.success();
                } // tx2
            } // while lines
            in.close();
            out.close();
            System.out.println("\r" + address[0] + "\t" + num_genes + "\t" + num_mRNAs + "\t" + num_tRNAs + "\t" + num_rRNAs + "\t");
            log_file.write("Genome "+address[0] + " : " + num_genes + " genes\t" + num_mRNAs + " mRNAs\t" + num_tRNAs + " tRNAs\t" + num_rRNAs + " rRNAs\n");
            log_file.write("----------------------------------------------------\n");
        }catch (IOException ioe) {
            System.out.println(ioe.getMessage());
            System.out.println("Could not open " + gbk_path + "!");
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
            if (f.ID.equals(id)){//.toLowerCase()
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
    private void set_protein_sequences(LinkedList<feature> gene_nodes, int genome, BufferedWriter log_file, String pangenome_path){
        IntPairComparator comp = new IntPairComparator();
        PriorityQueue<int[]> pq = new PriorityQueue(comp);
        Node cds_node, mrna_node, gene_node;
        Relationship start_edge;
        //IndexPointer start_ptr = new IndexPointer();
        int trsc, isoforms_num, gene_start_pos;
        boolean forward;
        feature f;
        String annotaion_id, protein, mRNA_id, gene_id;
        int[] address, begin_end;
        StringBuilder gene_builder = new StringBuilder();
        StringBuilder rna_builder = new StringBuilder();
        try (BufferedWriter out = new BufferedWriter(new FileWriter(pangenome_path + "/proteins/proteins_" + genome + ".fasta"))) {
            //System.out.print("Adding protein sequences...");
            while (!gene_nodes.isEmpty()){
                try (Transaction tx = graphDb.beginTx()) {
                    for (trsc = 0; !gene_nodes.isEmpty()&& trsc < 2 * MAX_TRANSACTION_SIZE; ++trsc){
                        f = gene_nodes.remove();
                        gene_node = f.node;
                        gene_id = f.ID;
                        address = (int[])gene_node.getProperty("address");
                        forward = gene_node.getProperty("strand").equals("+");
                        gene_start_pos = address[2];
                        /*extract gene sequence as appears in the sequence and connects it to its nodes
                        gene_builder.setLength(0);
                        start_edge = gene_node.getSingleRelationship(RelTypes.starts, Direction.OUTGOING);
                        start_ptr.node_id = start_edge.getEndNode().getId();
                        start_ptr.offset = (int)start_edge.getProperty("offset");
                        start_ptr.canonical = (boolean)start_edge.getProperty("forward");
                        extract_sequence(gene_builder, start_ptr, address, K);*/
                        address[2] -= 1;
                        address[3] -= 1;
                        gene_builder.setLength(0);
                        genome_scanner.get_sub_sequence(gene_builder, address, true);//(boolean)start_edge.getProperty("forward"));
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
                                        log_file.write("Gene ID = " + gene_id + " miss-annotated!\n");
                                        break;
                                    }
                                }
                                if (!forward)
                                    reverse_complement(rna_builder);
                                protein = translate(rna_builder);
                                ++num_proteins;
                                if (protein.length() > 0){
                                    //if (num_proteins % 11 == 1)
                                    //    System.out.print("\rAdding protein sequences... " + num_proteins);
                                    mRNA_id = (String)mrna_node.getProperty("id");
                                    annotaion_id = (String)mrna_node.getProperty("annotation_id");
                                    out.write(">" + mRNA_id + "\n"); // " " + (String)mrna_node.getProperty("attribute") +
                                    write_fasta(out, protein, 70);
                                    mrna_node.setProperty("sequence", rna_builder.toString());
                                    mrna_node.setProperty("protein", protein.toString());
                                    mrna_node.setProperty("protein_ID", mRNA_id);
                                    mrna_node.setProperty("protein_length", protein.length());
                                    if (protein.charAt(0) != 'M'  || protein.charAt(protein.length() - 1) != '*'){
                                        log_file.write("Protein ID = " + mRNA_id + " miss-annotated!\n");
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
            //System.out.print("\r                                                        \r"); //clear the line
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
        for (int i = 0; i < fields.length; ++i)
            if (fields[i].startsWith(property)){//.toLowerCase()
                if (fields[i].contains("="))
                    return fields[i].trim().split("=")[1];//.toLowerCase();
                else
                    return fields[i].trim().split("\\s+")[1];//.toLowerCase();
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
    public void retrieve_feature() {
        if (PATH_TO_THE_GENOME_NUMBERS_FILE == null){
            System.out.println("PATH_TO_THE_GENOME_NUMBERS_FILE is empty.");
            System.exit(1);
        }  
        if (! new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH).exists()) {
            System.out.println("No database found in " + PATH_TO_THE_PANGENOME_DATABASE);
            System.exit(1);
        }
        BufferedReader in;
        ResourceIterator<Node> feature_nodes;
        //Relationship rstart;
        Node feature;//start, 
        String ID;
        String line;
        boolean forward;
        int i, j, num_genomes, begin, end, genome = 0;
        int[] address, genome_numbers;
        StringBuilder feature_seq;
        Label feature_label = labels.get(FEATURE);
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        startTime = System.currentTimeMillis();
        genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
        indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, "sorted");
        K_SIZE = indexDb.get_K();
        genome_scanner = new SequenceScanner(genomeDb, 1, 1, K_SIZE, indexDb.get_pre_len());
        num_genomes = 0;
        feature_seq = new StringBuilder();
        BufferedWriter[] out = new BufferedWriter[genomeDb.num_genomes + 1];
        try {
            in = new BufferedReader(new FileReader(PATH_TO_THE_GENOME_NUMBERS_FILE));
            while (in.ready()) {
                line = in.readLine();
                if (line.equals("")) {
                    continue;
                }
                try{
                    genome = Integer.parseInt(line.trim());
                }catch(NumberFormatException e){
                    System.out.println(genome + "is not a valid genome number.");
                    continue;
                }
                num_genomes++;
            }
            in.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
        genome_numbers = new int[num_genomes];
        try (Transaction tx = graphDb.beginTx()) {
            try {
                in = new BufferedReader(new FileReader(PATH_TO_THE_GENOME_NUMBERS_FILE));
                // fill all the records in an array to be sorted    
                for (i = 0; in.ready();) {
                    line = in.readLine();
                    if (line.equals("")) {
                        continue;
                    }
                    try{
                        genome = genome_numbers[i] = Integer.parseInt(line.trim());
                    }catch(NumberFormatException e){
                        System.out.println(genome + " skipped.");
                        continue;
                    }
                    out[genome] = new BufferedWriter(new FileWriter(FEATURE + "s." + genome + ".fasta"));
                    ++i;
                }
                in.close();
                Arrays.sort(genome_numbers);
            } catch (IOException ioe) {
                System.out.println("Failed to read file names!");
                System.exit(1);
            }
            try {
                // for all the genes in the database    
                for (i = 0, j = 0, feature_nodes = graphDb.findNodes(feature_label); feature_nodes.hasNext();){
                    feature = feature_nodes.next();
                    genome = (int) feature.getProperty("genome");
                    ID = (String) feature.getProperty("id");
                    feature_seq.setLength(0);
                    if (Arrays.binarySearch(genome_numbers, genome) >= 0){
                        ++i;
                        if (FEATURE == "mRNA")
                            feature_seq.append(feature.getProperty("sequence"));
                        else{
                            //rstart = feature.getSingleRelationship(RelTypes.starts, Direction.OUTGOING);
                           // start = rstart.getEndNode();
                            address = (int[]) feature.getProperty("address");
                            begin = address[2];
                            end = address[3];
                            forward = feature.getProperty("strand").equals("+");
                            //extract_sequence(gene_seq, new IndexPointer(start.getId(), (boolean) rstart.getProperty("forward"), (int) rstart.getProperty("offset"),-1l), address, K);//
                            address[2] -= 1;
                            address[3] -= 1;
                            genome_scanner.get_sub_sequence(feature_seq, address, forward);
                            //genomeDb=new sequence_database(pangenome_path+GENOME_DATABASE_PATH);
                            //if(gene_seq.toString().equals(genomeDb.get_sequence(genome, sequence, begin-1, end-begin+1, strand))
                            //|| gene_seq.toString().equals(genomeDb.get_sequence(genome, sequence, begin-1, end-begin+1, !strand)) )//gene_seq.length() == end-begin+1)//
                        }
                        ++j;
                        out[genome].write(">" + ID + "\n");
                        //if (strand) {
                            write_fasta(out[genome], feature_seq.toString(), 70);
                        //} else {
                        //    reverse_complement(gene_seq);
                        //    write_fasta(out, gene_seq.toString(), 70);
                        //}
                    }
                    //if (i % (num_genes / 100 + 1) == 0) {
                    //    System.out.print((long) i * 100 / num_genes + 1 + "%\r");
                    //}
                }//for i
                System.out.println(j + " out of " + i + " " + FEATURE + "s found and retrieved successfully.");
                for (i = 0; i < genome_numbers.length; ++i)
                   out[genome_numbers[i]].close();
            } catch (IOException ioe) {
                System.out.println("Failed to read file names!");
                System.exit(1);
            }
            tx.success();
        }
        graphDb.shutdown();
        genomeDb.close();
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
        for(i = 0; i <= mRNA.length() - 3; i += 3)
            protein_builder.append(aminoacid_table[binary[mRNA.charAt(i)] * 16 
                                                 + binary[mRNA.charAt(i + 1)] * 4
                                                 + binary[mRNA.charAt(i + 2)]]);
        return protein_builder.toString();
    }    

    public void remove_annotaions() {
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
