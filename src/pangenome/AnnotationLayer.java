/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import alignment.Alignment;
import genome.SequenceDatabase;
import index.IndexPointer;
import protein.protein_builder;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.PriorityQueue;
import org.neo4j.graphdb.Direction;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import static pangenome.SequenceLayer.append_fwd;
import static pangenome.SequenceLayer.append_rev;

import static pantools.Pantools.GENOME_DATABASE_PATH;
import static pantools.Pantools.GRAPH_DATABASE_PATH;
import static pantools.Pantools.K;
import static pantools.Pantools.PATH;
import static pantools.Pantools.RelTypes;
import static pantools.Pantools.gene_label;
import static pantools.Pantools.genomeDb;
import static pantools.Pantools.graphDb;
import static pantools.Pantools.homolog_lable;
import static pantools.Pantools.mRNA_label;
import static pantools.Pantools.ncRNA_label;
import static pantools.Pantools.num_edges;
import static pantools.Pantools.num_nodes;
import static pantools.Pantools.ortholog_lable;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.pgRNA_label;
import static pantools.Pantools.sequence_label;
import static pantools.Pantools.startTime;
import static pantools.Pantools.tRNA_label;
import static pangenome.SequenceLayer.locate;
import static pangenome.SequenceLayer.write_fasta;
import static pangenome.SequenceLayer.get_outgoing_edge;
import static pantools.Pantools.db_trsc_limit;

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
    public class PairComparator implements Comparator<int[]> {
        @Override
        public int compare(int[] x, int[] y) {
            if (x[0] < y[0]) 
                return -1;
            else if (x[0] > y[0]) 
                return 1;
            else if (x[1] < y[1]) 
                return -1;
            else if (x[1] > y[1]) 
                return 1;
            else
                return 0;
        }
    }

    /**
     * Adds nodes of different genomic features to the pangenome.
     * All features should be annotated in the hierarchical order of gene, RNA, exon, CDS  
     * 
     * @param gff_paths_file 
     */
    public void annotate(String gff_paths_file) {
        int i, j, num_genes=0, num_mRNAs, num_tRNAs, num_ncRNAs, num_pgRNAs, num_exons, protein_num;
        int begin, end, rna_len, isoforms_num, gene_start_pos = 0, phase;
        boolean forward=false;
        num_exons = 0;
        String sequence_id, current_sequence_id=null, origin = null;
        Node db_node, gene_node = null, rna_node = null;
        protein_builder pb = new protein_builder();
        StringBuilder log = new StringBuilder();
        StringBuilder coding_RNA = new StringBuilder();
        String[] fields;
        StringBuilder gene_builder = new StringBuilder();
        String strand = null, gff_name, line=null, protein;
        List<Long>[][] genes_list;
        int[] address = new int[3];
        long[] genes_array;
        int[] pair;
        PairComparator comp = new PairComparator();
        PriorityQueue<int[]> pq = new PriorityQueue(comp);
        if (new File(PATH + GRAPH_DATABASE_PATH).exists()) {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH + GRAPH_DATABASE_PATH))
                    .setConfig("keep_logical_logs", "100M size").newGraphDatabase();
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
            try (BufferedReader gff_files = new BufferedReader(new FileReader(gff_paths_file))) {
                System.out.println("genome\tgenes\tmRNAs\ttRNAs\tncRNAs\tpgRNAs");
                for (address[0] = 1; gff_files.ready() && address[0] <= genomeDb.num_genomes; ++address[0]) // for each gff file
                {
                    protein_num = 0;
                    num_genes = num_mRNAs = num_tRNAs = num_ncRNAs = num_pgRNAs = 0;
                    rna_len = isoforms_num = num_exons = 0;
                    gff_name = gff_files.readLine();
                    if (gff_name.equals(""))
                        continue;
                    BufferedReader in = new BufferedReader(new FileReader(gff_name));
                    BufferedWriter out = new BufferedWriter(new FileWriter(gff_name + ".proteins.fasta"));
                // for each record of gff file
                    while (in.ready()) 
                    {
                        try (Transaction tx2 = graphDb.beginTx()) {
                            for (i = 0; i < db_trsc_limit && in.ready(); ++i) {
                                line = in.readLine();
                                if (line.equals("") || line.charAt(0) == '#') // if line is empty or a comment skip it
                                {
                                    continue;
                                }
                                fields = line.split("\\t");
                                if (fields.length < 8) {
                                    continue;
                                }
                                begin = Integer.parseInt(fields[3]) - 1;
                                end = Integer.parseInt(fields[4]) - 1;
                                strand = fields[6];
                                phase = !".".equals(fields[7]) ? Integer.parseInt(fields[7]) : 0 ;
                                forward = strand.equals("+");
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
                                    origin = address[0] + "_" + address[1];
                                    switch (fields[2]) {
                                        case "CDS":
                                        case "five_prime_UTR":
                                        case "three_prime_UTR":
                                        // if the RNA codes for a protein (is not a pseudogene)    
                                            if (fields[2].equals("CDS") && rna_node.hasLabel(mRNA_label)) {
                                                rna_len += end - begin + 1 - phase;
                                            // For forward strand features, phase is counted from the start field. 
                                            // For reverse strand features, phase is counted from the end field.
                                                if ( forward )
                                                    pq.add(new int[]{begin + 1 + phase, end + 1});
                                                else
                                                    pq.add(new int[]{begin + 1, end + 1 - phase});
                                            }
                                            if (!rna_node.hasProperty(fields[2])) 
                                                rna_node.setProperty(fields[2], fields[3] + "_" + fields[4]);
                                            else
                                                rna_node.setProperty(fields[2], (String) rna_node.getProperty(fields[2]) + " " + fields[3] + "_" + fields[4]);
                                            break;
                                        case "exon":
                                        case "pseudogenic_exon":
                                            ++num_exons;
                                            break;
                                        case "mRNA":
                                        case "tRNA":
                                        case "ncRNA":
                                        case "pseudogenic_transcript":
                                        // if is not the first RNA of the gene
                                            if (num_exons != 0) 
                                            {
                                            // if is a coding RNA
                                                if (rna_len != 0) 
                                                {
                                                    while (pq.size() != 0) {
                                                        pair = pq.remove();
                                                        coding_RNA.append(gene_builder.substring(pair[0] - gene_start_pos, pair[1] - gene_start_pos + 1));
                                                    }
                                                    protein = pb.translate(coding_RNA, forward);
                                                    if (protein.length() != 0) {
                                                        ++protein_num;
                                                        out.write(">" + origin + "_" + protein_num + "\n");
                                                        write_fasta(out, protein, 70);
                                                        rna_node.setProperty("protein", protein);
                                                        rna_node.setProperty("protein_number", origin + "_" + protein_num);
                                                        rna_node.setProperty("coding_length", rna_len);
                                                    }
                                                    coding_RNA.setLength(0);
                                                    rna_len = 0;
                                                }
                                                rna_node.setProperty("num_exons", num_exons);
                                                num_exons = 0;
                                            }
                                            if (fields[2].charAt(0) == 'm') // is a mRNA
                                            {
                                                ++num_mRNAs;
                                                rna_node = graphDb.createNode(mRNA_label);
                                                isoforms_num++;
                                            } else if (fields[2].charAt(0) == 't') // is a tRNA
                                            {
                                                ++num_tRNAs;
                                                rna_node = graphDb.createNode(tRNA_label);
                                            } else if (fields[2].charAt(0) == 'n')// is a ncRNA
                                            {
                                                ++num_ncRNAs;
                                                rna_node = graphDb.createNode(ncRNA_label);
                                            } else // is a pseudogenic_transcript
                                            {
                                                ++num_pgRNAs;
                                                rna_node = graphDb.createNode(pgRNA_label);
                                            }
                                            gene_node.createRelationshipTo(rna_node, RelTypes.has);
                                            rna_node.setProperty("origin", origin);
                                            rna_node.setProperty("begin", begin + 1);
                                            rna_node.setProperty("end", end + 1);
                                            rna_node.setProperty("strand", strand);
                                            rna_node.setProperty("attribute", fields[fields.length - 1]);
                                            break;
                                        case "gene":
                                        case "pseudogene":
                                        case "transposable_element_gene":
                                        // is not the first RNA of the gene
                                            if (num_exons != 0) 
                                            {
                                            // if is a coding RNA
                                                if (rna_len != 0) 
                                                {
                                                    // concatenates CDS in increasing order of their start point using a priority queue
                                                    while (pq.size() != 0) {
                                                        pair = pq.remove();
                                                        coding_RNA.append(gene_builder.substring(pair[0] - gene_start_pos, pair[1] - gene_start_pos + 1));
                                                    }
                                                    protein = pb.translate(coding_RNA, forward);
                                                    if (protein.length() != 0) {
                                                        ++protein_num;
                                                        out.write(">" + origin + "_" + protein_num + "\n");
                                                        write_fasta(out, protein, 70);
                                                        rna_node.setProperty("protein", protein);
                                                        rna_node.setProperty("protein_number", origin + "_" + protein_num);
                                                        rna_node.setProperty("coding_length", rna_len);
                                                    }
                                                    else 
                                                        System.out.println(line);
                                                    coding_RNA.setLength(0);
                                                    rna_len = 0;
                                                }
                                                rna_node.setProperty("num_exons", num_exons);
                                                num_exons = 0;
                                            }
                                            if (isoforms_num != 0) // is not the first gene
                                            {
                                                gene_node.setProperty("isoforms_num", isoforms_num);
                                                isoforms_num = 0;
                                            }
                                            ++num_genes;
                                            gene_start_pos = begin + 1;
                                            // initial a gene_node
                                            gene_node = graphDb.createNode(gene_label);
                                            gene_node.setProperty("type", fields[2]);
                                            gene_node.setProperty("genome", address[0]);
                                            gene_node.setProperty("origin", origin);
                                            gene_node.setProperty("begin", begin + 1);
                                            gene_node.setProperty("end", end + 1);
                                            gene_node.setProperty("length", end - begin + 1);
                                            gene_node.setProperty("strand", strand);
                                            gene_node.setProperty("attribute", fields[fields.length - 1]); // is needed in retrieve genes function
                                        // extract gene sequence as appears in the sequence and connects it to its nodes   
                                            connect_gene_to_nodes(gene_builder, address, begin, end,gene_node);
                                            gene_node.setProperty("sequence", gene_builder.toString());
                                        // adding gene_node id to the sequence node
                                            genes_list[address[0]][address[1]].add(gene_node.getId());
                                            break;
                                    } // switch
                                } else // if sequence not found
                                {
                                    log.append(sequence_id).append(" missed in genome ").append(address[0]).append("\n"); // usually organal genes
                                }
                                System.out.print("\r" + address[0] + "\t" + num_genes + "\t" + num_mRNAs + "\t" + num_tRNAs + "\t" + num_ncRNAs + "\t" + num_pgRNAs);
                            }// for trsc
                            tx2.success();
                        } // tx2
                    } // while lines
                    try (Transaction tx3 = graphDb.beginTx()) {
                    // for the last RNA: if is not the first RNA of the gene    
                        if (num_exons != 0) 
                        {
                        // if is a coding RNA    
                            if (rna_len != 0) 
                            {
                                while (pq.size() != 0) {
                                    pair = pq.remove();
                                    coding_RNA.append(gene_builder.substring(pair[0] - gene_start_pos, pair[1] - gene_start_pos + 1));
                                }
                                protein = pb.translate(coding_RNA, forward);
                                if (protein.length() != 0) {
                                    ++protein_num;
                                    out.write(">" + origin + " _" + protein_num + "\n");
                                    write_fasta(out, protein, 70);
                                    rna_node.setProperty("protein", protein);
                                    rna_node.setProperty("protein_number", origin + "_" + protein_num);
                                    rna_node.setProperty("coding_length", rna_len);
                                }
                                else 
                                    System.out.println(line);
                                coding_RNA.setLength(0);
                                rna_len = 0;
                            }
                            rna_node.setProperty("num_exons", num_exons);
                            num_exons = 0;
                        }
                        if (isoforms_num != 0) 
                        {
                            gene_node.setProperty("isoforms_num", isoforms_num);
                            isoforms_num = 0;
                        }
                        in.close();
                        out.close();
                        System.out.println("\r" + address[0] + "\t" + num_genes + "\t" + num_mRNAs + "\t" + num_tRNAs + "\t" + num_ncRNAs + "\t" + num_pgRNAs);
                        tx3.success();
                    }
                } // for genomes
                gff_files.close();
            } catch (IOException ioe) {
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
                tx4.success();
            }
            try {
                BufferedWriter out = new BufferedWriter(new FileWriter(PATH + "/annotation.log"));
                out.write(log.toString());
                out.close();
            } catch (IOException ioe) {
            }
            graphDb.shutdown();
            genomeDb.close();
            // delete the database transaction files
            File directory = new File(PATH + GRAPH_DATABASE_PATH);
            for (File f : directory.listFiles())
                if (f.getName().startsWith("neostore.transaction.db."))
                    f.delete();
        } else {
            System.out.println("No database found in " + PATH);
            System.exit(1);
        }
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
     * @param address Contains the genome and the sequence number of the gene
     * @param begin Start position of the gene
     * @param end Stop position of the gene
     * @param gene_node The gene node itself
     */
    public void connect_gene_to_nodes(StringBuilder gene_builder, int[] address, int begin, int end, Node gene_node) {
        Relationship rel;
        int loc, node_len, neighbor_len, seq_len = end - begin + 1, position;
        String rel_name;
        gene_builder.setLength(0);
        loc = begin;
        IndexPointer start_ptr = locate(address[0], address[1], begin);
        position = start_ptr.position;
        Node neighbor, node = graphDb.getNodeById(start_ptr.node_id);
        node_len = (int) node.getProperty("length");
        rel = gene_node.createRelationshipTo(node, RelTypes.visits);
        rel.setProperty("forward", start_ptr.canonical);
        rel.setProperty("position", position);
        gene_node.setProperty("start_node_id", node.getId());
        gene_node.setProperty("start_edge_id", rel.getId());
        if (start_ptr.canonical) {
            if (position + seq_len - 1 <= node_len - 1) {
                loc += append_fwd(gene_builder, (String) node.getProperty("sequence"), position, position + seq_len - 1);
            } else {
                loc += append_fwd(gene_builder, (String) node.getProperty("sequence"), position, node_len - 1);
            }
        } else {
            if (position - (seq_len - 1) >= 0) {
                loc += append_rev(gene_builder, (String) node.getProperty("sequence"), position - (seq_len - 1), position);
            } else {
                loc += append_rev(gene_builder, (String) node.getProperty("sequence"), 0, position);
            }
        }
        while (gene_builder.length() < seq_len ) {
            address[2] = loc - K + 1;
            rel = get_outgoing_edge(node, address);
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
            node_len = (int) node.getProperty("length");
            gene_node.createRelationshipTo(node, RelTypes.visits);
        } // while
        gene_node.setProperty("stop_node", node.getId());
    }

    /**
     * Groups the highly similar genes together
     */
    public void denovo_homology_annotation() {
        int i, step, gene_len=0, current_gene_len, num_groups=0, total_genes=0, group_size;
        int[] copy_number;
        String current_gene_seq, gene_seq;
        Node group_node, gene_node, node, current_gene_node;
        long gene_id, current_gene_id;
        ResourceIterator<Node> genes_iterator;
        if (new File(PATH + GRAPH_DATABASE_PATH).exists()) {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH + GRAPH_DATABASE_PATH))
                    .setConfig("keep_logical_logs", "100M size").newGraphDatabase();
            registerShutdownHook(graphDb);
            startTime = System.currentTimeMillis();
            genomeDb = new SequenceDatabase(PATH + GENOME_DATABASE_PATH);
            LinkedList<Node> homolog_gene_nodes = new LinkedList(); 
            PriorityQueue<Long> pq = new PriorityQueue();
            copy_number = new int[genomeDb.num_genomes+1];
            Alignment aligner = new Alignment();
            System.out.println("genes\tgroups");
            try (Transaction tx1 = graphDb.beginTx()) {
                genes_iterator = graphDb.findNodes(gene_label);
                while (genes_iterator.hasNext()) {
                try (Transaction tx2 = graphDb.beginTx()) {
                    for (step = 0; step < db_trsc_limit && genes_iterator.hasNext(); ++step) {
                        gene_node=genes_iterator.next();
                        gene_len=(int)gene_node.getProperty("length");
                        gene_seq=(String)gene_node.getProperty("sequence");
                        if (!gene_node.hasRelationship(RelTypes.contains, Direction.INCOMING)) { // To avoid having one gene in different groups
                            for (Relationship r1: gene_node.getRelationships(Direction.OUTGOING,RelTypes.visits)) {
                                node = r1.getEndNode();
                                for (Relationship r2: node.getRelationships(Direction.INCOMING,RelTypes.visits)) {
                                    if(r1 != r2) // avoid comparisons to the gene itself  
                                        pq.offer(r2.getStartNode().getId());
                                }
                            }
                            current_gene_id = -1l; // To be different from the first in the proprity queue
                            while (!pq.isEmpty()) { // for all the candidates with some shared node with the gene
                                gene_id = pq.remove();
                                if(gene_id != current_gene_id) {
                                    current_gene_id = gene_id;
                                    current_gene_node = graphDb.getNodeById(gene_id);
                                    if ( ! current_gene_node.hasRelationship(RelTypes.contains, Direction.INCOMING) // To avoid having one gene in different groups
                                            && ! have_overlap(gene_node,current_gene_node )
                                            ) {
                                        current_gene_len = (int)current_gene_node.getProperty("length");
                                    // If gene and the current mate are almost equally long 
                                        if(Math.abs(current_gene_len - gene_len) <= Math.max(current_gene_len, gene_len) / 10 ) {
                                            current_gene_seq=(String)current_gene_node.getProperty("sequence");
                                        // If genes are highly similar 
                                            if ( aligner.get_similarity(current_gene_seq, gene_seq) > 0.75 ) {
                                                homolog_gene_nodes.add(current_gene_node);
                                                copy_number[(int)current_gene_node.getProperty("genome")] += 1;
                                            }
                                        }
                                    }
                                }
                            }// while
                            group_size = homolog_gene_nodes.size()+1;
                            if( group_size > 1 ) {
                                group_node = graphDb.createNode(homolog_lable);
                                ++num_groups;
                                total_genes += group_size;
                                group_node.setProperty("num_members", group_size );
                            // Because the gene has not been added to the group itself    
                                group_node.createRelationshipTo(gene_node, RelTypes.contains);
                                copy_number[(int)gene_node.getProperty("genome")] += 1;
                                while (!homolog_gene_nodes.isEmpty()) 
                                    group_node.createRelationshipTo(homolog_gene_nodes.removeFirst(), RelTypes.contains);
                                for (i = 1; i <= genomeDb.num_genomes; ++i) {
                                    group_node.setProperty("copy_number_in_genome_" + i, copy_number[i]);
                                    copy_number[i] = 0;
                                }
                            }
                            else
                                total_genes += 1;
                        }
                        if (num_groups % 10  == 0)
                            System.out.print("\r" + total_genes + "\t" + num_groups);
                        }// for
                    tx2.success();
                    }// transaction 2
                } // while 
                tx1.success();
                System.out.println("\r" + total_genes + "\t" + num_groups);
            } // transaction 1
        } else {
            System.out.println("pangenome database not found!");
            System.exit(0);
        }
    }

    /**
     * Decides if two genes overlap.
     * 
     * @param gene1 The first gene
     * @param gene2 The second gene
     * @return Decision
     */
    private boolean have_overlap(Node gene1, Node gene2){
        int gene1_begin, gene1_end, gene2_begin, gene2_end;
        String gene1_origin, gene2_origin;
        gene1_origin=(String)gene1.getProperty("origin");
        gene1_begin=(int)gene1.getProperty("begin");
        gene1_end=(int)gene1.getProperty("end");        
        gene2_origin=(String)gene2.getProperty("origin");
        gene2_begin=(int)gene2.getProperty("begin");
        gene2_end=(int)gene2.getProperty("end");
        if (!gene1_origin.equals(gene2_origin))
            return false;
        else if (gene1_begin > gene2_end || gene1_end < gene2_begin)
            return false;
        else if (gene2_begin > gene1_end || gene2_end < gene1_begin)
            return false;
        else
            return true;
    }
    
    /**
     * Groups the gene families found by a orthogroup finder tool.
     * 
     * @param group_file A text file containing inferred orthogroups in orthoMCL format 
     */
    public void group_ortholog_proteins(String group_file) {
        if (new File(PATH + GRAPH_DATABASE_PATH).exists()) {
            Node group_node = null, mRNA_node;
            Relationship rel;
            String line, origin, genome;
            String protein_number, group_name = null;
            String[] fields;
            ResourceIterator<Node> mRNAs;
            List<String> p_numbers = new LinkedList();
            List<Long> mRNA_ids = new LinkedList();
            ListIterator<String> p_numbers_iterator;
            int i, j, total_proteins = 0, num_proteins, num_groups = 0;
            boolean found;
            System.out.println("Grouping protein_coding genes ...");
            try (BufferedReader in = new BufferedReader(new FileReader(group_file))) {
                graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH + GRAPH_DATABASE_PATH))
                        .setConfig("keep_logical_logs", "100M size").newGraphDatabase();
                registerShutdownHook(graphDb);
                startTime = System.currentTimeMillis();
                try (Transaction tx = graphDb.beginTx()) {
                    mRNAs = graphDb.findNodes(mRNA_label);
                    while (mRNAs.hasNext()) {
                        mRNA_node = mRNAs.next();
                        if (mRNA_node.hasProperty("protein_number")) {
                            mRNA_ids.add(mRNA_node.getId());
                            p_numbers.add((String) mRNA_node.getProperty("protein_number"));
                        }
                    }
                    mRNAs.close();
                    System.out.println("proteins\tgroups");
                    while (in.ready()) {
                        line = in.readLine();
                        if (line.equals("")) {
                            continue;
                        }
                        fields = line.split("\\s");
                        num_proteins = fields.length - 1;
                        total_proteins += num_proteins;
                        group_name = fields[0].substring(0, fields[0].length() - 1);
                        if (num_proteins > 1) {
                            group_node = graphDb.createNode(ortholog_lable);
                            ++num_groups;
                            group_node.setProperty("num_proteins", num_proteins);
                            group_node.setProperty("name", group_name);
                            for (i = 1; i <= num_proteins; ++i) {
                                protein_number = fields[i];
                                p_numbers_iterator = p_numbers.listIterator();
                                for (found = false, j = 0; !found && p_numbers_iterator.hasNext(); ++j) {
                                    if (p_numbers_iterator.next().equals(protein_number)) {
                                        mRNA_node = graphDb.getNodeById(mRNA_ids.get(j));
                                        origin = (String) mRNA_node.getProperty("origin");
                                        genome = origin.split("_")[0];
                                        rel = group_node.createRelationshipTo(mRNA_node, RelTypes.has);
                                        rel.setProperty("origin", "in_" + origin);
                                        rel.setProperty("genome", "in_" + genome);
                                        if (group_node.hasProperty("num_" + origin)) {
                                            group_node.setProperty("num_" + origin, (int) group_node.getProperty("num_" + origin) + 1);
                                        } else {
                                            group_node.setProperty("num_" + origin, 1);
                                        }
                                        if (group_node.hasProperty("num_" + genome)) {
                                            group_node.setProperty("num_" + genome, (int) group_node.getProperty("num_" + genome) + 1);
                                        } else {
                                            group_node.setProperty("num_" + genome, 1);
                                        }
                                        found = true;
                                    }
                                }
                                if (!found) {
                                    System.out.println(protein_number + " not found!");
                                }
                            }
                        }
                        if (num_groups % 100  == 0)
                            System.out.print("\r" + total_proteins + "\t" + num_groups);
                    }
                    tx.success();
                    System.out.println("\r" + total_proteins + "\t" + num_groups);
                }
                in.close();
            } catch (IOException ioe) {
                System.out.println("Failed to open " + group_file);
                System.exit(1);
            }
            File directory = new File(PATH + GRAPH_DATABASE_PATH);
            for (File f : directory.listFiles()) {
                if (f.getName().startsWith("neostore.transaction.db.")) {
                    f.delete();
                }
            }
            graphDb.shutdown();
        } else {
            System.out.println("No database found in " + PATH + GRAPH_DATABASE_PATH);
            System.exit(1);
        }
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
