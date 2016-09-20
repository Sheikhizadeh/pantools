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
import static pangenome.SequenceLayer.extract_sequence;
import static pangenome.SequenceLayer.write_fasta;
import static pangenome.SequenceLayer.b_search;
import static alignment.Alignment.*;
import static pantools.Pantools.db_trsc_limit;

/**
 *
 * @author sheik005
 */
public class AnnotationLayer {
    /*
     To add gene_nodes to the graph. 
     gff_paths : a text file containing paths to annotation files 
     */

    public class PairComparator implements Comparator<int[]> {

        @Override
        public int compare(int[] x, int[] y) {
            // Assume neither string is null. Real code should
            // probably be more robust
            // You could also just return x.length() - y.length(),
            // which would be more efficient.
            if (x[0] < y[0]) 
                return -1;
            else if (x[0] > y[0]) 
                return 1;
            else
                return 0;
        }
    }

    public void annotate(String gff_paths) {
        int i, j, num_genes, num_mRNAs, num_tRNAs, num_ncRNAs, num_pgRNAs, num_exons, protein_num;
        int begin, end, s = 0, g_number, rna_len, isoforms_num, gene_start_pos = 0, phase;
        boolean forward=false;
        num_exons = 0;
        String sequence_id, current_sequence_id=null, origin = null;
        Node db_node, start_node, gene_node = null, rna_node = null;
        IndexPointer start_ptr;
        protein_builder pb = new protein_builder();
        Relationship rel;
        StringBuilder log = new StringBuilder();
        StringBuilder coding_RNA = new StringBuilder();
        String[] fields;
        StringBuilder gene = new StringBuilder();
        String strand = null, gff_name, line=null, protein;
        List<Long>[][] genes_list;
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
                // Read graph information    
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
            try (BufferedReader gff_files = new BufferedReader(new FileReader(gff_paths))) {
                System.out.println("genome\tgenes\tmRNAs\ttRNAs\tncRNAs\tpgRNAs");
                for (g_number = 1; gff_files.ready() && g_number <= genomeDb.num_genomes; ++g_number) // for each gff file
                {
                    protein_num = 0;
                    num_genes = num_mRNAs = num_tRNAs = num_ncRNAs = num_pgRNAs = 0;
                    rna_len = isoforms_num = num_exons = 0;
                    gff_name = gff_files.readLine();
                    if (gff_name.equals(""))
                        continue;
                    BufferedReader in = new BufferedReader(new FileReader(gff_name));
                    //BufferedWriter out = new BufferedWriter(new FileWriter(gff_name + ".proteins.fasta"));
                    while (in.ready()) // for each record of gff file
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
                                   s = find_sequence(sequence_id, g_number);
                                   current_sequence_id = sequence_id;
                                }
                                if (s > 0) // if sequence found
                                {
                                    origin = g_number + "_" + s;
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
                                            if (num_exons != 0) // is not the first RNA of the gene
                                            {
                                                if (rna_len != 0) // if is a coding RNA
                                                {
                                                    while (pq.size() != 0) {
                                                        pair = pq.remove();
                                                        coding_RNA.append(gene.substring(pair[0] - gene_start_pos, pair[1] - gene_start_pos + 1));
                                                    }
                                                    protein = pb.translate(coding_RNA, forward);
                                                    if (protein.length() != 0) {
                                                        ++protein_num;
                                                        //out.write(">" + origin + "_" + protein_num + "\n");
                                                        //write_fasta(out, protein, 70);
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
                                            if (num_exons != 0) // is not the first RNA of the gene
                                            {
                                                if (rna_len != 0) // if is a coding RNA
                                                {
                                                    // concatenates CDS in increasing order of their start point using a priority queue
                                                    while (pq.size() != 0) {
                                                        pair = pq.remove();
                                                        coding_RNA.append(gene.substring(pair[0] - gene_start_pos, pair[1] - gene_start_pos + 1));
                                                    }
                                                    protein = pb.translate(coding_RNA, forward);
                                                    if (protein.length() != 0) {
                                                        ++protein_num;
                                                        //out.write(">" + origin + "_" + protein_num + "\n");
                                                        //write_fasta(out, protein, 70);
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
                                            gene_node.setProperty("genome", g_number);
                                            gene_node.setProperty("origin", origin);
                                            gene_node.setProperty("begin", begin + 1);
                                            gene_node.setProperty("end", end + 1);
                                            gene_node.setProperty("length", end - begin + 1);
                                            gene_node.setProperty("strand", strand);
                                            gene_node.setProperty("attribute", fields[fields.length - 1]); // is needed in retrieve genes function
                                            start_ptr = locate(origin, begin);
                                            //stop_ptr =locate(origin,end-K+1); // points to the start position of the last k_mer of the gene
                                            start_node = graphDb.getNodeById(start_ptr.node_id);
                                            start_node.setProperty("gene_starts", true);
                                            rel = gene_node.createRelationshipTo(start_node, RelTypes.begin);
                                            rel.setProperty("forward", start_ptr.canonical);
                                            rel.setProperty("pos", start_ptr.position);
                                            // There is no need to point to the stop node.
                                            //stop_node=graphDb.getNodeById(stop_ptr.node_id);
                                            //stop_node.setProperty("gene_stops", true);
                                            //rel=gene_node.createRelationshipTo(stop_node, RelTypes.end);
                                            //rel.setProperty("forward", stop_ptr.canonical);
                                            //rel.setProperty("pos", stop_ptr.position);
                                            gene.setLength(0);
                                            // extract gene as appears in the sequence    
                                            extract_sequence(gene, start_ptr, origin, begin, end);
                                            /*if (gene.toString().contains("N")) {
                                                gene_node.setProperty("broken", "true");
                                            }*/
                                            gene_node.setProperty("node_ids", annotate_nodes_of_gene(start_ptr, origin, begin, end,gene_node.getId()));
                                            gene_node.setProperty("sequence", gene.toString());
                                            // adding gene_node id to the sequence node
                                            genes_list[g_number][s].add(gene_node.getId());
                                            break;
                                    } // switch
                                } else // if sequence not found
                                {
                                    log.append(sequence_id).append(" missed in genome ").append(g_number).append("\n"); // usually organal genes
                                }
                                System.out.print("\r" + g_number + "\t" + num_genes + "\t" + num_mRNAs + "\t" + num_tRNAs + "\t" + num_ncRNAs + "\t" + num_pgRNAs);
                            }// for trsc
                            tx2.success();
                        } // tx2
                    } // while lines
                    try (Transaction tx3 = graphDb.beginTx()) {
                        if (num_exons != 0) // for the last RNA: if is not the first RNA of the gene
                        {
                            if (rna_len != 0) // if is a coding RNA
                            {
                                while (pq.size() != 0) {
                                    pair = pq.remove();
                                    coding_RNA.append(gene.substring(pair[0] - gene_start_pos, pair[1] - gene_start_pos + 1));
                                }
                                protein = pb.translate(coding_RNA, forward);
                                if (protein.length() != 0) {
                                    ++protein_num;
                                    //out.write(">" + origin + " _" + protein_num + "\n");
                                    //write_fasta(out, protein, 70);
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
                        in.close();
                        //out.close();
                        System.out.println("\r" + g_number + "\t" + num_genes + "\t" + num_mRNAs + "\t" + num_tRNAs + "\t" + num_ncRNAs + "\t" + num_pgRNAs);
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
                        graphDb.findNode(sequence_label, "origin", i + "_" + j).setProperty("genes", genes_array);
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
            File directory = new File(PATH + GRAPH_DATABASE_PATH);
            for (File f : directory.listFiles())
                if (f.getName().startsWith("neostore.transaction.db."))
                    f.delete();
        } else {
            System.out.println("No database found in " + PATH);
            System.exit(1);
        }
    }
    /*
     To return the number of the sequence from genome "g" whose name contains "prefix".
     */

    private int find_sequence(String prefix, int g) {
        boolean found = false;
        int s;
        for (found = false, s = 1; !found && s <= genomeDb.num_sequences[g]; ++s) {
            if (genomeDb.sequence_titles[g][s].contains(prefix)) {
                found = true;
            }
        }
        --s;
        if (found) {
            return s;
        } else {
            return -1;
        }
    }
    /*
     To annotate nodes of the gene with gene_id and returning them in the order they traversed
     */
    public Long[] annotate_nodes_of_gene(IndexPointer start_ptr, String origin, int begin, int end, long gene_id) {
        boolean found;
        Node neighbor, node=graphDb.getNodeById(start_ptr.node_id);
        String lable, name;
        String[] prp = {"F" + origin, "R" + origin};
        LinkedList<Long> node_ids_list = new LinkedList();
        int i, loc, node_len, neighbor_len, s_side, d_side, gene_len = end - begin + 1, length = 0, num_genes, position=start_ptr.position;
        long[] gene_ids, new_gene_ids;
        loc = begin;
        node_len = (int) node.getProperty("length");
        //System.out.println("gene_id "+gene_id);
        //System.out.println("forward:"+forward+" position:"+position+" sequence_number:"+sequence_number+" begin:"+begin+" end:"+end+" gene_len:"+gene_len+" K:"+K);
        if (start_ptr.canonical) {
            if (position + gene_len - 1 <= node_len - 1) {
                length += gene_len;
                loc += gene_len;
            } else {
                length += node_len - position;
                loc += node_len - position;
            }
        } else {
            if (position - (gene_len - 1) >= 0) {
                length += gene_len;
                loc += gene_len;
            } else {
                length += position + 1;
                loc += position + 1;
            }
        }
        found = true;
        if (!node.hasProperty("gene_ids")) {
            node.setProperty("gene_ids", new long[]{gene_id});
        } else {
            gene_ids = (long[]) node.getProperty("gene_ids");
            num_genes = gene_ids.length;
            if (gene_ids[num_genes - 1] != gene_id) // the node is not repeated on the path
            {
                new_gene_ids = new long[num_genes + 1];
                for (i = 0; i < num_genes; ++i) {
                    new_gene_ids[i] = gene_ids[i];
                }
                new_gene_ids[i] = gene_id;
                node.setProperty("gene_ids", new_gene_ids);
            }
        }
        //System.out.println("loc before start:"+loc);
        while (length < gene_len && found) {
            node_ids_list.add(node.getId());
            found = false;
            for (Relationship r : node.getRelationships(Direction.OUTGOING)) {
                neighbor = r.getEndNode();
                name = r.getType().name();
                d_side = (name.charAt(2) - 48) % 2;
                lable = (d_side == 0 ? prp[0] : prp[1]);
                if (neighbor.hasProperty(lable)) {
                    neighbor_len = (int) neighbor.getProperty("length");
                    if (b_search((int[]) neighbor.getProperty(lable), loc - K + 1) >= 0) {
                        found = true;
                        //System.out.println(neighbor.getId()+" "+(loc-K+1));
                        //System.out.println("length:"+length+" neighbor_len:"+neighbor_len+" gene_len:"+gene_len);
                        //System.out.println((length+neighbor_len-K+1)+" ? "+(gene_len));
                        if (length + neighbor_len - K + 1 > gene_len) {
                            //System.out.println((length+neighbor_len-K+1)+" here "+(gene_len));
                            //loc+=gene_len-length;
                            length = gene_len;
                        } else {
                            //System.out.println(" there ");
                            length += neighbor_len - K + 1;
                            loc += neighbor_len - K + 1;
                        }
                        //System.out.println("loc :"+loc+" length :"+length);
                        node = neighbor;
                        if (!node.hasProperty("gene_ids")) {
                            node.setProperty("gene_ids", new long[]{gene_id});
                        } else {
                            gene_ids = (long[]) node.getProperty("gene_ids");
                            num_genes = gene_ids.length;
                            if (gene_ids[num_genes - 1] != gene_id) // the node is not repeated on the path
                            {
                                new_gene_ids = new long[num_genes + 1];
                                for (i = 0; i < num_genes; ++i) {
                                    new_gene_ids[i] = gene_ids[i];
                                }
                                new_gene_ids[i] = gene_id;
                                node.setProperty("gene_ids", new_gene_ids);
                            }
                        }
                        break;
                    }
                    //else
                    //  System.out.println("else: "+neighbor.getId()+" "+((int[])neighbor.getProperty(lable))[0]+" "+(loc-K+1));
                }
            }
            if (!found) {
                for (Relationship r : node.getRelationships(Direction.INCOMING)) {
                    neighbor = r.getStartNode();
                    name = r.getType().name();
                    s_side = (name.charAt(2) - 48) / 2;
                    lable = (s_side == 0 ? prp[1] : prp[0]);
                    if (neighbor.hasProperty(lable)) {
                        neighbor_len = (int) neighbor.getProperty("length");
                        if (b_search((int[]) neighbor.getProperty(lable), loc - K + 1) >= 0) {
                            found = true;
                            if (length + neighbor_len - K + 1 > gene_len) {
                                //loc+=gene_len-length;
                                length = gene_len;
                            } else {
                                length += neighbor_len - K + 1;
                                loc += neighbor_len - K + 1;
                            }
                            node = neighbor;
                            if (!node.hasProperty("gene_ids")) {
                                node.setProperty("gene_ids", new long[]{gene_id});
                            } else {
                                gene_ids = (long[]) node.getProperty("gene_ids");
                                num_genes = gene_ids.length;
                                if (gene_ids[num_genes - 1] != gene_id) // the node is not repeated on the path
                                {
                                    new_gene_ids = new long[num_genes + 1];
                                    for (i = 0; i < num_genes; ++i) {
                                        new_gene_ids[i] = gene_ids[i];
                                    }
                                    new_gene_ids[i] = gene_id;
                                    node.setProperty("gene_ids", new_gene_ids);
                                }
                            }
                            break;
                        }
                    }
                }
            }
        }//for loc
        Long[] node_ids_array = new Long[node_ids_list.size()];
        return node_ids_array;
    }
    /*
    
    */
    public void denovo_homology_annotaion() {
        int i, j, gene_len=0, current_gene_len, num_group_members;
        long[] node_ids;
        long[] gene_ids;
        int copy_number[] = new int[genomeDb.num_genomes+1];
        String current_gene_seq, gene_seq;
        Node group_node, gene_node, node, current_gene_node;
        long gene_id, current_gene_id=-1;
        LinkedList<Node> homolog_gene_nodes = new LinkedList(); 
        PriorityQueue<Long> pq = new PriorityQueue();
        ResourceIterator<Node> genes_iterator;
        Alignment aligner = new Alignment();
        try (Transaction tx = graphDb.beginTx()) {
            genes_iterator=graphDb.findNodes(gene_label);
            while (genes_iterator.hasNext()) {
                num_group_members = 1;
                for (i = 1; i <= genomeDb.num_genomes; ++i) {
                    copy_number[i] = 0;
                }
                gene_node=genes_iterator.next();
                gene_len=(int)gene_node.getProperty("length");
                gene_seq=(String)gene_node.getProperty("sequence");
                if (!gene_node.hasRelationship(RelTypes.contains, Direction.INCOMING)) {
                    node_ids=(long[])gene_node.getProperty("node_ids");
                    for (i=0;i<node_ids.length;++i) {
                        node = graphDb.getNodeById(node_ids[i]);
                        gene_ids=(long[])node.getProperty("gene_ids");
                        for (j = 0 ; j < gene_ids.length ; ++j)
                            pq.offer(gene_ids[j]);
                    }
                    while (!pq.isEmpty()) {
                        gene_id = pq.remove();
                        if(gene_id != current_gene_id) {
                            current_gene_id = gene_id;
                            current_gene_node = graphDb.getNodeById(gene_id);
                            if ( ! current_gene_node.hasRelationship(RelTypes.contains, Direction.INCOMING) ) {
                                // If gene and the current mate are almost equally long and highly similar 
                                current_gene_len = (int)current_gene_node.getProperty("length");
                                if(Math.abs(current_gene_len - gene_len) <= Math.max(current_gene_len, gene_len) / 10 ) {
                                    current_gene_seq=(String)current_gene_node.getProperty("sequence");
                                    if (aligner.similarity(current_gene_seq, gene_seq) > match_score * 0.9) {
                                        num_group_members += 1;
                                        homolog_gene_nodes.add(current_gene_node);
                                        copy_number[(int)current_gene_node.getProperty("genome")] += 1;
                                    }
                                }
                            }
                        }
                    }// while
                    if( num_group_members > 1) {
                        group_node = graphDb.createNode(homolog_lable);
                        group_node.setProperty("num_members", num_group_members);
                        group_node.createRelationshipTo(gene_node, RelTypes.contains);
                        while (!homolog_gene_nodes.isEmpty()) 
                            group_node.createRelationshipTo(homolog_gene_nodes.removeFirst(), RelTypes.contains);
                        for (i = 1; i <= genomeDb.num_genomes; ++i) 
                            group_node.setProperty("copy_number_in_genome_" + i, copy_number[i]);
                    }
                }
            }
            tx.success();
        } 
    }
    
    /*
     To group ortholog/homolog/etc genes. 
     group_names : a text file with lines starting with a group name followed by a colon, followed by genes IDs seperated by one space.
     */
    public void group_ortholog_proteins(String group_names) {
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
            try (BufferedReader in = new BufferedReader(new FileReader(group_names))) {
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
                        System.out.print("\r" + total_proteins + "\t" + num_groups);
                    }
                    tx.success();
                    System.out.println("\r" + total_proteins + "\t" + num_groups);
                }
                in.close();
            } catch (IOException ioe) {
                System.out.println("Failed to open " + group_names);
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
    /*
     To register the action to be taken if the program halts unexpectedly
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
