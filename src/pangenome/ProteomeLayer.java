/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import alignment.ProteinAlignment;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import static java.lang.Integer.max;
import static java.lang.Integer.min;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Queue;
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
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import static org.neo4j.graphdb.factory.GraphDatabaseSettings.keep_logical_logs;
import static pangenome.GenomeLayer.getFolderSize;

import static pantools.Pantools.GRAPH_DATABASE_PATH;
import static pantools.Pantools.RelTypes;
import static pantools.Pantools.graphDb;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.startTime;
import static pantools.Pantools.MAX_TRANSACTION_SIZE;
import static pantools.Pantools.cores;
import static pantools.Pantools.executeCommand_for;
import static pantools.Pantools.heapSize;
import static pantools.Pantools.homology_group_lable;
import static pantools.Pantools.mRNA_label;

/**
 * Implements all the functionalities related to the annotation layer of the pangenome
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class ProteomeLayer {
    private double INTERSECTION;
    private double CONTRAST;
    private double INFLATION;
    private int K_SIZE;
    private int THRESHOLD;
    private int MAX_ALIGNMENT_LENGTH;
    private int MAX_INTERSECTIONS;
    private int MAX_KMER_FREQ;
    private int THREADS;
    private int MAX_KMERS_NUM;
    private AtomicInteger num_intersections;
    private AtomicInteger num_similarities;
    private AtomicInteger num_components;
    private int num_hexamers;
    private int num_genomes;
    private int num_proteins;
    private int num_groups;
    private LinkedList[] kmers_proteins_list;
    private int[] kmer_frequencies;
    private BlockingQueue<Node> proteins;
    private BlockingQueue<intersection> intersections;
    private BlockingQueue<intersection> similarities;
    private BlockingQueue<LinkedList> components; 
    private BlockingQueue<LinkedList> homology_groups_list; 
    private String pangenome_path;
    private Node pangenome_node;
    
    
    public ProteomeLayer(){
        K_SIZE = 6;
        INTERSECTION = 0.09;
        CONTRAST = 8;
        INFLATION = 9.6;
        THRESHOLD = 95;
        MAX_ALIGNMENT_LENGTH  = 1000;
        MAX_INTERSECTIONS  = 10000000;
        THREADS = cores;
        MAX_KMERS_NUM = (int)Math.round(Math.pow(20, K_SIZE));
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
            
    /**
     * Iterates over the proteins available in the database and put them in the 
     * size-limited Blocking-queue "proteins".
     * Only one thread should be called for this runnable.
     */
    public class Generate_proteins implements Runnable {
        public Generate_proteins() {
        }

        @Override
        public void run() {
            ResourceIterator<Node>  proteins_iterator;
            try(Transaction tx = graphDb.beginTx()){
                pangenome_node = graphDb.findNodes(pangenome_label).next();
                proteins_iterator = graphDb.findNodes(mRNA_label);
                while (proteins_iterator.hasNext())
                    proteins.offer(proteins_iterator.next());
                proteins_iterator.close();
                tx.success();
            }
        }
    }    
     
    /**
     * Takes proteins from the Blocking-queue "proteins" and stores the count
     * of k-mers of the proteome in the aaray "kmer_frequencies". K-mers are
     * represented by numbers in the base of 20.
     * Only one thread should be called for this runnable.
     */
    public class count_kmers implements Runnable {
        int num_proteins;
        public count_kmers(int num) {
            num_proteins = num;
        }

        @Override
        public void run() {
            int i = 0, c;
            Node protein_node;
            int protein_length, kmer_index;
            String protein;
            int[] code = new int[256];
            char[] aminoacids = new char[]
            {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
            for (i = 0; i < 20; ++i)
                code[aminoacids[i]] = i;
            kmer_frequencies = new int[MAX_KMERS_NUM];
            try{
                try (Transaction tx = graphDb.beginTx()) {
                    for (c = 0; c < num_proteins; ++c){
                        protein_node = proteins.take();
                        protein = (String)protein_node.getProperty("protein", "");
                        protein_length = protein.length();
                        if (protein_length > K_SIZE){
                            kmer_index = 0;
                            for (i = 0; i < K_SIZE; ++i)
                                kmer_index = kmer_index * 20 + code[protein.charAt(i)];
                            for (; i < protein_length; ++i){// for each kmer of the protein
                                if (kmer_frequencies[kmer_index] == 0)
                                    ++num_hexamers;
                                kmer_frequencies[kmer_index] += 1;
                                kmer_index = kmer_index % (MAX_KMERS_NUM / 20) * 20 + code[protein.charAt(i)];
                            }                            
                        }
                    }
                    tx.success();
                }
            } catch(InterruptedException e){
                System.err.println(e.getMessage());
            }
        }
    }    

    /**
     * Takes proteins from the Blocking-queue "proteins" and k-merizes the 
     * proteins. K-mers are represented by numbers in the base of 20.
     * For each k-mer a list of proteins containing that k-mer is stored in
     * the array "kmers_proteins_list" at the index represent by that k-mer.
     * Only one thread should be called for this runnable.
     */
    public class Kmerize_proteins implements Runnable {
        int num_proteins;
        public Kmerize_proteins(int num) {
            num_proteins = num;
        }

        @Override
        public void run() {
            int i = 0, c, chunk = num_proteins / 40;
            Node protein_node;
            int protein_length, kmer_index;
            String protein;
            long protein_id;
            int[] code = new int[256];
            char[] aminoacids = new char[]
            {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
            for (i = 0; i < 20; ++i)
                code[aminoacids[i]] = i;
            try{
                try (Transaction tx = graphDb.beginTx()) {
                    kmers_proteins_list = new LinkedList[MAX_KMERS_NUM];
                    for (c = 0; c < num_proteins; ++c){
                        protein_node = proteins.take();
                        protein = (String)protein_node.getProperty("protein", "");
                        protein_length = protein.length();
                        protein_id = protein_node.getId();
                        if (protein_length > K_SIZE){
                            kmer_index = 0;
                            for (i = 0; i < K_SIZE; ++i)
                                kmer_index = kmer_index * 20 + code[protein.charAt(i)];
                            for (; i < protein_length; ++i){// for each kmer of the protein
                            // ignore extremely rare and abundant k-mers    
                                if (kmer_frequencies[kmer_index] > 1 + num_genomes / 2 && kmer_frequencies[kmer_index] < MAX_KMER_FREQ){
                                    if (kmers_proteins_list[kmer_index] == null)
                                        kmers_proteins_list[kmer_index] = new LinkedList();
                                    kmers_proteins_list[kmer_index].add(protein_id);
                                }
                                kmer_index = kmer_index % (MAX_KMERS_NUM / 20) * 20 + code[protein.charAt(i)];
                            }                            
                        }
                        if (c % chunk == 0)
                            System.out.print("|");
                    }
                    tx.success();
                    kmer_frequencies = null;
                }
            } catch(InterruptedException e){
                System.err.println(e.getMessage());
            }
        }
    }    

    /**
     * Takes proteins from the Blocking-queue "proteins" and detects pairs of 
     * intersecting proteins and put them in a Blocking-queue "intersections".
     * Only one thread should be called for this runnable.
     */
    public class Find_intersections implements Runnable {
        int num_proteins;
        double frac = INTERSECTION;
        int max_intersection = MAX_INTERSECTIONS;
        public Find_intersections(int num) {
            num_proteins = num;
        }

        @Override
        public void run() {
            int i, chunk = num_proteins / 40;
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
                try{
                    for (p = 0; p < num_proteins; ++p) {
                        protein_node = proteins.take();
                        protein = (String)protein_node.getProperty("protein");
                        protein_length = protein.length();
                        if (protein_length > K_SIZE){
                            protein_id = protein_node.getId();
                            num_ids = 0;
                            kmer_index = 0;
                            for (i = 0; i < K_SIZE; ++i)
                                kmer_index = kmer_index * 20 + code[protein.charAt(i)];
                            for (; i < protein_length && num_ids < max_intersection; ++i){// for each kmer of the protein
                                if (kmers_proteins_list[kmer_index] != null){
                                    proteins_list = kmers_proteins_list[kmer_index].iterator();
                                    while(proteins_list.hasNext() && num_ids < max_intersection){
                                        crossing_protein_id = proteins_list.next();
                                    // Only crossing proteins with a higher ID
                                        if (crossing_protein_id > protein_id){ 
                                            crossing_protein_ids[num_ids++] = crossing_protein_id;
                                        }
                                    }
                                }
                                kmer_index = kmer_index % (MAX_KMERS_NUM / 20) * 20 + code[protein.charAt(i)];
                            }
                        // Sorts the crossing protein IDs to count the number of shared k-mers.    
                            Arrays.sort(crossing_protein_ids, 0, num_ids);
                            for (i = 0, counter = 0, crossing_protein_id = crossing_protein_ids[0]; i < num_ids ; ++i){
                                p_id = crossing_protein_ids[i];
                            // New run of protein IDs    
                                if (crossing_protein_id != p_id){
                                    if(counter > 1){
                                        crossing_protein_node = graphDb.getNodeById(crossing_protein_id);
                                        crossing_protein = (String)crossing_protein_node.getProperty("protein");
                                        shorter_len = Math.min(protein_length, crossing_protein.length());
                                        if (counter >= frac * (shorter_len - K_SIZE + 1)){
                                            intersections.offer(new intersection(protein_node, crossing_protein_node,0));
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
                                if (counter >= frac * (shorter_len - K_SIZE + 1)){
                                    intersections.offer(new intersection(protein_node, crossing_protein_node,0));
                                    num_intersections.getAndIncrement();
                                }
                            }
                            if (p % chunk == 0)
                                System.out.print("|");
                        }
                    }// for protein
                // Signify the end of intersections queue.    
                    for (i = 0; i < THREADS; ++i)
                        intersections.offer(new intersection(null, null,0));// end of queue
                } catch(InterruptedException e){
                    System.err.println(e.getMessage());
                }
                tx.success();
            }
        }
    } 
    
    /**
     * Takes the intersections from the Blocking-queue "intersections" and 
     * calculates the similarity score of the protein pair, if it is high enough
     * writes it in the intersection object and put the object in the Blocking-
     * queue "similarities".
     * Multiple threads can call this runnable.
     */
    public class Find_similarities implements Runnable {
        int threshold = THRESHOLD;
        int max_alg_len = MAX_ALIGNMENT_LENGTH;
        ProteinAlignment aligner;
        public Find_similarities() {
            aligner = new ProteinAlignment(-10,-1,max_alg_len);
        }

        @Override
        public void run() {
            Node protein_node1, protein_node2;
            String protein1, protein2;
            intersection ints;
            try{
                try(Transaction tx = graphDb.beginTx()){
                    ints = intersections.take();
                    while (ints.protein1 != null) {
                        protein_node1 = ints.protein1;
                        protein_node2 = ints.protein2;
                        protein1 = (String)protein_node1.getProperty("protein");
                        protein2 = (String)protein_node2.getProperty("protein");
                        ints.similarity = ((double)protein_similarity(aligner, protein1, protein2))/perfect_score(aligner, protein1, protein2)*100;
                        if (ints.similarity > threshold){
                            similarities.offer(ints);
                            num_similarities.getAndIncrement();
                        }
                        ints = intersections.take();
                    }
                // Signify the end of the similarities queue.   
                    similarities.offer(new intersection(null, null,0));
                    tx.success();
                }
            }catch(InterruptedException e){
                System.err.println(e.getMessage());
            }
        }
    }    

    /**
     * Takes the similarities from the Blocking-queue "similarities" and write 
     * it in the graph database.
     * Only one thread should call this runnable to avoid Deadlock exceptions.
     */
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
     * Build lists containing proteins of similarity components and put then in
     * the Blocking-queue "components".
     * Only one thread should call this runnable.
     */
    public class build_similarity_components implements Runnable {
        int num_proteins;
        String pangenome_path;
        public build_similarity_components(String path, int num) {
            pangenome_path = path;
            num_proteins = num;
        }

        @Override
        public void run(){
            int i;
            Node protein_node;
            LinkedList<Node> component = new LinkedList();
            try{
                for (i = 0; i< num_proteins; ++i) {
                    protein_node = proteins.take();
                    breadth_first_search(component, protein_node);
                    if (component.size() > 0){ 
                        components.offer(component);
                        num_components.getAndIncrement();
                        component = new LinkedList();
                    }
                } 
            // Signifies the end of the components queue for all the threads    
                for (i = 0; i < THREADS; ++i)
                   components.offer(new LinkedList());
            } catch(InterruptedException e){
                System.err.println(e.getMessage());
            }
        }
        
        /**
         * Puts all the proteins similar to a protein in one component
         * @param component An empty list to be filled with proteins similar 
         * to the query protein and the query protein itself
         * @param start_protein The query protein 
         */
        private void breadth_first_search(LinkedList<Node> component, Node start_protein){
            long start_id = start_protein.getId();
            Node crossing_protein;
            try (Transaction tx = graphDb.beginTx()) {
            // To avoid having one protein in different groups    
                if (!start_protein.hasProperty("component")) { 
                    start_protein.setProperty("component", start_id);
                    Queue<Node> homologs = new LinkedList();
                    homologs.add(start_protein);
                    // for all the candidates with some shared node with the protein    
                    while (!homologs.isEmpty()) {
                        start_protein = homologs.remove();
                        component.add(start_protein);
                        for (Relationship crossing_edge: start_protein.getRelationships(RelTypes.is_similar_to)) {
                            crossing_protein = crossing_edge.getOtherNode(start_protein);
                            if(!crossing_protein.hasProperty("component")){
                                crossing_protein.setProperty("component", start_id);
                                homologs.add(crossing_protein);
                            }
                        }
                    }// while
                }
                tx.success(); 
            }
        }  
    }      

    /**
     * Takes similarity components from the Blocking-queue "components" and 
     * splits them into homology groups represented by a list stored in the
     * Blocking-queue "homology_groups_list".
     * Multiple threads can call this runnable.
     */
    public class build_homology_groups implements Runnable {
        int num_proteins;
        String pangenome_path;
        double[][] phylogeny_distance;
        int[][] count;
        public build_homology_groups(String path, int num) {
            pangenome_path = path;
            num_proteins = num;
        }

        @Override
        public void run(){
            int i;
            LinkedList<Node> component;
            phylogeny_distance = new double[num_genomes + 1][];
            count = new int[num_genomes + 1][];
            for (i = 1; i < phylogeny_distance.length; ++i){
                phylogeny_distance[i] = new double[num_genomes + 1];
                count[i] = new int[num_genomes + 1];
            } 
            try{
                component = components.take();
                while (!component.isEmpty()) {// Not finished
                    if (component.size() == 1){
                        homology_groups_list.offer(component);
                    }
                    else
                        break_component(component, pangenome_path);
                    component = components.take();
                }
            } catch(InterruptedException e){
                System.err.println(e.getMessage());
            }
        }

        /**
         * Splits the similarity component into homology groups using MCL algorithm 
         * @param component The similarity component
         * @param pangenome_path The path to the current graph database
         */   
        void break_component(LinkedList<Node> component, String pangenome_path){
            int i, group_size, wating_time, time;
            double infl;
            LinkedList<Node> homology_group, singletons_group = new LinkedList();
            String graph_path, clusters_path, line, command;
            String[] fields;
            BufferedReader clusters_file;
            File tmp_file;
            group_size = component.size();
            graph_path = pangenome_path + "/" + component.getFirst().getId() + ".graph";
            clusters_path = pangenome_path + "/" + component.getFirst().getId() + ".clusters";
        // Prepare the input file for MCL
            write_similaity_matrix(component, graph_path);
        //  Estimate the run-time of MCL
            time = wating_time = 1 + (int)Math.round(group_size / 100000000.0 * group_size);
            for( infl = INFLATION; infl < 30; infl += 0.5){
                command = "mcl " + graph_path + " --abc -I " + infl + " -o " + clusters_path;
                if(executeCommand_for(command, wating_time))
                    break;
                if((tmp_file = new File(clusters_path)).exists())
                    tmp_file.delete();
                wating_time += time;
            }
            if (infl >= 30 ){
                System.err.println("Failed to split group ID = " + component.getFirst().getId());
                homology_groups_list.offer(component);
            } else {
                try (Transaction tx = graphDb.beginTx()) {
                    try{
                        clusters_file = new BufferedReader(new FileReader(clusters_path));
                    // For each line of the MCL output    
                        while (clusters_file.ready()){
                            line = clusters_file.readLine().trim();
                            if (line.equals("")) // if line is empty
                                continue;
                            fields = line.split("\\s");
                            if (fields.length > 1){
                                homology_group = new LinkedList();
                                for (i = 0; i < fields.length; ++i)
                                    homology_group.add(graphDb.getNodeById(Long.parseLong(fields[i])));
                                homology_groups_list.offer(homology_group);
                            } else { // if is a singleton group
                                singletons_group.add(graphDb.getNodeById(Long.parseLong(fields[0])));
                            }
                        }
                    // put all singletons of the component in one homology groups    
                        if (singletons_group.size() > 0){
                                homology_groups_list.offer(singletons_group);
                        }
                        clusters_file.close();
                        new File(clusters_path).delete();
                        new File(graph_path).delete();
                    }catch (IOException ex){
                        System.out.print(ex.getMessage());
                    }         
                    tx.success();
                }
            }
        }   

        /**
         * Given a homology group complete all the pairwise similarities between its members and 
         * writes the weighted edges of the graph in SIMILARITY_GRAPH_FILE_NAME to be used by MCL clustering algorithm.
         * @param homology_group_node The homology group
         */
        private void write_similaity_matrix(LinkedList<Node> component, String graph_path){
            ListIterator<Node> itr;
            Node protein1_node, protein2_node;
            int genome1, genome2;
            double similarity;
            try (PrintWriter graph = new PrintWriter(graph_path)){
                try (Transaction tx = graphDb.beginTx()) {
                    calculate_phylogeny_distances(component);
                    for (itr = component.listIterator(); itr.hasNext(); ){
                        protein1_node = itr.next();
                        genome1 = (int)protein1_node.getProperty("genome");
                        for (Relationship homology_edge: protein1_node.getRelationships(RelTypes.is_similar_to, Direction.OUTGOING)){
                            protein2_node = homology_edge.getEndNode();
                            similarity = (double)homology_edge.getProperty("similarity");
                            similarity -= THRESHOLD;
                            genome2 = (int)protein2_node.getProperty("genome");
                            similarity += phylogeny_distance[genome1][genome2];
                            graph.write(protein1_node.getId()+" "+protein2_node.getId()+" "+ Math.pow(similarity, CONTRAST) + "\n");
                        }
                    }
                    tx.success();
                }
                graph.close();
            } catch (IOException ex){
            }
        }

        private void calculate_phylogeny_distances(LinkedList<Node> component){
            int i, j, g1, g2;
            Node p1, p2;
            double similarity;
            ListIterator<Node> proteins_itr = component.listIterator();
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

    }      
    
    /**
     * Takes the homology groups from the Blocking-queue "homology_groups_list"
     * and created homology groups in the graph database and writes the groups
     * in "pantools_homologs.txt".
     * Only one thread can call this runnable.
     */
    public class write_homology_groups implements Runnable {
        String pangenome_path;
        int num_proteins;
        public write_homology_groups(String path, int num) {
            pangenome_path = path;
            num_proteins = num;
        }

        @Override
        public void run(){
            int i, trs, p = 0, chunk = num_proteins / 40;
            int[] copy_number = new int[num_genomes + 1];
            LinkedList<Node> homology_group;
            Node homology_node, protein_node;
            BufferedWriter groups_file;
            try{
                try{
                    groups_file = new BufferedWriter(new FileWriter(pangenome_path + "/pantools_homologs.txt"));
                    while (p < num_proteins) {
                        try(Transaction tx = graphDb.beginTx()){
                            for (trs = 0; trs < MAX_TRANSACTION_SIZE && p < num_proteins; ++trs){
                                homology_group = homology_groups_list.take();
                                ++num_groups;
                                homology_node = graphDb.createNode(homology_group_lable); 
                                homology_node.setProperty("num_members", homology_group.size());
                                groups_file.write(Long.toString(homology_node.getId()) + ":");
                                while (!homology_group.isEmpty()){
                                    protein_node = homology_group.remove();
                                    homology_node.createRelationshipTo(protein_node, RelTypes.has_homolog);
                                    ++copy_number[(int)protein_node.getProperty("genome")];
                                    groups_file.write(" " + ((String)protein_node.getProperty("protein_ID")).replace(' ', '_'));
                                    if (p % chunk == 0)
                                        System.out.print("|");
                                    ++p;
                                }
                                groups_file.write("\n");
                                homology_node.setProperty("copy_number", copy_number);
                                for (i=0; i < copy_number.length; ++i)
                                    copy_number[i] = 0;
                            }
                            tx.success();
                        }
                    }
                    groups_file.close();
                }catch (IOException ex){
                    System.out.print(ex.getMessage());
                }                
            } catch(InterruptedException e){
                System.err.println(e.getMessage());
            }
        }
    }      

    /**
     * Constructs the proteome layer from a set of proteins.
     * @param protein_paths_file A test file containing paths to the proteome FASTA files
     * @param pangenome_path Path to the pan-genome graph database
     */
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
     * Groups the similar proteins into homology groups
     * @param args The command line arguments, args[1] Path to the database folder
     */
    public void group(String[] args) {
        pangenome_path = args[1];
        int i, n, d, p;
        double x;
        startTime = System.currentTimeMillis();
        for (i = 2; i < args.length; ++i){
            switch (args[i]){
                case "-i":
                    x = Double.parseDouble(args[i + 1]);
                    if (x >= 0.001 && x <= 0.1)
                        INTERSECTION = x;
                    ++i;
                    break;
                case "-t":
                    n = Integer.parseInt(args[i + 1]);
                    if (n > 0 && n < 100)
                        THRESHOLD = n;
                    ++i;
                    break;
                case "-m":
                    x = Double.parseDouble(args[i + 1]);
                    if (x > 1 && x < 19)
                        INFLATION = x;
                    ++i;
                    break;
                case "-c":
                    x = Double.parseDouble(args[i + 1]);
                    if (x > 0 && x < 10)
                        CONTRAST = x;
                    ++i;
                    break;
                case "-d":
                    d = Integer.parseInt(args[i + 1]);
                    d = d < 1 ? 1 : d;
                    d = d > 8 ? 8 : d;
                    INTERSECTION = new double[] {0, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02}[d];
                    THRESHOLD = new int[]   {0, 95, 85, 75, 65, 55, 45, 35, 25 }[d];
                    INFLATION = new double[]{0, 9.6, 8.4, 7.2, 6.0, 4.8, 3.6, 2.4, 1.2}[d];
                    CONTRAST = new double[] {0, 8, 7, 6, 5, 4, 3, 2, 1 }[d];
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
        System.out.println("Intersection rate = " + INTERSECTION);
        System.out.println("Threshold = " + THRESHOLD);
        System.out.println("MCL inflation = " + INFLATION);
        System.out.println("Contrast = " + CONTRAST);
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(pangenome_path + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        proteins = new LinkedBlockingQueue<>((int)(heapSize/50));
        intersections = new LinkedBlockingQueue<>((int)(heapSize/200));
        similarities = new LinkedBlockingQueue<>((int)(heapSize/200));
        num_intersections = new AtomicInteger(0);
        num_similarities = new AtomicInteger(0);
        num_components = new AtomicInteger(0);
        try(Transaction tx = graphDb.beginTx()){
            pangenome_node = graphDb.findNodes(pangenome_label).next();
            num_genomes = (int)pangenome_node.getProperty("num_genomes");
            num_proteins = (int)pangenome_node.getProperty("num_proteins");
            tx.success();
        }

        MAX_KMER_FREQ = num_proteins / 1000 + 50 * num_genomes; //  because of probability p and copy number of 2
        num_hexamers = 0;
        System.out.println("\nKmerizing proteins :");
        System.out.print("0 ......................................... 100\n  "); 
        try{
            ExecutorService es = Executors.newFixedThreadPool(2);
            es.execute(new Generate_proteins());
            es.execute(new count_kmers(num_proteins));
            es.shutdown();
            es.awaitTermination(10, TimeUnit.DAYS);        
        } catch (InterruptedException e){
            
        }
        try{
            ExecutorService es = Executors.newFixedThreadPool(2);
            es.execute(new Generate_proteins());
            es.execute(new Kmerize_proteins(num_proteins));
            es.shutdown();
            es.awaitTermination(10, TimeUnit.DAYS);        
        } catch (InterruptedException e){
            
        }
        
        /*try(Transaction tx = graphDb.beginTx()){
            Node node;
            ResourceIterator<Node> itr = graphDb.findNodes(mRNA_label);
            while(itr.hasNext()){
                node = itr.next();
                node.removeProperty("component");
                for (Relationship r: node.getRelationships(RelTypes.has_homolog, Direction.INCOMING))
                    r.delete();
            }
            tx.success();
        }*/
        System.out.println("\n\nFinding intersections and similarities:");
        System.out.print("0 ......................................... 100\n  "); 
        try{
            ExecutorService es = Executors.newFixedThreadPool(THREADS);
            es.execute(new Generate_proteins());
            es.execute(new Find_intersections(num_proteins));
            for(i = 1; i <= THREADS - 2; i++)
                es.execute(new Find_similarities());;
            es.execute(new Write_similarities());
            es.shutdown();
            es.awaitTermination(10, TimeUnit.DAYS);        
        } catch (InterruptedException e){
            
        }
        kmers_proteins_list = null;   
        intersections = null;
        similarities = null;
        System.gc();

        System.out.println("\n\nBuilding homology groups : ");
        System.out.print("0 ......................................... 100\n  "); 
        components = new LinkedBlockingQueue<>((int)(heapSize/200));
        homology_groups_list = new LinkedBlockingQueue<>((int)(heapSize/200));
        try{
            ExecutorService es = Executors.newFixedThreadPool(2);
            es.execute(new Generate_proteins());
            es.execute(new build_similarity_components(pangenome_path, num_proteins));
            for(i = 1; i <= THREADS - 2; i++)
                es.execute(new build_homology_groups(pangenome_path, num_proteins));
            es.execute(new write_homology_groups(pangenome_path, num_proteins));
            es.shutdown();
            es.awaitTermination(10, TimeUnit.DAYS);        
        } catch (InterruptedException e){
            
        }
        System.out.println("\n\nHexamers = " + num_hexamers);
        System.out.println("Proteins = " + num_proteins);
        System.out.println("Intersections = " + num_intersections.intValue());
        System.out.println("Similarities = " + num_similarities.intValue());
        System.out.println("Components = " + num_components);
        System.out.println("Groups = " + num_groups);
        System.out.println("Database size = " + getFolderSize(new File(pangenome_path + GRAPH_DATABASE_PATH)) + " MB");        

        File directory = new File(pangenome_path + GRAPH_DATABASE_PATH);
        for (File f : directory.listFiles()) {
            if (f.getName().startsWith("neostore.transaction.db.")) {
                f.delete();
            }
        }
    }
        
    /**
     * The similarity score of the shorter protein with itself  
     * @param aligner The protein aligner object
     * @param p1 The first protein
     * @param p2 The second protein
     * @return 
     */
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
