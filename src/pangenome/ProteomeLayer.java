/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import alignment.LocalSequenceAlignment;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
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
import java.util.logging.Level;
import java.util.logging.Logger;
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
import static pantools.Pantools.startTime;
import static pantools.Pantools.MAX_TRANSACTION_SIZE;
import static pantools.Pantools.PATH_TO_THE_PANGENOME_DATABASE;
import static pantools.Pantools.MIN_PROTEIN_IDENTITY;
import static pantools.Pantools.INTERSECTION_RATE;
import static pantools.Pantools.CONTRAST;
import static pantools.Pantools.GAP_EXT;
import static pantools.Pantools.GAP_OPEN;
import static pantools.Pantools.MCL_INFLATION;
import static pantools.Pantools.PATH_TO_THE_PROTEOMES_FILE;
import static pantools.Pantools.THREADS;
import static pantools.Pantools.executeCommand_for;
import static pantools.Pantools.heapSize;
import static pantools.Pantools.homology_group_label;
import static pantools.Pantools.mRNA_label;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.genomeDb;
import static pantools.Pantools.genomeSc;
import static pantools.Pantools.graphDb;
import static pantools.Pantools.indexDb;

/**
 * Implements all the functionalities related to the annotation layer of the pangenome
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class ProteomeLayer {
    private int PEPTIDE_SIZE = 6;
    private int MAX_INTERSECTIONS;
    private int MAX_KMER_FREQ;
    private int MAX_KMERS_NUM;
    private AtomicInteger num_intersections;
    private AtomicInteger num_similarities;
    private AtomicInteger num_components;
    private AtomicInteger similarity_bars;
    private int num_hexamers;
    private int num_genomes;
    private int num_proteins;
    private int num_groups;
    private long[][] kmers_proteins_list;
    private int[] kmer_frequencies;
    private BlockingQueue<Node> proteins;
    private BlockingQueue<intersection> intersections;
    private BlockingQueue<intersection> similarities;
    private BlockingQueue<LinkedList> components; 
    private BlockingQueue<LinkedList> homology_groups_list; 
    private Node pangenome_node;
    
    public ProteomeLayer(){
        MAX_INTERSECTIONS  = 10000000;
        MAX_KMERS_NUM = (int)Math.round(Math.pow(21, PEPTIDE_SIZE));
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
                proteins.clear();
                pangenome_node = graphDb.findNodes(pangenome_label).next();
                proteins_iterator = graphDb.findNodes(mRNA_label);
                while (proteins_iterator.hasNext())
                    try {
                        proteins.put(proteins_iterator.next());
                    } catch (InterruptedException ex) {
                        Logger.getLogger(ProteomeLayer.class.getName()).log(Level.SEVERE, null, ex);
                    }
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
            {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'};
            for (i = 0; i < 21; ++i)
                code[aminoacids[i]] = i;
            kmer_frequencies = new int[MAX_KMERS_NUM];
            try{
                try (Transaction tx = graphDb.beginTx()) {
                    for (c = 0; c < num_proteins; ++c){
                        protein_node = proteins.take();
                        if (protein_node.hasProperty("protein") ){
                            protein = (String)protein_node.getProperty("protein", "");
                            protein_length = protein.length();
                            if (protein_length > PEPTIDE_SIZE){
                                kmer_index = 0;
                                for (i = 0; i < PEPTIDE_SIZE; ++i)
                                    kmer_index = kmer_index * 21 + code[protein.charAt(i)];
                                for (; i < protein_length; ++i){// for each kmer of the protein
                                    if (kmer_frequencies[kmer_index] == 0)
                                        ++num_hexamers;
                                    kmer_frequencies[kmer_index] += 1;
                                    kmer_index = kmer_index % (MAX_KMERS_NUM / 21) * 21 + code[protein.charAt(i)];
                                }                            
                            }
                        }
                    }
                    tx.success();
                }
                kmers_proteins_list = new long[MAX_KMERS_NUM][];
                for (i = 0; i < MAX_KMERS_NUM; ++i){
                    if (kmer_frequencies[i] > 1 && kmer_frequencies[i] < MAX_KMER_FREQ)// + num_genomes / 2
                        kmers_proteins_list[i] = new long[kmer_frequencies[i]];
                    else
                        kmers_proteins_list[i] = null;
                    kmer_frequencies[i] = 0;
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
            int i = 0, c, chunk = num_proteins > 40 ? num_proteins / 40 : 1;
            Node protein_node;
            int protein_length, kmer_index;
            String protein;
            long protein_id;
            int[] code = new int[256];
            char[] aminoacids = new char[]
            {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'};
            for (i = 0; i < 21; ++i)
                code[aminoacids[i]] = i;
            try{
                try (Transaction tx = graphDb.beginTx()) {
                    for (c = 0; c < num_proteins; ++c){
                        protein_node = proteins.take();
                        if (protein_node.hasProperty("protein") ){
                            protein = (String)protein_node.getProperty("protein", "");
                            protein_length = protein.length();
                            protein_id = protein_node.getId();
                            if (protein_length > PEPTIDE_SIZE){
                                kmer_index = 0;
                                for (i = 0; i < PEPTIDE_SIZE; ++i)
                                    kmer_index = kmer_index * 21 + code[protein.charAt(i)];
                                for (; i < protein_length; ++i){// for each kmer of the protein
                                // ignore extremely rare and abundant k-mers    
                                    if (kmers_proteins_list[kmer_index] != null){
                                        kmers_proteins_list[kmer_index][kmer_frequencies[kmer_index]] = protein_id;
                                        ++kmer_frequencies[kmer_index];
                                    }
                                    kmer_index = kmer_index % (MAX_KMERS_NUM / 21) * 21 + code[protein.charAt(i)];
                                }                            
                            }
                            if (c % chunk == 0)
                                System.out.print("|");
                        }
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
        double frac = INTERSECTION_RATE;
        int max_intersection = MAX_INTERSECTIONS;
        public Find_intersections(int num) {
            num_proteins = num;
        }

        @Override
        public void run() {
            int i, j, len, chunk = num_proteins > 40 ? num_proteins / 40 : 1;
            int p, counter, num_ids, kmer_index, num_ins = 0;
            long[] crossing_protein_ids= new long[max_intersection];
            int[] code = new int[256];
            char[] aminoacids = new char[]
            {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'};
            for (i = 0; i < 21; ++i)
                code[aminoacids[i]] = i;
            Node protein_node, crossing_protein_node;
            long protein_id;
            int protein_length, shorter_len;
            String protein, crossing_protein;
            long crossing_protein_id, p_id;
            
            try (Transaction tx = graphDb.beginTx()) {
                try{
                    for (p = 0; p < num_proteins; ++p) {
                        protein_node = proteins.take();
                        if (protein_node.hasProperty("protein") ){
                            protein = (String)protein_node.getProperty("protein");
                            protein_length = protein.length();
                            if (protein_length > PEPTIDE_SIZE){
                                protein_id = protein_node.getId();
                                num_ids = 0;
                                kmer_index = 0;
                                for (i = 0; i < PEPTIDE_SIZE; ++i)
                                    kmer_index = kmer_index * 21 + code[protein.charAt(i)];
                                for (; i < protein_length && num_ids < max_intersection; ++i){// for each kmer of the protein
                                    if (kmers_proteins_list[kmer_index] != null){
                                        len = kmers_proteins_list[kmer_index].length;
                                        for (j = 0; j < len && num_ids < max_intersection; ++j){
                                            crossing_protein_id = kmers_proteins_list[kmer_index][j];
                                        // Only crossing proteins with a higher ID
                                            if (crossing_protein_id > protein_id){ 
                                                crossing_protein_ids[num_ids++] = crossing_protein_id;
                                            }
                                        }
                                    }
                                    kmer_index = kmer_index % (MAX_KMERS_NUM / 21) * 21 + code[protein.charAt(i)];
                                }
                            // Sorts the crossing protein IDs to count the number of shared k-mers.    
                                Arrays.sort(crossing_protein_ids, 0, num_ids);
                                for (i = 0, counter = 1, crossing_protein_id = crossing_protein_ids[0]; i < num_ids ; ++i){
                                    p_id = crossing_protein_ids[i];
                                // New run of protein IDs    
                                    if (crossing_protein_id != p_id){
                                        if(counter > 1){
                                            crossing_protein_node = graphDb.getNodeById(crossing_protein_id);
                                            crossing_protein = (String)crossing_protein_node.getProperty("protein");
                                            shorter_len = Math.min(protein_length, crossing_protein.length());
                                            if (counter >= frac * (shorter_len - PEPTIDE_SIZE + 1)){
                                                intersections.put(new intersection(protein_node, crossing_protein_node,0));
                                                ++num_ins;
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
                                    if (counter >= frac * (shorter_len - PEPTIDE_SIZE + 1)){
                                        intersections.put(new intersection(protein_node, crossing_protein_node,0));
                                        ++num_ins;
                                    }
                                }
                                if (p % chunk == 0)
                                    System.out.print("|");
                            }
                        }
                    }// for protein
                    num_intersections.getAndAdd(num_ins);
                    System.out.println("\nIntersections = " + num_ins + "\n");
                    System.out.println("Calculating similarities:");
                    System.out.print("0 ......................................... 100\n  "); 
                // Signify the end of intersections queue.    
                    for (i = 0; i < THREADS; ++i)
                        intersections.put(new intersection(null, null,-1));// end of queue
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
        int MAX_ALIGNMENT_LENGTH = 1000;
        int m, n, threshold = MIN_PROTEIN_IDENTITY;
        int processed = 0;
        StringBuilder query;
        StringBuilder subject;
        LocalSequenceAlignment aligner;
        public Find_similarities() {
            aligner = new LocalSequenceAlignment(GAP_OPEN, GAP_EXT,MAX_ALIGNMENT_LENGTH, 0, 'P');
            query = new StringBuilder();
            subject = new StringBuilder();
        }

        @Override
        public void run() {
            Node protein_node1, protein_node2;
            String protein1, protein2;
            intersection ints;
            int num_ints = 0;
            boolean all_intersections_found = false;
            try{
                try(Transaction tx = graphDb.beginTx()){
                    ints = intersections.take();
                    while (ints.protein1 != null) {
                        ++processed;
                        protein_node1 = ints.protein1;
                        protein_node2 = ints.protein2;
                        protein1 = (String)protein_node1.getProperty("protein");
                        protein2 = (String)protein_node2.getProperty("protein");
                        m = protein1.length();
                        n = protein2.length();
                        if (m == n)
                            ints.similarity = aligner.get_match_percentage(protein1, protein2);
                        else 
                            if (m < n)
                            ints.similarity = similarity_percantage(protein1, protein2);
                        else
                            ints.similarity = similarity_percantage(protein2, protein1);
                        if (ints.similarity > threshold){
                            similarities.put(ints);
                            num_similarities.getAndIncrement();
                        }
                        ints = intersections.take();
                        if (!all_intersections_found){
                            num_ints = num_intersections.intValue();
                            if (num_ints > 0){
                                all_intersections_found = true;
                                num_ints -= processed;
                            }
                            if (num_ints < 40)
                                num_ints = 40;
                        }
                        if (all_intersections_found && processed % (num_ints / 40) == 0){
                            System.out.print("|");
                            similarity_bars.getAndIncrement();
                        }
                    }
                // Signify the end of the similarities queue.   
                    similarities.put(new intersection(null, null,0));
                    tx.success();
                }
            }catch(InterruptedException e){
                System.err.println(e.getMessage());
            }
        }
        /**
         * Given two proteins calculates the normalized similarity score between them which is less or equal to 1.
         * Proteins longer than MAX_LENGTH will be broken in smaller parts to be compared correspondingly.  
         * @param p1 The first protein
         * @param p2 The second protein
         * @return The normalized similarity score which is less or equal to 1
         */
        double similarity_percantage(String p1, String p2){
            int m, n,i, parts_num = 1, part_len1, part_len2;
            long score = 0, p_score = 0;
            m = p1.length();
            n = p2.length();
            if (n > MAX_ALIGNMENT_LENGTH){
                parts_num = (n / MAX_ALIGNMENT_LENGTH) + (n % MAX_ALIGNMENT_LENGTH == 0 ? 0 : 1);
                part_len1 = m / parts_num;
                part_len2 = n / parts_num;
                for (i = 0; i < parts_num; ++i){
                    query.setLength(0);
                    subject.setLength(0);
                    query.append(p1.substring(i * part_len1, Math.min(m, (i + 1) * part_len1)));
                    subject.append(p2.substring(i * part_len2, Math.min(n, (i + 1) * part_len2)));
                    aligner.align(query, subject);
                    score += aligner.get_similarity();
                    p_score += aligner.get_match_score(query, query);//5 * query.length();
                }
            } else {
                query.setLength(0);
                subject.setLength(0);
                query.append(p1);
                subject.append(p2);
                aligner.align(query, subject);
                score = aligner.get_similarity();
                p_score = aligner.get_match_score(query, query);//5 * query.length();
            }
            return score * 100.0 / p_score;
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
                        components.put(component);
                        num_components.getAndIncrement();
                        component = new LinkedList();
                    }
                } 
            // Signifies the end of the components queue for all the threads    
                for (i = 0; i < THREADS; ++i)
                   components.put(new LinkedList());
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
                        homology_groups_list.put(component);
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
        void break_component(LinkedList<Node> component, String pangenome_path) throws InterruptedException{
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
            for( infl = MCL_INFLATION; infl < 30; infl += 0.5){
                command = "mcl " + graph_path + " --abc -I " + infl + " -o " + clusters_path;
                if(executeCommand_for(command, wating_time))
                    break;
                if((tmp_file = new File(clusters_path)).exists())
                    tmp_file.delete();
                wating_time += time;
            }
            if (infl >= 30 ){
                System.err.println("Failed to split group ID = " + component.getFirst().getId());
                homology_groups_list.put(component);
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
                                homology_groups_list.put(homology_group);
                            } else { // if is a singleton group
                                singletons_group.add(graphDb.getNodeById(Long.parseLong(fields[0])));
                            }
                        }
                    // put all singletons of the component in one homology groups    
                        if (singletons_group.size() > 0){
                                homology_groups_list.put(singletons_group);
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
                            similarity -= MIN_PROTEIN_IDENTITY;
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
                        phylogeny_distance[i][j] = MIN_PROTEIN_IDENTITY;
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
            int i, trs, p = 0, chunk = num_proteins > 40 ? num_proteins / 40 : 1;
            int[] copy_number = new int[num_genomes + 1];
            LinkedList<Node> homology_group;
            Node homology_node, protein_node;
            BufferedWriter homology_file;
            try{
                try{
                    homology_file = new BufferedWriter(new FileWriter(pangenome_path + "/pantools_homology_groups.txt"));
                    while (p < num_proteins) {
                        try(Transaction tx = graphDb.beginTx()){
                            for (trs = 0; trs < MAX_TRANSACTION_SIZE && p < num_proteins; ++trs){
                                homology_group = homology_groups_list.take();
                                ++num_groups;
                                homology_node = graphDb.createNode(homology_group_label); 
                                homology_node.setProperty("num_members", homology_group.size());
                                homology_file.write(Long.toString(homology_node.getId()) + ":");
                                while (!homology_group.isEmpty()){
                                    protein_node = homology_group.remove();
                                    homology_node.createRelationshipTo(protein_node, RelTypes.has_homolog);
                                    ++copy_number[(int)protein_node.getProperty("genome")];
                                    if (protein_node.hasProperty("protein_ID") )
                                        homology_file.write(" " + ((String)protein_node.getProperty("protein_ID")).replace(' ', '_'));
                                    else
                                        homology_file.write(" " + protein_node.getId());
                                    if (p % chunk == 0)
                                        System.out.print("|");
                                    ++p;
                                }
                                homology_file.write("\n");
                                homology_node.setProperty("copy_number", copy_number);
                                for (i=0; i < copy_number.length; ++i)
                                    copy_number[i] = 0;
                            }
                            tx.success();
                        }
                    }
                    homology_file.close();
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
    public void initialize_panproteome(){
        String file_path, file_type, line, protein_ID;
        StringBuilder protein = new StringBuilder();
        Node protein_node = null, panproteome;
        int trsc, num_proteins = 0, genome;
        String[] fields;
    // If a database folder is already exist in the specified path, removes all the content of it.    
        File theDir = new File(PATH_TO_THE_PANGENOME_DATABASE);
        if (theDir.exists()) {
            try {
                FileUtils.deleteRecursively(new File(PATH_TO_THE_PANGENOME_DATABASE));
            } catch (IOException ioe) {
                System.out.println("Failed to delete the database " + PATH_TO_THE_PANGENOME_DATABASE);
                System.exit(1);  
            }
        } else {
            try {
                theDir.mkdir();
            } catch (SecurityException se) {
                System.out.println("Failed to create directory " + PATH_TO_THE_PANGENOME_DATABASE);
                System.exit(1);
            }
        }
        if (PATH_TO_THE_PROTEOMES_FILE == null){
            System.out.println("PATH_TO_THE_PROTEOMES_FILE is empty.");
            System.exit(1);
        }  
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        startTime = System.currentTimeMillis();
        try(BufferedReader protein_paths = new BufferedReader(new FileReader(PATH_TO_THE_PROTEOMES_FILE))) {
            for (genome = 1; protein_paths.ready(); ++genome){
                file_path = protein_paths.readLine().trim();
                if (file_path.equals("")) // if line is empty
                    continue;
                fields = file_path.split("\\.");
                file_type = fields[fields.length - 1].toLowerCase();
                if (file_type.equals("fasta") || file_type.equals("faa")){
                    BufferedReader in = new BufferedReader(new FileReader(file_path));
                // skip lines till get to the first id line   
                    do{
                        line = in.readLine().trim();
                    } while(line.equals(""));
                    protein_ID = line.substring(1);
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
                // For the last protein    
                    try (Transaction tx = graphDb.beginTx()) {
                        protein_node = graphDb.createNode(mRNA_label);
                        protein_node.setProperty("protein_ID", protein_ID);
                        protein_node.setProperty("protein", protein.toString());
                        protein_node.setProperty("protein_length", protein.length());
                        protein_node.setProperty("genome",genome);
                        protein.setLength(0);
                        tx.success();
                    }
                } else {
                    System.out.println(file_path + " does not have a valid extention (fasta, faa)");
                    System.exit(1);
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
    public void group() {
        startTime = System.currentTimeMillis();
        if (! new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH).exists()) {
            System.out.println("No database found in " + PATH_TO_THE_PANGENOME_DATABASE);
            System.exit(1);
        }
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        proteins = new LinkedBlockingQueue<>();
        intersections = new LinkedBlockingQueue<>((int)(heapSize/2/100));
        similarities = new LinkedBlockingQueue<>((int)(heapSize/2/100));
        num_intersections = new AtomicInteger(0);
        num_similarities = new AtomicInteger(0);
        num_components = new AtomicInteger(0);
        similarity_bars = new AtomicInteger(0);
        try(Transaction tx = graphDb.beginTx()){
            pangenome_node = graphDb.findNodes(pangenome_label).next();
            num_genomes = (int)pangenome_node.getProperty("num_genomes");
            num_proteins = (int)pangenome_node.getProperty("num_proteins");
            tx.success();
        }
        if (THREADS < 3)
            THREADS = 3;
        MAX_KMER_FREQ = num_proteins / 1000 + 50 * num_genomes; //  because of probability p and copy number of 50
        num_hexamers = 0;
        System.out.println("Proteins = " + num_proteins);
        try{
            ExecutorService es = Executors.newFixedThreadPool(2);
            es.execute(new Generate_proteins());
            es.execute(new count_kmers(num_proteins));
            es.shutdown();
            es.awaitTermination(10, TimeUnit.DAYS);        
        } catch (InterruptedException e){
            
        }
        System.out.println("Kmerizing proteins :");
        System.out.print("0 ......................................... 100\n  "); 
        try{
            ExecutorService es = Executors.newFixedThreadPool(2);
            es.execute(new Generate_proteins());
            es.execute(new Kmerize_proteins(num_proteins));
            es.shutdown();
            es.awaitTermination(10, TimeUnit.DAYS);        
        } catch (InterruptedException e){
            
        }
        System.out.println("\nKmers = " + num_hexamers + "\n");
        //System.exit(1);
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
        System.out.println("Finding intersections:");
        System.out.print("0 ......................................... 100\n  "); 
        try{
            ExecutorService es = Executors.newFixedThreadPool(THREADS);
            es.execute(new Generate_proteins());
            es.execute(new Find_intersections(num_proteins));
            for(int i = 1; i <= THREADS - 2; i++)
                es.execute(new Find_similarities());;
            es.execute(new Write_similarities());
            es.shutdown();
            es.awaitTermination(10, TimeUnit.DAYS);        
        } catch (InterruptedException e){
            
        }
        while (similarity_bars.intValue() < 41){
            System.out.print("|");
            similarity_bars.getAndIncrement();
        }
        System.out.println("\nSimilarities = " + num_similarities.intValue() + "\n");
    //Free a considerable space    
        for (int i = 0; i < MAX_KMERS_NUM; ++i)
            kmers_proteins_list[i] = null;
        kmers_proteins_list = null;
        kmer_frequencies = null;
        intersections = null;
        similarities = null;
        System.gc();

        System.out.println("Building homology groups : ");
        System.out.print("0 ......................................... 100\n  "); 
        components = new LinkedBlockingQueue<>((int)(heapSize/200));
        homology_groups_list = new LinkedBlockingQueue<>((int)(heapSize/200));
        try{
            ExecutorService es = Executors.newFixedThreadPool(2);
            es.execute(new Generate_proteins());
            es.execute(new build_similarity_components(PATH_TO_THE_PANGENOME_DATABASE, num_proteins));
            for(int i = 1; i <= THREADS - 2; i++)
                es.execute(new build_homology_groups(PATH_TO_THE_PANGENOME_DATABASE, num_proteins));
            es.execute(new write_homology_groups(PATH_TO_THE_PANGENOME_DATABASE, num_proteins));
            es.shutdown();
            es.awaitTermination(10, TimeUnit.DAYS);        
        } catch (InterruptedException e){
            
        }
        System.out.println("\nComponents = " + num_components);
        System.out.println("Groups = " + num_groups + "\n");
        System.out.println("Database size = " + getFolderSize(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH)) + " MB");        

        File directory = new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH);
        for (File f : directory.listFiles()) {
            if (f.getName().startsWith("neostore.transaction.db.")) {
                f.delete();
            }
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
