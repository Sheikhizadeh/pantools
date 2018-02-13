/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.neo4j.graphdb.Label;
import static org.neo4j.graphdb.Label.label;
import static pantools.Pantools.print_peak_memory;
import static pantools.Pantools.startTime;
import static pantools.Pantools.PATH_TO_THE_PANGENOME_DATABASE;
import static pantools.Pantools.PATH_TO_THE_ANNOTATIONS_FILE;
import static pantools.Pantools.PATH_TO_THE_GENOME_NUMBERS_FILE;
import static pantools.Pantools.labels;

/**
 *
 * @author sheik005
 */
public class AnnotationLayerTest {
    private static String test_directory;
    
    public AnnotationLayerTest() {
        test_directory = System.getProperty("user.home") + "/pantools/example/";
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of retrieve_genes method, of class AnnotationLayer.
     */
    @Test
    public void testRetrieve_genes() {
        System.out.println("retrieve_features");
        labels = new HashMap<String,Label>();
        String[] label_strings = new String[]{
        "pangenome", "genome","sequence","nucleotide","degenerate",
        "annotation","variation","gene","coding_gene", "mRNA", 
        "tRNA", "rRNA", "CDS", "exon", "intron", "feature", 
        "broken_protein", "homology_group", "low_complexity"};        
        for (int i = 0; i < label_strings.length; ++i)
            labels.put(label_strings[i], label(label_strings[i]));
        AnnotationLayer instance = new AnnotationLayer();
        PATH_TO_THE_ANNOTATIONS_FILE = test_directory + "sample_annotations_path.txt";
        PATH_TO_THE_GENOME_NUMBERS_FILE = test_directory + "sample_genome_numbers.txt";
        PATH_TO_THE_PANGENOME_DATABASE = System.getProperty("user.home") + "/test/";
        instance.add_annotaions();
        instance.retrieve_feature();
        if( are_the_same("genes.1.fasta", test_directory + "sample_genes.1.fasta") &&
            are_the_same("genes.2.fasta", test_directory + "sample_genes.2.fasta"))
            System.out.println("Genes extracted properly.");
        else
            System.out.println("Test failed.");
        System.out.println("Total time : " + (System.currentTimeMillis() - startTime) / 1000 + "." + (System.currentTimeMillis() - startTime) % 1000 + " seconds");
        print_peak_memory();
        System.out.println("-----------------------------------------------------------------------");
    }
    
    public boolean are_the_same(String file_name1, String file_name2) {
        BufferedReader in1, in2;
        try {
            System.out.println("Comparing " + file_name1 + " and " + file_name2 + "...");
            in1 = new BufferedReader(new FileReader(file_name1));
            in2 = new BufferedReader(new FileReader(file_name2));
            do
                if( !in1.readLine().toUpperCase().trim().equals(in2.readLine().toUpperCase().trim()) )
                        return false;
            while (in1.ready() && in2.ready()) ;
            in1.close();
            in2.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            return false;
        }
        return true;
    }
    
}
