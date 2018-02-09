/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static pantools.Pantools.print_peak_memory;
import static pantools.Pantools.startTime;
import static pantools.Pantools.PATH_TO_THE_PANGENOME_DATABASE;
import static pantools.Pantools.PATH_TO_THE_ANNOTATIONS_FILE;
import static pantools.Pantools.PATH_TO_THE_GENE_RECORDS;

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
        System.out.println("retrieve_genes");
        AnnotationLayer instance = new AnnotationLayer();
        PATH_TO_THE_ANNOTATIONS_FILE = test_directory + "sample_annotations_path.txt";
        PATH_TO_THE_GENE_RECORDS = test_directory + "sample_annotation_records.txt";
        PATH_TO_THE_PANGENOME_DATABASE = System.getProperty("user.home") + "/test/";
        instance.add_annotaions();
        instance.retrieve_genes();
        //String [] fields = gene_records_file.split("\\/");//sample_annotation_records.txt
        if( are_the_same(PATH_TO_THE_PANGENOME_DATABASE + "/sample_annotation_records.txt.fasta", test_directory + "sample_annotation_records.fasta") )
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
            in1 = new BufferedReader(new FileReader(file_name1));
            in2 = new BufferedReader(new FileReader(file_name2));
            System.out.println("Comparing " + file_name1 + " and " + file_name2 + "...");
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
