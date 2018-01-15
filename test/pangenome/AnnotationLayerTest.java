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
import static org.junit.Assert.*;
import static pantools.Pantools.print_peak_memory;
import static pantools.Pantools.startTime;

/**
 *
 * @author sheik005
 */
public class AnnotationLayerTest {
    private static String test_directory;
    private static String database_path;
    private static String gff_paths_file;
    private static String gene_records_file;
    private static String genes_file;
    
    public AnnotationLayerTest() {
        test_directory = System.getProperty("user.home") + "/pantools/example/";
        database_path = System.getProperty("user.home") + "/test/";
        gff_paths_file = "sample_annotations_path.txt";
        gene_records_file = "sample_annotation_records.txt";
        genes_file = "sample_annotation_records.fasta";
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
        instance.add_annotaions(test_directory + gff_paths_file, database_path);
        instance.retrieve_genes(test_directory + gene_records_file, database_path);
        String [] fields = gene_records_file.split("\\/");
        if( are_the_same(database_path + "/" + fields[fields.length - 1] + ".fasta", test_directory + genes_file) )
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
