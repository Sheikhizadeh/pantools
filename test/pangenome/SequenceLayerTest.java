/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import genome.SequenceDatabase;
import index.IndexPointer;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.FixMethodOrder;
import org.junit.runners.MethodSorters;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import pantools.Pantools;
import static pantools.Pantools.print_peak_memory;
import static pantools.Pantools.startTime;

@FixMethodOrder(MethodSorters.NAME_ASCENDING)
public class SequenceLayerTest {
    private static String input_path;
    private static String output_path;
    private static String region_records_file;
    private static String regions_file;
   
    public SequenceLayerTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        input_path = "/home/sheik005/yeast";
        output_path = "/dev/shm/sheik005/yeast";
        region_records_file = "/03.regions";
        regions_file = "/regions_03.fasta";
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
     * Test of build method, of class SequenceLayer.
     
    @Test
    public void test1_Build() {
        System.out.println("Testing build:");
        SequenceLayer instance = new SequenceLayer();
        instance.build(input_path + "/03.txt", output_path + "/03");
        System.out.println("Total time : " + (System.currentTimeMillis() - startTime) / 1000 + "." + (System.currentTimeMillis() - startTime) % 1000 + " seconds");
        print_peak_memory();
        System.out.println("-----------------------------------------------------------------------");
    }

    /**
     * Test of add method, of class SequenceLayer.
     
    @Test
    public void test2_Add() {
        System.out.println("Testing add:");
        SequenceLayer instance = new SequenceLayer();
        instance.build(input_path + "/01.txt", output_path + "/03_cumulative");
        instance.add(input_path + "/2and3.txt", output_path + "/03_cumulative");
        System.out.println("Total time : " + (System.currentTimeMillis() - startTime) / 1000 + "." + (System.currentTimeMillis() - startTime) % 1000 + " seconds");
        print_peak_memory();
        System.out.println("-----------------------------------------------------------------------");
    }
  

    /**
     * Test of retrieve_regions method, of class SequenceLayer.
     */
    @Test
    public void test4_Retrieve_regions() {
        System.out.println("Testing retrieve_regions:");
        GenomeLayer instance = new GenomeLayer();
        instance.retrieve_regions(input_path + region_records_file, output_path + "/03");
        if( are_the_same(input_path + region_records_file + ".fasta", regions_file) )
            System.out.println("Regions extracted properly.");
        else
            System.out.println("Test failed.");
        print_peak_memory();
        System.out.println("-----------------------------------------------------------------------");
    }

    /**
     * Test of reconstruct_genomes method, of class SequenceLayer.
     */     
    @Test
    public void test5_Reconstruct_genomes() {
        System.out.println("reconstruct_genomes");
        GenomeLayer instance = new GenomeLayer();
        instance.retrieve_genomes("all",output_path + "/03");
        String genome_file;
        try {
            BufferedReader in = new BufferedReader(new FileReader(input_path + "/03.txt"));
            for (int i = 1 ; in.ready() ; ++i) {
                genome_file = in.readLine().trim();
            if( are_the_same(output_path + "/03/genome_" + i + ".fasta", genome_file) )
                System.out.println("Genome " + i +" extracted properly.");
            else
                System.out.println("Genome " + i +" extraction failed.");
            }
            in.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
        print_peak_memory();
        System.out.println("-----------------------------------------------------------------------");
    }
    
    public boolean are_the_same(String file_name1, String file_name2) {
        BufferedReader in1, in2;
        try {
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
        }
        return true;
    }

    
    
    /**
     * Test of retrieve_genes method, of class SequenceLayer.
     
    @Test
    public void testRetrieve_genes() {
        System.out.println("retrieve_genes");
        SequenceLayer instance = new SequenceLayer();
        //instance.retrieve_genes(annotation_records_file);
        // TODO review the generated test code and remove the default call to fail.
        //fail("The test case is a prototype.");
    }
    */
}
