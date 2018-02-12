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
import org.junit.FixMethodOrder;
import org.junit.runners.MethodSorters;
import static pantools.Pantools.print_peak_memory;
import static pantools.Pantools.startTime;
import static pantools.Pantools.K_SIZE;
import static pantools.Pantools.PATH_TO_THE_GENOMES_FILE;
import static pantools.Pantools.PATH_TO_THE_GENOME_NUMBERS_FILE;
import static pantools.Pantools.PATH_TO_THE_PANGENOME_DATABASE;
import static pantools.Pantools.PATH_TO_THE_REGIONS_FILE;

@FixMethodOrder(MethodSorters.NAME_ASCENDING)
public class GenomeLayerTest {
    private static String test_directory;
   
    public GenomeLayerTest() {
        test_directory = System.getProperty("user.home") + "/pantools/example/";
        K_SIZE = -1;
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
     * Test of build method, of class GenomeLayer.
     */
    @Test
    public void test1_Build() {
        System.out.println("Testing build:");
        GenomeLayer instance = new GenomeLayer();
        PATH_TO_THE_GENOMES_FILE = test_directory + "sample_genomes_path.txt";
        PATH_TO_THE_PANGENOME_DATABASE = System.getProperty("user.home") + "/test/";
        instance.initialize_pangenome();
        System.out.println("Total time : " + (System.currentTimeMillis() - startTime) / 1000 + "." + (System.currentTimeMillis() - startTime) % 1000 + " seconds");
        print_peak_memory();
        System.out.println("-----------------------------------------------------------------------");
    }

    /**
     * Test of add method, of class GenomeLayer.
     */
    @Test
    public void test2_Add() {
        System.out.println("Testing add:");
        GenomeLayer instance = new GenomeLayer();
        PATH_TO_THE_GENOMES_FILE = test_directory + "sample_genomes_path_1.txt";
        PATH_TO_THE_PANGENOME_DATABASE = System.getProperty("user.home") + "/cumulative_test/";
        instance.initialize_pangenome();
        PATH_TO_THE_GENOMES_FILE = test_directory + "sample_genomes_path_2.txt";
        instance.add_genomes();
        System.out.println("Total time : " + (System.currentTimeMillis() - startTime) / 1000 + "." + (System.currentTimeMillis() - startTime) % 1000 + " seconds");
        print_peak_memory();
        System.out.println("-----------------------------------------------------------------------");
    }
  

    /**
     * Test of retrieve_regions method, of class GenomeLayer.
     */
    @Test
    public void test4_Retrieve_regions() {
        System.out.println("Testing retrieve_regions:");
        GenomeLayer instance = new GenomeLayer();
        PATH_TO_THE_REGIONS_FILE = test_directory + "sample_genomic_regions.txt";
        PATH_TO_THE_PANGENOME_DATABASE = System.getProperty("user.home") + "/test/";
        instance.retrieve_regions();
        String[] fields = PATH_TO_THE_REGIONS_FILE.split("\\/");
        if( are_the_same(PATH_TO_THE_PANGENOME_DATABASE + "/" + fields[fields.length - 1] + ".fasta", test_directory + "sample_genomic_regions.fasta") )
            System.out.println("Regions extracted properly.");
        else
            System.out.println("Test failed.");
        System.out.println("Total time : " + (System.currentTimeMillis() - startTime) / 1000 + "." + (System.currentTimeMillis() - startTime) % 1000 + " seconds");
        print_peak_memory();
        System.out.println("-----------------------------------------------------------------------");
    }

    /**
     * Test of reconstruct_genomes method, of class GenomeLayer.
     */     
    @Test
    public void test5_Reconstruct_genomes() {
        System.out.println("reconstruct_genomes");
        GenomeLayer instance = new GenomeLayer();
        PATH_TO_THE_GENOME_NUMBERS_FILE = test_directory + "sample_genome_numbers.txt";
        PATH_TO_THE_PANGENOME_DATABASE = System.getProperty("user.home") + "/test/";
        instance.retrieve_genomes();
        String genome_file;
        try {
            PATH_TO_THE_GENOMES_FILE = test_directory + "sample_genomes_path.txt";
            BufferedReader in = new BufferedReader(new FileReader(PATH_TO_THE_GENOMES_FILE));
            for (int i = 1 ; in.ready() ; ++i) {
                genome_file = in.readLine().trim();
                if (genome_file.equals("")) {
                    continue;
                }
                if( are_the_same(PATH_TO_THE_PANGENOME_DATABASE + "Genome_" + i + ".fasta", genome_file) )
                    System.out.println("Genome " + i +" extracted properly.");
                else
                    System.out.println("Genome " + i +" extraction failed.");
            }
            in.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
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
