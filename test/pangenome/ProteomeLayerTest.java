/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

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
import static pantools.Pantools.labels;

/**
 *
 * @author sheik005
 */
public class ProteomeLayerTest {
    
    public ProteomeLayerTest() {
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
     * Test of group method, of class ProteomeLayer.
     */
    @Test
    public void testGroup() {
        System.out.println("group");
        labels = new HashMap<String,Label>();
        String[] label_strings = new String[]{
        "pangenome", "genome","sequence","nucleotide","degenerate",
        "annotation","variation","gene","coding_gene", "mRNA", 
        "tRNA", "rRNA", "CDS", "exon", "intron", "feature", 
        "broken_protein", "homology_group", "low_complexity"};        
        for (int i = 0; i < label_strings.length; ++i)
            labels.put(label_strings[i], label(label_strings[i]));
        PATH_TO_THE_PANGENOME_DATABASE = System.getProperty("user.home") + "/test";
        String[] args = new String[]{"", PATH_TO_THE_PANGENOME_DATABASE, "-d", "8"};
        ProteomeLayer instance = new ProteomeLayer();
        instance.group();
        System.out.println("Total time : " + (System.currentTimeMillis() - startTime) / 1000 + "." + (System.currentTimeMillis() - startTime) % 1000 + " seconds");
        print_peak_memory();
        System.out.println("-----------------------------------------------------------------------");
    }

}
