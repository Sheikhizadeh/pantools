/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static pantools.Pantools.print_peak_memory;
import static pantools.Pantools.startTime;
import static pantools.Pantools.PATH_TO_THE_PANGENOME_DATABASE;

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
        PATH_TO_THE_PANGENOME_DATABASE = System.getProperty("user.home") + "/test";
        String[] args = new String[]{"", PATH_TO_THE_PANGENOME_DATABASE, "-d", "8"};
        ProteomeLayer instance = new ProteomeLayer();
        instance.group(args);
        System.out.println("Total time : " + (System.currentTimeMillis() - startTime) / 1000 + "." + (System.currentTimeMillis() - startTime) % 1000 + " seconds");
        print_peak_memory();
        System.out.println("-----------------------------------------------------------------------");
    }

}
