/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package alignment;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import static pantools.Pantools.GAP_EXT;
import static pantools.Pantools.GAP_OPEN;

/**
 *
 * @author sheik005
 */
public class LocalSequenceAlignmentTest {
    
    public LocalSequenceAlignmentTest() {
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

    @Test
    public void test(){
        int i, j, m, n;
        LocalSequenceAlignment instance = new LocalSequenceAlignment(GAP_OPEN, GAP_EXT, 1000, 0, 'N');
        StringBuilder seq1, seq2;
        seq2 = new StringBuilder("AGTAGCAAATTAAGAGAACAAAA");
        seq1 = new StringBuilder(  "AGGCAAATTAAGAGAACAA");
        instance.align(seq1, seq2);
    }
    
}
