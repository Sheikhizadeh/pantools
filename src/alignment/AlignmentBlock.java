/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package alignment;

/**
 * Implements the structure for storing the result of pairwise nucleotide alignment. 
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class AlignmentBlock {

    public StringBuilder one;
    public StringBuilder two;
    public double score;

    /**
     * The constructor.
     */
    public AlignmentBlock() {
        one = new StringBuilder();
        two = new StringBuilder();
        score = 0;
    }
    
    /**
     * clears the structure.
     */
    public void reset() {
        one.setLength(0);
        two.setLength(0);
        score = 0;
    }
    
    /**
     * Prints the alignment
     */
    public void print() {
        System.out.println("Similarity: " + score);
        System.out.println(one.reverse());
        System.out.println(two.reverse());
    }
}
