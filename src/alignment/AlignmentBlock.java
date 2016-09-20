/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package alignment;

/**
 *
 * @author sheik005
 */
public class AlignmentBlock {

    public StringBuilder one;
    public StringBuilder two;
    public int score;

    public AlignmentBlock(int length) {
        one = new StringBuilder(length);
        two = new StringBuilder(length);
        score = 0;
    }

    public void reset() {
        one.setLength(0);
        two.setLength(0);
        score = 0;
    }
}
