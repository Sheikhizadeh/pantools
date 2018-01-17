package alignment;

import java.util.LinkedList;
import java.util.Stack;

/**
 * Implements required functionalities for nucleotide pairwise alignment 
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class SequenceAlignment {

    private int match[][];
    private long matrix[][];
    private long up[][];
    private long left[][];
    private StringBuilder seq1;
    private StringBuilder seq2;
    private StringBuilder one;
    private StringBuilder two;
    StringBuilder cigar;
    Stack<Character> operation_stack;
    Stack<Integer> count_stack;
    private  long score;
    private int equals;
    public static int GAP_OPEN;
    public static int GAP_EXT;
    public static int MAX_LENGTH;
    public static int GAPS_INTERVAL;
    
    /**
     * The constructor of the class
     * @param gap_open
     * @param gap_ext
     * @param max_length
     */
    public SequenceAlignment(int o, int e, int m, int g) {
        GAP_OPEN = o;
        GAP_EXT = e;
        MAX_LENGTH = m;
        GAPS_INTERVAL = g;
        one = new StringBuilder();
        two = new StringBuilder();
    // initialize matrixes
        matrix = new long[MAX_LENGTH + 1][MAX_LENGTH + 1];
        up = new long[MAX_LENGTH + 1][MAX_LENGTH + 1];
        left = new long[MAX_LENGTH + 1][MAX_LENGTH + 1];
        cigar = new StringBuilder();
        operation_stack = new Stack();
        count_stack = new Stack();
        equals = 0;
        for (int i = 0; i <= MAX_LENGTH; i++) {
            for (int j = 0; j <= MAX_LENGTH; j++) {
                up[i][j] = Integer.MIN_VALUE;
                left[i][j] = Integer.MIN_VALUE;
                matrix[i][j] = (i == 0 || j == 0 ? 0 : Integer.MIN_VALUE);//GAP_OPEN + (i + j) * GAP_EXT
            }
        }    
    // initialize similarity matrix NUCC.4.4
        match = new int[256][256];

        match['A']['A'] = 5;
        match['A']['T'] = -4;
        match['A']['G'] = -4;
        match['A']['C'] = -4;
        match['A']['S'] = -4;
        match['A']['W'] = 1;
        match['A']['R'] = 1;
        match['A']['Y'] = -4;
        match['A']['K'] = -4;
        match['A']['M'] = 1;
        match['A']['B'] = -4;
        match['A']['V'] = -1;
        match['A']['H'] = -1;
        match['A']['D'] = -1;
        match['A']['N'] = -2;
        match['T']['A'] = -4;
        match['T']['T'] = 5;
        match['T']['G'] = -4;
        match['T']['C'] = -4;
        match['T']['S'] = -4;
        match['T']['W'] = 1;
        match['T']['R'] = -4;
        match['T']['Y'] = 1;
        match['T']['K'] = 1;
        match['T']['M'] = -4;
        match['T']['B'] = -1;
        match['T']['V'] = -4;
        match['T']['H'] = -1;
        match['T']['D'] = -1;
        match['T']['N'] = -2;
        match['G']['A'] = -4;
        match['G']['T'] = -4;
        match['G']['G'] = 5;
        match['G']['C'] = -4;
        match['G']['S'] = 1;
        match['G']['W'] = -4;
        match['G']['R'] = 1;
        match['G']['Y'] = -4;
        match['G']['K'] = 1;
        match['G']['M'] = -4;
        match['G']['B'] = -1;
        match['G']['V'] = -1;
        match['G']['H'] = -4;
        match['G']['D'] = -1;
        match['G']['N'] = -2;
        match['C']['A'] = -4;
        match['C']['T'] = -4;
        match['C']['G'] = -4;
        match['C']['C'] = 5;
        match['C']['S'] = 1;
        match['C']['W'] = -4;
        match['C']['R'] = -4;
        match['C']['Y'] = 1;
        match['C']['K'] = -4;
        match['C']['M'] = 1;
        match['C']['B'] = -1;
        match['C']['V'] = -1;
        match['C']['H'] = -1;
        match['C']['D'] = -4;
        match['C']['N'] = -2;
        match['S']['A'] = -4;
        match['S']['T'] = -4;
        match['S']['G'] = 1;
        match['S']['C'] = 1;
        match['S']['S'] = -1;
        match['S']['W'] = -4;
        match['S']['R'] = -2;
        match['S']['Y'] = -2;
        match['S']['K'] = -2;
        match['S']['M'] = -2;
        match['S']['B'] = -1;
        match['S']['V'] = -1;
        match['S']['H'] = -3;
        match['S']['D'] = -3;
        match['S']['N'] = -1;
        match['W']['A'] = 1;
        match['W']['T'] = 1;
        match['W']['G'] = -4;
        match['W']['C'] = -4;
        match['W']['S'] = -4;
        match['W']['W'] = -1;
        match['W']['R'] = -2;
        match['W']['Y'] = -2;
        match['W']['K'] = -2;
        match['W']['M'] = -2;
        match['W']['B'] = -3;
        match['W']['V'] = -3;
        match['W']['H'] = -1;
        match['W']['D'] = -1;
        match['W']['N'] = -1;
        match['R']['A'] = 1;
        match['R']['T'] = -4;
        match['R']['G'] = 1;
        match['R']['C'] = -4;
        match['R']['S'] = -2;
        match['R']['W'] = -2;
        match['R']['R'] = -1;
        match['R']['Y'] = -4;
        match['R']['K'] = -2;
        match['R']['M'] = -2;
        match['R']['B'] = -3;
        match['R']['V'] = -1;
        match['R']['H'] = -3;
        match['R']['D'] = -1;
        match['R']['N'] = -1;
        match['Y']['A'] = -4;
        match['Y']['T'] = 1;
        match['Y']['G'] = -4;
        match['Y']['C'] = 1;
        match['Y']['S'] = -2;
        match['Y']['W'] = -2;
        match['Y']['R'] = -4;
        match['Y']['Y'] = -1;
        match['Y']['K'] = -2;
        match['Y']['M'] = -2;
        match['Y']['B'] = -1;
        match['Y']['V'] = -3;
        match['Y']['H'] = -1;
        match['Y']['D'] = -3;
        match['Y']['N'] = -1;
        match['K']['A'] = -4;
        match['K']['T'] = 1;
        match['K']['G'] = 1;
        match['K']['C'] = -4;
        match['K']['S'] = -2;
        match['K']['W'] = -2;
        match['K']['R'] = -2;
        match['K']['Y'] = -2;
        match['K']['K'] = -1;
        match['K']['M'] = -4;
        match['K']['B'] = -1;
        match['K']['V'] = -3;
        match['K']['H'] = -3;
        match['K']['D'] = -1;
        match['K']['N'] = -1;
        match['M']['A'] = 1;
        match['M']['T'] = -4;
        match['M']['G'] = -4;
        match['M']['C'] = 1;
        match['M']['S'] = -2;
        match['M']['W'] = -2;
        match['M']['R'] = -2;
        match['M']['Y'] = -2;
        match['M']['K'] = -4;
        match['M']['M'] = -1;
        match['M']['B'] = -3;
        match['M']['V'] = -1;
        match['M']['H'] = -1;
        match['M']['D'] = -3;
        match['M']['N'] = -1;
        match['B']['A'] = -4;
        match['B']['T'] = -1;
        match['B']['G'] = -1;
        match['B']['C'] = -1;
        match['B']['S'] = -1;
        match['B']['W'] = -3;
        match['B']['R'] = -3;
        match['B']['Y'] = -1;
        match['B']['K'] = -1;
        match['B']['M'] = -3;
        match['B']['B'] = -1;
        match['B']['V'] = -2;
        match['B']['H'] = -2;
        match['B']['D'] = -2;
        match['B']['N'] = -1;
        match['V']['A'] = -1;
        match['V']['T'] = -4;
        match['V']['G'] = -1;
        match['V']['C'] = -1;
        match['V']['S'] = -1;
        match['V']['W'] = -3;
        match['V']['R'] = -1;
        match['V']['Y'] = -3;
        match['V']['K'] = -3;
        match['V']['M'] = -1;
        match['V']['B'] = -2;
        match['V']['V'] = -1;
        match['V']['H'] = -2;
        match['V']['D'] = -2;
        match['V']['N'] = -1;
        match['H']['A'] = -1;
        match['H']['T'] = -1;
        match['H']['G'] = -4;
        match['H']['C'] = -1;
        match['H']['S'] = -3;
        match['H']['W'] = -1;
        match['H']['R'] = -3;
        match['H']['Y'] = -1;
        match['H']['K'] = -3;
        match['H']['M'] = -1;
        match['H']['B'] = -2;
        match['H']['V'] = -2;
        match['H']['H'] = -1;
        match['H']['D'] = -2;
        match['H']['N'] = -1;
        match['D']['A'] = -1;
        match['D']['T'] = -1;
        match['D']['G'] = -1;
        match['D']['C'] = -4;
        match['D']['S'] = -3;
        match['D']['W'] = -1;
        match['D']['R'] = -1;
        match['D']['Y'] = -3;
        match['D']['K'] = -1;
        match['D']['M'] = -3;
        match['D']['B'] = -2;
        match['D']['V'] = -2;
        match['D']['H'] = -2;
        match['D']['D'] = -1;
        match['D']['N'] = -1;
        match['N']['A'] = -2;
        match['N']['T'] = -2;
        match['N']['G'] = -2;
        match['N']['C'] = -2;
        match['N']['S'] = -1;
        match['N']['W'] = -1;
        match['N']['R'] = -1;
        match['N']['Y'] = -1;
        match['N']['K'] = -1;
        match['N']['M'] = -1;
        match['N']['B'] = -1;
        match['N']['V'] = -1;
        match['N']['H'] = -1;
        match['N']['D'] = -1;
        match['N']['N'] = -1;
        /*seq1 = new StringBuilder("ACGTAACGTT");
        seq2 = new StringBuilder("ACGTAACGTT");
        this.align(seq1, seq2);
        System.out.println(this.get_alignment());
        System.out.println(this.get_cigar());
        System.out.println(this.get_score());
        System.exit(0);*/
    }

    /**
     * Calculates the similarity score of two nucleotide sequences (score is not greater than 1).
     * @param s1 First sequence
     * @param s2 Second sequence
     * @return The similarity score
     */
    public long align(StringBuilder s1, StringBuilder s2, int anchor_pos, int K) {
        int i, j;
        int m, n;
        seq1 = s1;
        seq2 = s2;
        m = seq1.length();
        n = seq2.length();
        for (i = 1; i <= m; i++) {
                //System.out.print(seq1.charAt(i-1) + " ");
            for (j = i; j <= i + GAPS_INTERVAL; j++) {
                up[i][j] = Math.max( up[i-1][j] + GAP_EXT , matrix[i-1][j]+GAP_OPEN + GAP_EXT);
                left[i][j] = Math.max( left[i][j-1] + GAP_EXT , matrix[i][j-1]+GAP_OPEN + GAP_EXT);
                matrix[i][j] = Math.max( match[seq1.charAt(i-1)][seq2.charAt(j-1)] + matrix[i - 1][j - 1] , Math.max( up[i][j] , left[i][j]) );
                //System.out.print(String.format("%7d", matrix[i][j] ));
            }
            //System.out.println();
        }
        return score = matrix[m][n];
    }

    /**
     * Generates the alignment of two nucleotide sequences.
     * 
     * @param s1 First sequence
     * @param s2 Second sequence
     * @return The aligned sequences
     */
    public String get_alignment() {
        int i, j;
        one.setLength(0);
        two.setLength(0);
        i = seq1.length();
        j = seq2.length();
        while (i > 0 && j > 0) {
            if (matrix[i][j] == up[i][j]) {
                one.append( seq1.charAt(i-1) );
                two.append( '-' );
                i = i - 1;
            } else if (matrix[i][j] == left[i][j]) {
                one.append( '-' );
                two.append( seq2.charAt(j-1) );
                j = j - 1;
            } else {
                one.append( seq1.charAt(i-1) );
                two.append( seq2.charAt(j-1) );
                i = i - 1;
                j = j - 1;
            }
        } 
        for (;i > 0; --i){
            one.append( seq1.charAt(i-1) );
            two.append( '-' );
        }
        for (;j > 0; --j){
            one.append( '-' );
            two.append( seq2.charAt(j-1) );
        }
        return one.reverse() + "\n" + two.reverse();
    }

    public int calculate_cigar(int anchor_pos) {
        int i, j, genomic_offset = 0, move_counts = 1, count;
        char curr_move, prev_move, operation;
        cigar.setLength(0);
        i = seq1.length();
        j = seq2.length();
        if (matrix[i][j] == up[i][j]){ 
            prev_move = 'I';
            i = i - 1;
        }
        else if (matrix[i][j] == left[i][j]){
            prev_move = 'D';
            j = j - 1;
        } else {
            if (seq1.charAt(i-1) == seq2.charAt(j-1))
                prev_move = '=';
            else
                prev_move = 'M';               
            i = i - 1;
            j = j - 1;
        }
        while (i > 0 && j > 0) {
            if (matrix[i][j] == up[i][j]) {
                curr_move = 'I';
                i = i - 1;
            } else if (matrix[i][j] == left[i][j]) {
                curr_move = 'D';
                j = j - 1;
            } else {
               if (seq1.charAt(i-1) == seq2.charAt(j-1))
                    curr_move = '=';
                else
                    curr_move = 'M';               
                i = i - 1;
                j = j - 1;
            }
            if (prev_move == curr_move)
                ++move_counts;
            else{
                operation_stack.push(prev_move);
                count_stack.push(move_counts);
                if (operation_stack.size() == 1 && prev_move == 'D'){ //Remove the deletion at the end of the read
                    operation_stack.pop();
                    score += (5 + count_stack.pop());// make up for the penalties
                } 
                if (prev_move == '=')
                    equals += move_counts;
                move_counts = 1;
            }
            prev_move = curr_move;
        } 
        operation_stack.push(prev_move);
        count_stack.push(move_counts);
        if (prev_move == '=')
            equals += move_counts;
        if (i > 0){
            operation_stack.push('I');
            count_stack.push(i);
        }
        // Avoid D at the start of cigar
        /*else if (j > 0){
            operation_stack.push('D');
            count_stack.push(j);
        }*/
        //System.out.println(anchor_pos);
        int loc = 0;
        while (!operation_stack.isEmpty()){
            operation = operation_stack.pop();
            count = count_stack.pop();
            cigar.append(count).append(operation);
            //System.out.println(loc + "<=?"+anchor_pos+ " && " +" "+count+operation);
            if (operation == 'D' && loc <= anchor_pos){
                genomic_offset += count;
                //System.out.println(genomic_offset);
                loc += count;
            }
            if (operation == 'I' && loc < anchor_pos){
                genomic_offset -= count;
                //System.out.println(genomic_offset);
            } 
            if (operation == '=' || operation == 'M')
                loc += count;
        }
            //System.out.println(">" + loc + " "+ genomic_offset);
        return genomic_offset;
    }

    public StringBuilder get_cigar(){
        return cigar;
    }
    

    public long get_score(){
        return score;
    }

    public int get_equals(){
        int i, j;
        equals = 0;
        i = seq1.length();
        j = seq2.length();
        while (i > 0 && j > 0) {
            if (matrix[i][j] == up[i][j]) {
                i = i - 1;
            } else if (matrix[i][j] == left[i][j]) {
                j = j - 1;
            } else {
                if (seq1.charAt(i - 1) == seq2.charAt(j - 1))
                    ++equals;
                i = i - 1;
                j = j - 1;
            }
        } 
        return equals;
    }

}
