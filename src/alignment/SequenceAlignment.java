package alignment;

/**
 * Implements required functionalities for nucleotide pairwise alignment 
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class SequenceAlignment {

    private double match[][];
    private double matrix[][];
    private double up[][];
    private double left[][];
    public int alignment_length;
    public static int MAX_LENGTH = 9000;
    public static int MATCH_SCORE = 4;
    public static int MISMATCH_SCORE = -2;
    public static int GAP_OPEN = -2;
    public static int GAP_EXT = -1;
    public static double THRESHOLD = 0.75;
    
    /**
     * The constructor of the class
     */
    public SequenceAlignment() {
    // initialize matrixes
        matrix = new double[MAX_LENGTH+1][MAX_LENGTH+1];
        up = new double[MAX_LENGTH+1][MAX_LENGTH+1];
        left = new double[MAX_LENGTH+1][MAX_LENGTH+1];
    // initialize match matrix
        match = new double[256][256];
        for (int i = 0; i < 256; i++) {
            for (int j = 0; j < 256; j++) {
                if (i == j) 
                    match[i][j] = MATCH_SCORE;
                else 
                    match[i][j] = MISMATCH_SCORE; 
            }
        }
    // Degenerate bases are supported    
        match['A']['D'] = MATCH_SCORE;
        match['A']['H'] = MATCH_SCORE;
        match['A']['M'] = MATCH_SCORE;
        match['A']['N'] = MATCH_SCORE;
        match['A']['R'] = MATCH_SCORE;
        match['A']['V'] = MATCH_SCORE;
        match['A']['W'] = MATCH_SCORE;
        match['C']['B'] = MATCH_SCORE;
        match['C']['H'] = MATCH_SCORE;
        match['C']['M'] = MATCH_SCORE;
        match['C']['N'] = MATCH_SCORE;
        match['C']['S'] = MATCH_SCORE;
        match['C']['V'] = MATCH_SCORE;
        match['C']['Y'] = MATCH_SCORE;
        match['G']['B'] = MATCH_SCORE;
        match['G']['D'] = MATCH_SCORE;
        match['G']['K'] = MATCH_SCORE;
        match['G']['N'] = MATCH_SCORE;
        match['G']['R'] = MATCH_SCORE;
        match['G']['S'] = MATCH_SCORE;
        match['G']['V'] = MATCH_SCORE;
        match['T']['B'] = MATCH_SCORE;
        match['T']['D'] = MATCH_SCORE;
        match['T']['H'] = MATCH_SCORE;
        match['T']['K'] = MATCH_SCORE;
        match['T']['N'] = MATCH_SCORE;
        match['T']['W'] = MATCH_SCORE;
        match['T']['Y'] = MATCH_SCORE;
        match['D']['A'] = MATCH_SCORE;
        match['H']['A'] = MATCH_SCORE;
        match['M']['A'] = MATCH_SCORE;
        match['N']['A'] = MATCH_SCORE;
        match['R']['A'] = MATCH_SCORE;
        match['V']['A'] = MATCH_SCORE;
        match['W']['A'] = MATCH_SCORE;
        match['B']['C'] = MATCH_SCORE;
        match['H']['C'] = MATCH_SCORE;
        match['M']['C'] = MATCH_SCORE;
        match['N']['C'] = MATCH_SCORE;
        match['S']['C'] = MATCH_SCORE;
        match['V']['C'] = MATCH_SCORE;
        match['Y']['C'] = MATCH_SCORE;
        match['B']['G'] = MATCH_SCORE;
        match['D']['G'] = MATCH_SCORE;
        match['K']['G'] = MATCH_SCORE;
        match['N']['G'] = MATCH_SCORE;
        match['R']['G'] = MATCH_SCORE;
        match['S']['G'] = MATCH_SCORE;
        match['V']['G'] = MATCH_SCORE;
        match['B']['T'] = MATCH_SCORE;
        match['D']['T'] = MATCH_SCORE;
        match['H']['T'] = MATCH_SCORE;
        match['K']['T'] = MATCH_SCORE;
        match['N']['T'] = MATCH_SCORE;
        match['W']['T'] = MATCH_SCORE;
        match['Y']['T'] = MATCH_SCORE;
    }

    /**
     * Generates the alignment of two nucleotide sequences.
     * 
     * @param s1 First sequence
     * @param s2 Second sequence
     * @return The aligned sequences
     */
    public AlignmentBlock get_alignment(String s1, String s2) {
        int i, j;
        int m = Math.min(s1.length() - 1, MAX_LENGTH), n = Math.min(s2.length() - 1, MAX_LENGTH);
        AlignmentBlock alignment = new AlignmentBlock();
        alignment.score = get_similarity(s1, s2);
        i = m;
        j = n;
        while (i > 0 && j > 0) {
            if (matrix[i][j] == up[i][j]) {
                alignment.one.append( s1.charAt(i) );
                alignment.two.append( '-' );
                i = i - 1;
            } else if (matrix[i][j] == left[i][j]) {
                alignment.one.append( '-' );
                alignment.two.append( s2.charAt(j) );
                j = j - 1;
            } else {
                alignment.one.append( s1.charAt(i) );
                alignment.two.append( s2.charAt(j) );
                i = i - 1;
                j = j - 1;
            }
        }
        return alignment;
    }

    /**
     * Calculates the similarity score of two nucleotide sequences (score is not greater than 1).
     * @param s1 First sequence
     * @param s2 Second sequence
     * @return ThE similarity score
     */
    public double get_similarity(String s1, String s2) {
        int i, j;
        int m = Math.min(s1.length() - 1, MAX_LENGTH), n = Math.min(s2.length() - 1, MAX_LENGTH);
        for (i = 1; i <= m; i++) {
            left[i][0] = -1000;
            matrix[i][0] = GAP_OPEN + i*GAP_EXT;
        }        
        for (j = 1; j <= n; j++) {
            up[0][j] = -1000;
            matrix[0][j] = GAP_OPEN + j*GAP_EXT;
        }   
        for (i = 1; i <= m; i++) {
            for (j = 1; j <= n; j++) {
                up[i][j] = Math.max( up[i-1][j] + GAP_EXT , matrix[i-1][j]+GAP_OPEN + GAP_EXT);
                left[i][j] = Math.max( left[i][j-1] + GAP_EXT , matrix[i][j-1]+GAP_OPEN + GAP_EXT);
                matrix[i][j] = Math.max( match[s1.charAt(i)][s2.charAt(j)] + matrix[i - 1][j - 1] , Math.max( up[i][j] , left[i][j]) );
            }
        }
        return matrix[m][n] / (Math.max(m, n) * MATCH_SCORE);
    }
}
