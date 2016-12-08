package alignment;

/**
 * Implements required functionalities for nucleotide pairwise alignment 
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class ProteinAlignment {

    private double match[][];
    private double matrix[][];
    private double up[][];
    private double left[][];
    public int alignment_length;
    public int MAX_LENGTH;
    public double GAP_OPEN;
    public double GAP_EXT;
    
    /**
     * The constructor of the class
     * @param gap_open
     * @param gap_ext
     * @param max_length
     */
    public ProteinAlignment(double gap_open, double gap_ext, int max_length) {
        MAX_LENGTH = max_length;
        GAP_OPEN = gap_open;
        GAP_EXT = gap_ext;
    // initialize matrixes
        matrix = new double[MAX_LENGTH+1][MAX_LENGTH+1];
        up = new double[MAX_LENGTH+1][MAX_LENGTH+1];
        left = new double[MAX_LENGTH+1][MAX_LENGTH+1];
    // initialize similarity matrix BLOSUM62
        match = new double[256][256];
        match['A']['A'] = 4;
        match['A']['R'] = -1;
        match['A']['N'] = -2;
        match['A']['D'] = -2;
        match['A']['C'] = 0;
        match['A']['Q'] = -1;
        match['A']['E'] = -1;
        match['A']['G'] = 0;
        match['A']['H'] = -2;
        match['A']['I'] = -1;
        match['A']['L'] = -1;
        match['A']['K'] = -1;
        match['A']['M'] = -1;
        match['A']['F'] = -2;
        match['A']['P'] = -1;
        match['A']['S'] = 1;
        match['A']['T'] = 0;
        match['A']['W'] = -3;
        match['A']['Y'] = -2;
        match['A']['V'] = 0;
        match['A']['B'] = -2;
        match['A']['Z'] = -1;
        match['A']['X'] = 0;
        match['A']['*'] = -4;


        match['R']['A'] = -1;
        match['R']['R'] = 5;
        match['R']['N'] = 0;
        match['R']['D'] = -2;
        match['R']['C'] = -3;
        match['R']['Q'] = 1;
        match['R']['E'] = 0;
        match['R']['G'] = -2;
        match['R']['H'] = 0;
        match['R']['I'] = -3;
        match['R']['L'] = -2;
        match['R']['K'] = 2;
        match['R']['M'] = -1;
        match['R']['F'] = -3;
        match['R']['P'] = -2;
        match['R']['S'] = -1;
        match['R']['T'] = -1;
        match['R']['W'] = -3;
        match['R']['Y'] = -2;
        match['R']['V'] = -3;
        match['R']['B'] = -1;
        match['R']['Z'] = 0;
        match['R']['X'] = -1;
        match['R']['*'] = -4;


        match['N']['A'] = -2;
        match['N']['R'] = 0;
        match['N']['N'] = 6;
        match['N']['D'] = 1;
        match['N']['C'] = -3;
        match['N']['Q'] = 0;
        match['N']['E'] = 0;
        match['N']['G'] = 0;
        match['N']['H'] = 1;
        match['N']['I'] = -3;
        match['N']['L'] = -3;
        match['N']['K'] = 0;
        match['N']['M'] = -2;
        match['N']['F'] = -3;
        match['N']['P'] = -2;
        match['N']['S'] = 1;
        match['N']['T'] = 0;
        match['N']['W'] = -4;
        match['N']['Y'] = -2;
        match['N']['V'] = -3;
        match['N']['B'] = 3;
        match['N']['Z'] = 0;
        match['N']['X'] = -1;
        match['N']['*'] = -4;


        match['D']['A'] = -2;
        match['D']['R'] = -2;
        match['D']['N'] = 1;
        match['D']['D'] = 6;
        match['D']['C'] = -3;
        match['D']['Q'] = 0;
        match['D']['E'] = 2;
        match['D']['G'] = -1;
        match['D']['H'] = -1;
        match['D']['I'] = -3;
        match['D']['L'] = -4;
        match['D']['K'] = -1;
        match['D']['M'] = -3;
        match['D']['F'] = -3;
        match['D']['P'] = -1;
        match['D']['S'] = 0;
        match['D']['T'] = -1;
        match['D']['W'] = -4;
        match['D']['Y'] = -3;
        match['D']['V'] = -3;
        match['D']['B'] = 4;
        match['D']['Z'] = 1;
        match['D']['X'] = -1;
        match['D']['*'] = -4;


        match['C']['A'] = 0;
        match['C']['R'] = -3;
        match['C']['N'] = -3;
        match['C']['D'] = -3;
        match['C']['C'] = 9;
        match['C']['Q'] = -3;
        match['C']['E'] = -4;
        match['C']['G'] = -3;
        match['C']['H'] = -3;
        match['C']['I'] = -1;
        match['C']['L'] = -1;
        match['C']['K'] = -3;
        match['C']['M'] = -1;
        match['C']['F'] = -2;
        match['C']['P'] = -3;
        match['C']['S'] = -1;
        match['C']['T'] = -1;
        match['C']['W'] = -2;
        match['C']['Y'] = -2;
        match['C']['V'] = -1;
        match['C']['B'] = -3;
        match['C']['Z'] = -3;
        match['C']['X'] = -2;
        match['C']['*'] = -4;


        match['Q']['A'] = -1;
        match['Q']['R'] = 1;
        match['Q']['N'] = 0;
        match['Q']['D'] = 0;
        match['Q']['C'] = -3;
        match['Q']['Q'] = 5;
        match['Q']['E'] = 2;
        match['Q']['G'] = -2;
        match['Q']['H'] = 0;
        match['Q']['I'] = -3;
        match['Q']['L'] = -2;
        match['Q']['K'] = 1;
        match['Q']['M'] = 0;
        match['Q']['F'] = -3;
        match['Q']['P'] = -1;
        match['Q']['S'] = 0;
        match['Q']['T'] = -1;
        match['Q']['W'] = -2;
        match['Q']['Y'] = -1;
        match['Q']['V'] = -2;
        match['Q']['B'] = 0;
        match['Q']['Z'] = 3;
        match['Q']['X'] = -1;
        match['Q']['*'] = -4;


        match['E']['A'] = -1;
        match['E']['R'] = 0;
        match['E']['N'] = 0;
        match['E']['D'] = 2;
        match['E']['C'] = -4;
        match['E']['Q'] = 2;
        match['E']['E'] = 5;
        match['E']['G'] = -2;
        match['E']['H'] = 0;
        match['E']['I'] = -3;
        match['E']['L'] = -3;
        match['E']['K'] = 1;
        match['E']['M'] = -2;
        match['E']['F'] = -3;
        match['E']['P'] = -1;
        match['E']['S'] = 0;
        match['E']['T'] = -1;
        match['E']['W'] = -3;
        match['E']['Y'] = -2;
        match['E']['V'] = -2;
        match['E']['B'] = 1;
        match['E']['Z'] = 4;
        match['E']['X'] = -1;
        match['E']['*'] = -4;


        match['G']['A'] = 0;
        match['G']['R'] = -2;
        match['G']['N'] = 0;
        match['G']['D'] = -1;
        match['G']['C'] = -3;
        match['G']['Q'] = -2;
        match['G']['E'] = -2;
        match['G']['G'] = 6;
        match['G']['H'] = -2;
        match['G']['I'] = -4;
        match['G']['L'] = -4;
        match['G']['K'] = -2;
        match['G']['M'] = -3;
        match['G']['F'] = -3;
        match['G']['P'] = -2;
        match['G']['S'] = 0;
        match['G']['T'] = -2;
        match['G']['W'] = -2;
        match['G']['Y'] = -3;
        match['G']['V'] = -3;
        match['G']['B'] = -1;
        match['G']['Z'] = -2;
        match['G']['X'] = -1;
        match['G']['*'] = -4;


        match['H']['A'] = -2;
        match['H']['R'] = 0;
        match['H']['N'] = 1;
        match['H']['D'] = -1;
        match['H']['C'] = -3;
        match['H']['Q'] = 0;
        match['H']['E'] = 0;
        match['H']['G'] = -2;
        match['H']['H'] = 8;
        match['H']['I'] = -3;
        match['H']['L'] = -3;
        match['H']['K'] = -1;
        match['H']['M'] = -2;
        match['H']['F'] = -1;
        match['H']['P'] = -2;
        match['H']['S'] = -1;
        match['H']['T'] = -2;
        match['H']['W'] = -2;
        match['H']['Y'] = 2;
        match['H']['V'] = -3;
        match['H']['B'] = 0;
        match['H']['Z'] = 0;
        match['H']['X'] = -1;
        match['H']['*'] = -4;


        match['I']['A'] = -1;
        match['I']['R'] = -3;
        match['I']['N'] = -3;
        match['I']['D'] = -3;
        match['I']['C'] = -1;
        match['I']['Q'] = -3;
        match['I']['E'] = -3;
        match['I']['G'] = -4;
        match['I']['H'] = -3;
        match['I']['I'] = 4;
        match['I']['L'] = 2;
        match['I']['K'] = -3;
        match['I']['M'] = 1;
        match['I']['F'] = 0;
        match['I']['P'] = -3;
        match['I']['S'] = -2;
        match['I']['T'] = -1;
        match['I']['W'] = -3;
        match['I']['Y'] = -1;
        match['I']['V'] = 3;
        match['I']['B'] = -3;
        match['I']['Z'] = -3;
        match['I']['X'] = -1;
        match['I']['*'] = -4;


        match['L']['A'] = -1;
        match['L']['R'] = -2;
        match['L']['N'] = -3;
        match['L']['D'] = -4;
        match['L']['C'] = -1;
        match['L']['Q'] = -2;
        match['L']['E'] = -3;
        match['L']['G'] = -4;
        match['L']['H'] = -3;
        match['L']['I'] = 2;
        match['L']['L'] = 4;
        match['L']['K'] = -2;
        match['L']['M'] = 2;
        match['L']['F'] = 0;
        match['L']['P'] = -3;
        match['L']['S'] = -2;
        match['L']['T'] = -1;
        match['L']['W'] = -2;
        match['L']['Y'] = -1;
        match['L']['V'] = 1;
        match['L']['B'] = -4;
        match['L']['Z'] = -3;
        match['L']['X'] = -1;
        match['L']['*'] = -4;


        match['K']['A'] = -1;
        match['K']['R'] = 2;
        match['K']['N'] = 0;
        match['K']['D'] = -1;
        match['K']['C'] = -3;
        match['K']['Q'] = 1;
        match['K']['E'] = 1;
        match['K']['G'] = -2;
        match['K']['H'] = -1;
        match['K']['I'] = -3;
        match['K']['L'] = -2;
        match['K']['K'] = 5;
        match['K']['M'] = -1;
        match['K']['F'] = -3;
        match['K']['P'] = -1;
        match['K']['S'] = 0;
        match['K']['T'] = -1;
        match['K']['W'] = -3;
        match['K']['Y'] = -2;
        match['K']['V'] = -2;
        match['K']['B'] = 0;
        match['K']['Z'] = 1;
        match['K']['X'] = -1;
        match['K']['*'] = -4;


        match['M']['A'] = -1;
        match['M']['R'] = -1;
        match['M']['N'] = -2;
        match['M']['D'] = -3;
        match['M']['C'] = -1;
        match['M']['Q'] = 0;
        match['M']['E'] = -2;
        match['M']['G'] = -3;
        match['M']['H'] = -2;
        match['M']['I'] = 1;
        match['M']['L'] = 2;
        match['M']['K'] = -1;
        match['M']['M'] = 5;
        match['M']['F'] = 0;
        match['M']['P'] = -2;
        match['M']['S'] = -1;
        match['M']['T'] = -1;
        match['M']['W'] = -1;
        match['M']['Y'] = -1;
        match['M']['V'] = 1;
        match['M']['B'] = -3;
        match['M']['Z'] = -1;
        match['M']['X'] = -1;
        match['M']['*'] = -4;


        match['F']['A'] = -2;
        match['F']['R'] = -3;
        match['F']['N'] = -3;
        match['F']['D'] = -3;
        match['F']['C'] = -2;
        match['F']['Q'] = -3;
        match['F']['E'] = -3;
        match['F']['G'] = -3;
        match['F']['H'] = -1;
        match['F']['I'] = 0;
        match['F']['L'] = 0;
        match['F']['K'] = -3;
        match['F']['M'] = 0;
        match['F']['F'] = 6;
        match['F']['P'] = -4;
        match['F']['S'] = -2;
        match['F']['T'] = -2;
        match['F']['W'] = 1;
        match['F']['Y'] = 3;
        match['F']['V'] = -1;
        match['F']['B'] = -3;
        match['F']['Z'] = -3;
        match['F']['X'] = -1;
        match['F']['*'] = -4;


        match['P']['A'] = -1;
        match['P']['R'] = -2;
        match['P']['N'] = -2;
        match['P']['D'] = -1;
        match['P']['C'] = -3;
        match['P']['Q'] = -1;
        match['P']['E'] = -1;
        match['P']['G'] = -2;
        match['P']['H'] = -2;
        match['P']['I'] = -3;
        match['P']['L'] = -3;
        match['P']['K'] = -1;
        match['P']['M'] = -2;
        match['P']['F'] = -4;
        match['P']['P'] = 7;
        match['P']['S'] = -1;
        match['P']['T'] = -1;
        match['P']['W'] = -4;
        match['P']['Y'] = -3;
        match['P']['V'] = -2;
        match['P']['B'] = -2;
        match['P']['Z'] = -1;
        match['P']['X'] = -2;
        match['P']['*'] = -4;


        match['S']['A'] = 1;
        match['S']['R'] = -1;
        match['S']['N'] = 1;
        match['S']['D'] = 0;
        match['S']['C'] = -1;
        match['S']['Q'] = 0;
        match['S']['E'] = 0;
        match['S']['G'] = 0;
        match['S']['H'] = -1;
        match['S']['I'] = -2;
        match['S']['L'] = -2;
        match['S']['K'] = 0;
        match['S']['M'] = -1;
        match['S']['F'] = -2;
        match['S']['P'] = -1;
        match['S']['S'] = 4;
        match['S']['T'] = 1;
        match['S']['W'] = -3;
        match['S']['Y'] = -2;
        match['S']['V'] = -2;
        match['S']['B'] = 0;
        match['S']['Z'] = 0;
        match['S']['X'] = 0;
        match['S']['*'] = -4;


        match['T']['A'] = 0;
        match['T']['R'] = -1;
        match['T']['N'] = 0;
        match['T']['D'] = -1;
        match['T']['C'] = -1;
        match['T']['Q'] = -1;
        match['T']['E'] = -1;
        match['T']['G'] = -2;
        match['T']['H'] = -2;
        match['T']['I'] = -1;
        match['T']['L'] = -1;
        match['T']['K'] = -1;
        match['T']['M'] = -1;
        match['T']['F'] = -2;
        match['T']['P'] = -1;
        match['T']['S'] = 1;
        match['T']['T'] = 5;
        match['T']['W'] = -2;
        match['T']['Y'] = -2;
        match['T']['V'] = 0;
        match['T']['B'] = -1;
        match['T']['Z'] = -1;
        match['T']['X'] = 0;
        match['T']['*'] = -4;


        match['W']['A'] = -3;
        match['W']['R'] = -3;
        match['W']['N'] = -4;
        match['W']['D'] = -4;
        match['W']['C'] = -2;
        match['W']['Q'] = -2;
        match['W']['E'] = -3;
        match['W']['G'] = -2;
        match['W']['H'] = -2;
        match['W']['I'] = -3;
        match['W']['L'] = -2;
        match['W']['K'] = -3;
        match['W']['M'] = -1;
        match['W']['F'] = 1;
        match['W']['P'] = -4;
        match['W']['S'] = -3;
        match['W']['T'] = -2;
        match['W']['W'] = 11;
        match['W']['Y'] = 2;
        match['W']['V'] = -3;
        match['W']['B'] = -4;
        match['W']['Z'] = -3;
        match['W']['X'] = -2;
        match['W']['*'] = -4;


        match['Y']['A'] = -2;
        match['Y']['R'] = -2;
        match['Y']['N'] = -2;
        match['Y']['D'] = -3;
        match['Y']['C'] = -2;
        match['Y']['Q'] = -1;
        match['Y']['E'] = -2;
        match['Y']['G'] = -3;
        match['Y']['H'] = 2;
        match['Y']['I'] = -1;
        match['Y']['L'] = -1;
        match['Y']['K'] = -2;
        match['Y']['M'] = -1;
        match['Y']['F'] = 3;
        match['Y']['P'] = -3;
        match['Y']['S'] = -2;
        match['Y']['T'] = -2;
        match['Y']['W'] = 2;
        match['Y']['Y'] = 7;
        match['Y']['V'] = -1;
        match['Y']['B'] = -3;
        match['Y']['Z'] = -2;
        match['Y']['X'] = -1;
        match['Y']['*'] = -4;


        match['V']['A'] = 0;
        match['V']['R'] = -3;
        match['V']['N'] = -3;
        match['V']['D'] = -3;
        match['V']['C'] = -1;
        match['V']['Q'] = -2;
        match['V']['E'] = -2;
        match['V']['G'] = -3;
        match['V']['H'] = -3;
        match['V']['I'] = 3;
        match['V']['L'] = 1;
        match['V']['K'] = -2;
        match['V']['M'] = 1;
        match['V']['F'] = -1;
        match['V']['P'] = -2;
        match['V']['S'] = -2;
        match['V']['T'] = 0;
        match['V']['W'] = -3;
        match['V']['Y'] = -1;
        match['V']['V'] = 4;
        match['V']['B'] = -3;
        match['V']['Z'] = -2;
        match['V']['X'] = -1;
        match['V']['*'] = -4;


        match['B']['A'] = -2;
        match['B']['R'] = -1;
        match['B']['N'] = 3;
        match['B']['D'] = 4;
        match['B']['C'] = -3;
        match['B']['Q'] = 0;
        match['B']['E'] = 1;
        match['B']['G'] = -1;
        match['B']['H'] = 0;
        match['B']['I'] = -3;
        match['B']['L'] = -4;
        match['B']['K'] = 0;
        match['B']['M'] = -3;
        match['B']['F'] = -3;
        match['B']['P'] = -2;
        match['B']['S'] = 0;
        match['B']['T'] = -1;
        match['B']['W'] = -4;
        match['B']['Y'] = -3;
        match['B']['V'] = -3;
        match['B']['B'] = 4;
        match['B']['Z'] = 1;
        match['B']['X'] = -1;
        match['B']['*'] = -4;


        match['Z']['A'] = -1;
        match['Z']['R'] = 0;
        match['Z']['N'] = 0;
        match['Z']['D'] = 1;
        match['Z']['C'] = -3;
        match['Z']['Q'] = 3;
        match['Z']['E'] = 4;
        match['Z']['G'] = -2;
        match['Z']['H'] = 0;
        match['Z']['I'] = -3;
        match['Z']['L'] = -3;
        match['Z']['K'] = 1;
        match['Z']['M'] = -1;
        match['Z']['F'] = -3;
        match['Z']['P'] = -1;
        match['Z']['S'] = 0;
        match['Z']['T'] = -1;
        match['Z']['W'] = -3;
        match['Z']['Y'] = -2;
        match['Z']['V'] = -2;
        match['Z']['B'] = 1;
        match['Z']['Z'] = 4;
        match['Z']['X'] = -1;
        match['Z']['*'] = -4;


        match['X']['A'] = 0;
        match['X']['R'] = -1;
        match['X']['N'] = -1;
        match['X']['D'] = -1;
        match['X']['C'] = -2;
        match['X']['Q'] = -1;
        match['X']['E'] = -1;
        match['X']['G'] = -1;
        match['X']['H'] = -1;
        match['X']['I'] = -1;
        match['X']['L'] = -1;
        match['X']['K'] = -1;
        match['X']['M'] = -1;
        match['X']['F'] = -1;
        match['X']['P'] = -2;
        match['X']['S'] = 0;
        match['X']['T'] = 0;
        match['X']['W'] = -2;
        match['X']['Y'] = -1;
        match['X']['V'] = -1;
        match['X']['B'] = -1;
        match['X']['Z'] = -1;
        match['X']['X'] = -1;
        match['X']['*'] = -4;


        match['*']['A'] = -4;
        match['*']['R'] = -4;
        match['*']['N'] = -4;
        match['*']['D'] = -4;
        match['*']['C'] = -4;
        match['*']['Q'] = -4;
        match['*']['E'] = -4;
        match['*']['G'] = -4;
        match['*']['H'] = -4;
        match['*']['I'] = -4;
        match['*']['L'] = -4;
        match['*']['K'] = -4;
        match['*']['M'] = -4;
        match['*']['F'] = -4;
        match['*']['P'] = -4;
        match['*']['S'] = -4;
        match['*']['T'] = -4;
        match['*']['W'] = -4;
        match['*']['Y'] = -4;
        match['*']['V'] = -4;
        match['*']['B'] = -4;
        match['*']['Z'] = -4;
        match['*']['X'] = -4;
        match['*']['*'] = 1;
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
        int m = Math.min(s1.length(), MAX_LENGTH), n = Math.min(s2.length(), MAX_LENGTH);
        AlignmentBlock alignment = new AlignmentBlock();
        alignment.score = get_similarity(s1, s2);
        i = m;
        j = n;
        while (i > 0 && j > 0) {
            if (matrix[i][j] == up[i][j]) {
                alignment.one.append( s1.charAt(i-1) );
                alignment.two.append( '-' );
                i = i - 1;
            } else if (matrix[i][j] == left[i][j]) {
                alignment.one.append( '-' );
                alignment.two.append( s2.charAt(j-1) );
                j = j - 1;
            } else {
                alignment.one.append( s1.charAt(i-1) );
                alignment.two.append( s2.charAt(j-1) );
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
     * @return The similarity score
     */
    public double get_similarity(String s1, String s2) {
        int i, j;
        double match_score = 0, max_score = 0, match_length = 0;
        int m = Math.min(s1.length(), MAX_LENGTH), n = Math.min(s2.length(), MAX_LENGTH);
        /*if ( m > n )
            for (i=1; i<=m; ++i)
                max_score += match[s1.charAt(i-1)][s1.charAt(i-1)];
        else
            for (j=1; j<=n; ++j)
                max_score += match[s2.charAt(j-1)][s2.charAt(j-1)];
        */    
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
                matrix[i][j] = Math.max( match[s1.charAt(i-1)][s2.charAt(j-1)] + matrix[i - 1][j - 1] , Math.max( up[i][j] , left[i][j]) );
            }
        }
        i = m;
        j = n;
        while (i > 0 && j > 0) {
            if (matrix[i][j] == up[i][j]) {
                i = i - 1;
            } else if (matrix[i][j] == left[i][j]) {
                j = j - 1;
            } else {
                match_score += match[s1.charAt(i-1)][s2.charAt(j-1)];
                max_score += (match[s1.charAt(i-1)][s1.charAt(i-1)] + match[s2.charAt(j-1)][s2.charAt(j-1)])/2;
                ++match_length;
                i = i - 1;
                j = j - 1;
            }
        }      
        return (match_score / max_score) * (match_length / Math.min(m, n));
        //return matrix[m][n] / max_score;
    }
}
