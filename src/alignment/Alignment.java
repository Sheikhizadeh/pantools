/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package alignment;

//Harry Hull
//Implements the following Smith-Waterman algorithm http://en.wikipedia.org/wiki/Smith_waterman
//	Affine Gap algorithm taken from:
//  http://en.wikipedia.org/wiki/Gap_penalty#Affine_gap_penalty
//	gap = o + (l-1)*e;
//	o:	gap opening penalty  (o < 0)
//	l:	length of the gap
//	e:	gap extension penalty (o < e < 0)
public class Alignment {

    //private String one, two;
    private int match[][];
    private int matrix[][];
    private char path[][];
    private int maxlength;
    private int gap;
    private int l;
    private StringBuilder alignone;
    private StringBuilder aligntwo;
    private Alignment aligner;
    private Alignment alignm;
    public static int match_score=2;
    public static int gap_open=-4;
    public static int gap_extn=-1;
    
    public Alignment() {
        // initialize matrixes
        l = 0;
        // initialize match matrix
        match = new int[256][256];
        for (int i = 0; i < 256; i++) {
            for (int j = 0; j < 256; j++) {
                if (i == j) {
                    match[i][j] = match_score;
                } else {
                    match[i][j] = gap_extn; // mismatch is penalized like a gap extension
                }
            }
        }
        match['A']['D'] = match_score;
        match['A']['H'] = match_score;
        match['A']['M'] = match_score;
        match['A']['N'] = match_score;
        match['A']['R'] = match_score;
        match['A']['V'] = match_score;
        match['A']['W'] = match_score;
        match['C']['B'] = match_score;
        match['C']['H'] = match_score;
        match['C']['M'] = match_score;
        match['C']['N'] = match_score;
        match['C']['S'] = match_score;
        match['C']['V'] = match_score;
        match['C']['Y'] = match_score;
        match['G']['B'] = match_score;
        match['G']['D'] = match_score;
        match['G']['K'] = match_score;
        match['G']['N'] = match_score;
        match['G']['R'] = match_score;
        match['G']['S'] = match_score;
        match['G']['V'] = match_score;
        match['T']['B'] = match_score;
        match['T']['D'] = match_score;
        match['T']['H'] = match_score;
        match['T']['K'] = match_score;
        match['T']['N'] = match_score;
        match['T']['W'] = match_score;
        match['T']['Y'] = match_score;
        match['D']['A'] = match_score;
        match['H']['A'] = match_score;
        match['M']['A'] = match_score;
        match['N']['A'] = match_score;
        match['R']['A'] = match_score;
        match['V']['A'] = match_score;
        match['W']['A'] = match_score;
        match['B']['C'] = match_score;
        match['H']['C'] = match_score;
        match['M']['C'] = match_score;
        match['N']['C'] = match_score;
        match['S']['C'] = match_score;
        match['V']['C'] = match_score;
        match['Y']['C'] = match_score;
        match['B']['G'] = match_score;
        match['D']['G'] = match_score;
        match['K']['G'] = match_score;
        match['N']['G'] = match_score;
        match['R']['G'] = match_score;
        match['S']['G'] = match_score;
        match['V']['G'] = match_score;
        match['B']['T'] = match_score;
        match['D']['T'] = match_score;
        match['H']['T'] = match_score;
        match['K']['T'] = match_score;
        match['N']['T'] = match_score;
        match['W']['T'] = match_score;
        match['Y']['T'] = match_score;
    }

    // reset the matrixes
    public void reset(int max) {
        alignone.setLength(0);
        aligntwo.setLength(0);
        for (int i = 0; i <= max; i++) {
            for (int j = 0; j <= max; j++) {
                matrix[i][j] = 0;
                path[i][j] = 0;
            }
        }
    }

    // returns the alignment score
    public void computeSmithWaterman(AlignmentBlock a) {
        int i, j, align, insert, delete;
        int m = a.one.length() - 1, n = a.two.length() - 1, max = Math.max(m, n);
        matrix = new int[m + 1][n + 1];
        path = new char[m + 1][n + 1];
        alignone = new StringBuilder(2 * max);
        aligntwo = new StringBuilder(2 * max);
        //reset(2*max);
        //System.out.println("aligning:"+m+" "+n);
        //System.out.println(a.one.toString());
        //System.out.println(a.two.toString());
        for (i = 1; i <= m; i++) {
            for (j = 1; j <= n; j++) {
                gap = gap_open + (l - 1) * gap_extn;
                align = match[a.one.charAt(i)][a.two.charAt(j)] + matrix[i - 1][j - 1];
                insert = gap + matrix[i][j - 1];
                delete = gap + matrix[i - 1][j];
                if (align >= insert && align >= delete) // are equal
                {
                    //System.out.println("match");
                    l = 0; // reset l
                    matrix[i][j] = align;
                    path[i][j] = 'a';
                } else if (insert >= align && insert >= delete) {
                    //System.out.println("up");
                    l++;
                    matrix[i][j] = insert;
                    path[i][j] = 'i';
                } else {
                    //System.out.println("left");
                    l++;
                    matrix[i][j] = delete;
                    path[i][j] = 'd';
                }
                /*if(matrix[i][j]>longest)
                 {
                 longest=matrix[i][j];
                 il=i;
                 jl=j;
                 }*/
            }
        }
        //printMatrix(one,two);

        // Backtrack to reconstruct the path
        i = m;
        j = n;
        while (i > 0 && j > 0) {
            // diag case
            //System.out.println(path[i][j]);
            //System.out.println(one.charAt(i));
            //System.out.println(two.charAt(j));
            if (path[i][j] == 'a') {
                alignone.append(a.one.charAt(i));
                aligntwo.append(a.two.charAt(j)); // or one.charAt(i) 
                i = i - 1;
                j = j - 1;
                // left case
            } else if (path[i][j] == 'i') {
                alignone.append('-');
                aligntwo.append(a.two.charAt(j));
                j = j - 1;
                // up case
            } else {
                alignone.append(a.one.charAt(i));
                aligntwo.append('-');
                i = i - 1;
            }
        }
        for (; i > 0; --i) {
            alignone.append(a.one.charAt(i));
            aligntwo.append('-');
        }
        for (; j > 0; --j) {
            alignone.append('-');
            aligntwo.append(a.two.charAt(j));
        }
        //System.out.println("alignone "+alignone.reverse().toString());
        //System.out.println("aligntwo "+aligntwo.reverse().toString());
        a.reset();
        for (i = alignone.length() - 1; i >= 0; --i) {
            a.one.append(alignone.charAt(i));
            a.two.append(aligntwo.charAt(i));
        }
        // print alignment
        //System.out.println(alignone.toString() + "\n" + aligntwo.toString());
        //return longest;
    }

    public double similarity(String s1, String s2) {
        int i, j, align, insert, delete, longest;
        int m = s1.length() - 1, n = s2.length() - 1, max = Math.max(m, n);
        matrix = new int[m + 1][n + 1];
        longest = -max * match_score;
        for (i = 1; i <= m; i++) {
            for (j = 1; j <= n; j++) {
                gap = gap_open + (l - 1) * gap_extn;
                align = match[s1.charAt(i)][s2.charAt(j)] + matrix[i - 1][j - 1];
                insert = gap + matrix[i][j - 1];
                delete = gap + matrix[i - 1][j];
                if (align >= insert && align >= delete) // are equal
                {
                    //System.out.println("match");
                    l = 0; // reset l
                    matrix[i][j] = align;
                } else if (insert >= align && insert >= delete) {
                    //System.out.println("up");
                    l++;
                    matrix[i][j] = insert;
                } else {
                    //System.out.println("left");
                    l++;
                    matrix[i][j] = delete;
                }
                if (matrix[i][j] > longest) {
                    longest = matrix[i][j];
                }
            }
        }
        return (double) longest / max;
    }
    /*
    
     */

    public void printMatrix(String one, String two) {
        for (int i = 0; i < one.length(); i++) {
            if (i == 0) {
                for (int z = 0; z < two.length(); z++) {
                    if (z == 0) {
                        System.out.print("   ");
                    }
                    System.out.print(two.charAt(z) + "  ");

                    if (z == two.length() - 1) {
                        System.out.println();
                    }
                }
            }

            for (int j = 0; j < two.length(); j++) {
                if (j == 0) {
                    System.out.print(one.charAt(i) + "  ");
                }
                System.out.print(matrix[i][j] + "  ");
            }
            System.out.println();
        }
        System.out.println();
    }

    /*
     To compute identity between two sequences
     */
    int identity(StringBuilder s1, StringBuilder s2) {
        int i, sum;
        for (i = sum = 0; i < s1.length(); ++i) {
            if (s1.charAt(i) == s2.charAt(i)) {
                ++sum;
            }
        }
        return sum * 100 / s1.length();
    }
}
