/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package protein;

/**
 * Implements some functionality to work with proteins.
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class protein_builder {   
    char[] table;
    int[] binary;
    StringBuilder sequence;
    public protein_builder()
    {
        table = new char[]
          {'K','N','K','N',
           'T','T','T','T',
           'R','S','R','S',
           'I','I','M','I',
           'Q','H','Q','H',
           'P','P','P','P',
           'R','R','R','R',
           'L','L','L','L',
           'E','D','E','D',
           'A','A','A','A',
           'G','G','G','G',
           'V','V','V','V',
           '*','Y','*','Y',
           'S','S','S','S',
           '*','C','W','C',
           'L','F','L','F'
          };
        binary=new int[256];
        binary['A'] = 0; 
        binary['C'] = 1; 
        binary['G'] = 2; 
        binary['T'] = 3; 
    // All degenerate bases will be replaces by 'A' in the translation  
        sequence=new StringBuilder();
    }
    
    /**
     * Translates a coding nucleotide sequence to a protein sequence.
     * 
     * @param mRNA The sequence of the codinf RNA
     * @return The protein sequence
     */
    public StringBuilder translate(StringBuilder mRNA)
    {
        sequence.setLength(0);
        int i;   
        for(i=0;i<=mRNA.length()-3;i+=3)
            sequence.append(table[binary[mRNA.charAt(i)]*16+binary[mRNA.charAt(i+1)]*4+binary[mRNA.charAt(i+2)]]);
        return sequence;
    }
} 
