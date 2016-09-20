/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author sheik005
 */
package protein;

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
    // degenerate bases will be replaces by 'A'    
        sequence=new StringBuilder();
    }
    public String translate(StringBuilder mRNA, boolean forward)
    {
        sequence.setLength(0);
        int i;   
        //if (mRNA.length() % 3 == 0)
            if ( forward )
            {
                for(i=0;i<=mRNA.length()-6;i+=3)
                    sequence.append(table[binary[mRNA.charAt(i)]*16+binary[mRNA.charAt(i+1)]*4+binary[mRNA.charAt(i+2)]]);
            }
            else 
            {
                for(i=mRNA.length()-1;i>=5;i-=3)
                    sequence.append(table[(3-binary[mRNA.charAt(i)])*16+(3-binary[mRNA.charAt(i-1)])*4+(3-binary[mRNA.charAt(i-2)])]);            
            }
        //System.out.println("proteinn:\n"+sequence.toString());
    return sequence.toString();
    }

} 
