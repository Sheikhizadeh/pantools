/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package index;

/**
 * Implements the data structure for a kmer. 
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class kmer {
    private int K; 
    private int prefix_length;   
    private int suffix_length;
    private int prefix_mask;   // A binary code with just prefix_length 1's on the right
    private int shift;         // Number of left shifts to reach to the position of left most base in the kmer 
    private int prefix;
    private byte[] suffix;
    private boolean canonical;
    
    /**
     * The constructor
     * 
     * @param k Size of K
     * @param p_len The length of the prefix of the kmer
     * @param s_len The length of the suffix of the kmer
     */
    public kmer(int k, int p_len, int s_len)
    {
        K = k;
        prefix_length = p_len;
        suffix_length = s_len;
        prefix_mask = (1 << (2*prefix_length)) - 1;
        shift = 2 * (prefix_length - 1);
        suffix=new byte[suffix_length / 4];
    }

    /**
     * Copy constructor
     * 
     * @param k_mer The kmer object we want to copy. 
     */
    public kmer(kmer k_mer)
    {
        K = k_mer.K;
        prefix_length = k_mer.prefix_length;
        suffix_length = k_mer.suffix_length;
        prefix_mask = k_mer.prefix_mask;
        shift = k_mer.shift;
        prefix = k_mer.prefix;
        suffix = k_mer.suffix;
        canonical = k_mer.canonical;
    }
    
    /**
     * Gives the prefix of the k-mer
     * @return The prefix of the k-mer in the form of an integer
     */
    public int get_prefix(){
        return prefix;
    }
    
    /**
     * Gives the suffix of the k-mer
     * @return The suffix of the k-mer in the form of a byte array
     */
    public byte[] get_suffix(){
        return suffix;
    }

    /**
     * Determines if the k-mer is canonical
     * @return True or False
     */
    public boolean get_canonical(){
        return canonical;
    }

    /**
     * Sets the prefix to a new value
     * 
     * @param p The new prefix for the k-mer  
     */
    public void set_prefix(int p){
        prefix = p;
    }    
    
    /**
     * Recalculates the suffix and, in turn, the prefix based on the new suffix length
     * 
     * @param s_len The new suffix length for the k-mer  
     */
    public void set_suffix(int s_len){
        int i, j, old_suffix_length = suffix_length;
        byte[] old_suffix = suffix;
        suffix = new byte[s_len / 4];
        
        for(i = old_suffix_length/4 - 1, j = s_len/4 - 1; i >= 0 && j >= 0; --i, --j)
            suffix[j] = old_suffix[i];        
        
        suffix_length = s_len;
        
        for (;i >= 0; --i){
            prefix = (prefix << 8) | (old_suffix[i] & 0x00FF);
        }
        
        for (;j >= 0; --j){
            suffix[j] = (byte)(prefix & 0x0FF);
            prefix = prefix >> 8;
        }
        
        prefix_length = K - suffix_length;
        prefix_mask=(1<<(2*prefix_length))-1;
        shift=2*(prefix_length-1);
    }
    
    /**
     * Sets the canonical value of the kmer
     * 
     * @param c The canonical value  
     */
    public void set_canonical(boolean c){
        canonical = c;
    }
    
    /**
     * Clears the content of the kmer
     */
    public void reset()
    {
        prefix=0;
        for(int i=0;i<suffix.length;++i)
            suffix[i]=0;
        canonical=false;
    }
    
    /**
     * Compares the suffices of two kmers
     * @param k_mer The second kmer
     * @return -1 if the first suffix is smaller, 0 if they are equal and 1 if the second suffix is smaller
     */
    public int compare_suffix(kmer k_mer)
    {
        int i;
        for(i=0;i<suffix.length;++i)
            if(suffix[i]!=k_mer.suffix[i])
                break;
        if(i==suffix.length)
            return 0;
        else if((suffix[i] & 0x0ff) < (k_mer.suffix[i] & 0x0ff) )
            return -1;
        else
            return 1;
    }
    
    /**
     * Compares two kmers
     * @param k_mer The second kmer
     * @return -1 if the first kmer is smaller, 0 if they are equal and 1 if the second kmer is smaller
     */
    public int compare(kmer k_mer)
    {
        if(prefix<k_mer.prefix)
            return -1;
        else if(prefix>k_mer.prefix)
            return 1;
        else
            return compare_suffix(k_mer);

    }

    /**
     * Gives the next forward kmer
     * @param base_code The binary code of the base at right end of the kmer
     */
    public void next_fwd_kmer(int base_code)
    {
        int i;
        prefix=((prefix<<2) & prefix_mask) | ((suffix[0]>>6) & 0x03 );
        for(i=0;i<suffix.length-1;++i)
           suffix[i]=(byte)((suffix[i]<<2) | (( suffix[i+1]>>6) & 0x03)); 
        suffix[i]=(byte)((suffix[i]<<2) | base_code); 
    }

    /**
     * Gives the next reverse kmer
     * @param base_code The binary code of the base at the left end of the kmer
     */
    public void next_rev_kmer(int base_code)
    {
        int i;
        for(i=suffix.length-1;i>0;--i)
            suffix[i]=(byte)(((suffix[i]>>2) & 0x03f) | ((suffix[i-1] & 0x03 )<<6)); 
        suffix[i]=(byte)(((suffix[i]>>2) & 0x03f) | ((prefix & 0x03 )<<6)); 
        prefix=(prefix>>2) | (base_code<<shift); 
    }

    /**
     * Represents a kmer in the form of a string.
     * @return The string representing the kmer
     */
    @Override
    public String toString()
    {
        char[] sym=new char[]{ 'A', 'C', 'G' , 'T', 'M','R','W','S','Y','K','V','H','D','B','N'};
        StringBuilder seq=new StringBuilder(K);
        build_sequence(sym, prefix, prefix_length,seq);
        for(int i=0; i < suffix.length; ++i)
            build_sequence(sym, suffix[i], 4, seq);
        return seq.toString();
    }  

    /**
     * A recursive function for decoding the binary kmer to a string kmer  
     * @param sym Array of symbols
     * @param word The binary string to be decoded
     * @param k Size of the K
     * @param seq The string builder containing the decoded kmer
     */
    private void build_sequence(char[] sym, int word,int k, StringBuilder seq )
    {
        if(k > 0)
        {
            build_sequence(sym, word >> 2,k-1,seq);
            seq.append(sym[word & 0x00000003]);
        }
    } 
}