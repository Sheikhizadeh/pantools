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
    public int K;
    public int pre_len;
    public int suf_len;
    public long value; 
    public boolean canonical;
    
    /**
     * A constructor
     * 
     * @param k Size of K
     * @param p The length of the prefix of the kmer
     * @param s The length of the suffix of the kmer
     */
    public kmer(int k, int p, int s)
    {
        K=k;
        pre_len=p;
        suf_len=s;
        value = 0;
        canonical=false;
    }

    /**
     * Copy constructor
     * 
     * @param k_mer The kmer object we want to copy. 
     */
    public kmer(kmer k_mer)
    {
        K=k_mer.K;
        pre_len=k_mer.pre_len;
        suf_len=k_mer.suf_len;
        value=k_mer.value;
        canonical=k_mer.canonical;
    }

    /**
     * Clears the content of the kmer
     */
    public void reset()
    {
        value = 0;
        canonical=false;
    }
    public void adjust(int s_len){
        suf_len = s_len;
        pre_len = K - s_len;
    }
    /**
     * Compares the suffices of two kmers
     * @param k_mer The second kmer
     * @return -1 if the first suffix is smaller, 0 if they are equal and 1 if the second suffix is smaller
     */
    public int compare_suffix(byte[] suf)
    {
        int i;
        long given_suffix = 0, my_suffix = get_suffix();
        for(i = 0; i < suf_len / 4; ++i){
            given_suffix = (given_suffix << 8) | (0x00000000000000ff & suf[i]);
        }
        if(given_suffix == my_suffix)
            return 0;
        else if(my_suffix < given_suffix)
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
        if(value < k_mer.value)
            return -1;
        else if(value > k_mer.value)
            return 1;
        else
            return 0;

    }

    /**
     * Gives the next forward kmer
     * @param base_code The binary code of the nucleotide at the end of the kmer
     */
    public void next_up_kmer(long base_code)
    {
        value = ((value & ((1L << (2 * (K-1))) - 1) ) << 2) | base_code;
    }

    /**
     * Gives the next reverse kmer
     * @param base_code The binary code of the nucleotide at the end of the kmer
     */
    public void next_down_kmer(long base_code)
    {
        value = (value >> 2) | (base_code << (2 * (K - 1)));
    }
    
    public long get_prefix(){
        return value >> (2 * suf_len);
    }
    
    public long get_suffix(){
        return value & ((1L << (2 * suf_len)) - 1);
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
        append_word(sym, value,K,seq);
        return seq.toString();
    }  

    /**
     * A recursive function for decoding the binary kmer to a string kmer  
     * @param sym Array of symbols
     * @param word The binary string to be decoded
     * @param k Size of the K
     * @param seq The string builder containing the decoded kmer
     */
    private void append_word(char[] sym, long value,int k, StringBuilder seq )
    {
        if(k>0)
        {
            append_word(sym, value>>2,k-1,seq);
            seq.append(sym[(int)(value & 3l)]);
        }
    } 
}
