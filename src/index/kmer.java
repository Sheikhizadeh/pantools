/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package index;

import org.bouncycastle.util.Arrays;

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
    private int fwd_prefix;
    private byte[] fwd_suffix;
    private int rev_prefix;
    private byte[] rev_suffix;
    private boolean canonical;
    
    /**
     * The constructor
     * 
     * @param k Size of K
     * @param p_len The length of the prefix of the kmer
     * @param s_len The length of the suffix of the kmer
     */
    public kmer(int k, int p_len)
    {
        K = k;
        prefix_length = p_len;
        suffix_length = K - p_len;
        prefix_mask = (1 << (2*prefix_length)) - 1;
        shift = 2 * (prefix_length - 1);
        fwd_suffix=new byte[suffix_length / 4];
        rev_suffix=new byte[suffix_length / 4];
        canonical = false;
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
        fwd_prefix = k_mer.fwd_prefix;
        fwd_suffix = Arrays.clone(k_mer.fwd_suffix);
        rev_prefix = k_mer.rev_prefix;
        rev_suffix = Arrays.clone(k_mer.rev_suffix);
        canonical = k_mer.canonical;
    }
    
    /**
     * Gives the prefix of the k-mer
     * @return The prefix of the k-mer in the form of an integer
     */
    public int get_canonical_prefix(){
        return canonical ? fwd_prefix : rev_prefix;
    }
    
    /**
     * Gives the suffix of the k-mer
     * @return The suffix of the k-mer in the form of a byte array
     */
    public byte[] get_canonical_suffix(){
        return canonical ? fwd_suffix : rev_suffix;
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
    public void set_fwd_prefix(int p){
        fwd_prefix = p;
    }    
    
    /**
     * Sets the prefix to a new value
     * 
     * @param p The new prefix for the k-mer  
     */
    public void set_fwd_suffix(byte[] s){
        fwd_suffix = s;
    }    

    /**
     * Sets the prefix to a new value
     * 
     * @param p The new prefix for the k-mer  
     */
    public void set_canonical(boolean c){
        canonical = c;
    } 

    /**
     * Sets the prefix to a new value
     * 
     * @param p The new prefix for the k-mer  
     */
    public int get_fwd_prefix(){
        return fwd_prefix;
    }    
    
    /**
     * Sets the prefix to a new value
     * 
     * @param p The new prefix for the k-mer  
     */
    public byte[] get_fwd_suffix(){
        return fwd_suffix;
    }    
    
    /**
     * Recalculates the suffix and, in turn, the prefix based on the new suffix length
     * 
     * @param s_len The new suffix length for the k-mer  
     */
    public static void adjust_fwd_kmer(kmer kmer2, kmer kmer1){
        int i, j, k, suffix_length1, suffix_length2;
        suffix_length1 = kmer1.suffix_length;
        suffix_length2 = kmer2.suffix_length;
        if (suffix_length1 != suffix_length2){
            for(i = suffix_length1/4 - 1, j = suffix_length2/4 - 1; i >= 0 && j >= 0; --i, --j)
                kmer2.fwd_suffix[j] = kmer1.fwd_suffix[i];  
            if (i >= 0){
                for (kmer2.fwd_prefix = kmer1.fwd_prefix, k = 0; k <= i; ++k)
                    kmer2.fwd_prefix = (kmer2.fwd_prefix << 8) | (kmer1.fwd_suffix[k] & 0x00FF);
            }
            if (j >= 0){
                kmer2.fwd_prefix = kmer1.fwd_prefix;
                for (;j >= 0; --j){
                    kmer2.fwd_suffix[j] = (byte)(kmer2.fwd_prefix & 0x0FF);
                    kmer2.fwd_prefix = kmer2.fwd_prefix >> 8;
                }
            }
        } else{
            //kmer2.set_fwd_suffix(Arrays.clone(kmer1.get_fwd_suffix()));
            //kmer2.set_fwd_prefix(kmer1.get_fwd_prefix());
            kmer2 = kmer1;
        }
        kmer2.canonical = true; // to be found by find()
    }
    
    /**
     * Clears the content of the kmer
     */
    public void reset()
    {
        fwd_prefix=0;
        for(int i = 0;i < fwd_suffix.length; ++i)
            fwd_suffix[i] = 0;
        rev_prefix = 0;
        for(int i = 0; i < rev_suffix.length; ++i)
            rev_suffix[i] = 0;
        canonical=false;
    }
    
    /**
     * Compares two kmers
     * @param k_mer The second kmer
     * @return -1 if the first kmer is smaller, 0 if they are equal and 1 if the second kmer is smaller
     */
    public void is_canonical()
    {
        if(fwd_prefix<rev_prefix)
            canonical = true;
        else if(fwd_prefix>rev_prefix)
            canonical = false;
        else
            canonical = compare_suffix(fwd_suffix, rev_suffix) <= 0;

    }

    /**
     * Compares the suffices of two kmers
     * @param k_mer The second kmer
     * @return -1 if the first suffix is smaller, 0 if they are equal and 1 if the second suffix is smaller
     */
    public static int compare_suffix(byte[] suf1, byte[] suf2){
        int i;
        for(i=0;i<suf1.length;++i)
            if(suf1[i]!=suf2[i])
                break;
        if(i==suf1.length)
            return 0;
        else if((suf1[i] & 0x0ff) < (suf2[i] & 0x0ff) )
            return -1;
        else
            return 1;        
    }
    
    /**
     * Gives the next forward kmer
     * @param base_code The binary code of the base at right end of the kmer
     */
    public void next_kmer(int base_code)
    {
        int i;
        fwd_prefix=((fwd_prefix<<2) & prefix_mask) | ((fwd_suffix[0]>>6) & 0x03 );
        for(i=0;i<fwd_suffix.length-1;++i)
           fwd_suffix[i]=(byte)((fwd_suffix[i]<<2) | (( fwd_suffix[i+1]>>6) & 0x03)); 
        fwd_suffix[i]=(byte)((fwd_suffix[i]<<2) | base_code); 
        base_code = 3 - base_code;
        for(i=rev_suffix.length-1;i>0;--i)
            rev_suffix[i]=(byte)(((rev_suffix[i]>>2) & 0x03f) | ((rev_suffix[i-1] & 0x03 )<<6)); 
        rev_suffix[i]=(byte)(((rev_suffix[i]>>2) & 0x03f) | ((rev_prefix & 0x03 )<<6)); 
        rev_prefix=(rev_prefix>>2) | (base_code<<shift); 
        is_canonical(); 
    }

    /**
     * Gives the next forward kmer
     * @param base_code The binary code of the base at right end of the kmer
     */
    public void prev_kmer(int base_code)
    {
        int i;
        for(i=fwd_suffix.length-1;i>0;--i){
            fwd_suffix[i]=(byte)(((fwd_suffix[i]>>2) & 0x03f) | ((fwd_suffix[i-1] & 0x03 )<<6)); 
        }
        fwd_suffix[i]=(byte)(((fwd_suffix[i]>>2) & 0x03f) | ((fwd_prefix & 0x03 )<<6)); 
        fwd_prefix=(fwd_prefix>>2) | (base_code<<shift); 
        base_code = 3 - base_code;
        rev_prefix=((rev_prefix<<2) & prefix_mask) | ((rev_suffix[0]>>6) & 0x03 );
        for(i=0;i<rev_suffix.length-1;++i)
           rev_suffix[i]=(byte)((rev_suffix[i]<<2) | (( rev_suffix[i+1]>>6) & 0x03)); 
        rev_suffix[i]=(byte)((rev_suffix[i]<<2) | base_code); 
        is_canonical(); 
    }    

    /**
     * Represents a kmer in the form of a string.
     * @return The string representing the kmer
     */
    @Override
    public String toString()
    {
        char[] sym=new char[]{ 'A', 'C', 'G' , 'T', 'M','R','W','S','Y','K','V','H','D','B','N'};
        StringBuilder fwd_seq=new StringBuilder(K);
        StringBuilder rev_seq=new StringBuilder(K);
        build_sequence(sym, fwd_prefix, prefix_length,fwd_seq);
        for(int i=0; i < fwd_suffix.length; ++i)
            build_sequence(sym, fwd_suffix[i], 4, fwd_seq);
        build_sequence(sym, rev_prefix, prefix_length,rev_seq);
        for(int i=0; i < rev_suffix.length; ++i)
            build_sequence(sym, rev_suffix[i], 4, rev_seq);
        return fwd_seq.append(" ").append(rev_seq).toString();
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