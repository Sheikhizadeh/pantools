/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pantools;
/**
 *
 * @author sheik005
 */
public class kmer {
    public int K;
    public int pre_len;
    public int suf_len;
    public int mask;
    public int shift;
    public int prefix;
    public byte[] suffix;
    public boolean canonical;
    public kmer(int k, int p, int s)
    {
        K=k;
        pre_len=p;
        suf_len=s;
        mask=(1<<(2*pre_len))-1;
        shift=2*(pre_len-1);
        suffix=new byte[suf_len/4];
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public kmer(kmer k_mer)
    {
        K=k_mer.K;
        pre_len=k_mer.pre_len;
        suf_len=k_mer.suf_len;
        mask=k_mer.mask;
        shift=k_mer.shift;
        prefix=k_mer.prefix;
        suffix=k_mer.suffix;
        canonical=k_mer.canonical;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void reset()
    {
        prefix=0;
        for(int i=0;i<suffix.length;++i)
            suffix[i]=0;
        canonical=false;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
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
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public int compare(kmer k_mer)
    {
        if(prefix<k_mer.prefix)
            return -1;
        else if(prefix>k_mer.prefix)
            return 1;
        else
            return compare_suffix(k_mer);

    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void next_up_kmer(int code)
    {
        int i;
        prefix=((prefix<<2) & mask) | ((suffix[0]>>6) & 0x03 );
        for(i=0;i<suffix.length-1;++i)
           suffix[i]=(byte)((suffix[i]<<2) | (( suffix[i+1]>>6) & 0x03)); 
        suffix[i]=(byte)((suffix[i]<<2) | code); 
        //if(prefix<0)
          //  System.out.println("Up: "+toString());
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void next_down_kmer(int code)
    {
        int i;
        for(i=suffix.length-1;i>0;--i)
            suffix[i]=(byte)(((suffix[i]>>2) & 0x03f) | ((suffix[i-1] & 0x03 )<<6)); 
        suffix[i]=(byte)(((suffix[i]>>2) & 0x03f) | ((prefix & 0x03 )<<6)); 
        prefix=(prefix>>2) | (code<<shift); 
        //if(prefix<0)
          //  System.out.println("Do: "+toString());
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    @Override
    public String toString()
    {
        char[] sym=new char[]{ 'A', 'C', 'G' , 'T', 'M','R','W','S','Y','K','V','H','D','B','N'};
        StringBuilder seq=new StringBuilder(K);
        append_word(sym, prefix, pre_len,seq);
        for(int i=0;i<suffix.length;++i)
            append_word(sym, suffix[i],4,seq);
        return seq.toString();
    }  
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    private void append_word(char[] sym, int word,int k, StringBuilder seq )
    {
        if(k>0)
        {
            append_word(sym, word>>2,k-1,seq);
            seq.append(sym[word & 0x00000003]);
        }
    } 
}
