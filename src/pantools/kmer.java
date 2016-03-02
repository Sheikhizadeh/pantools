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
    public byte[] bytes;
    public kmer(int d)
    {
        bytes=new byte[d];
        for(int i=0;i<d;++i)
            bytes[i]=0;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public kmer(byte[] b)
    {
        bytes=b;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void reset(int d)
    {
        for(int i=0;i<d;++i)
            bytes[i]=0;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public int compare(kmer k_mer, int d)
    {
        int i;
        for(i=0;i<d;++i)
            if(bytes[i]!=k_mer.bytes[i])
                break;
        if(i==d)
            return 0;
        else if((bytes[i] & 0x000000ff) < (k_mer.bytes[i] & 0x000000ff) )
            return -1;
        else
            return 1;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void next_up_kmer(int code, int dim, int K)
    {
        int i;
        int frame,left_code;
        frame=K%4==0?255:(1<<(K%4*2))-1;
        left_code=3<<6;
        for(i=0;i<dim;++i)
        {
            bytes[i]=(byte)(bytes[i]<<2);  
            if(i==0)
                bytes[0]=(byte)(bytes[0] & frame); 
            if(i==dim-1)
                bytes[i]=(byte)(bytes[i] | code); 
            else
                bytes[i]=(byte)(bytes[i] | ((( bytes[i+1] & left_code )>>6) & 0x00000003)); 
        }
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void next_down_kmer(int code, int dim, int K)
    {
        int i;
        int shift;
        shift=(K%4==0?6:(K%4-1)*2);
        for(i=dim-1;i>0;--i)
        {
            bytes[i]=(byte)((bytes[i]>>2) & 0x0000003f);  
            bytes[i]=(byte)(bytes[i] | ((bytes[i-1] & 3 )<<6)); 
        }
        bytes[0]=(byte)((bytes[0]>>2) & 0x0000003f);  
        bytes[0]=(byte)(bytes[0] | (code<<shift)); 
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    void print_kmer(int dim,int K)
    {
        print_byte(bytes[0] & 0x000000ff,K%4==0?4:K%4);
        for(int i=1;i<dim;++i)
            print_byte(bytes[i],4);
        System.out.print("\n");
    }  
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    void print_byte(int word,int k)
    {
        if(k>0)
        {
            print_byte(word>>2,k-1);
            System.out.print(Pantools.sym[word & 0x00000003]);
        }
    } 
}
