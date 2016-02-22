/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pantools;

import pantools.*;

/**
 *
 * @author siavashs
 */
public class pointer {
    public int[] kmer;
    public long node_id;
    public int position;
    public long next_index;
    public byte format;
    public pointer(int d)
    {
        kmer=new int[d];
        for(int i=0;i<d;++i)
            kmer[i]=0;
        node_id=-1;
        next_index=-1L;
        format=-1;
    }
    public pointer(int d, long id,byte s,int p,int next)
    {
        node_id=id;
        position=p;
        format=s;
        next_index=next;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public int compare_kmer(pointer p, int d)
    {
        int i;
        for(i=0;i<d;++i)
            if(kmer[i]!=p.kmer[i])
                break;
        if(i==d)
            return 0;
        else if((kmer[i] & 0x00000000ffffffffL) < (p.kmer[i] & 0x00000000ffffffffL) )
            return -1;
        else
            return 1;
    }    
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,    
    public pointer revcom(int dim, int K)
    {
        int i;
        pointer rc=new pointer(dim);
        int frame=K%16==0?-1:(1<<(K%16*2))-1;
        for(i=0;i<dim;++i)
        {
            rc.kmer[i]=~kmer[i];  
        }
        rc.kmer[0]&=frame;
        return rc;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void next_up_kmer(int code, int dim, int K)
    {
        int i;
        int frame,left_code;
        frame=K%16==0?-1:(1<<(K%16*2))-1;
        left_code=3<<30;
        for(i=0;i<dim;++i)
        {
            kmer[i]=kmer[i]<<2;  
            if(i==0)
                kmer[0]=kmer[0] & frame; 
            if(i==dim-1)
                kmer[i]=kmer[i] | code; 
            else
                kmer[i]=kmer[i] | (( kmer[i+1] & left_code )>>>30); 
        }
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void next_down_kmer(int code, int dim, int K)
    {
        int i;
        int shift;
        shift=K%16==0?30:(K%16-1)*2;
        for(i=dim-1;i>0;--i)
        {
            kmer[i]=kmer[i]>>>2;  
            kmer[i]=kmer[i] | ((kmer[i-1] & 3 )<<30); 
        }
        kmer[i]=kmer[i]>>>2;  
        kmer[i]=kmer[i] | (code<<shift); 
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    void print_kmer(int dim,int K)
    {
        print_word(kmer[0],K%16);
        for(int i=1;i<dim;++i)
            print_word(kmer[i],16);
        System.out.print("\n");
    }  
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    void print_word(int word,int k)
    {
        if(k>0)
        {
            print_word(word>>>2,k-1);
            System.out.print(Pantools.sym[word & 3]);
        }
    } 
    
}
