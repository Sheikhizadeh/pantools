package pantools;


import pantools.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author siavashs
 */
public class parallel_sort {
    pointer[] pointers;
    pointer[] tmp_pointers;
    int len;
    int dim;
    int cores;
    parallel_sort(pointer[] a, int n, int d, int c)
    {
        pointers=a;
        tmp_pointers=new pointer[n];
        len=n;
        dim=d;
        cores=c;
    }
    void sort()
    {
        mergesort(0,len-1);
        tmp_pointers=null;
    }
    void psort()
    {
        int size;
        int chunks=cores;
        int i;
        int[] lows,highs;
        ExecutorService e = Executors.newFixedThreadPool(cores);
        lows=new int[cores+1];
        highs=new int[cores];
        for(i=0,lows[0]=0;i<cores;++i)
        {
            size=(int)Math.ceil(len/chunks);
            highs[i]=lows[i]+size-1;
            lows[i+1]=highs[i]+1;
            chunks--;
            len-=size;
            e.execute(new Runnable() {
                int low;
                int high;
                @Override                            
                public void run() {
                    mergesort(low,high);                              
                }
                public Runnable init(int l, int h) {
                    this.low=l;
                    this.high=h;
                    return(this);
                }
            }.init(lows[i],highs[i])); 
        }
        e.shutdown();
        try {
            e.awaitTermination(1, TimeUnit.HOURS);
        }
        catch( InterruptedException ex ){
            System.out.println("Timeout exception for parallel tasks!");
            System.exit(1);
        }
        for(i=0;i<cores-1;++i)
            merge(lows[0],highs[i],highs[i+1]);
    }
//,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void merge(int low, int mid, int high)
    {
        int i,j,k;
        for(i=k=low,j=mid+1;i<=mid && j<=high;++k)
            if(pointers[i].compare_kmer(pointers[j],dim)==-1)
                tmp_pointers[k]=pointers[i++];
            else
                tmp_pointers[k]=pointers[j++];
        while(i<=mid)
                tmp_pointers[k++]=pointers[i++];
        while(j<=high)
                tmp_pointers[k++]=pointers[j++];
        for(i=low;i<=high;++i)
                pointers[i]=tmp_pointers[i];
    }
//,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void mergesort(int low,int high)
    {
        if(low<high)
        {
            int mid=(low+high)/2;
            mergesort(low,mid);
            mergesort(mid+1,high);
            merge(low,mid,high);
        }
    }
}
