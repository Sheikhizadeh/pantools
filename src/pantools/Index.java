/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pantools;

import pantools.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.FileOutputStream;
import java.io.FileInputStream;
import org.neo4j.graphdb.DynamicLabel;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import static pantools.Pantools.Genomes;
import static pantools.Pantools.sequence_length;
import static pantools.Pantools.binary;
import static pantools.Pantools.graphDb;

/**
 *
 * @author siavashs
 */
public class Index {
    char sym[]={ 'A', 'C', 'G' , 'T' , 'N'};
    public pointer[][] pointers;
    public long totallength;
    public int[] length;
    public int parts;
    public int dim;
    private String path;
    private String name;
    private int genome,sequence,K;
    private int position,seq_len;
    private int f_bin,r_bin;
    private pointer f_kmer,r_kmer,kmer;
    private boolean canonical;
    private boolean finish;


    public Index(String file, String p, String n, int K_size)throws IOException
    {
        int mod,cores=Runtime.getRuntime().availableProcessors();
        String line;
        parallel_sort ps;
        path=p;
        name=n;
        K=K_size;
        dim=(int)Math.ceil(K/16.0);
        int j,k,mini;
        long i=0;
        int[] inx;
        if(K<10)
        {
            generate_short_kmers();
        }
        else
        {
            System.out.println("Running KMC2...");
            String output=Pantools.executeCommand("kmc -r -k"+K+" -t"+cores+" -ci1 -fm "+file+" "+path+"/tmp "+path+name);
            String[] fields=output.split(" +");
            for(j=0;!fields[j].equals("unique");++j);
            totallength=Long.parseLong(fields[j+3].trim());//
            System.out.println("Indexing "+totallength+" kmers...");
            Pantools.executeCommand("kmc_dump -r -cx0 "+path+"/tmp "+path+name);
            try(BufferedReader in = new BufferedReader(new FileReader(path+name))){
            parts=(int)(totallength/2100000000)+1;
            inx=new int[parts];
            pointers=new pointer[parts][];
            length=new int[parts];
            for(k=0;k<parts;++k)
            {
                length[k]=k<parts-1?2100000000:(int)(totallength%2100000000);
                pointers[k]=new pointer[length[k]];
            }
            mod=(16-K%16)%16;
            for(i=0,line=in.readLine();i<totallength;line=in.readLine(),++i)
            {
                pointers[(int)(i/2100000000)][(int)(i%2100000000)]=new pointer(dim);
                for(j=mod;j<K+mod;++j)
                    pointers[(int)(i/2100000000)][(int)(i%2100000000)].kmer[j/16]=(pointers[(int)(i/2100000000)][(int)(i%2100000000)].kmer[j/16]<<2 | binary[line.charAt(j-mod)]); 
                if(i%(totallength/100+1)==0) 
                    System.out.print((long)i*100/totallength+1+"%\r");
            }
            in.close();
            Pantools.executeCommand("rm "+path+"/tmp.txt.kmc_pre");
            Pantools.executeCommand("rm "+path+"/tmp.txt.kmc_suf");
            Pantools.executeCommand("rm "+path+name);
            System.out.println("Sorting index...");
            for(k=0;k<parts;++k) // sort each part
            {
                ps=new parallel_sort(pointers[k],length[k],dim,cores);
                if(length[k]<10000)
                   ps.sort();
                else
                   ps.psort();
            }
            if(parts>1)
            {
                pointer[][] tmp_pointers=new pointer[parts][];
                for(k=0;k<parts;++k)
                    tmp_pointers[k]=new pointer[length[k]];
                for(i=0;i<totallength;++i) // sort array of arrays
                {
                    for(mini=0;mini<parts && inx[mini]==length[mini];++mini)
                       {}
                    for(k=0;k<parts;++k)
                        if(inx[k]<length[k] && pointers[k][inx[k]].compare_kmer(pointers[mini][inx[mini]],dim)==-1)
                            mini=k;
                    tmp_pointers[(int)(i/2100000000)][(int)(i%2100000000)]=pointers[mini][inx[mini]];
                    pointers[mini][inx[mini]]=null;
                    ++inx[mini];
                    if(inx[mini]==length[mini])
                        pointers[mini]=null;
                }
                pointers=tmp_pointers;
            }
            }
            catch(IOException ioe){
               System.out.println("Failed to open "+path+name);  
               System.exit(1);
            } 
        }
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public Index(String index_file, long n, int K_size)throws IOException
    {
        totallength=n;
        K=K_size;
        dim=(int)Math.ceil(K/16.0);
        read(index_file,totallength);
        System.out.println(totallength+" kmers loaded.");
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public Index(Index old_index, Index new_index, GraphDatabaseService graphDb, int K_size)
    {
        int k,comp;
        long i,j,inx,new_kmer_num,old_kmer_num;
        pointer[][] tmp_pointers;
        new_kmer_num=new_index.totallength;
        old_kmer_num=old_index.totallength;
        totallength=new_kmer_num+old_kmer_num;
        K=K_size;
        dim=(int)Math.ceil(K/16.0);
        parts=(int)(totallength/2100000000)+1;
        tmp_pointers=new pointer[parts][];
        length=new int[parts];
        for(k=0;k<parts;++k)
        {
            length[k]=k<parts-1?2100000000:(int)(totallength%2100000000);
            tmp_pointers[k]=new pointer[length[k]];
        }
        System.out.println("merging indexes...");
        for(i=j=inx=0;i<new_kmer_num && j<old_kmer_num;inx++)
        {
            comp=new_index.pointers[(int)(i/2100000000)][(int)(i%2100000000)].compare_kmer(old_index.pointers[(int)(j/2100000000)][(int)(j%2100000000)],dim);
            if(comp==-1)
            {
                tmp_pointers[(int)(inx/2100000000)][(int)(inx%2100000000)]=new_index.pointers[(int)(i/2100000000)][(int)(i%2100000000)];
                ++i;
            }
            else if(comp==1)
            {
                tmp_pointers[(int)(inx/2100000000)][(int)(inx%2100000000)]=old_index.pointers[(int)(j/2100000000)][(int)(j%2100000000)];
                ++j;
            }
            else
            {
                tmp_pointers[(int)(inx/2100000000)][(int)(inx%2100000000)]=old_index.pointers[(int)(j/2100000000)][(int)(j%2100000000)];
                ++i;
                ++j;
            }
        }
        while(i<new_kmer_num)
            {
                tmp_pointers[(int)(inx/2100000000)][(int)(inx%2100000000)]=new_index.pointers[(int)(i/2100000000)][(int)(i%2100000000)];
                ++i;
                ++inx;
            }
        while(j<old_kmer_num)
            {
                tmp_pointers[(int)(inx/2100000000)][(int)(inx%2100000000)]=old_index.pointers[(int)(j/2100000000)][(int)(j%2100000000)];
                ++j;
                ++inx;
            }
        totallength=inx;
        parts=(int)(totallength/2100000000)+1;
        pointers=new pointer[parts][];
        length=new int[parts];
        for(k=0;k<parts;++k)
        {
            length[k]=k<parts-1?2100000000:(int)(totallength%2100000000);
            pointers[k]=new pointer[length[k]];
        }
        for(inx=0;inx<totallength;++inx)
            pointers[(int)(inx/2100000000)][(int)(inx%2100000000)]=tmp_pointers[(int)(inx/2100000000)][(int)(inx%2100000000)];
        for(k=0;k<parts;++k)
            tmp_pointers[k]=null;
        System.out.println(totallength+" kmers indexed.");
        tmp_pointers=null;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    long find(pointer k_mer)
    {
        int low, high,i=0;
        for(i=0;i<parts && k_mer.compare_kmer(pointers[i][length[i]-1], dim)==1;++i);
        low=0;
        high=length[i]-1; 
        int mid;
        int comp;
        while(low<=high)
        {
            mid=(low+high)/2;
            comp=k_mer.compare_kmer(pointers[i][mid],dim);
            if(comp==-1)
                high=mid-1;
            else if(comp==1)
                low=mid+1;
            else
                return i*2100000000+mid;
        }
        System.out.println("Unable to find kmer!");
        return -1; // not found
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    void jump()
    {
        int j;
        f_bin=binary[Genomes[genome][sequence][position]];
        r_bin=3-f_bin;
        do{
            while(f_bin==4 && position<seq_len)
            {
                ++position;
                f_bin=binary[Genomes[genome][sequence][position]];
                r_bin=3-f_bin;
            }
            f_kmer=new pointer(dim);
            r_kmer=new pointer(dim);
            f_kmer.next_up_kmer(f_bin,dim,K);
            r_kmer.next_down_kmer(r_bin,dim,K);
            for(j=0; j<K-1 && position<seq_len;++j)
            {
                ++position;
                f_bin=binary[Genomes[genome][sequence][position]];
                r_bin=3-f_bin;
                if(f_bin==4 )
                    break;
                f_kmer.next_up_kmer(f_bin,dim,K);
                r_kmer.next_down_kmer(r_bin,dim,K);
            }
        }while(f_bin==4 && position<seq_len);
        canonical=f_kmer.compare_kmer(r_kmer,dim)==-1;
        kmer=canonical?f_kmer:r_kmer;
        if(position==seq_len)
            finish=true;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    void next_kmer()
    {
        ++position;
        f_bin=binary[Genomes[genome][sequence][position]];
        r_bin=3-f_bin;
        f_kmer.next_up_kmer(f_bin,dim,K);
        r_kmer.next_down_kmer(r_bin,dim,K);
        canonical=f_kmer.compare_kmer(r_kmer,dim)==-1;
        kmer=canonical?f_kmer:r_kmer;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void initial_kmers()
    {
        int i;
        f_kmer=new pointer(dim);
        r_kmer=new pointer(dim);
        for(i=0;i<K && position<seq_len;++i)
        {
            if(binary[Genomes[genome][sequence][position+1]]== 4 )
            {
                ++position;
                jump();
                break;
            }
            next_kmer();
        }   
        canonical=f_kmer.compare_kmer(r_kmer,dim)==-1;
        kmer=canonical?f_kmer:r_kmer;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    private void generate_short_kmers()
    {
        int i,j;
        boolean kmers[]=new boolean[1<<(2*K)];
        for(i=0;i<kmers.length;++i)
            kmers[i]=false;
        totallength=0;
        for(genome=Pantools.old_num_genomes+1;genome<=Pantools.num_genomes;++genome) 
            for(sequence=1;sequence<=Pantools.num_sequences[genome];++sequence) 
            {
                position=-1;
                seq_len=Pantools.sequence_length[genome][sequence]-1;
                initial_kmers();
                finish=false;
                while(!finish)
                {
                    if(!kmers[kmer.kmer[0]])
                    {
                        ++totallength;
                        kmers[kmer.kmer[0]]=true;
                    }
                    if(position==seq_len)
                        break;
                    if(binary[Genomes[genome][sequence][position+1]]==4)
                    {
                        ++position;
                        jump();
                        kmer=canonical?f_kmer:r_kmer;
                    }
                    else
                        next_kmer();
                }//while
            }
        parts=1;
        pointers=new pointer[parts][];
        length=new int[parts];
        length[0]=(int)totallength;
        pointers[0]=new pointer[length[0]];
        for(i=j=0;i<kmers.length;++i)
            if(kmers[i])
            {
                pointers[0][j]=new pointer(dim);
                pointers[0][j++].kmer[0]=i;
            }
        kmers=null;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public int byteArrayToInt(byte[] b) 
    {
        return   b[3] & 0xFF |
                (b[2] & 0xFF) << 8 |
                (b[1] & 0xFF) << 16 |
                (b[0] & 0xFF) << 24;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void intToByteArray(byte[] data, int a)
    {
        data[0]=(byte) ((a >> 24) & 0xFF);
        data[1]=(byte) ((a >> 16) & 0xFF);
        data[2]=(byte) ((a >> 8) & 0xFF);
        data[3]=(byte) (a & 0xFF);
    } 
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public long byteArrayToLong(byte[] b) 
    {
        return   b[7] & 0xFF |
                (b[6] & 0xFF) << 8  |
                (b[5] & 0xFF) << 16 |
                (b[4] & 0xFF) << 24 |
                (b[3] & 0xFF) << 32 |
                (b[2] & 0xFF) << 40 |
                (b[1] & 0xFF) << 48 |
                (b[0] & 0xFF) << 56;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void longtToByteArray(byte[] data, long a)
    {
        data[0]=(byte) ((a >> 56) & 0xFF);
        data[1]=(byte) ((a >> 48) & 0xFF);   
        data[2]=(byte) ((a >> 40) & 0xFF);   
        data[3]=(byte) ((a >> 32) & 0xFF);   
        data[4]=(byte) ((a >> 24) & 0xFF);
        data[5]=(byte) ((a >> 16) & 0xFF);   
        data[6]=(byte) ((a >> 8) & 0xFF);  
        data[7]=(byte) (a & 0xFF);
    } 
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void write(String index_name) throws IOException
    {
        System.out.println("Writing index...");
        pointer p;
        int j;
        long i;
        byte[] data=new byte[8];
        try{
            FileOutputStream index_out = new FileOutputStream(index_name);
            for(i = 0; i < this.totallength; ++i)
            {
                p=this.pointers[(int)(i/2100000000)][(int)(i%2100000000)];
                for(j=0;j<dim;++j)
                {
                    intToByteArray(data,this.pointers[(int)(i/2100000000)][(int)(i%2100000000)].kmer[j]);
                    index_out.write(data,0,4);
                }
                longtToByteArray(data,p.node_id);
                index_out.write(data,0,8);
                index_out.write(p.format);
                intToByteArray(data,p.position);
                index_out.write(data,0,4);
                longtToByteArray(data,p.next_index);
                index_out.write(data,0,8);
                this.pointers[(int)(i/2100000000)][(int)(i%2100000000)]=null;
            }
            index_out.close();
        }
        catch(IOException ioe){
           System.out.println("Failed to read file names!");  
        }        
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void read(String index_name, long len) throws IOException
    {
        System.out.println("Reading available index...");
        int j,k;
        long i;
        byte[] int_record = new byte[4];
        byte[] long_record = new byte[8];
        totallength=len;
        parts=(int)(totallength/2100000000)+1;
        pointers=new pointer[parts][];
        length=new int[parts];
        for(k=0;k<parts;++k)
        {
            length[k]=k<parts-1?2100000000:(int)(totallength%2100000000);
            pointers[k]=new pointer[length[k]];
        }
        try (FileInputStream index_in = new FileInputStream(index_name)) {
            for(i = 0; i < len; ++i)
            {
                pointers[(int)(i/2100000000)][(int)(i%2100000000)]=new pointer(dim);
                for(j=0;j<dim;++j)
                {
                    index_in.read(int_record);
                    pointers[(int)(i/2100000000)][(int)(i%2100000000)].kmer[j]= byteArrayToInt(int_record);
                }
                index_in.read(long_record);
                pointers[(int)(i/2100000000)][(int)(i%2100000000)].node_id=byteArrayToLong(long_record);
                pointers[(int)(i/2100000000)][(int)(i%2100000000)].format=(byte)index_in.read();
                index_in.read(int_record);
                pointers[(int)(i/2100000000)][(int)(i%2100000000)].position=byteArrayToInt(int_record);
                index_in.read(long_record);
                pointers[(int)(i/2100000000)][(int)(i%2100000000)].next_index=byteArrayToLong(long_record);
                if(i%(len/100+1)==0) 
                    System.out.print((long)i*100/len+1+"%\r");
            }
        }
        catch(IOException ioe){
           System.out.println("Failed to read file names!");  
        }     
    }  
}


