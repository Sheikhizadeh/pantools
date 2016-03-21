/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pantools;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import org.neo4j.graphdb.DynamicLabel;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 *
 * @author siavashs
 */
public class index_database {
    private kmer key;
    private long[] prefix_ptr;
    private int ptr_parts_num;
    private long[] ptr_parts_size;
    private int suf_parts_num;
    private long[] suf_parts_size;
    private int max_byte=2100000000;
    private int header_pos;
    
    private int K;
    private int mode;
    private int ctr_size;
    private int pre_len;
    private int min_count;
    private int max_count;
    private long kmers_num;

    private int suf_rec_size;
    private int suf_len;
    private final int ptr_len=21;
    private RandomAccessFile suf_file;
    private RandomAccessFile pre_file;  
    private MappedByteBuffer[] suf_buff;
    private RandomAccessFile ptr_file;
    private MappedByteBuffer[] ptr_buff;
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public index_database(String path,String file, int K_size,genome_database genomeDb)
    {
        int cores=Runtime.getRuntime().availableProcessors()/2+1;
        int k;
        long p;
        byte[] minus_one=null;
        try{
            if(file!=null) // we are not loading an available index
            {
                Files.createDirectory(Paths.get(path));
                System.out.println("Running KMC2...                      ");
                executeCommand("kmc -r -k"+K_size+" -t"+cores+" -ci1 -fm "+(genomeDb.num_genomes>1?"@"+file.trim():genomeDb.genome_names[1])+" "+path+"/kmers "+path);
                String output=executeCommand("kmc_tools sort "+path+"/kmers "+path+"/sorted");
                if(output.startsWith("This database contains sorted k-mers already!"))
                {
                    Files.copy(Paths.get(path+"/kmers.kmc_pre"), Paths.get(path+"/sorted.kmc_pre"));
                    Files.copy(Paths.get(path+"/kmers.kmc_suf"), Paths.get(path+"/sorted.kmc_suf"));
                }
                minus_one=new byte[max_byte];
                Arrays.fill( minus_one, (byte)(-1) );
            }
            pre_file = new RandomAccessFile(path+"/sorted.kmc_pre","r");
            pre_file.seek(pre_file.length()-8);
            header_pos=read_int(pre_file);
            pre_file.seek(pre_file.length()-8-header_pos);
        // read the index properties    
            K=read_int(pre_file);
            mode=read_int(pre_file);
            ctr_size=read_int(pre_file);
            pre_len=read_int(pre_file);
            min_count=read_int(pre_file);
            max_count=read_int(pre_file);
            kmers_num=read_long(pre_file);
            suf_len=K-pre_len;
            key=new kmer(K,pre_len,suf_len);
            System.out.println("Indexing "+kmers_num+" kmers...                    ");
        // load the prefix file into the memory    
            pre_file.seek(4);
            int q,len=1<<(2*pre_len);
            prefix_ptr=new long[len];
            MappedByteBuffer pre_buff;
            for(q=0,p=0;p<8;++p)
            {
                pre_buff=pre_file.getChannel().map(FileChannel.MapMode.READ_ONLY, 4+p*len, len);
                for(k=0;k<len/8;++k,++q)
                {
                    prefix_ptr[q]=read_long(pre_buff);
                }
            }
            pre_file.close();
            pre_buff=null;
        // mapping suffix file into the memory
            suf_rec_size=ctr_size+suf_len/4;
            max_byte=max_byte/suf_rec_size*suf_rec_size;
            suf_parts_num=(int)((kmers_num*suf_rec_size)%max_byte==0?(kmers_num*suf_rec_size)/max_byte:(kmers_num*suf_rec_size)/max_byte+1);
            suf_parts_size=new long[suf_parts_num];
            suf_file = new RandomAccessFile(path+"/sorted.kmc_suf","r");
            suf_buff= new MappedByteBuffer[suf_parts_num];
            for(k=0;k<suf_parts_num;++k)
            {
                suf_parts_size[k]=(int)(k==suf_parts_num-1?(kmers_num*suf_rec_size)%max_byte:max_byte);
                suf_buff[k]=suf_file.getChannel().map(FileChannel.MapMode.READ_ONLY, 4+k*suf_parts_size[0], suf_parts_size[k]);
            }
        // mapping pointers file into the memory
            max_byte=max_byte/ptr_len*ptr_len;
            ptr_parts_num=(int)((kmers_num*ptr_len)%max_byte==0?(kmers_num*ptr_len)/max_byte:(kmers_num*ptr_len)/max_byte+1);
            ptr_parts_size=new long[ptr_parts_num];
            ptr_file = new RandomAccessFile(path+"/pointers.db", "rw");
            ptr_buff= new MappedByteBuffer[ptr_parts_num];
            for(k=0;k<ptr_parts_num;++k)
            {
                ptr_parts_size[k]=(int)(k==ptr_parts_num-1?(kmers_num*ptr_len)%max_byte:max_byte);
                ptr_buff[k]=ptr_file.getChannel().map(FileChannel.MapMode.READ_WRITE, k*ptr_parts_size[0], ptr_parts_size[k]);
                if(file!=null) // we are not lading an available index
                    ptr_buff[k].put(minus_one,0,(int)ptr_parts_size[k]);
            }
        }catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }         
        minus_one=null;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public index_database(String path,String file,genome_database genomeDb, GraphDatabaseService graphDb)
    {
        int cores=Runtime.getRuntime().availableProcessors()/2+1;
        int k,p;
        long c_index,p_index,l;
        Node node;
        ResourceIterator<Node> nodes;
    // move current index files to directory old_index
        try{
            Files.createDirectory(Paths.get(path+"/old_index"));
            Files.move(Paths.get(path+"/kmers.kmc_pre"),  Paths.get(path+"/old_index/kmers.kmc_pre"));
            Files.move(Paths.get(path+"/kmers.kmc_suf"),  Paths.get(path+"/old_index/kmers.kmc_suf"));
            Files.move(Paths.get(path+"/sorted.kmc_pre"), Paths.get(path+"/old_index/sorted.kmc_pre"));
            Files.move(Paths.get(path+"/sorted.kmc_suf"), Paths.get(path+"/old_index/sorted.kmc_suf"));
            Files.move(Paths.get(path+"/pointers.db"),    Paths.get(path+"/old_index/pointers.db"));
        // load old_index
            index_database old_index=new index_database(path+"/old_index",null,0,genomeDb);
            K=old_index.K;
        // make new index for new genomes
            System.out.println("Running KMC2...                      ");
           executeCommand("kmc -r -k"+K+" -t"+cores+" -ci1 -fm "+(genomeDb.num_genomes-genomeDb.previous_num_genomes>1?"@"+file.trim():genomeDb.genome_names[genomeDb.previous_num_genomes+1])+" "+path+"/new_kmers "+path);
        // merge two indeces    
            executeCommand("kmc_tools union "+path+"/old_index/kmers "+path+"/new_kmers "+path+"/sorted");
            pre_file = new RandomAccessFile(path+"/sorted.kmc_pre","r");
            pre_file.seek(pre_file.length()-8);
            header_pos=read_int(pre_file);
            pre_file.seek(pre_file.length()-8-header_pos);
        // read the merged index properties    
            K=read_int(pre_file);
            mode=read_int(pre_file);
            ctr_size=read_int(pre_file);
            pre_len=read_int(pre_file);
            min_count=read_int(pre_file);
            max_count=read_int(pre_file);
            kmers_num=read_long(pre_file);
            System.out.println("Indexing "+kmers_num+" kmers...                    ");
        // load the prefix file into the memory    
            pre_file.seek(4);
            int q,len=1<<(2*pre_len);
            prefix_ptr=new long[len];
            MappedByteBuffer pre_buff;
            for(q=0,p=0;p<8;++p)
            {
                pre_buff=pre_file.getChannel().map(FileChannel.MapMode.READ_ONLY, 4+p*len, len);
                for(k=0;k<len/8;++k,++q)
                {
                    prefix_ptr[q]=read_long(pre_buff);
                    //System.out.println(q+" "+prefix_ptr[q]);
                }
            }
            pre_file.close();
            pre_buff=null;
        // mapping suffix file into the memory    
            suf_len=K-pre_len;
            suf_rec_size=ctr_size+suf_len/4;
            max_byte=max_byte/suf_rec_size*suf_rec_size;
            suf_parts_num=(int)((kmers_num*suf_rec_size)%max_byte==0?(kmers_num*suf_rec_size)/max_byte:(kmers_num*suf_rec_size)/max_byte+1);
            suf_parts_size=new long[suf_parts_num];
            suf_file = new RandomAccessFile(path+"/sorted.kmc_suf","r");
            suf_buff= new MappedByteBuffer[suf_parts_num];
            key=new kmer(K,pre_len,suf_len);
            for(k=0;k<suf_parts_num;++k)
            {
                suf_parts_size[k]=(int)(k==suf_parts_num-1?(kmers_num*suf_rec_size)%max_byte:max_byte);
                suf_buff[k]=suf_file.getChannel().map(FileChannel.MapMode.READ_ONLY, 4+k*suf_parts_size[0], suf_parts_size[k]);
            }
        // mapping pointers file into the memory    
            max_byte=max_byte/ptr_len*ptr_len;
            ptr_parts_num=(int)((kmers_num*ptr_len)%max_byte==0?(kmers_num*ptr_len)/max_byte:(kmers_num*ptr_len)/max_byte+1);
            ptr_parts_size=new long[ptr_parts_num];
            ptr_file = new RandomAccessFile(path+"/pointers.db", "rw");
            ptr_buff= new MappedByteBuffer[ptr_parts_num];
            byte[] minus_one=new byte[max_byte];
            Arrays.fill( minus_one, (byte)(-1) );
            for(k=0;k<ptr_parts_num;++k)
            {
                ptr_parts_size[k]=(int)(k==ptr_parts_num-1?(kmers_num*ptr_len)%max_byte:max_byte);
                ptr_buff[k]=ptr_file.getChannel().map(FileChannel.MapMode.READ_WRITE, k*ptr_parts_size[0], ptr_parts_size[k]);
                ptr_buff[k].put(minus_one,0,(int)ptr_parts_size[k]);
            }
            minus_one=null;

        // adjusting available pointers
            try(Transaction tx = graphDb.beginTx()){
            nodes = graphDb.findNodes( DynamicLabel.label( "node" ));
            pointer ptr=new pointer();
            byte[] byte_ptr=new byte[ptr_len];
            for(;nodes.hasNext();)
            {
                node=nodes.next();
                l=(long)node.getProperty("first_kmer");
                p_index=find(old_index.get_kmer(l));
                old_index.get_pointer(byte_ptr,ptr,l);
                put_pointer(ptr,p_index);
                node.setProperty("first_kmer", p_index);
                for(l=old_index.get_next_index(l);l!=-1L;l=old_index.get_next_index(l))
                {
                    c_index=find(old_index.get_kmer(l));
                    old_index.get_pointer(byte_ptr,ptr,l);
                    put_pointer(ptr,c_index);
                    put_next_index(c_index,p_index);
                    p_index=c_index;
                }
                put_next_index(-1L,p_index);
                node.setProperty("last_kmer", p_index);
            }
            tx.success();} 
            old_index.close();
            Files.delete(Paths.get(path+"/old_index/kmers.kmc_pre"));
            Files.delete(Paths.get(path+"/old_index/kmers.kmc_suf"));
            Files.delete(Paths.get(path+"/old_index/sorted.kmc_pre"));
            Files.delete(Paths.get(path+"/old_index/sorted.kmc_suf"));
            Files.delete(Paths.get(path+"/old_index/pointers.db"));
            Files.delete(Paths.get(path+"/new_kmers.kmc_pre"));
            Files.delete(Paths.get(path+"/new_kmers.kmc_suf"));
        }catch (IOException e) {
            System.out.println(e.getMessage()+"\nFailed to make index!");
            System.exit(1);
        }        
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public long length()
    {
        return kmers_num;
    } 
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public int get_pre_len()
    {
        return pre_len;
    }    
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public int get_suf_len()
    {
        return suf_len;
    } 
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void close()throws IOException
    {
        ptr_file.close();
        suf_file.close();
        for(int k=0;k<ptr_parts_num;++k)
            ptr_buff[k]=null;
        for(int k=0;k<suf_parts_num;++k)
            suf_buff[k]=null;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    int read_int(RandomAccessFile file)throws IOException
    {
        int number=0;
        for (int i=0;i<4;++i)
            number+=(file.readByte() & 0xFF)<<(8*i);
        return number;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    long read_long(RandomAccessFile file)throws IOException
    {
        long number=0;
        for (int i=0;i<8;++i)
            number+=((long)(file.readByte() & 0xFF))<<(8*i);
        return number;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    long read_long(MappedByteBuffer buff)throws IOException
    {
        long number=0;
        for (int i=0;i<8;++i)
            number+=((long)(buff.get() & 0xFF))<<(8*i);
        return number;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public kmer get_kmer(long i)
    {
        kmer k_mer=new kmer(K,pre_len,suf_len);
        boolean found=false;
        int low=0, high=(1<<(2*pre_len))-1, mid=0;
        while(low<=high && !found)
        {
            mid=(low+high)/2;
            if(i<prefix_ptr[mid])
                high=mid-1;
            else if(i>prefix_ptr[mid])
                low=mid+1;
            else
                found=true;
        }
        if(!found)
            mid=high;
        else
            for(;mid<prefix_ptr.length-1 && prefix_ptr[mid+1]==prefix_ptr[mid];++mid);
        k_mer.prefix=mid;
        suf_buff[(int)(i*suf_rec_size/suf_parts_size[0])].position((int)(i*suf_rec_size%suf_parts_size[0]));
        suf_buff[(int)(i*suf_rec_size/suf_parts_size[0])].get(k_mer.suffix,0,suf_len/4);
        return k_mer;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void get_pointer(byte[] b, pointer p, long i)
    {
        p.node_id=   ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])].getLong((int)(i*ptr_len%ptr_parts_size[0]));
        p.canonical=    ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])].get((int)(i*ptr_len%ptr_parts_size[0]+8))==0;
        p.position=  ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])]. getInt((int)(i*ptr_len%ptr_parts_size[0]+9));
        p.next_index=ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])].getLong((int)(i*ptr_len%ptr_parts_size[0]+13));
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void put_pointer(pointer p, long i)
    {
        ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])].putLong((int)(i*ptr_len%ptr_parts_size[0]),   p.node_id);
        ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])].    put((int)(i*ptr_len%ptr_parts_size[0]+8), (byte)(p.canonical?0:1));
        ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])]. putInt((int)(i*ptr_len%ptr_parts_size[0]+9), p.position);
        ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])].putLong((int)(i*ptr_len%ptr_parts_size[0]+13),p.next_index);
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public long get_node_id(long i)
    {
        return ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])].getLong((int)(i*ptr_len%ptr_parts_size[0]));
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void put_node_id(long node_id, long i)
    {
        ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])].putLong((int)(i*ptr_len%ptr_parts_size[0]), node_id);
    }     
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public boolean get_canonical(long i)
    {
        return ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])].get((int)(i*ptr_len%ptr_parts_size[0]+8))==0;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void put_canonical(boolean side, long i)
    {
        ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])].put((int)(i*ptr_len%ptr_parts_size[0]+8),(byte)(side?0:1));
    }     
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public int get_position(long i)
    {
        return ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])].getInt((int)(i*ptr_len%ptr_parts_size[0]+9));
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void put_position(int pos, long i)
    {
        ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])].putInt((int)(i*ptr_len%ptr_parts_size[0]+9),pos);
    }  
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public long get_next_index(long i)
    {
        return ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])].getLong((int)(i*ptr_len%ptr_parts_size[0]+13));
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public void put_next_index(long next_index, long i)
    {
        ptr_buff[(int)(i*ptr_len/ptr_parts_size[0])].putLong((int)(i*ptr_len%ptr_parts_size[0]+13),next_index);
    }       
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    long find(kmer k_mer)
    {
        long low, mid, high;
        int j,comp;
        low=prefix_ptr[k_mer.prefix];
        if(k_mer.prefix==prefix_ptr.length-1)
            high=kmers_num-1;
        else
            high=prefix_ptr[k_mer.prefix+1]-1;
        while(low<=high)
        {
            mid=(low+high)/2;
            j=(int)(mid*suf_rec_size/suf_parts_size[0]);
            suf_buff[j].position((int)(mid*suf_rec_size%suf_parts_size[0]));
            suf_buff[j].get(key.suffix, 0, suf_len/4);
            comp=k_mer.compare_suffix(key);
            if(comp==-1)
                high=mid-1;
            else if(comp==1)
                low=mid+1;
            else
                return mid;
        }
        System.out.print("Failed to find kmer: ");
        System.out.print(k_mer.toString());
        System.exit(0);
        return -1; // not found
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public int b2i(byte[] b, int i) 
    {
        return   (b[i+3] & 0xFF) |
                ((b[i+2] & 0xFF) << 8) |
                ((b[i+1] & 0xFF) << 16) |
                ((b[i] & 0xFF) << 24);
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public byte[] i2b(int a)
    {
        byte[] data=new byte[4];
        data[0]=(byte) ((a >> 24) & 0xFF);
        data[1]=(byte) ((a >> 16) & 0xFF);
        data[2]=(byte) ((a >> 8) & 0xFF);
        data[3]=(byte) (a & 0xFF);
        return data;
    } 
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public long b2l(byte[] b, int i) 
    {
        return   (b[i+7] & 0xFF) |
                ((b[i+6] & 0xFF) << 8)  |
                ((b[i+5] & 0xFF) << 16) |
                ((b[i+4] & 0xFF) << 24) |
                ((b[i+3] & 0xFF) << 32) |
                ((b[i+2] & 0xFF) << 40) |
                ((b[i+1] & 0xFF) << 48) |
                ((b[i] & 0xFF) << 56);
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    public byte[] l2b(long a)
    {
        byte[] data=new byte[8];
        data[0]=(byte) ((a >> 56) & 0xFF);
        data[1]=(byte) ((a >> 48) & 0xFF);   
        data[2]=(byte) ((a >> 40) & 0xFF);   
        data[3]=(byte) ((a >> 32) & 0xFF);   
        data[4]=(byte) ((a >> 24) & 0xFF);
        data[5]=(byte) ((a >> 16) & 0xFF);   
        data[6]=(byte) ((a >> 8) & 0xFF);  
        data[7]=(byte) (a & 0xFF);
        return data;
    } 
    public static String executeCommand(String command) {
        StringBuilder exe_output=new StringBuilder();
        String line = "";			
        Process p;
        try {
                p = Runtime.getRuntime().exec(command);
                p.waitFor();
                BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
                while ((line = reader.readLine())!= null) {
                        exe_output.append(line + "\n");
                }

        } catch (Exception e) {
                e.printStackTrace();
        }
        return exe_output.toString();
    }    
}


