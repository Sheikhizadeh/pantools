/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pantools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import static pantools.Pantools.PATH;
import static pantools.Pantools.complement;
import static pantools.Pantools.executeCommand;

/**
 *
 * @author sheik005
 */
public class genome_bank {
    public long size;
    public int num_genomes;
    public int num_sequences[];        // Number of sequences in each genome    
    public String[][] sequence_titles; // Name of sequences for each genome
    public String[] genome_names;     // Paths to the genome FASTA files
    public long genome_length[];    // Length of sequences for each genome
    public long sequence_length[][];    // Length of sequences for each genome
    public long sequence_start[][];    // Length of sequences for each genome    
    private RandomAccessFile genomes_file;
    private MappedByteBuffer[] genomes_buff;
    private byte[] binary;
    private byte[] complement;
    private int parts_num;
    private int[] parts_size;
    public int previous_num_genomes;
    
    public genome_bank(String file) throws IOException
    {
        int g,s,i;
        previous_num_genomes=0;
        complement=new byte[124];
        binary=new byte[124];
        for(i=0;i<124;++i)
            binary[i]=(byte)(i-62);
        binary[65] = binary[97]= 0; // A , a
        binary[67] = binary[99]= 1; // C , c
        binary[71] = binary[103]= 2;// G , g
        binary[84] = binary[116]= 3;// T , t
        complement[0]=3; // A T
        complement[3]=0; // T A
        complement[1]=2; // C G
        complement[2]=1; // G C
        complement[20]=27; // R Y
        complement[27]=20; // Y R
        complement[21]=21; // S S
        complement[25]=25; // W W
        complement[13]=15; // K M
        complement[15]=13; // M K
        complement[4]=24; // B V
        complement[24]=4; // V B
        complement[6]=10; // D H
        complement[10]=6; // H D
        complement[16]=16; // N N     
        File f = new File(PATH+"/genomes.db");
        if(f.exists())  
        {
            genomes_file=new RandomAccessFile(PATH+"/genomes.db","rw");
            genomes_file.seek(f.length());
            BufferedReader in;
            try{
                in = new BufferedReader(new FileReader(PATH+"/genomes.info")); 
                size=Long.valueOf(in.readLine().split(":")[1]);
                previous_num_genomes=Integer.parseInt(in.readLine().split(":")[1]);
                num_genomes=previous_num_genomes+Integer.parseInt(executeCommand("wc -l "+file).trim().split("\\s")[0]);
                genome_names=new String[num_genomes+1];
                genome_length=new long[num_genomes+1];
                sequence_titles=new String[num_genomes+1][];
                sequence_length=new long[num_genomes+1][];
                sequence_start=new long[num_genomes+1][];
                num_sequences=new int[num_genomes+1];
                for(g=1;g<=previous_num_genomes;++g)
                {
                    genome_length[g]=Long.valueOf(in.readLine().split(":")[1]);
                    num_sequences[g]=Integer.parseInt(in.readLine().split(":")[1]);
                    for(s=1;s<=num_sequences[g];++s)
                    {
                        sequence_titles[g][s]=in.readLine().split(":")[1];
                        sequence_length[g][s]=Long.valueOf(in.readLine().split(":")[1]);
                        sequence_start[g][s]=Long.valueOf(in.readLine().split(":")[1]);
                    }
                }
                in.close();
            }
            catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            }            
        }
        else
        {
            genomes_file=new RandomAccessFile(PATH+"/genomes.db","rw");
            size=0;
            num_genomes=Integer.parseInt(executeCommand("wc -l "+file).trim().split("\\s")[0]);
            genome_names=new String[num_genomes+1];
            genome_length=new long[num_genomes+1];
            sequence_titles=new String[num_genomes+1][];
            sequence_length=new long[num_genomes+1][];
            sequence_start=new long[num_genomes+1][];
            num_sequences=new int[num_genomes+1];
        }
        try{
            BufferedReader in = new BufferedReader(new FileReader(file));
            for(g=previous_num_genomes+1;g<=num_genomes;++g)
            {
                genome_names[g]=in.readLine();
                num_sequences[g]=Integer.parseInt(executeCommand("grep -c '>' "+genome_names[g]).trim());
                sequence_titles[g]=new String[num_sequences[g]+1];
                sequence_length[g]=new long[num_sequences[g]+1];
                sequence_start[g]=new long[num_sequences[g]+1];
            }
            load_genomes();
        }
        catch(FileNotFoundException ioe){
           System.out.println("Failed to read file names: "+ioe);  
           System.exit(1);
        } 
        try{
            BufferedWriter out = new BufferedWriter(new FileWriter(PATH+"/genomes.info"));
            out.write("number_of_bases:"+size+"\n");
            out.write("number_of_genomes:"+num_genomes+"\n");
            for(g=1;g<=num_genomes;++g)
            {
                out.write("genome length:"+genome_length[g]+"\n");
                out.write("number_of_sequences:"+num_sequences[g]+"\n");
                for(s=1;s<=num_sequences[g];++s)
                {
                    out.write("sequence title:"+sequence_titles[g][s]+"\n");
                    out.write("sequence length:"+sequence_length[g][s]+"\n");
                    out.write("sequence start:"+sequence_start[g][s]+"\n");
                }
            }
            out.close();
        }
        catch(FileNotFoundException ioe){
           System.out.println("Failed to read file names: "+ioe);  
           System.exit(1);
        }        
    }
    /*
    This function loads genomes'sequences to Genome array. 
    */      
    public void load_genomes()  throws IOException
    {
        String line;
        int i=0,j,g,s,len,k,max_byte=2100000000,l;
        BufferedReader in;
        System.out.println("Reading "+num_genomes+" genomes...");
        for(g=previous_num_genomes+1;g<=num_genomes;++g) 
        {
            try{
                in = new BufferedReader(new FileReader(genome_names[g]));
                s=0;
                len=0;
                while(in.ready())
                {
                    line=in.readLine();
                    l=line.length();
                    if (line==null) 
                        continue;
                    if(line.charAt(0)=='>')
                    {
                        ++s;
                        sequence_titles[g][s]=line.substring(1);
                        sequence_start[g][s]=size;
                        len=0;
                    }
                    else
                    {
                        len += l;
                        sequence_length[g][s]+=l;
                        genome_length[g]+=l;
                        size+=l;
                    }
                }
                in.close();
            }
            catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            }
        } 
        parts_num=(int)(size%max_byte==0?size/max_byte:size/max_byte+1);
        parts_size=new int[parts_num];
        genomes_buff= new MappedByteBuffer[parts_num];
        for(k=0;k<parts_num;++k)
        {
            parts_size[k]=(int)(k==parts_num-1?size%max_byte:max_byte);
            genomes_buff[k]=genomes_file.getChannel().map(FileChannel.MapMode.READ_WRITE, k*parts_size[0], parts_size[k]);
        } 
        for(g=previous_num_genomes+1;g<=num_genomes;++g) 
        {
            try{
                in = new BufferedReader(new FileReader(genome_names[g]));
                s=0;
                while(in.ready())
                {
                    line=in.readLine().toUpperCase();
                    if (line==null) 
                        continue;
                    if(line.charAt(0)=='>')
                    {
                        i=0;
                        ++s;
                        System.out.print("Sequence "+s+"/"+num_sequences[g]+" of genome "+g+"         \r");
                    }
                    else
                    {
                        for(j=0;j<line.length();++j,++i)
                            put(g,s,i,binary[line.charAt(j)]);
                    }
                }
                in.close();
            }
            catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            }
        } 

    }  
    ///////////////////////////////////////////////////////////////
    public byte get(int g, int s, int p)
    {
        long position=sequence_start[g][s]+p;
        return genomes_buff[(int)(position/parts_size[0])].get((int)(position%parts_size[0]));
    }
    ///////////////////////////////////////////////////////////////
    public void put(int g, int s, int p, byte b)
    {
        long position=sequence_start[g][s]+p;
        genomes_buff[(int)(position/parts_size[0])].put((int)(position%parts_size[0]),b);
    }
    /*
    This function determines the equality of two subsequences of s1 and s2 of length len starting at start1 and start2, respectively.
    forward determines the direction of comparion.
    */
    public boolean compare(int[] a1, int[] a2, int len, boolean forward)
    {
        int i;
        boolean equal=true;
        if(forward)
            for(i=0;i<len && equal;++i)
                if(get(a1[0],a1[1],a1[2]+i)!=get(a2[0],a2[1],a2[2]+i))
                    equal=false;
        else
            for(i=0;i<len && equal;++i)
                if(get(a1[0],a1[1],a1[2]+i)!=complement[get(a2[0],a2[1],a2[2]+len-i-1)])
                    equal=false;
        return equal;
    }
}
