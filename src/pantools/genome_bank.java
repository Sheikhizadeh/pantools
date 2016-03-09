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
import static pantools.Pantools.executeCommand;
import static pantools.Pantools.sym;

/**
 *
 * @author sheik005
 */
public class genome_bank {
    public long num_bytes;
    public int num_genomes;
    public int num_sequences[];        // Number of sequences in each genome    
    public String[][] sequence_titles; // Name of sequences for each genome
    public String[] genome_names;     // Paths to the genome FASTA files
    public long genome_length[];    // Length of sequences for each genome
    public long sequence_length[][];    // Length of sequences for each genome
    public long sequence_start[][];    // Length of sequences for each genome    
    private RandomAccessFile genomes_file;
    private MappedByteBuffer[] genomes_buff;
    private int[] binary;
    private int[] complement;
    private int parts_num;
    private long[] parts_size;
    public int previous_num_genomes;
    
    public genome_bank(String file) throws IOException
    {
        int g,s,i;
        previous_num_genomes=0;
        complement=new int[]{3,2,1,0,9,8,6,7,5,4,13,12,11,10,14};
        binary=new int[256];
        for(i=0;i<256;++i)
            binary[i]=14;
        binary['A'] = 0; 
        binary['C'] = 1; 
        binary['G'] = 2; 
        binary['T'] = 3; 
        binary['M'] = 4;  
        binary['R'] = 5;  
        binary['W'] = 6;  
        binary['S'] = 7;  
        binary['Y'] = 8;  
        binary['K'] = 9;  
        binary['V'] = 10; 
        binary['H'] = 11; 
        binary['D'] = 12; 
        binary['B'] = 13; 
        binary['N'] = 14; 
 
        File f = new File(PATH+"/genomes.db");
        if(f.exists())  
        {
            BufferedReader in;
            try{
                in = new BufferedReader(new FileReader(PATH+"/genomes.info")); 
                num_bytes=Long.valueOf(in.readLine().split(":")[1]);
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
                    sequence_titles[g]=new String[num_sequences[g]+1];
                    sequence_length[g]=new long[num_sequences[g]+1];
                    sequence_start[g]=new long[num_sequences[g]+1];
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
            num_bytes=0;
            previous_num_genomes=0;
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
            BufferedWriter out = new BufferedWriter(new FileWriter(PATH+"/genomes.info"));
            out.write("number_of_bytes:"+num_bytes+"\n");
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
        char carry;
        boolean havecarry;
        long size=0,byte_number=num_bytes;
        int j,g,s,len,k,max_byte=2100000000;
        BufferedReader in;
        System.out.println("Reading "+(num_genomes-previous_num_genomes)+" genomes...");
        for(g=previous_num_genomes+1;g<=num_genomes;++g) 
        {
            try{
                in = new BufferedReader(new FileReader(genome_names[g]));
                s=0;
                while(in.ready())
                {
                    line=in.readLine();
                    if (line==null) 
                        continue;
                    if(line.charAt(0)=='>')
                    {
                        ++s;
                        sequence_titles[g][s]=line.substring(1);
                        if(size%2==1)
                            ++size;
                        sequence_start[g][s]=num_bytes+size/2;
                    }
                    else
                    {
                        len=line.length();
                        sequence_length[g][s]+=len;
                        genome_length[g]+=len;
                        size+=len;
                    }
                }
                if(size%2==1)
                    ++size;                
                in.close();
            }
            catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            }
        } 
        num_bytes+=size/2;
        genomes_file=new RandomAccessFile(PATH+"/genomes.db","rw");
        parts_num=(int)(num_bytes%max_byte==0?num_bytes/max_byte:num_bytes/max_byte+1);
        parts_size=new long[parts_num];
        genomes_buff= new MappedByteBuffer[parts_num];
        for(k=0;k<parts_num;++k)
        {
            parts_size[k]=(int)(k==parts_num-1?num_bytes%max_byte:max_byte);
            genomes_buff[k]=genomes_file.getChannel().map(FileChannel.MapMode.READ_WRITE, k*parts_size[0], parts_size[k]);
        } 
        genomes_buff[(int)(byte_number/parts_size[0])].position((int)(byte_number%parts_size[0]));
        for(g=previous_num_genomes+1;g<=num_genomes;++g) 
        {
            try{
                in = new BufferedReader(new FileReader(genome_names[g]));
                carry=' ';
                havecarry=false;
                s=0;
                while(in.ready())
                {
                    line=in.readLine();
                    if (line==null) 
                        continue;
                    if(line.charAt(0)!='>' && havecarry)
                        line=carry+line;
                    if(line.charAt(0)=='>')
                    {
                        if(havecarry)
                        {
                            genomes_buff[(int)(byte_number/parts_size[0])].put((byte)(binary[carry]<<4 ));
                            ++byte_number;                            
                        }
                        havecarry=false;
                        ++s;
                        System.out.print("Sequence "+s+"/"+num_sequences[g]+" of genome "+g+"         \r");
                    }
                    else
                    {
                        line=line.toUpperCase();
                        len=line.length();
                        havecarry=(len%2==1);
                        if(havecarry)
                        {
                            carry=line.charAt(len-1);
                            --len;
                        }
                        for(j=0;j<len;j+=2,++byte_number)
                            genomes_buff[(int)(byte_number/parts_size[0])].put((byte)((binary[line.charAt(j)]<<4) | binary[line.charAt(j+1)]));
                    }
                }
                if(havecarry)
                {
                    genomes_buff[(int)(byte_number/parts_size[0])].put((byte)(binary[carry]<<4 ));
                    ++byte_number;                            
                }       
                havecarry=false;
                in.close();
            }
            catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            }
        } 
    }  
    ///////////////////////////////////////////////////////////////
    public int get(int g, int s, int p)
    {
        byte b;
        long position=sequence_start[g][s]+p/2;
        b=genomes_buff[(int)(position/parts_size[0])].get((int)(position%parts_size[0]));
        if(p%2==0)
            return (b >> 4) & 0x0f;
        else
            return (b & 0x0f);
    }
    ///////////////////////////////////////////////////////////////
    public String get_sequence(int[] a,int l)
    {
        StringBuilder seq=new StringBuilder();
        for(int i=0;i<l;++i)
            seq.append(sym[get(a[0],a[1],a[2]+i)]);
        return seq.toString();
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
