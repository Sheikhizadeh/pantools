/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pantools;
import java.io.InputStreamReader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.io.UnsupportedEncodingException;
import java.util.Random;
import java.lang.management.MemoryUsage;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryPoolMXBean; 

import org.neo4j.graphdb.DynamicLabel;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.io.fs.FileUtils;
import org.neo4j.graphdb.Direction;

/**
 *
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen University, Netherlands 
 */
public class Pantools {

    /*
    All possible relationship types between two nodes in the graph.
    */
    private static enum RelTypes implements RelationshipType
    {
        AA0,AA1,AA2,AA3,
        AC0,AC1,AC2,AC3,
        AG0,AG1,AG2,AG3,
        AT0,AT1,AT2,AT3,        
        AN0,AN1,AN2,AN3,

        CA0,CA1,CA2,CA3,
        CC0,CC1,CC2,CC3,
        CG0,CG1,CG2,CG3,
        CT0,CT1,CT2,CT3,        
        CN0,CN1,CN2,CN3,
        
        GA0,GA1,GA2,GA3,
        GC0,GC1,GC2,GC3,
        GG0,GG1,GG2,GG3,
        GT0,GT1,GT2,GT3,        
        GN0,GN1,GN2,GN3,
        
        TA0,TA1,TA2,TA3,
        TC0,TC1,TC2,TC3,
        TG0,TG1,TG2,TG3,
        TT0,TT1,TT2,TT3,        
        TN0,TN1,TN2,TN3,
        
        NA0,NA1,NA2,NA3,
        NC0,NC1,NC2,NC3,
        NG0,NG1,NG2,NG3,
        NT0,NT1,NT2,NT3,  
        NN0,NN1,NN2,NN3,
        
        RA0,RA1,RA2,RA3,
        RC0,RC1,RC2,RC3,
        RG0,RG1,RG2,RG3,
        RT0,RT1,RT2,RT3, 
        RN0,RN1,RN2,RN3,
        
        PA0,PA1,PA2,PA3,
        PC0,PC1,PC2,PC3,
        PG0,PG1,PG2,PG3,
        PT0,PT1,PT2,PT3, 
        PN0,PN1,PN2,PN3,
        
        begin,end,
        has
    } 
    private class location
    {
        Node node;
        int format;
        int position;
        public location(Node n, int s, int p)
        {
            node=n;
            format=s;
            position=p;
        }
        
    }
    /*
    Some shared attributes
    */    
    private static String help_comment="Requirements:\n" +
    "\n" +
    "1- KMC: add KMC/bin to the path.\n" +
    "2- neo4j-community-2.3.1-unix\n" +
    "3- Java Virtual Machine jdk1.7 : add the corresponding path to java executable\n" +
    "\n" +
    "To run the program:\n" +
    "java  [-server] [-XX:+UseConcMarkSweepGC]  [-Xmx??g] -jar ./pantools/dist/pantools.jar command arguments\n" +
    "\n" +
    "\n" +
    "List of commands:\n" +
    "\n" +
    "build:    To build a pan-genome out of a set of genomes\n" +
    "          java -jar pantools.jar build K database_path fasta_names_file\n" +
    "\n" +
    "annotate: To add annotations to a pan-genome\n" +
    "          java -jar pantools.jar annotate database_path gff_names_file fasta_names_file\n" +
    "\n" +
    "add:      To add new genomes to an available pan-genome\n" +
    "          java -jar pantools.jar add database_path fasta_names_file\n" +
    "\n" +
    "retrieve_genes:    To extract sequence of genes according to the gene IDs\n" +
    "          java -jar pantools.jar retrieve_genes database_path gene_ids_file\n" +
    "\n" +
    "retrieve_regions:  To extract region sequence\n" +
    "          java -jar pantools.jar retrieve_regions database_path regions_file\n" +
    "        \n" +
    "group:    To group some genes by adding group nodes pointing to them\n" +
    "          java -jar pantools.jar group database_path group_file_file\n" +
    "\n" +
    "compare:     To compare two pan-genomes\n" +
    "          java -jar pantools.jar compare database_1_path database_2_path\n" +
    "\n" +
    "List of arguments:\n" +
    "\n" +
    "K :  Size of k\n" +
    "\n" +
    "fasta_names_file: A text file containing paths to FASTA files; each in one line\n" +
    "\n" +
    "gff_names_file  : A text file containing paths to GFF files corresponding to the genomes in fasta_names_file; each in one line   \n" +
    "\n" +
    "gene_ids_file   : A text file containing genes IDs; each in one line\n" +
    "\n" +
    "regions_file    : A text file containing genome_number, sequence_number, begin, end of a region in each line seperated by one space \n" +
    "\n" +
    "group_file_file : A text file with each line starting with a group name followed by a colon, followed by genes IDs seperated by one space";

    private static String PATH;
    private static String DATABASE="/database.db";
    private static String INDEX="/index.txt";
    private static String KMERS="/kmers.txt";
    public static GraphDatabaseService graphDb;
    public static Index index;
    public final int trsc_limit=1000;    //   The number of transactions to be committed in batch
    public final int anchor_distance=100; // Distance between two anchor nodes
    
    /*
    There are following types of nodes:
    - pangenome   
    - genome
    - sequence
    - node
    - degenerate
    - gene
    - group
    */    
    private static Label pangenome_label= DynamicLabel.label( "pangenome" );
    private static Label genome_label= DynamicLabel.label( "genome" );
    private static Label sequence_label= DynamicLabel.label( "sequence" );
    private static Label node_label = DynamicLabel.label( "node" );
    private static Label degenerate_label = DynamicLabel.label( "degenerate" );
    private static Label gene_label = DynamicLabel.label( "gene" );
    private static Label group_lable = DynamicLabel.label( "group" );


    private static long startTime;
    private static long phaseTime;
    private static int K;
    private static int dim;
    private static int seq_nodes;
    private static int gene_nodes;
    private static int num_edges;
    private static int num_bases;
    public static byte[] binary;
    public static byte[] complement;

    public static byte [][][] Genomes; // Each genome could contain several sequences in the FASTA file
    private static String[][] sequence_names; // Name of sequences for each genome
    private static String[] genome_paths;     // Paths to the genome FASTA files
    public static int sequence_length[][];    // Length of sequences for each genome
    public static int num_sequences[];        // Number of sequences in each genome
    public static char sym[]={ 'A', 'C', 'G' , 'T', 'N'};
    
    public static int num_genomes;
    public static int old_num_genomes;
    private pointer f_kmer,r_kmer,kmer;
    private long seq_len;
    private int genome,sequence,position;
    private static long curr_index;
    private static byte curr_side;
    private Node curr_node;
    private Node db_node;
    private Node new_node;
    private Node split_node;
    private Node degenerate_node;
    private pointer pointer;
    
    private static StringBuilder exe_output;
    private int f_bin,r_bin;
    private boolean canonical,annotate_canonical;
    private boolean finish=false;
    private StringBuilder occ=new StringBuilder();
    private StringBuilder s=new StringBuilder();    
    /**
     * @param args the command line arguments
     * build : To build a pan-genome out of a set of genomes
     * add   : To add new genomes to an available pan-genome
     * genes : To extract gene sequence according to their IDs
     * group : To add group nodes pointing to genes 
     * comp  : To compare two pan-genomes
     * @throws java.io.IOException
     */
    /*
    The main function
    */
    public static void main(String[] args) throws IOException {
        if(args.length < 2 || args[1].equals("--help") || args[1].equals("-h")){
            System.out.println(help_comment);
            System.exit(1);
        } 
        int i;
        Pantools program=new Pantools();
        binary=new byte[124];
        for(i=0;i<124;++i)
            binary[i]=4;
        binary[65] = binary[97]= 0; // A , a
        binary[67] = binary[99]= 1; // C , c
        binary[71] = binary[103]= 2;// G , g
        binary[84] = binary[116]= 3;// T , t
        complement=new byte[124];
        complement[65]=84; // A T
        complement[84]=65; // T A
        complement[67]=71; // C G
        complement[71]=67; // G C
        complement[82]=89; // R Y
        complement[89]=82; // Y R
        complement[83]=83; // S S
        complement[87]=87; // W W
        complement[75]=77; // K M
        complement[77]=75; // M K
        complement[66]=86; // B V
        complement[86]=66; // V B
        complement[68]=72; // D H
        complement[72]=68; // H D
        complement[78]=78; // N N

        System.out.println("------------------------------- PanGea -------------------------------");        
        switch (args[0]) {
            case "build":
                K=Integer.parseInt(args[1]);
                PATH=args[2];
                program.build(args[3]);
                break;
            case "add":
                PATH=args[1];
                program.add(args[2]);
                break;
            case "annotate":
                PATH=args[1];
                program.annotate(args[2],args[3]);
                break;
            case "group":
                PATH=args[1];
                program.group(args[2]);
                break;
            case "compare":
                if(program.compare_pan(args[1],args[2]))
                    System.out.println("Databases are equal.");
                else
                    System.out.println("Databases are different");
                break;
            case "retrieve_genes":
                PATH=args[1];
                program.genes(args[2]);
                break;
            case "retrieve_regions":
                PATH=args[1];
                program.regions(args[2]);
                break;
            default:
                System.out.println(help_comment);
                System.exit(1);            
        }
        System.out.println("-------------------------------- End --------------------------------");        
    }
    /*
    This function builds a pan-genome database out of a set of input genomes.
    fasta_names is a text file containing paths to FASTA files; each in a new line
    */  
    private void build(String fasta_names) throws IOException
    {
        int i;
        //long[] genome_ids;
        File theDir = new File(PATH);
        if (theDir.exists()) 
            FileUtils.deleteRecursively( new File( PATH+DATABASE ) ); // deletes old files
        else
            {
                try{
                    theDir.mkdir();
                } 
                catch(SecurityException se){
                   System.out.println("Failed to create "+PATH);  
                   System.exit(1);
                }  
            }
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+DATABASE)
                .setConfig( "keep_logical_logs","100M size").newGraphDatabase();
        registerShutdownHook( graphDb );            
        startTime = System.currentTimeMillis();
        dim=(int)Math.ceil(K/16.0);
        seq_nodes=0;
        num_edges=0;
        try(Transaction tx = graphDb.beginTx()){
            db_node=graphDb.createNode(pangenome_label);
            db_node.setProperty("k_mer_size", K);
        tx.success();}
        old_num_genomes=0;
        num_genomes=Integer.parseInt(executeCommand("wc -l "+fasta_names).trim().split("\\s")[0]);
        genome_paths=new String[num_genomes+1];
        Genomes=new byte[num_genomes+1][][];
        sequence_names=new String[num_genomes+1][];
        sequence_length=new int[num_genomes+1][];
        num_sequences=new int[num_genomes+1];
        try{
            BufferedReader in = new BufferedReader(new FileReader(fasta_names));
            for(i=1;i<=num_genomes;++i)
            {
                genome_paths[i]=in.readLine().trim();
                num_sequences[i]=Integer.parseInt(executeCommand("grep -c '>' "+genome_paths[i]).trim());
                sequence_names[i]=new String[num_sequences[i]+1];
                Genomes[i]=new byte[num_sequences[i]+1][];
                sequence_length[i]=new int[num_sequences[i]+1];
            }
        }
        catch(IOException ioe){
           System.out.println("Failed to read file names: "+ioe);  
           System.exit(1);
        } 
        load_genomes();
        index=new Index(num_genomes>1?"@"+fasta_names.trim():genome_paths[1],PATH,KMERS,K);
        construct_pangenome();
        finalise();
        try(Transaction tx = graphDb.beginTx()){
            db_node.setProperty("num_k_mers", index.totallength);
            db_node.setProperty("num_nodes",seq_nodes);
            db_node.setProperty("num_edges",num_edges);
            db_node.setProperty("num_genomes",num_genomes);
            db_node.setProperty("num_genes",0);
            db_node.setProperty("num_bases",num_bases);
            tx.success();}
        try(Transaction tx = graphDb.beginTx()){
            graphDb.schema().indexFor( gene_label ).on( "id" ).create();
            tx.success();}
        graphDb.shutdown();  
        System.out.println("Total time : "+(System.currentTimeMillis()-startTime)/1000+"."+(System.currentTimeMillis()-startTime)%1000+" seconds");
        print_peak_memory();
        executeCommand("rm "+PATH+"/database.db/neostore.transaction.db*");
        System.out.print(Pantools.executeCommand("du --apparent-size -shc -m "+PATH+DATABASE));
        System.out.print(Pantools.executeCommand("du --apparent-size -shc -m "+PATH+INDEX));
    }
    /*
    This function adds new genomes to an available pan-genome database.
    fasta_names is a text file containing paths to FASTA files; each in a new line
    */      
    private void add(String fasta_names) throws IOException
    {
        int i,j;
        Index old_index;
        Index new_index;
        long old_kmer_num=0,c_index,p_index,l;
        ResourceIterator<Node> nodes;
        if (new File(PATH+DATABASE).exists()) 
        {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+DATABASE)
                    .setConfig( "keep_logical_logs","100M size").newGraphDatabase();      
            registerShutdownHook( graphDb );
            startTime = System.currentTimeMillis();
            try(Transaction tx = graphDb.beginTx()){
                db_node=graphDb.findNodes(pangenome_label).next();
                if(db_node==null)
                {
                    System.out.println("Can not locate database node!");
                    System.exit(1);
                }
                K = (int)db_node.getProperty("k_mer_size");
                dim=(int)Math.ceil(K/16.0);
                seq_nodes=(int)db_node.getProperty("num_nodes");
                num_edges=(int)db_node.getProperty("num_edges");
                old_kmer_num=(long)db_node.getProperty("num_k_mers");
                old_num_genomes=(int)db_node.getProperty("num_genomes");
                gene_nodes=(int)db_node.getProperty("num_genes");
                num_genomes=old_num_genomes+Integer.parseInt(executeCommand("wc -l "+fasta_names).trim().split("\\s")[0]);
                genome_paths=new String[num_genomes+1];
                Genomes=new byte[num_genomes+1][][];
                sequence_names=new String[num_genomes+1][];
                sequence_length=new int[num_genomes+1][];
                num_sequences=new int[num_genomes+1];
                String[] old_genome_paths=new String[old_num_genomes+1];
                nodes = graphDb.findNodes( genome_label );
                while(nodes.hasNext())
                {
                    curr_node=nodes.next();
                    old_genome_paths[(int)curr_node.getProperty("genome_number")]=(String)curr_node.getProperty("genome_path");
                }                
                try{
                    BufferedReader in = new BufferedReader(new FileReader(fasta_names));
                    for(i=1;i<=num_genomes;++i)
                    {
                        genome_paths[i]= i<=old_num_genomes?old_genome_paths[i]:in.readLine();
                        num_sequences[i]=Integer.parseInt(executeCommand("grep -c '>' "+genome_paths[i]).trim());
                        sequence_names[i]=new String[num_sequences[i]+1];
                        Genomes[i]=new byte[num_sequences[i]+1][];
                        sequence_length[i]=new int[num_sequences[i]+1];
                    }
                    in.close();
                }
                catch(IOException ioe){
                   System.out.println("Failed to read file names!");  
                   System.exit(1);
                } 
                load_genomes();
                old_index=new Index(PATH+INDEX,old_kmer_num,K);
                new_index=new Index("@"+fasta_names.trim(),PATH,KMERS,K);
                index=new Index(old_index,new_index, graphDb, K);
                System.out.println("Updating pointers...");
                    nodes = graphDb.findNodes( node_label );
                    for(j=1;nodes.hasNext();++j)
                    {
                        curr_node=nodes.next();
                        l=(long)curr_node.getProperty("first_kmer");
                        p_index=index.find(old_index.pointers[(int)(l/2100000000)][(int)(l%2100000000)]);
                        curr_node.setProperty("first_kmer", p_index);
                        for(l=old_index.pointers[(int)(l/2100000000)][(int)(l%2100000000)].next_index;l!=-1L;l=old_index.pointers[(int)(l/2100000000)][(int)(l%2100000000)].next_index)
                        {
                            c_index=index.find(old_index.pointers[(int)(l/2100000000)][(int)(l%2100000000)]);
                            index.pointers[(int)(p_index/2100000000)][(int)(p_index%2100000000)].next_index=c_index;
                            p_index=c_index;
                        }
                        index.pointers[(int)(p_index/2100000000)][(int)(p_index%2100000000)].next_index=-1L;
                        if(j%(seq_nodes/100+1)==0) 
                            System.out.print(j*100/seq_nodes+1+"%\r");
                    }
                old_index=null;
                new_index=null;
                Runtime.getRuntime().gc();
            tx.success();} 
            construct_pangenome();
            finalise();
            try(Transaction tx = graphDb.beginTx())
            {
                db_node.setProperty("k_mer_size", K);
                db_node.setProperty("num_k_mers", index.totallength);
                db_node.setProperty("num_nodes",seq_nodes);
                db_node.setProperty("num_edges",num_edges);
                db_node.setProperty("num_genomes",num_genomes);
                db_node.setProperty("num_genes",gene_nodes);
                db_node.setProperty("num_bases",num_bases);
                tx.success();
            }
            graphDb.shutdown(); 
            System.out.println("Total time : "+(System.currentTimeMillis()-startTime)/1000+"."+(System.currentTimeMillis()-startTime)%1000+" seconds");
            print_peak_memory();
            executeCommand("rm "+PATH+"/database.db/neostore.transaction.db*");
            System.out.print(Pantools.executeCommand("du --apparent-size -shc -m "+PATH+DATABASE));
            System.out.print(Pantools.executeCommand("du --apparent-size -shc -m "+PATH+INDEX));
        }
        else
        {
            System.out.println("No database found in "+PATH); 
            System.exit(1);
        }  
    }
    /*
    This function adds gene_nodes to the graph. 
    GFF file should be in the same path as the FASTA file with the same name but .gff extention. 
    */      
    private void annotate(String gff_names, String fasta_names) throws IOException
    {
        int i;
        long kmer_num=0;
        int begin,end,g,s,line_len;
        String seq,id,name;
        boolean start_canonical,stop_canonical;
        Node gene_node,genome_node;
        pointer start_point,stop_point;
        Relationship rel;
        String[] fields;
        String strand,fasta_name,gff_name;
        if(new File(PATH+DATABASE).exists()) 
        {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+DATABASE)
                    .setConfig( "keep_logical_logs","100M size").newGraphDatabase();        
            registerShutdownHook( graphDb );
            try(Transaction tx = graphDb.beginTx()){
                db_node=graphDb.findNodes(pangenome_label).next();
                if(db_node==null)
                {
                    System.out.println("Can not locate database node!");
                    System.exit(1);
                }
                K = (int)db_node.getProperty("k_mer_size");
                dim=(int)Math.ceil(K/16.0);
                seq_nodes=(int)db_node.getProperty("num_nodes");
                num_edges=(int)db_node.getProperty("num_edges");
                num_genomes=(int)db_node.getProperty("num_genomes");
                kmer_num=(long)db_node.getProperty("num_k_mers");
                gene_nodes=(int)db_node.getProperty("num_genes");
                genome_paths=new String[num_genomes+1];
                Genomes=new byte[num_genomes+1][][];
                sequence_names=new String[num_genomes+1][];
                sequence_length=new int[num_genomes+1][];
                num_sequences=new int[num_genomes+1];
                try{
                    BufferedReader in = new BufferedReader(new FileReader(fasta_names));
                    for(i=1;i<=num_genomes;++i)
                    {
                        genome_paths[i]=in.readLine();
                        num_sequences[i]=Integer.parseInt(executeCommand("grep -c '>' "+genome_paths[i]).trim());
                        sequence_names[i]=new String[num_sequences[i]+1];
                        Genomes[i]=new byte[num_sequences[i]+1][];
                        sequence_length[i]=new int[num_sequences[i]+1];
                    }
                    in.close();
                }
                catch(IOException ioe){
                   System.out.println("Failed to read file names!");  
                   System.exit(1);
                }
                load_genomes();
                index=new Index(PATH+INDEX,kmer_num,K);
                startTime = System.currentTimeMillis();
                try(BufferedReader gff_files = new BufferedReader(new FileReader(gff_names))){
                BufferedReader fasta_files = new BufferedReader(new FileReader(fasta_names));
                 while(gff_files.ready())
                    {
                        gff_name=gff_files.readLine().trim();
                        fasta_name=fasta_files.readLine().trim();
                        genome_node=graphDb.findNode(genome_label, "genome_path", fasta_name);
                        if(genome_node.hasProperty("annotated"))
                        {
                            System.out.println("genome "+fasta_name+" already annotated.");
                            continue;
                        }
                        g=(int)genome_node.getProperty("genome_number");
                        BufferedReader in = new BufferedReader(new FileReader(gff_name));
                        System.out.println("Annotating genome "+g);
                        while(in.ready())
                        {
                            fields=in.readLine().split("\\t");
                            line_len=fields.length;
                            for(i=0;i<line_len && !fields[i].equals("gene");++i)
                                {}
                            if(i<line_len && fields[i].equals("gene"))
                            {
                                seq=fields[0];
                                s=find_sequence(seq, g);
                                if(s>0) // if sequence found
                                {
                                    begin=Integer.parseInt(fields[i+1])-1;
                                    end=Integer.parseInt(fields[i+2])-1;
                                    strand=fields[i+4];
                                    id=fields[i+6].contains("ID=")?fields[i+6].split("ID=")[1].split(";")[0]:"";
                                    name=fields[i+6].contains("Name=")?fields[i+6].split("Name=")[1].split(";")[0]:"";
                                    start_point=find_node(g,s,begin);
                                    start_canonical=annotate_canonical;
                                    stop_point=find_node(g,s,end-K+1);
                                    stop_canonical=annotate_canonical;
                                    if(start_point!=null && stop_point!=null && end-begin+1>=K) // genes should not be shorter than K
                                    {
                                        ++gene_nodes;
                                        gene_node=graphDb.createNode(gene_label);
                                        gene_node.setProperty("genome_name", genome_paths[g]);
                                        gene_node.setProperty("genome_number", g);
                                        gene_node.setProperty("sequence_name", seq);
                                        gene_node.setProperty("sequence_number", s);
                                        gene_node.setProperty("begin", begin+1);
                                        gene_node.setProperty("end", end+1);
                                        gene_node.setProperty("length", end-begin+1);
                                        gene_node.setProperty("strand", strand);
                                        gene_node.setProperty("id","genome"+g+"|"+seq+"|"+id+"|"+name);
                                        rel=gene_node.createRelationshipTo(graphDb.getNodeById(start_point.node_id), RelTypes.begin);
                                        rel.setProperty("side", start_canonical ?start_point.format:1-start_point.format);
                                        rel.setProperty("pos", start_point.position);
                                        rel=gene_node.createRelationshipTo(graphDb.getNodeById( stop_point.node_id), RelTypes.end);
                                        rel.setProperty("side", stop_canonical ?stop_point.format:1-stop_point.format);
                                        rel.setProperty("pos", stop_point.position+K-1);
                                    }
                                    else
                                        System.out.println(id+" not properly called!"); 
                                }
                                else 
                                    System.out.println(seq +" missed in "+genome_paths[g]);
                            }
                            genome_node.setProperty("annotated","Yes");
                        }
                        in.close();
                    }
                    gff_files.close();
                }
                catch(IOException ioe)
                { 
                }
                db_node.setProperty("num_genes",gene_nodes);
                tx.success();
            }
            graphDb.shutdown(); 
            System.out.println("Number of gene nodes :\t"+gene_nodes);
            System.out.println("Total time : "+(System.currentTimeMillis()-startTime)/1000+"."+(System.currentTimeMillis()-startTime)%1000+" seconds");
            print_peak_memory();
            executeCommand("rm "+PATH+"/database.db/neostore.transaction.db*");
        }
        else
        {
            System.out.println("No database found in "+PATH); 
            System.exit(1);
        }  
    }
    /*
    This function adds group_nodes to the graph. These nodes link ortholog/homolog/etc genes. 
    group_file is text file with each line starting with a group name followed by a colon, 
    followed by genes IDs seperated by one space.
    */     
    private void group(String group_file)
    {
        if (new File(PATH+DATABASE).exists()) 
        {
            Node group_node,gene;
            String line;
            String[] fields;
            ResourceIterator<Node> genes;
            int i,j,g;
            boolean has;
            System.out.println("Adding groups of "+group_file+" to "+PATH+DATABASE+"...");
            try(BufferedReader in = new BufferedReader(new FileReader(group_file))){
                graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+DATABASE)
                        .setConfig( "keep_logical_logs","100M size").newGraphDatabase();  
                registerShutdownHook( graphDb );
                startTime = System.currentTimeMillis();
                try(Transaction tx = graphDb.beginTx())
                {
                    for(i=0;in.ready();++i)
                    {
                        line=in.readLine(); 
                        group_node=graphDb.createNode(group_lable);
                        fields=line.split(":");
                        group_node.setProperty("group_name", fields[0]);
                        fields=fields[1].trim().split(" ");
                        for(g=0,j=0;j<fields.length;++j)
                        {
                            genes=graphDb.findNodes(gene_label,"id",fields[j]);//.split("\\|")[1]);
                            for(;genes.hasNext();++g)
                            {
                                gene=genes.next();
                                has=false;
                                for(Relationship r: group_node.getRelationships(RelTypes.has, Direction.OUTGOING))
                                    if(r.getEndNode().equals(gene))
                                    {
                                        has=true;
                                        break;
                                    }
                                if(!has)
                                    group_node.createRelationshipTo(gene,RelTypes.has);
                            }
                            genes.close();
                        }
                        group_node.setProperty("num_genes", g);
                    }
                    tx.success();
                    System.out.println("\n"+i+" groups added to the pan-genome");
                }
                in.close();
            }
            catch(IOException ioe){
               System.out.println("Failed to open "+group_file);  
               System.exit(1);
            }
            System.out.println("Total time : "+(System.currentTimeMillis()-startTime)/1000+"."+(System.currentTimeMillis()-startTime)%1000+" seconds");
            print_peak_memory();
            graphDb.shutdown();
        }
        else
        {
            System.out.println("No database found in "+PATH+DATABASE); 
            System.exit(1);
        }
    }
    /*
    This function retrieves genes' sequence from the graph
    Input:
        file: A text file containing ID of genes whose sequence should be retrieved, One ID per line
    Output:
        A FASTA file of gene sequences.
    */ 
    private void genes(String file) throws IOException
    {
        if (new File(PATH+DATABASE).exists()) 
        {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+DATABASE)
                    .setConfig( "keep_logical_logs","100M size").newGraphDatabase(); 
            registerShutdownHook( graphDb );
            startTime = System.currentTimeMillis();
            try(Transaction tx = graphDb.beginTx()){
                K=(int)graphDb.findNodes(pangenome_label).next().getProperty("k_mer_size");
            tx.success();}
            Relationship rstart,rstop;
            Node start,stop,gene;
            String prp,gene_id,strand,id;
            int begin,end;
            String[] occ;
            StringBuilder gene_seq;
            try(Transaction tx = graphDb.beginTx()){
            try (BufferedReader in = new BufferedReader(new FileReader(file))) {
                BufferedWriter out = new BufferedWriter(new FileWriter(file+".fasta"));
                for(gene_id=in.readLine().trim();gene_id!=null;gene_id=in.readLine())
                {
                    gene=graphDb.findNode(gene_label,"id",gene_id);
                    if(gene != null)
                    {
                        System.out.println("Retrieving "+ gene.getProperty("id").toString());
                        rstart=gene.getSingleRelationship(RelTypes.begin, Direction.OUTGOING);
                        rstop=gene.getSingleRelationship(RelTypes.end, Direction.OUTGOING);
                        start=rstart.getEndNode();
                        stop=rstop.getEndNode();
                        prp=gene.getProperty("genome_number").toString()+"_"+gene.getProperty("sequence_number").toString();
                        begin=(int)gene.getProperty("begin");
                        end=(int)gene.getProperty("end");
                        id=gene.getProperty("id").toString();
                        strand=gene.getProperty("strand").toString();
                        out.write(">"+id+" begin:"+begin+" end:"+end+" strand:"+strand+" length:");
                        gene_seq=sequence(start,stop,(int)rstart.getProperty("side"),(int)rstart.getProperty("pos"),prp, begin-1, end-1);
                        if(gene_seq.length()!=end-begin+1)
                        {
                            System.out.println("Failed to assemble "+gene.getProperty("id").toString());
                            System.out.println(gene_seq.length()+" != "+(end-begin+1));
                            //System.out.println(gene_seq.toString());
                        }
                        out.write(gene_seq.length()+"\n");
                        if(strand.equals("+"))
                            write_fasta(out,gene_seq.toString(),70);
                        else
                            write_fasta(out,reverse_complement(gene_seq.toString()),70);
                        gene_seq.setLength(0);
                    }
                    else
                    {}//System.out.println(gene_id+" not annotated!");
                }//for gene_id
                in.close();
                out.close();
            }
            catch(IOException ioe){
               System.out.println("Failed to read file names!");  
               System.exit(1);
            }
            tx.success();}
            System.out.println("Total time : "+(System.currentTimeMillis()-startTime)/1000+"."+(System.currentTimeMillis()-startTime)%1000+" seconds");
            print_peak_memory();
            graphDb.shutdown(); 
        }
        else
        {
            System.out.println("No database found in "+PATH); 
            System.exit(1);
        }  
    }
    /*
    This function retrieve a genomic region from the graph
    Input:
        from: The start index of the region
        to:   The stop index of the region
        g:    The genome number from which the region should be extracted 
        s:    The sequence number from which the region should be extracted 
    Output:
        sequence of the regions in FASTA format
    */ 
    private void regions(String file) throws IOException
    {
        if (new File(PATH+DATABASE).exists()) 
        {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+DATABASE)
                    .setConfig( "keep_logical_logs","100M size").newGraphDatabase(); 
            registerShutdownHook( graphDb );
            try(Transaction tx = graphDb.beginTx()){
                K=(int)graphDb.findNodes(pangenome_label).next().getProperty("k_mer_size");
            tx.success();}
            String[] fields;
            String line;
            StringBuilder seq;
            long[] anchor_nodes;
            int[] anchor_positions;
            location start,stop;
            Node seq_node;
            int c,g,s,i,j,low,high,mid,from, to, num_regions;
            String anchor_sides;
            num_regions=Integer.parseInt(Pantools.executeCommand("wc -l "+file).trim().split("\\s")[0]);
            startTime = System.currentTimeMillis();
            try(Transaction tx = graphDb.beginTx()){
            try (BufferedReader in = new BufferedReader(new FileReader(file))) {
                BufferedWriter out = new BufferedWriter(new FileWriter(file+".fasta"));
                for(c=0,line=in.readLine().trim();line!=null;line=in.readLine())
                {
                   fields=line.trim().split("\\s");
                    g= Integer.parseInt(fields[0]);
                    s= Integer.parseInt(fields[1]);
                    from= Integer.parseInt(fields[2]);
                    to= Integer.parseInt(fields[3]);
                    seq_node=graphDb.findNode(sequence_label, "sequence_number", g+"_"+s);
                    anchor_nodes=(long [])seq_node.getProperty("anchor_nodes");
                    anchor_positions=(int [])seq_node.getProperty("anchor_positions");
                    anchor_sides=(String)seq_node.getProperty("anchor_sides");
                    for(low=0,high=anchor_sides.length()-1,mid=(low+high)/2;low<high;mid=(low+high)/2)
                    {
                        if(from<anchor_positions[mid])
                            high=mid-1;
                        else if(from>anchor_positions[mid])
                            low=mid+1;
                        else
                            break;
                    }
                    if(from>=anchor_positions[mid])
                        i=mid;
                    else
                        i=mid-1;
                    for(low=0,high=anchor_sides.length()-1,mid=(low+high)/2;low<high;mid=(low+high)/2)
                    {
                        if(to<anchor_positions[mid])
                            high=mid-1;
                        else if(to>anchor_positions[mid])
                            low=mid+1;
                        else
                            break;
                    }
                    if(to>=anchor_positions[mid])
                        j=mid;
                    else
                        j=mid-1;
                    start=locate(graphDb.getNodeById(anchor_nodes[i]),anchor_positions[i],g+"_"+s,from);
                    stop =locate(graphDb.getNodeById(anchor_nodes[j]),anchor_positions[j],g+"_"+s,to);
                    if(start.format==0)
                        seq=sequence(start.node,stop.node,start.format, from-start.position,g+"_"+s, from, to);
                    else
                        seq=sequence(start.node,stop.node,start.format, (int)start.node.getProperty("length")-(from-start.position)-K,g+"_"+s,from, to);
                    if(seq.length()!=to-from+1)
                    {
                        System.out.println("Failed to assemble "+line);
                        System.out.println(seq.length()+" != "+(to-from+1));
                    }                    
                    out.write(">genome:"+g+" sequence:"+s+" from:"+from+" to:"+to+" length:"+seq.length()+"\n");
                    write_fasta(out,seq.toString(),70);
                    seq.setLength(0);
                    ++c;
                    if(c%(num_regions/100+1)==0) 
                        System.out.print((long)c*100/num_regions+1+"%\r");
                }
                in.close();
                out.close();
            }
            catch(IOException ioe){
                System.out.println("Failed to read file names!");  
                System.exit(1);
            }
            tx.success();}
            System.out.println("Total time : "+(System.currentTimeMillis()-startTime)/1000+"."+(System.currentTimeMillis()-startTime)%1000+" seconds");
            print_peak_memory();
            graphDb.shutdown(); 
        }
    }
    /*
    This function returns the sequence starting at (start_node, start_side, start_pos) and stopping at stop_node which belongs to 
    genoeme and sequence specified by "coordinates" and is of length end-begin+1
    */
    private StringBuilder sequence(Node start_node, Node stop_node, int start_side, int start_pos, String coordinates, int begin, int end) throws UnsupportedEncodingException
    {
        boolean found;
        Node node=start_node,neighbor;
        String lable,name;
        int loc,node_len,neighbor_len,s_side,d_side,seq_len=end-begin+1;
        String[] occ;
        StringBuilder seq=new StringBuilder();
        if(start_node.equals(stop_node) && start_side==0 && start_pos+end-begin<(int)start_node.getProperty("length")) // if gene is inside the node not a circle which starts and ends in the node
            append(seq,(byte[])node.getProperty("sequence"),start_pos,start_pos+end-begin+1,true);
        else if(start_node.equals(stop_node) && start_side==1 && start_pos+K-1-end+begin>=0)
            append(seq,(byte[])node.getProperty("sequence"),start_pos+K-1-end+begin,start_pos+K,false);
        else
        {
            node_len=(int)node.getProperty("length");
            if(start_side==0)
            {
                if(start_pos+end-begin<(int)start_node.getProperty("length"))
                    append(seq,(byte[])node.getProperty("sequence"),start_pos,start_pos+end-begin+1,true);
                else
                    append(seq,(byte[])node.getProperty("sequence"),start_pos,0,true);
                loc=begin-start_pos;
            }
            else
            {
                if(start_pos+K-1-end+begin>=0)
                    append(seq,(byte[])node.getProperty("sequence"),start_pos+K-1-end+begin,start_pos+K,false);
                else
                    append(seq,(byte[])node.getProperty("sequence"),0,start_pos+K,false);
                loc=begin-node_len+K+start_pos;
            }
            found=true;
            while(seq.length()<seq_len && found) 
                {
                    found=false;
                    for(Relationship r:node.getRelationships(Direction.OUTGOING))
                    {
                        neighbor=r.getEndNode();
                        name=r.getType().name();
                        d_side=(name.charAt(2)-48)%2;
                        lable=(d_side==0?"F"+coordinates:"R"+coordinates);
                        if(neighbor.hasProperty(lable) )
                        {
                            neighbor_len=(int)neighbor.getProperty("length");
                            if(b_search((int[])neighbor.getProperty(lable), loc+node_len-K+1)>=0)
                            {
                                found=true;
                                loc=loc+node_len-K+1;
                                if(d_side==0)
                                    if(seq.length()+neighbor_len-K+1>seq_len)
                                        append(seq,(byte[])neighbor.getProperty("sequence"),K-1,seq_len-seq.length()+K-1,true);
                                    else
                                        append(seq,(byte[])neighbor.getProperty("sequence"),K-1,0,true);
                                else
                                    if(seq.length()+neighbor_len-K+1>seq_len)
                                        append(seq,(byte[])neighbor.getProperty("sequence"),neighbor_len-K,neighbor_len-K-seq_len+seq.length()-1,false);
                                    else
                                        append(seq,(byte[])neighbor.getProperty("sequence"),neighbor_len-K,0,false);
                                node=neighbor;
                                node_len=(int)node.getProperty("length");
                                break;
                            }
                        }
                    }
                    if(!found)
                    for(Relationship r:node.getRelationships(Direction.INCOMING))
                    {
                        neighbor=r.getStartNode();
                        name=r.getType().name();
                        s_side=(name.charAt(2)-48)/2;
                        lable=(s_side==0?"R"+coordinates:"F"+coordinates);
                        if(neighbor.hasProperty(lable) )
                        {
                            neighbor_len=(int)neighbor.getProperty("length");
                            if(b_search((int[])neighbor.getProperty(lable), loc+node_len-K+1)>=0)
                            {
                                found=true;
                                loc=loc+node_len-K+1;
                                if(s_side==1)
                                    if(seq.length()+neighbor_len-K+1>seq_len)
                                        append(seq,(byte[])neighbor.getProperty("sequence"),K-1,seq_len-seq.length()+K-1,true);
                                    else
                                        append(seq,(byte[])neighbor.getProperty("sequence"),K-1,0,true);
                                else
                                    if(seq.length()+neighbor_len-K+1>seq_len)
                                        append(seq,(byte[])neighbor.getProperty("sequence"),neighbor_len-K,neighbor_len-K-seq_len+seq.length()-1,false);
                                    else
                                        append(seq,(byte[])neighbor.getProperty("sequence"),neighbor_len-K,0,false);
                                node=neighbor;
                                node_len=(int)node.getProperty("length");
                                break;
                            }
                        }
                    }
                }//for loc
        }
        return seq;
    }
    /*
    This function appends characters of "s" in range [from..to] to "seq". forward determines the appending direction. 
    */
    private void append(StringBuilder seq, byte[] s, int from, int to, boolean forward)
    {
        --to;
        if(to==-1)
            to=s.length-1;
        if(forward)
            for(int i=from;i<=to;++i) 
                seq.append((char)s[i]);
        else
            for(int i=to;i>=from;--i) 
                seq.append((char)(complement[s[i]]));            
    }
    /*
    This function returns the "location" (node, side, position) where a region located at "position" in "genoeme" and "sequence"
    specified by "coordinates" occurs. It starts the search from "start" node.
    */
    private location locate(Node start, int loc, String coordinates, int position)
    {
        int node_len,side=0;
        Node node,neighbor;
        String lable;
        boolean found;
        try(Transaction tx = graphDb.beginTx()){
        node=start;
        node_len=(int)node.getProperty("length");
        found=true;
        while(loc+node_len-K+1<position && found) // while have not reached to target
        {
            found=false;
            for(Relationship r:node.getRelationships(Direction.OUTGOING)) // check all the outgoing edges
            {
                neighbor=r.getEndNode();
                side=(r.getType().name().charAt(2)-48)%2;
                lable=(side==0?"F"+coordinates:"R"+coordinates);
                if(neighbor.hasProperty(lable) )
                {
                    if(b_search((int[])neighbor.getProperty(lable), loc+node_len-K+1)>=0)
                    {
                        found=true;
                        loc=loc+node_len-K+1;
                        node=neighbor;
                        node_len=(int)node.getProperty("length");
                        break;
                    }
                }
            }
            if(!found)
            for(Relationship r:node.getRelationships(Direction.INCOMING)) // check all the incoming edges
            {
                neighbor=r.getStartNode();
                side=(r.getType().name().charAt(2)-48)/2;
                lable=(side==0?"R"+coordinates:"F"+coordinates);
                if(neighbor.hasProperty(lable) )
                {
                    if(b_search((int[])neighbor.getProperty(lable), loc+node_len-K+1)>=0)
                    {
                        found=true;
                        loc=loc+node_len-K+1;
                        node=neighbor;
                        node_len=(int)node.getProperty("length");
                        break;
                    }
                }
            }
        }
        tx.success();}
        return new location(node,side,loc);
    }
    /*
    This function finds the node which contains the K-mer occuring at genome g, 
    sequence s at position pos and returns a pointer object to this node.
    */     
    pointer find_node(int g, int s, int pos)
    {
        ResourceIterator<Node> degenerate_nodes;
        long inx;
        if(pos+K>sequence_length[g][s])
            return null;
        boolean ambiguous=false,canonical;
        for(int i=0;i<K && !ambiguous;++i)
            if(binary[Genomes[g][s][pos+i]]==4)
                ambiguous=true;
        if(!ambiguous)
        {
            inx=index.find(make_kmer(g,s,pos));
            return index.pointers[(int)(inx/2100000000)][(int)(inx%2100000000)];
        }
        else // degenerate node found
        {
            Node degenerate_node=null;
            int loc, lc, len;
            for(loc=pos;binary[Genomes[g][s][loc]]!=4;)
                ++loc;
            lc=loc;
            do{
                loc=lc;
                for(lc=loc-1; lc>=loc-K+1;--lc)
                    if(lc<0 || binary[Genomes[g][s][lc]]==4)
                        break;
            }while(lc>=0 && binary[Genomes[g][s][lc]]==4);
            loc=lc+1;
            boolean found=false;
            byte[] seq;
            canonical=true;
            try(Transaction tx = graphDb.beginTx()){
            degenerate_nodes = graphDb.findNodes(degenerate_label,"genome",g);
            while(degenerate_nodes.hasNext() && !found)
            {
                degenerate_node=degenerate_nodes.next();
                len=(int)degenerate_node.getProperty("length");
                seq=(byte[])degenerate_node.getProperty("sequence");
                if(loc+len<=sequence_length[g][s] && are_equal(seq,0,Genomes[g][s],loc,len,true))
                    found=true;
            }
            degenerate_nodes.close();
            tx.success();} 
            if(found)
                return new pointer(dim,degenerate_node.getId(),(byte)0,(pos-loc),0);   
            else
            {
                System.out.println("Gap node not found!");
                return null;
            }
        }
    }

    /*
    This function returns the number of the sequence from genome "g" whose name starts with "prefix".
    */    
    int find_sequence(String prefix, int g)
    {
        boolean found=false;
        int s;
        for(found=false,s=1;!found && s<=num_sequences[g];++s)
            if(sequence_names[g][s].contains(prefix))
                found=true;
        --s;
        if(found)
            return s;
        else 
            return -1;
    }
    
    /*
    This function returns the reverse complement of "s".
    */    
    String reverse_complement(String s)
    {
        StringBuilder rv=new StringBuilder();
        for(int i=s.length();i>=0;--i)
            switch(s.charAt(i))
            {
                case 0: 
                    rv.append(3);
                    break;
                case 1: 
                    rv.append(2);
                    break;
                case 3: 
                    rv.append(1);
                    break;
                case 4: 
                    rv.append(0);
                    break;                    
                case 'R': 
                    rv.append('Y');
                    break;
                case 'Y': 
                    rv.append('R');
                    break;
                case 'K': 
                    rv.append('M');
                    break;
                case 'M': 
                    rv.append('K');
                    break;  
                case 'B': 
                    rv.append('V');
                    break;
                case 'V': 
                    rv.append('B');
                    break;
                case 'D': 
                    rv.append('H');
                    break;
                case 'H': 
                    rv.append('D');
                    break;  
                default:
                    rv.append(s.charAt(i));
            }
        return rv.toString();
    }

    /*
    This function loads genomes'sequences to Genome array. 
    */      
    private void load_genomes()
    {
        String line;
        int i=0,j,g,s,len;
        BufferedReader in;
        for(g=1;g<=num_genomes;++g) 
        {
            System.out.print("Loading genome "+g+"/"+num_genomes+"\r");
            try{
                in = new BufferedReader(new FileReader(genome_paths[g]));
                s=0;
                len=0;
                while(in.ready())
                {
                    line=in.readLine();
                    if (line.equals("")) 
                        continue;
                    if(line.charAt(0)=='>')
                    {
                        ++s;
                        sequence_names[g][s]=line.substring(1);
                        if(s>1)
                            sequence_length[g][s-1]=len;
                        len=0;
                    }
                    else
                        len += line.length();
                }
                in.close();
                sequence_length[g][s]=len;
                in = new BufferedReader(new FileReader(genome_paths[g]));
                s=0;
                while(in.ready())
                {
                    line=in.readLine().toUpperCase();
                    if (line.equals("")) 
                        continue;
                    if(line.charAt(0)=='>')
                    {
                        i=0;
                        ++s;
                        Genomes[g][s]=new byte[sequence_length[g][s]+1];
                    }
                    else
                    {
                        for(j=0;j<line.length();++j,++i)
                            Genomes[g][s][i] = (byte)line.charAt(j);
                    }
                }
                in.close();
            }
            catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            }
        }        
        System.out.println();
    }

    /*
    This function creates an edge between "src" and "des" nodes. 
    Type of the edge is determined by "s_side", "d_side" and transition bases.
    */     
    void connect(Node src, Node des,int s_side, int d_side, int stop)//stop=0,1
    {
        RelationshipType f_rt,r_rt;
        Relationship r;
        int s_gen,s_seq,s_loc,d_gen,d_seq,d_loc,s_base,d_base=0,s_len,d_len;
        int[] s_add,d_add;
        String rt_name;
        boolean has_fwd, has_rev;
        if(src.hasLabel(sequence_label))
            s_base=stop+5;
        else
        {
            s_len= (int)src.getProperty("length");
            s_add=(int[])src.getProperty("address");
            s_gen=s_add[0];
            s_seq=s_add[1];
            s_loc=s_add[2];
            s_base=s_side==0?binary[Genomes[s_gen][s_seq][s_loc+s_len-K]]:(3-binary[Genomes[s_gen][s_seq][s_loc+K-1]]);
            if(s_base==-1)
                s_base=4;
        }
        d_len= (int)des.getProperty("length") ;
        d_add=(int[])des.getProperty("address");
        d_gen=d_add[0];
        d_seq=d_add[1];
        d_loc=d_add[2];
        d_base=d_side==0?binary[Genomes[d_gen][d_seq][d_loc+K-1]]:(3-binary[Genomes[d_gen][d_seq][d_loc+d_len-K]]);
        if(d_base==-1)
            d_base=4;
        f_rt=RelTypes.values()[s_base*20+d_base*4+s_side*2+d_side];
        if(src.hasLabel(degenerate_label) || des.hasLabel(degenerate_label) || src.hasLabel(sequence_label))
        {
                src.createRelationshipTo( des, f_rt);
                ++num_edges;                
        }
        else
        {
            rt_name=f_rt.name();
            r_rt=RelTypes.values()[(3-binary[rt_name.charAt(1)])*20+(3-binary[rt_name.charAt(0)])*4+(rt_name.charAt(2)=='0'?3:(rt_name.charAt(2)-48)%3)];
            has_fwd=src.hasRelationship(f_rt, Direction.OUTGOING);
            has_rev=des.hasRelationship(r_rt, Direction.OUTGOING);
            if( !has_fwd && !has_rev)
                {
                    src.createRelationshipTo( des, f_rt);
                    ++num_edges;
                }
        }
    } 
    /*
    This function splits 'node' at position 'pos' by creating a 'split_node' as a part has been seperated from 'node'   
    */ 
    void split(Node node, int pos)
    {
        int l=(int)node.getProperty("length");
        int i,k=0,s_id,s_base,d_base,s_side=0,d_side=0,gen,seq,loc;
        long j,split_first_kmer,node_last_kmer=0;
        int[] address;
        Node neighbor;
        Relationship rel;
        address=(int[])node.getProperty("address");
        gen=address[0];
        seq=address[1];
        loc=address[2];
        ++seq_nodes;
        split_node = graphDb.createNode(node_label);
        split_node.setProperty("address",new int[]{gen,seq,loc+pos});
        split_node.setProperty("length",l-pos);
        String prp;
        for(Relationship r:node.getRelationships(RelTypes.begin, Direction.INCOMING))// To update genes whose beginning lies in split_node 
            {
                if((int)r.getProperty("pos")>=pos)
                {
                    rel=r.getStartNode().createRelationshipTo(split_node, RelTypes.begin);
                    rel.setProperty("side", r.getProperty("side"));
                    rel.setProperty("pos", (int)r.getProperty("pos")-pos);
                    r.delete();
                }
            }
        for(Relationship r:node.getRelationships(RelTypes.end, Direction.INCOMING)) // To update genes whose ending lies in split_node
            {
                if((int)r.getProperty("pos")-K+1>=pos)
                {
                    rel=r.getStartNode().createRelationshipTo(split_node, RelTypes.end);
                    rel.setProperty("side", r.getProperty("side"));
                    rel.setProperty("pos", (int)r.getProperty("pos")-pos);
                    r.delete();
                }
            }
        for(i=1;i<=genome;++i) // initial occurence properties for the split_node
            for(j=1;j<=num_sequences[i];++j)
            {
                prp="F"+i+"_"+j;
                if(node.hasProperty(prp) )
                {
                    int[] node_positions=(int[])node.getProperty(prp);
                    int len=node_positions.length;
                    int[] split_positions=new int[len];
                    for(k=0;k<len;++k)
                        split_positions[k]=node_positions[k]+pos;
                    split_node.setProperty(prp,split_positions);
                }
                prp="R"+i+"_"+j;
                if(node.hasProperty(prp) )
                {
                    int[] node_positions=(int[])node.getProperty(prp);
                    int len=node_positions.length;
                    int[] split_positions=new int[len];
                    for(k=0;k<len;++k)
                    {
                        split_positions[k]=node_positions[k];
                        node_positions[k]=node_positions[k]+l-pos-K+1;
                    }
                    split_node.setProperty(prp,split_positions);
                    node.setProperty(prp,node_positions);
                }
            }
        for(split_first_kmer=(long)node.getProperty("first_kmer"),i=0;i<pos;++i)// split the k-mer chain
        {
            node_last_kmer=split_first_kmer;
            split_first_kmer=index.pointers[(int)(node_last_kmer/2100000000)][(int)(node_last_kmer%2100000000)].next_index;
        }
        index.pointers[(int)(node_last_kmer/2100000000)][(int)(node_last_kmer%2100000000)].next_index=-1; 
        split_node.setProperty("first_kmer",split_first_kmer);
        split_node.setProperty("last_kmer",node.getProperty("last_kmer"));
        s_id=(int)split_node.getId();
        for(i=0,j=split_first_kmer;j!=-1;j=index.pointers[(int)(j/2100000000)][(int)(j%2100000000)].next_index,++i) // update kmer coordinates
        {
            index.pointers[(int)(j/2100000000)][(int)(j%2100000000)].node_id=s_id;
            index.pointers[(int)(j/2100000000)][(int)(j%2100000000)].position=i;
        }     
        s_base=binary[Genomes[gen][seq][loc+l-K]];
        for(d_base=0;d_base<=4;++d_base ) // update incoming and outgoing edges
            for(d_side=0;d_side<=1;++d_side)
                for(Relationship r:node.getRelationships(RelTypes.values()[s_base*20+d_base*4+d_side] , Direction.OUTGOING))
                {
                    neighbor=r.getEndNode();
                    if(neighbor.equals(node))
                        neighbor=d_side==0?node:split_node;
                    split_node.createRelationshipTo(neighbor,RelTypes.values()[s_base*20+d_base*4+d_side]);
                    r.delete();
                }
        d_base=3-binary[Genomes[gen][seq][loc+l-K]];
        for(s_base=0;s_base<=6;++s_base )
            for(s_side=0;s_side<=1;++s_side)
                for(Relationship r:node.getRelationships(RelTypes.values()[s_base*20+d_base*4+s_side*2+1] , Direction.INCOMING))
                {
                    neighbor=r.getStartNode();
                    if(neighbor.equals(node))
                        neighbor=node;
                    neighbor.createRelationshipTo(split_node,RelTypes.values()[s_base*20+d_base*4+s_side*2+1]);
                    r.delete();
                }
        node.setProperty("last_kmer",node_last_kmer);
        node.setProperty("length",pos+K-1);  
        connect(node,split_node,0,0,0);
    }
    /*
    This function extends a new node till reach to a previously visited K-mer
    */
    void extend(Node node)
    {
        int begin;
        long id=node.getId(),last_kmer;
        boolean broke=false;
        while(position<seq_len)
        {
            if(binary[Genomes[genome][sequence][position+1]]==4)
            {
                ++position;
                begin=position-K+1;
                jump();
                create_degenerate(begin);
                curr_node=degenerate_node;
                curr_side=0;
                break;
            }
            next_kmer();
            curr_index=index.find(kmer);
            pointer=index.pointers[(int)(curr_index/2100000000)][(int)(curr_index%2100000000)];
            if(pointer.node_id == -1)
            {
                node.setProperty("length",(int)node.getProperty("length")+1);
                pointer.node_id=id;
                pointer.format=(byte)(canonical?0:1);
                pointer.position=(int)node.getProperty("length")-K;
                last_kmer=(long)node.getProperty("last_kmer");
                index.pointers[(int)(last_kmer/2100000000)][(int)(last_kmer%2100000000)].next_index=curr_index; //adding kmer at the end of kmer_chain of current node
                node.setProperty("last_kmer",curr_index);
                if(position%(seq_len/100+1)==0) 
                    System.out.print((long)position*100/seq_len+1+"%\r");
            }
            else
            {
                broke=true;
                break;
            }
        }
        if(!broke && position==seq_len )
            finish=true;
    }
    /*
    This function creates a new node
    */
    void create()
    {
        ++seq_nodes;
        new_node=graphDb.createNode(node_label);
        new_node.setProperty("F"+genome+"_"+sequence,new int[]{position-K+1});
        new_node.setProperty("address",new int[] {genome,sequence,position-K+1});
        new_node.setProperty("length",K);
        new_node.setProperty("last_kmer",curr_index);
        new_node.setProperty("first_kmer",curr_index);
        pointer.node_id=new_node.getId();
        pointer.format=(byte)(canonical?0:1);
        pointer.position=0;
        connect(curr_node,new_node,curr_side,0,0);
        curr_node=new_node;
        curr_side=0;
    }
    /*
    This function enters a node and follow it in reverse direction till reaching to a split position or end of the node
    */
    void follow_forward()
    {
        int l,p=1,pos, begin,g,s,loc;
        pos=(int)pointer.position;
        Node node;
        int[] address;
        boolean degenerated=false;
        node=graphDb.getNodeById(pointer.node_id);
        String prp;
        if(pos>0 )
        {
            split(node,pos);
            if(curr_node.equals(node) && curr_side==0)
                   curr_node=split_node;
            node=split_node;
        }
        connect(curr_node,node,curr_side,0,0);                
        curr_node=node;
        curr_side=0;
        l=(int)curr_node.getProperty("length")-K;
        address=(int[])curr_node.getProperty("address");
        g=address[0];
        s=address[1];
        loc=address[2];
        p=position-K+1;
        for(pos=0;pos<=l && position<=seq_len && Genomes[g][s][loc+pos+K-1]==Genomes[genome][sequence][position]; ++pos)
        {
            ++position;
            if(position<=seq_len && binary[Genomes[genome][sequence][position]]==4)
            {
                begin=position-K+1;
                jump();
                if(pos+1<=l)
                    split(curr_node,pos+1);
                degenerated=true;
                create_degenerate(begin);
                break;
            }
        }
        if(position==seq_len+1 )
            finish=true;
        else if(!degenerated)
        {
            position-=K;
            initial_kmers();
        }
        if( !degenerated && pos<=l )
        {
            split(curr_node,pos);
        }
        prp="F"+genome+"_"+sequence;
        if(curr_node.hasProperty(prp))
        {
            int[] positions=(int[])curr_node.getProperty(prp);
            int i,len= positions.length;
            int[] new_positions=new int[len+1];
            for(i=0;i<len;++i)
                new_positions[i]=positions[i];
            new_positions[i]=p;
            curr_node.setProperty(prp, new_positions);
        }
        else
            curr_node.setProperty(prp,new int[]{p});
        curr_index=index.find(kmer);
        if(curr_index!=-1)
            pointer=index.pointers[(int)(curr_index/2100000000)][(int)(curr_index%2100000000)];
        if(degenerated)
        {
            curr_node=degenerate_node;
            curr_side=0;
        }
    }
    /*
    This function enters a node and follow it in forward direction till reaching to a split position or end of the node
    */
    void follow_reverse()
    {
        int p=0,pos, begin,g,s,loc;
        int[] address;
        pos=(int)pointer.position;
        Node node;
        boolean degenerated=false;
        node=graphDb.getNodeById(pointer.node_id);
        String prp;
        if(pos<(int)node.getProperty("length")-K) 
        {
            split(node,pos+1);
            if(curr_node.equals(node) && curr_side==0)
               curr_node=split_node;
        }
        connect(curr_node,node,curr_side,1,0);                
        curr_node=node;
        curr_side=1;
        p=position-K+1;
        address=(int[])curr_node.getProperty("address");
        g=address[0];
        s=address[1];
        loc=address[2];
        for(pos=(int)node.getProperty("length")-K;pos>=0 && position<=seq_len && Genomes[g][s][loc+pos]==complement[Genomes[genome][sequence][position]]; --pos) 
        {
            ++position;
            if(position<=seq_len && binary[Genomes[genome][sequence][position]]==4)
            {
                begin=position-K+1;
                jump();
                if(pos>0)
                {
                    split(curr_node,pos);
                    curr_node=split_node;
                }
                create_degenerate(begin);
                degenerated=true;
                break;
            }
        }
        if(position==seq_len+1 )
                finish=true;
        else if(!degenerated)
        {
            position-=K;
            initial_kmers();
        }
        if(!degenerated && pos>=0)
        {
            split(curr_node,pos+1);
            curr_node=split_node;
        }
        prp="R"+genome+"_"+sequence;
        if(curr_node.hasProperty(prp))
        {
            int[] positions=(int[])curr_node.getProperty(prp);
            int i,len= positions.length;
            int[] new_positions=new int[len+1];
            for(i=0;i<len;++i)
                new_positions[i]=positions[i];
            new_positions[i]=p;
            curr_node.setProperty(prp, new_positions);
        }
        else
            curr_node.setProperty(prp,new int[]{p});
        curr_index=index.find(kmer);
        if(curr_index!=-1)
            pointer=index.pointers[(int)(curr_index/2100000000)][(int)(curr_index%2100000000)];
        if(degenerated)
        {
            curr_node=degenerate_node;
            curr_side=0;
        }
    }
    /*
    This function jumps over an ambiguous region.
    position points to the first position which degenerate starts, after jumping it points to the last base of the first K-mer after the ambiguous region 
    */
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
    /*
    This function creates a degenerate node starting at "begin" ending at position-1
    and connects it to previous node and makes the first K-mer after the ambiguous region
    */
    void create_degenerate(int begin)
    {
        ++seq_nodes;
        degenerate_node=graphDb.createNode(degenerate_label);
        degenerate_node.setProperty("F"+genome+"_"+sequence,new int[]{begin});
        degenerate_node.setProperty("address",new int[] {genome,sequence,begin});
        degenerate_node.setProperty("length",position-begin);
        connect(curr_node,degenerate_node,curr_side,0,0);
        curr_index=index.find(canonical?f_kmer:r_kmer);
        if(curr_index!=-1)
            pointer=index.pointers[(int)(curr_index/2100000000)][(int)(curr_index%2100000000)];
    }
    /*
    This function produce the next "f_kmer", "r_kmer" and the canonical "kmer" from the current ones.
    */ 
    void next_kmer()
    {
        ++position;
        f_bin=binary[Genomes[genome][sequence][position]];
        r_bin=3-f_bin;
        f_kmer.next_up_kmer(f_bin,dim,K);
        r_kmer.next_down_kmer(r_bin,dim,K);
        if(position%(seq_len/100+1)==0) 
            System.out.print((long)position*100/seq_len+1+"%\r");
        canonical=f_kmer.compare_kmer(r_kmer,dim)==-1;
        kmer=canonical?f_kmer:r_kmer;
    }
    /*
    This function initialize the first K-mer of the "genome", "sequence" at "position". 
    It might jump over the degenerate regions creating a degenerate_node 
    */ 
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
                create_degenerate(0);
                curr_node=degenerate_node;
                curr_side=0;
                break;
            }
            next_kmer();
        }   
        canonical=f_kmer.compare_kmer(r_kmer,dim)==-1;
        kmer=canonical?f_kmer:r_kmer;
    }
    /*
    This function gives the K-mer in "genome", "sequence" at "position". 
    f_kmer: Forward sequence of the K-mer
    r_kmer: Reverse sequence of the K-mer
    kmer  : canonical k-mer
    canonical: a boolean determining if the K-mer is canonical  
    */ 
    pointer make_kmer(int genome, int sequence, int position)
    {
        int j,f_bin,r_bin;
        pointer f_kmer=new pointer(dim);
        pointer r_kmer=new pointer(dim);
        for(j=0;j<K;++j)
        {
            f_bin=binary[Genomes[genome][sequence][position+j]];
            r_bin=3-f_bin;
            f_kmer.next_up_kmer(f_bin,dim,K);
            r_kmer.next_down_kmer(r_bin,dim,K);
        }   
        annotate_canonical=f_kmer.compare_kmer(r_kmer,dim)==-1;
        return annotate_canonical?f_kmer:r_kmer;
    }
    /*
    The main function constructs the pangenome of input sequences.
    */ 
    void construct_pangenome() 
    {
        int i;
        Node genome_node,sequence_node;
        long[] sequence_ids;
        phaseTime = System.currentTimeMillis();
        for(genome=old_num_genomes+1;genome<=num_genomes;++genome) 
        {
            try(Transaction tx = graphDb.beginTx()){
                genome_node=graphDb.createNode(genome_label);
                genome_node.setProperty("genome_number", genome);
                genome_node.setProperty("genome_path",genome_paths[genome]);
                genome_node.setProperty("num_sequences",num_sequences[genome]);
                db_node.createRelationshipTo(genome_node, RelTypes.has);
                System.out.println("Processing genome "+genome+" :");
            tx.success();}
            sequence_ids=new long[num_sequences[genome]+1];
            for(sequence=1;sequence<=num_sequences[genome];++sequence) 
            {
                System.out.println("sequence "+sequence+"/"+num_sequences[genome]+" of genome "+genome+"\tlength="+sequence_length[genome][sequence]);
                try(Transaction tx = graphDb.beginTx()){
                    sequence_node=curr_node=graphDb.createNode(sequence_label);
                    sequence_node.setProperty("sequence_number",genome+"_"+sequence);
                    sequence_node.setProperty("sequence_name",sequence_names[genome][sequence]);
                    sequence_node.setProperty("sequence_length",sequence_length[genome][sequence]);
                    genome_node.createRelationshipTo(sequence_node, RelTypes.has);
                    sequence_ids[sequence]=curr_node.getId();
                    curr_side=0;
                    position=-1;
                    seq_len=sequence_length[genome][sequence]-1;
                    initial_kmers();
                tx.success();}
                curr_index=index.find(kmer);
                pointer=index.pointers[(int)(curr_index/2100000000)][(int)(curr_index%2100000000)];
                finish=false;
                while(!finish)
                {
                    try(Transaction tx = graphDb.beginTx()){
                        for(i=0;i<trsc_limit && !finish;++i)
                        {
                            if (pointer.node_id == -1) // kmer is new
                            {
                                create();
                                extend(curr_node);
                            }
                            else if( canonical ^ (pointer.format==0) )// if sides don't agree
                                follow_reverse();
                            else
                                follow_forward();
                        }
                    tx.success();}
                }//while position<n
                try(Transaction tx = graphDb.beginTx()){
                    connect(sequence_node,curr_node,0,1-curr_side,1);// to point to the last k-mer of the sequence located in the other strand
                tx.success();}
            }//sequences
            try(Transaction tx = graphDb.beginTx()){
                genome_node.setProperty("sequence_ids", sequence_ids);
            tx.success();}
            System.out.println("Running time : " + (System.currentTimeMillis() - phaseTime) / 1000 + " seconds");
        }//genomes
        Index_genomes();
    }
    /*
    This function adds list of "anchor_nodes", "anchor_sides" and "anchor_positions" to each sequence_node.
    These properties facilitate extracting genomic regions. 
    */ 
    void Index_genomes() 
    {
        int i,m,loc,n,node_len,s_side,d_side,count;
        long[] anchor_nodes;
        int[] anchor_positions;
        Node node,neighbor,seq_node;
        StringBuilder nds=new StringBuilder();
        StringBuilder pos=new StringBuilder();
        StringBuilder sds=new StringBuilder();
        String[] ids_list,posis_list;
        String lable, prp;
        boolean found;
        for(genome=1;genome<=num_genomes;++genome) 
        {
            for(sequence=1;sequence<=num_sequences[genome];++sequence) 
            {
                //System.out.println("Indexing sequence "+genome+"_"+sequence+"...");
                try(Transaction tx = graphDb.beginTx()){
                n=sequence_length[genome][sequence]-1;
                seq_node=node=graphDb.findNode(sequence_label, "sequence_number", genome+"_"+sequence);
                prp=genome+"_"+sequence;
                node_len=K-1;
                found=true;
                count=0;
                for(loc=0;loc+node_len-K+1<=n && found;)
                {
                    found=false;
                    for(Relationship r:node.getRelationships(Direction.OUTGOING))
                    {
                        neighbor=r.getEndNode();
                        d_side=(r.getType().name().charAt(2)-48)%2;
                        lable=(d_side==0?"F"+prp:"R"+prp);
                        if(neighbor.hasProperty(lable) )
                        {
                            if(b_search((int[])neighbor.getProperty(lable), loc+node_len-K+1)>=0)
                            {
                                found=true;
                                if(count % anchor_distance==0)
                                {
                                    nds.append(neighbor.getId()).append(" ");
                                    sds.append(d_side);
                                    pos.append(loc+node_len-K+1).append(" ");
                                }
                                count++;
                                loc=loc+node_len-K+1;
                                node=neighbor;
                                node_len=(int)node.getProperty("length");
                                break;
                            }
                        }
                    }
                    if(!found)
                        for(Relationship r:node.getRelationships(Direction.INCOMING))
                        {
                            neighbor=r.getStartNode();
                            s_side=(r.getType().name().charAt(2)-48)/2;
                            lable=(s_side==0?"R"+prp:"F"+prp);
                            if(neighbor.hasProperty(lable) )
                            {
                                if(b_search((int[])neighbor.getProperty(lable), loc+node_len-K+1)>=0)
                                {
                                    found=true;
                                    if(count % anchor_distance==0)
                                    {
                                        nds.append(neighbor.getId()).append(" ");
                                        sds.append(s_side);
                                        pos.append(loc+node_len-K+1).append(" ");
                                    }
                                    count++;
                                    loc=loc+node_len-K+1;
                                    node=neighbor;
                                    node_len=(int)node.getProperty("length");
                                    break;
                                }
                            }
                        }
                    if(loc%(n/100+1)==0) 
                        System.out.print((long)loc*100/n+1+"%    \r");
                }
                m=sds.length();
                ids_list=nds.toString().split("\\s");
                posis_list=pos.toString().split("\\s");
                anchor_nodes=new long[m];
                anchor_positions=new int[m];
                for(i=0;i<m;++i)
                {
                    anchor_nodes[i]=Long.valueOf(ids_list[i]);
                    anchor_positions[i]=Integer.valueOf(posis_list[i]);
                }
                seq_node.setProperty("anchor_nodes", anchor_nodes);
                seq_node.setProperty("anchor_positions", anchor_positions);
                seq_node.setProperty("anchor_sides", sds.toString());
                nds.setLength(0);
                pos.setLength(0);
                sds.setLength(0);
                tx.success();}
            }//sequences
        }//genomes
    }
    /*
    This function adds "sequence" property to the nodes and writes the K-mer index on disk
    */ 
    void finalise() 
    {
        int i,len;
        int[] address;
        byte[] sequence;
        num_bases=0;
        ResourceIterator<Node> nodes;
        Node node;
        System.out.println("Finalizing the construction...");
        try(Transaction tx = graphDb.beginTx()){
        nodes = graphDb.findNodes( node_label );
        tx.success();}
        while(nodes.hasNext())
            try(Transaction tx = graphDb.beginTx()){
                for(i=0;i<trsc_limit && nodes.hasNext();++i)
                {
                    node=nodes.next();
                    address=(int[])node.getProperty("address");
                    len=(int)node.getProperty("length");
                    num_bases+=len;
                    sequence=new byte[len];
                    for(i=0;i<len;++i)
                        sequence[i]=(byte)Genomes[address[0]][address[1]][address[2]+i];
                    node.setProperty("sequence",sequence);
                }
            tx.success();}
        nodes.close();
        try(Transaction tx = graphDb.beginTx()){
        nodes = graphDb.findNodes( degenerate_label );
        tx.success();}
        while(nodes.hasNext())
            try(Transaction tx = graphDb.beginTx()){
                for(i=0;i<trsc_limit && nodes.hasNext();++i)
                {
                    node=nodes.next();
                    address=(int[])node.getProperty("address");
                    len=(int)node.getProperty("length");
                    num_bases+=len;
                    sequence=new byte[len];
                    for(i=0;i<len;++i)
                        sequence[i]=(byte)Genomes[address[0]][address[1]][address[2]+i];
                    node.setProperty("sequence",sequence);
                }
            tx.success();}
        nodes.close();
        try
        {
            index.write(PATH+INDEX);
        }
        catch(IOException ioe)
        {
           System.out.println("Failed to read file names!");  
        } 
        System.out.println("Number of genomes: "+num_genomes); 
        System.out.println("Number of kmers:   "+index.totallength);
        System.out.println("Number of nodes:   "+seq_nodes);
        System.out.println("Number of edges:   "+num_edges);
        System.out.println("Number of bases:   "+num_bases);
    }
    /*
    This function registers the action to be taken if the program halts unexpectedly
    */ 
    private static void registerShutdownHook( final GraphDatabaseService graphDb )
    {
        Runtime.getRuntime().addShutdownHook( new Thread()
        {
            @Override
            public void run()
            {
                graphDb.shutdown();
            }
        } );
    }    
    /*
    This function executes command "cmd" in unix shell
    */ 
    public static String executeCommand(String cmd) {
        String[] command = {"/bin/sh", "-c", cmd};
        exe_output=new StringBuilder();
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
    /*
    This function compares two graph databases located in paths "path1" and "path2" and return the boolean result
    */    
    private boolean compare_pan(String path1, String path2)
    {
        int sq_num1=0,sq_num2=0,ed_num1,ed_num2,k1,k2,ng1,ng2,len1,len2; 
        long k,knum1=0,knum2;
        Index index1=null,index2=null;
        Node db_node1,db_node2;
        if (new File(path1+DATABASE).exists() && new File(path2+DATABASE).exists()) 
        {
            GraphDatabaseService g1 = new GraphDatabaseFactory().newEmbeddedDatabase(path1+DATABASE );
            GraphDatabaseService g2 = new GraphDatabaseFactory().newEmbeddedDatabase(path2+DATABASE );
            registerShutdownHook( g1 );        
            registerShutdownHook( g2 );
            try(Transaction tx1 = g1.beginTx();Transaction tx2 = g2.beginTx();)
            {
                db_node1=g1.findNodes(pangenome_label).next();
                db_node2=g2.findNodes(pangenome_label).next();
                if(db_node1==null)
                {
                    System.out.println("Can not locate database node for "+path1);
                    System.exit(1);
                }
                if(db_node2==null)
                {
                    System.out.println("Can not locate database node for "+path2);
                    System.exit(1);
                }                
                k1=(int)db_node1.getProperty("k_mer_size");
                k2=(int)db_node2.getProperty("k_mer_size");
                sq_num1=(int)db_node1.getProperty("num_nodes");
                sq_num2=(int)db_node2.getProperty("num_nodes");
                ed_num1=(int)db_node1.getProperty("num_edges");
                ed_num2=(int)db_node2.getProperty("num_edges");
                ng1=(int)db_node1.getProperty("num_genomes");
                ng2=(int)db_node2.getProperty("num_genomes");
                knum1=(long)db_node1.getProperty("num_k_mers");
                knum2=(long)db_node2.getProperty("num_k_mers");
                if(k1!=k2 || sq_num1!=sq_num2 || sq_num1!=sq_num2 || ed_num1!=ed_num2 || ng1!=ng2 || knum1!=knum2)
                    return false;
                try{
                    index1=new Index(path1+INDEX,knum1,k1);
                    index2=new Index(path2+INDEX,knum2,k2);
                }
                catch(IOException ioe){
                    System.out.println("Failed to open a file!");  
                    System.exit(1);
                }  
                System.out.println("comparing pangenomes...");        
                Node n1,n2;
                byte[] s1,s2;
                    for(k=0;k<knum1;++k)
                    {
                        if(index1.pointers[(int)(k/2100000000)][(int)(k%2100000000)].compare_kmer(index2.pointers[(int)(k/2100000000)][(int)(k%2100000000)],K)!=0)
                            return false;
                        if(index1.pointers[(int)(k/2100000000)][(int)(k%2100000000)].next_index==-1 && index2.pointers[(int)(k/2100000000)][(int)(k%2100000000)].next_index==-1)
                        {
                            n1=g1.getNodeById(index1.pointers[(int)(k/2100000000)][(int)(k%2100000000)].node_id);
                            s1=(byte[])n1.getProperty("sequence");
                            len1=(int)n1.getProperty("length");
                            n2=g2.getNodeById(index2.pointers[(int)(k/2100000000)][(int)(k%2100000000)].node_id);
                            s2=(byte[])n2.getProperty("sequence");
                            len2=(int)n2.getProperty("length");
                            if(len1!=len2 || (!are_equal(s1,0,s2,0,len1,true) && !are_equal(s1,0,s2,0,len1,false)) )
                                return false;
                            if(n1.getDegree(Direction.BOTH)!=n2.getDegree(Direction.BOTH))
                                return false;
                        }
                        if(k%(knum1/100+1)==0) 
                                System.out.print(k*100/knum1+1+"%\r");
                    }
                    System.out.println();
                tx1.success();         
                tx2.success();
            }         
            g1.shutdown();
            g2.shutdown();
            return true;
        }
        System.out.println("Error: a pair of databases are needed to be compared!");
        return false;
    }
    /*
    This function splits the "sequence" to lines of length "length" and append the lines to "fasta_file"
    */
    private void write_fasta(BufferedWriter fasta_file, String seq, int length) throws IOException
    {
        int i;
        for(i=1;i<=seq.length();++i)
        {
            fasta_file.write(seq.charAt(i-1));
            if(i%length==0)
                fasta_file.write("\n");
        }
        fasta_file.write("\n");
    }
    /*
    This function determines the equality of two subsequences of s1 and s2 of length len starting at start1 and start2, respectively.
    forward determines the direction of comparion.
    */
    private boolean are_equal(byte[] s1, int start1,byte[] s2, int start2, int len, boolean forward)
    {
        int i;
        boolean equal=true;
        if(forward)
            for(i=0;i<len && equal;++i)
                if(s1[start1+i]!=s2[start2+i])
                    equal=false;
        else
            for(i=0;i<len && equal;++i)
                if(s1[start1+i]!=complement[s2[start2+len-i-1]])
                    equal=false;
        return equal;
    }
    /*
    This function provides a simple binary search functionality.
    */
    private int b_search(int[] array, int key)
    {
        int low=0,mid,high=array.length-1;
        while(low<=high)
        {
            mid=(low+high)/2;
            if(key<array[mid])
                high=mid-1;
            else if(key>array[mid])
                low=mid+1;
            else return mid;
        }
        return -1;
    }
    /*
    This function calculates and prints peak of memory usage of the program in mega bytes.
    */
    private void print_peak_memory()
    {
        long memoryUsage=0;
    try {
            for (MemoryPoolMXBean pool : ManagementFactory.getMemoryPoolMXBeans()) 
            {
                MemoryUsage peak = pool.getPeakUsage();
                memoryUsage += peak.getUsed();
            }
            System.out.println("Peak memory : "+ memoryUsage/1024/1024+" MB");
        } 
    catch (Throwable t) 
        {
            System.err.println("Exception in agent: " + t);
        }
    }
}
