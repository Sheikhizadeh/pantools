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
import java.util.Arrays;

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
    "retrieve_genes:    To extract sequence of some annotated genes\n" +
    "          java -jar pantools.jar retrieve_genes database_path annotation_records\n" +
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
    "annotation_records : A text file containing gene annotation records to be retrieved\n" +
    "\n" +
    "regions_file    : A text file containing genome_number, sequence_number, begin, end of a region in each line seperated by one space \n" +
    "\n" +
    "group_file_file : A text file with each line starting with a group name followed by a colon, followed by genes IDs seperated by one space";

    public static String PATH;
    public static String DATABASE="/graph.db";
    public static GraphDatabaseService graphDb;
    public static Index index;
    public static genome_bank genomes;
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
    public static int[] binary;
    public static int[] complement;

    public static char sym[]={ 'A', 'C', 'G' , 'T', 'M','R','W','S','Y','K','V','H','D','B','N'};
    
    private kmer f_kmer,r_kmer,k_mer;
    private long seq_len;
    private int genome,sequence,position;
    private static long curr_index;
    private static byte curr_side;
    private Node curr_node;
    private Node db_node;
    private Node new_node;
    private Node split_node;
    private Node degenerate_node;
    private pointer pointer=new pointer();
    private byte[] byte_pointer=new byte[21];
    
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
        binary=new int[256];
        complement=new int[]{3,2,1,0,9,8,6,7,5,4,13,12,11,10,14};
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
        System.out.println("------------------------------- PanTools -------------------------------");        
        switch (args[0]) {
            case "build":
                K=Integer.parseInt(args[1]);
                if(K<10)
                {
                    System.out.println("K should be greater than 9!");
                    System.exit(1);
                }
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
    private void build(String file) throws IOException
    {
        int i;
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
        dim=(int)Math.ceil(K/4.0);
        seq_nodes=0;
        num_edges=0;
        try(Transaction tx = graphDb.beginTx()){
            db_node=graphDb.createNode(pangenome_label);
            db_node.setProperty("k_mer_size", K);
        tx.success();}
        genomes=new genome_bank(file);
        index=new Index(PATH,file,K);
        construct_pangenome();
        System.out.println("Number of kmers:   "+index.length());
        System.out.println("Number of nodes:   "+seq_nodes);
        System.out.println("Number of edges:   "+num_edges);
        System.out.println("Number of bases:   "+num_bases);
        try(Transaction tx = graphDb.beginTx()){
            db_node.setProperty("num_k_mers", index.length());
            db_node.setProperty("num_nodes",seq_nodes);
            db_node.setProperty("num_edges",num_edges);
            db_node.setProperty("num_genomes",genomes.num_genomes);
            db_node.setProperty("num_genes",0);
            db_node.setProperty("num_bases",num_bases);
            tx.success();}
        graphDb.shutdown();  
        System.out.println("Total time : "+(System.currentTimeMillis()-startTime)/1000+"."+(System.currentTimeMillis()-startTime)%1000+" seconds");
        print_peak_memory();
        executeCommand("rm "+PATH+"/database.db/neostore.transaction.db*");
        System.out.print(Pantools.executeCommand("du --apparent-size -shc -m "+PATH+DATABASE));
        System.out.print(Pantools.executeCommand("du --apparent-size -shc -m "+PATH+"/index*"));
    }
    /*
    This function adds new genomes to an available pan-genome database.
    fasta_names is a text file containing paths to FASTA files; each in a new line
    */      
    private void add(String file) throws IOException
    {
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
                dim=(int)Math.ceil(K/4.0);
                seq_nodes=(int)db_node.getProperty("num_nodes");
                num_edges=(int)db_node.getProperty("num_edges");
                gene_nodes=(int)db_node.getProperty("num_genes");
                genomes=new genome_bank(file);
                index=new Index(PATH,file, graphDb);
            tx.success();} 
            construct_pangenome();
            System.out.println("Number of kmers:   "+index.length());
            System.out.println("Number of nodes:   "+seq_nodes);
            System.out.println("Number of edges:   "+num_edges);
            System.out.println("Number of bases:   "+num_bases);
            try(Transaction tx = graphDb.beginTx())
            {
                db_node.setProperty("k_mer_size", K);
                db_node.setProperty("num_k_mers", index.length());
                db_node.setProperty("num_nodes",seq_nodes);
                db_node.setProperty("num_edges",num_edges);
                db_node.setProperty("num_genomes",genomes.num_genomes);
                db_node.setProperty("num_genes",gene_nodes);
                db_node.setProperty("num_bases",num_bases);
                tx.success();
            }
            graphDb.shutdown(); 
            System.out.println("Total time : "+(System.currentTimeMillis()-startTime)/1000+"."+(System.currentTimeMillis()-startTime)%1000+" seconds");
            print_peak_memory();
            executeCommand("rm "+PATH+"/database.db/neostore.transaction.db*");
            System.out.print(Pantools.executeCommand("du --apparent-size -shc -m "+PATH+DATABASE));
            System.out.print(Pantools.executeCommand("du --apparent-size -shc -m "+PATH+"/index*"));
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
    private void annotate(String gff_names, String file) throws IOException
    {
        int i;
        int begin,end,g,s,line_len;
        String seq,id,name;
        boolean start_canonical,stop_canonical;
        Node gene_node,genome_node;
        pointer start_point,stop_point;
        Relationship rel;
        String[] fields;
        String strand,fasta_name,gff_name,line;
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
                dim=(int)Math.ceil(K/4.0);
                seq_nodes=(int)db_node.getProperty("num_nodes");
                num_edges=(int)db_node.getProperty("num_edges");
                gene_nodes=(int)db_node.getProperty("num_genes");
                genomes=new genome_bank(file);
                startTime = System.currentTimeMillis();
                try(BufferedReader gff_files = new BufferedReader(new FileReader(gff_names))){
                BufferedReader fasta_files = new BufferedReader(new FileReader(file));
                 while(gff_files.ready())
                    {
                        gff_name=gff_files.readLine();
                        fasta_name=fasta_files.readLine();
                        genome_node=graphDb.findNode(genome_label, "name", fasta_name);
                        if(genome_node.hasProperty("annotated"))
                        {
                            System.out.println("genome "+fasta_name+" already annotated.");
                            continue;
                        }
                        g=(int)genome_node.getProperty("number");
                        BufferedReader in = new BufferedReader(new FileReader(gff_name));
                        System.out.println("Annotating genome "+g);
                        while(in.ready())
                        {
                            line=in.readLine();
                            fields=line.split("\\t");
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
                                    start_point=find_node(g,s,begin);
                                    start_canonical=annotate_canonical;
                                    stop_point=find_node(g,s,end-K+1);
                                    stop_canonical=annotate_canonical;
                                    if(start_point!=null && stop_point!=null && end-begin+1>=K) // genes should not be shorter than K
                                    {
                                        ++gene_nodes;
                                        gene_node=graphDb.createNode(gene_label);
                                        gene_node.setProperty("genome_number", g);
                                        gene_node.setProperty("sequence_number", s);
                                        gene_node.setProperty("begin", begin+1);
                                        gene_node.setProperty("end", end+1);
                                        gene_node.setProperty("length", end-begin+1);
                                        gene_node.setProperty("strand", strand);
                                        gene_node.setProperty("record",line);
                                        rel=gene_node.createRelationshipTo(graphDb.getNodeById(start_point.node_id), RelTypes.begin);
                                        rel.setProperty("side", start_canonical ?start_point.format:1-start_point.format);
                                        rel.setProperty("pos", start_point.position);
                                        rel=gene_node.createRelationshipTo(graphDb.getNodeById( stop_point.node_id), RelTypes.end);
                                        rel.setProperty("side", stop_canonical ?stop_point.format:1-stop_point.format);
                                        rel.setProperty("pos", stop_point.position+K-1);
                                    }
                                    else
                                        System.out.println("["+line+"] not properly called!"); 
                                }
                                else 
                                    System.out.println(seq +" missed in genome "+g);
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
                System.out.println("K="+K);
            tx.success();}
            ResourceIterator<Node> gene_nodes;
            Relationship rstart,rstop;
            Node start,stop,gene=null;
            String prp,gene_record,strand,record;
            int i,j,begin,end,num_genes,genome,sequence;
            String[] occ;
            StringBuilder gene_seq;
            num_genes=Integer.parseInt(executeCommand("wc -l "+file).trim().split("\\s")[0]);
            String[] records=new String[num_genes];
            try(Transaction tx = graphDb.beginTx()){
            try (BufferedReader in = new BufferedReader(new FileReader(file))) {
                BufferedWriter out = new BufferedWriter(new FileWriter(file+".fasta"));
                for(i=0;i<num_genes;++i)
                {
                    gene_record=in.readLine();
                    records[i]=gene_record;
                }
                Arrays.sort(records);
                for(i=1,j=0,gene_nodes=graphDb.findNodes(gene_label);gene_nodes.hasNext();++i)
                {
                    gene=gene_nodes.next();
                    record=(String)gene.getProperty("record");
                    if(Arrays.binarySearch(records, record) >=0)
                    {
                        rstart=gene.getSingleRelationship(RelTypes.begin, Direction.OUTGOING);
                        rstop=gene.getSingleRelationship(RelTypes.end, Direction.OUTGOING);
                        start=rstart.getEndNode();
                        stop=rstop.getEndNode();
                        genome=(int)gene.getProperty("genome_number");
                        sequence=(int)gene.getProperty("sequence_number");
                        prp=genome+"_"+sequence;
                        begin=(int)gene.getProperty("begin");
                        end=(int)gene.getProperty("end");
                        strand=gene.getProperty("strand").toString();
                        gene_seq=sequence(start,stop,(int)rstart.getProperty("side"),(int)rstart.getProperty("pos"),prp, begin-1, end-1);
                        if(gene_seq.length() == end-begin+1)
                        {
                            ++j;
                            out.write(">"+record+"\n");
                            if(strand.equals("+"))
                                write_fasta(out,gene_seq.toString(),70);
                            else
                                write_fasta(out,reverse_complement(gene_seq.toString()),70);
                        }
                        else
                        {
                            System.out.println("Failed to assemble "+i+"\'th gene.");
                            //System.out.println(gene_seq.length()+" != "+(end-begin+1));
                        }
                        gene_seq.setLength(0);
                    }
                    else
                    {}//System.out.println(gene_record+" not annotated!");
                    if(i%(num_genes/100+1)==0) 
                        System.out.print((long)i*100/num_genes+1+"%\r");
                }//for i
                System.out.println(j+" genes retrieved successfully.");
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
                for(c=0,line=in.readLine();line!=null;line=in.readLine())
                {
                   fields=line.split("\\s");
                    g= Integer.parseInt(fields[0]);
                    s= Integer.parseInt(fields[1]);
                    from= Integer.parseInt(fields[2]);
                    to= Integer.parseInt(fields[3]);
                    seq_node=graphDb.findNode(sequence_label, "number", g+"_"+s);
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
            append_fwd(seq,(byte[])node.getProperty("sequence"),start_pos,start_pos+end-begin);
        else if(start_node.equals(stop_node) && start_side==1 && start_pos+K-1-(end-begin)>=0)
            append_rev(seq,(byte[])node.getProperty("sequence"),start_pos+K-1-(end-begin),start_pos+K-1);
        else
        {
            node_len=(int)node.getProperty("length");
            if(start_side==0)
            {
                if(start_pos+end-begin<(int)start_node.getProperty("length"))
                    append_fwd(seq,(byte[])node.getProperty("sequence"),start_pos,start_pos+end-begin);
                else
                    append_fwd(seq,(byte[])node.getProperty("sequence"),start_pos,(int)node.getProperty("length")-1);
                loc=begin-start_pos;
            }
            else
            {
            System.out.println("else");
                if(start_pos+K-1-(end-begin)>=0)
                    append_rev(seq,(byte[])node.getProperty("sequence"),start_pos+K-1-(end-begin),start_pos+K-1);
                else
                    append_rev(seq,(byte[])node.getProperty("sequence"),0,start_pos+K-1);
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
                                        append_fwd(seq,(byte[])neighbor.getProperty("sequence"),K-1,seq_len-seq.length()+K-2);
                                    else
                                        append_fwd(seq,(byte[])neighbor.getProperty("sequence"),K-1,neighbor_len-1);
                                else
                                    if(seq.length()+neighbor_len-K+1>seq_len)
                                        append_rev(seq,(byte[])neighbor.getProperty("sequence"),neighbor_len-K-(seq_len-seq.length())+1,neighbor_len-K);
                                    else
                                        append_rev(seq,(byte[])neighbor.getProperty("sequence"),0,neighbor_len-K);
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
                                        append_fwd(seq,(byte[])neighbor.getProperty("sequence"),K-1,seq_len-seq.length()+K-2);
                                    else
                                        append_fwd(seq,(byte[])neighbor.getProperty("sequence"),K-1,neighbor_len-1);
                                else
                                    if(seq.length()+neighbor_len-K+1>seq_len)
                                        append_rev(seq,(byte[])neighbor.getProperty("sequence"),neighbor_len-K-(seq_len-seq.length())+1,neighbor_len-K);
                                    else
                                        append_rev(seq,(byte[])neighbor.getProperty("sequence"),0,neighbor_len-K);
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
    This function appends s[from..to] to "seq". 
    */
    private void append_fwd(StringBuilder seq, byte[] s, int from, int to)
    {
        for(int i=from;i<=to;++i) 
            seq.append(sym[s[i]]);
        
    }
    /*
    This function appends reverse_complement of s[from..to] to "seq". 
    */
    private void append_rev(StringBuilder seq, byte[] s, int from, int to)
    {
        for(int i=to;i>=from;--i) 
            seq.append(sym[complement[s[i]]]);            
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
    pointer find_node(int g, int s, int pos) throws IOException
    {
        ResourceIterator<Node> degenerate_nodes;
        long inx;
        if(pos+K>genomes.sequence_length[g][s])
            return null;
        boolean ambiguous=false,canonical;
        for(int i=0;i<K && !ambiguous;++i)
            if(genomes.get(g,s,pos+i)>3)
                ambiguous=true;
        if(!ambiguous)
        {
            pointer p=new pointer();
            byte[] bp=new byte[21];
            inx=index.find(make_kmer(g,s,pos));
            index.get_pointer(bp,p,inx);
            return p;
        }
        else // degenerate node found
        {
            Node degenerate_node=null;
            int loc, lc, len;
            for(loc=pos;genomes.get(g,s,loc)<=3;)
                ++loc;
            lc=loc;
            do{
                loc=lc;
                for(lc=loc-1; lc>=loc-K+1;--lc)
                    if(lc<0 || genomes.get(g,s,lc)>3)
                        break;
            }while(lc>=0 && genomes.get(g,s,lc)>3);
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
                if(loc+len<=genomes.sequence_length[g][s] && genomes.compare((int[])degenerate_node.getProperty("address"),new int[]{g,s,loc},len,true))
                    found=true;
            }
            degenerate_nodes.close();
            tx.success();} 
            if(found)
                return new pointer(degenerate_node.getId(),(byte)0,(pos-loc),-1L);//,0);   
            else
            {
                System.out.println("Gap node not found!");
                return null;
            }
        }
    }
    /*
    This function returns the number of the sequence from genome "g" whose name contains "prefix".
    */    
    int find_sequence(String prefix, int g)
    {
        boolean found=false;
        int s;
        for(found=false,s=1;!found && s<=genomes.num_sequences[g];++s)
            if(genomes.sequence_titles[g][s].contains(prefix))
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
        for(int i=s.length()-1;i>=0;--i)
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
            s_base=s_side==0?genomes.get(s_gen,s_seq,s_loc+s_len-K):(3-genomes.get(s_gen,s_seq,s_loc+K-1));
            if(s_base<0 || s_base>3)
                s_base=4;
        }
        d_len= (int)des.getProperty("length") ;
        d_add=(int[])des.getProperty("address");
        d_gen=d_add[0];
        d_seq=d_add[1];
        d_loc=d_add[2];
        d_base=d_side==0?genomes.get(d_gen,d_seq,d_loc+K-1):(3-genomes.get(d_gen,d_seq,d_loc+d_len-K));
        if(d_base<0 || d_base>3)
            d_base=4;
        //System.out.println(s_base);
        //System.out.println(d_base);
        //System.out.println(s_side);
        //System.out.println(s_side);
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
        //System.out.println("split "+position);
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
            for(j=1;j<=genomes.num_sequences[i];++j)
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
            split_first_kmer=index.get_next_index(node_last_kmer);
        }
        index.put_next_index(-1L, node_last_kmer); 
        split_node.setProperty("first_kmer",split_first_kmer);
        split_node.setProperty("last_kmer",node.getProperty("last_kmer"));
        s_id=(int)split_node.getId();
        for(i=0,j=split_first_kmer;j!=-1L;j=index.get_next_index(j),++i) // update kmer coordinates
        {
            index.put_node_id(s_id, j);
            index.put_position(i, j);
        }     
        s_base=genomes.get(gen,seq,loc+l-K);
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
        d_base=3-genomes.get(gen,seq,loc+l-K);
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
    void extend(Node node) throws IOException
    {
        //System.out.println("extend "+position);
        int begin;
        long id=node.getId();
        boolean broke=false;
        while(position<seq_len)
        {
            //System.out.println("extend "+position);
            if(genomes.get(genome,sequence,position+1)>3)
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
            curr_index=index.find(k_mer.bytes);
            index.get_pointer(byte_pointer,pointer,curr_index);
            if(pointer.node_id == -1L)
            {
                index.put_next_index(curr_index,(long)node.getProperty("last_kmer"));
                node.setProperty("length",(int)node.getProperty("length")+1);
                pointer.node_id=id;
                pointer.format=(byte)(canonical?0:1);
                pointer.position=(int)node.getProperty("length")-K;
                pointer.next_index=-1L;
                index.put_pointer(pointer, curr_index);
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
        //System.out.println("create "+position);
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
        pointer.next_index=-1L;
        index.put_pointer(pointer, curr_index);
        connect(curr_node,new_node,curr_side,0,0);
        curr_node=new_node;
        curr_side=0;
    }
    /*
    This function enters a node and follow it in reverse direction till reaching to a split position or end of the node
    */
    void follow_forward() throws IOException
    {
        //System.out.println("forward "+position);
        int l,p=1,pos, begin,g,s,loc;
        Node node;
        int[] address;
        boolean degenerated=false;
        pos=pointer.position;
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
        for(pos=0;pos<=l && position<=seq_len && genomes.get(g, s, loc+pos+K-1)==genomes.get(genome,sequence,position); ++pos)
        {
            ++position;
            if(position<=seq_len && genomes.get(genome,sequence,position)>3)
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
        curr_index=index.find(k_mer.bytes);
        if(curr_index!=-1L)
            index.get_pointer(byte_pointer,pointer,curr_index);
        if(degenerated)
        {
            curr_node=degenerate_node;
            curr_side=0;
        }
    }
    /*
    This function enters a node and follow it in forward direction till reaching to a split position or end of the node
    */
    void follow_reverse()throws IOException
    {
        //System.out.println("reverse "+position);
        int p=0,pos, begin,g,s,loc;
        int[] address;
        Node node;
        boolean degenerated=false;
        pos=pointer.position;
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
        for(pos=(int)node.getProperty("length")-K;pos>=0 && position<=seq_len && genomes.get(g,s,loc+pos)==complement[genomes.get(genome,sequence,position)]; --pos) 
        {
            ++position;
            if(position<=seq_len && genomes.get(genome,sequence,position)>3)
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
        curr_index=index.find(k_mer.bytes);
        if(curr_index!=-1L)
            index.get_pointer(byte_pointer,pointer,curr_index);
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
        f_bin=genomes.get(genome,sequence,position);
        r_bin=3-f_bin;
        do{
            while(f_bin>3 && position<seq_len)
            {
                ++position;
                f_bin=genomes.get(genome,sequence,position);
                r_bin=3-f_bin;
            }
            f_kmer=new kmer(dim);
            r_kmer=new kmer(dim);
            f_kmer.next_up_kmer(f_bin,dim,K);
            r_kmer.next_down_kmer(r_bin,dim,K);
            for(j=0; j<K-1 && position<seq_len;++j)
            {
                ++position;
                f_bin=genomes.get(genome,sequence,position);
                r_bin=3-f_bin;
                if(f_bin>3 )
                    break;
                f_kmer.next_up_kmer(f_bin,dim,K);
                r_kmer.next_down_kmer(r_bin,dim,K);
            }
        }while(f_bin>3 && position<seq_len);
        canonical=f_kmer.compare(r_kmer,dim)==-1;
        k_mer=canonical?f_kmer:r_kmer;
        if(position==seq_len)
            finish=true;
    }
    /*
    This function creates a degenerate node starting at "begin" ending at position-1
    and connects it to previous node and makes the first K-mer after the ambiguous region
    */
    void create_degenerate(int begin) throws IOException
    {
        //System.out.println("create_degenerate "+position);
        ++seq_nodes;
        degenerate_node=graphDb.createNode(degenerate_label);
        degenerate_node.setProperty("F"+genome+"_"+sequence,new int[]{begin});
        degenerate_node.setProperty("address",new int[] {genome,sequence,begin});
        degenerate_node.setProperty("length",position-begin);
        connect(curr_node,degenerate_node,curr_side,0,0);
        curr_index=index.find(canonical?f_kmer.bytes:r_kmer.bytes);
        if(curr_index!=-1L)
            index.get_pointer(byte_pointer,pointer,curr_index);
    }
    /*
    This function produce the next "f_kmer", "r_kmer" and the canonical "kmer" from the current ones.
    */ 
    void next_kmer()
    {
        ++position;
        f_bin=genomes.get(genome,sequence,position);
        r_bin=3-f_bin;
        f_kmer.next_up_kmer(f_bin,dim,K);
        r_kmer.next_down_kmer(r_bin,dim,K);
        if(position%(seq_len/100+1)==0) 
            System.out.print((long)position*100/seq_len+1+"%\r");
        canonical=f_kmer.compare(r_kmer,dim)==-1;
        k_mer=canonical?f_kmer:r_kmer;
    }
    /*
    This function initialize the first K-mer of the "genome", "sequence" at "position". 
    It might jump over the degenerate regions creating a degenerate_node 
    */ 
    public void initial_kmers() throws IOException
    {
        int i;
        f_kmer=new kmer(dim);
        r_kmer=new kmer(dim);
        for(i=0;i<K && position<seq_len;++i)
        {
            if(genomes.get(genome,sequence,position+1)>3 )
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
        canonical=f_kmer.compare(r_kmer,dim)==-1;
        k_mer=canonical?f_kmer:r_kmer;
    }
    /*
    This function gives the K-mer in "genome", "sequence" at "position". 
    f_kmer: Forward sequence of the K-mer
    r_kmer: Reverse sequence of the K-mer
    kmer  : canonical k-mer
    canonical: a boolean determining if the K-mer is canonical  
    */ 
    byte[] make_kmer(int genome, int sequence, int position)
    {
        int j,f_bin,r_bin;
        kmer f_kmer=new kmer(dim);
        kmer r_kmer=new kmer(dim);
        for(j=0;j<K;++j)
        {
            f_bin=genomes.get(genome,sequence,position+j);
            r_bin=3-f_bin;
            f_kmer.next_up_kmer(f_bin,dim,K);
            r_kmer.next_down_kmer(r_bin,dim,K);
        }   
        annotate_canonical=f_kmer.compare(r_kmer,dim)==-1;
        return annotate_canonical?f_kmer.bytes:r_kmer.bytes;
    }
    /*
    The main function constructs the pangenome of input sequences.
    */ 
    void construct_pangenome() throws IOException
    {
        int i;
        Node genome_node,sequence_node;
        long[] sequence_ids;
        phaseTime = System.currentTimeMillis();
        for(genome=genomes.previous_num_genomes+1;genome<=genomes.num_genomes;++genome) 
        {
            try(Transaction tx = graphDb.beginTx()){
                genome_node=graphDb.createNode(genome_label);
                genome_node.setProperty("number", genome);
                genome_node.setProperty("num_sequences",genomes.num_sequences[genome]);
                db_node.createRelationshipTo(genome_node, RelTypes.has);
                System.out.println("Processing genome "+genome+" :             ");
            tx.success();}
            sequence_ids=new long[genomes.num_sequences[genome]+1];
            for(sequence=1;sequence<=genomes.num_sequences[genome];++sequence) 
            {
                System.out.println("sequence "+sequence+"/"+genomes.num_sequences[genome]+" of genome "+genome+"\tlength="+genomes.sequence_length[genome][sequence]);
                try(Transaction tx = graphDb.beginTx()){
                    sequence_node=curr_node=graphDb.createNode(sequence_label);
                    sequence_node.setProperty("number",genome+"_"+sequence);
                    sequence_node.setProperty("sequence_name",genomes.sequence_titles[genome][sequence]);
                    sequence_node.setProperty("sequence_length",genomes.sequence_length[genome][sequence]);
                    genome_node.createRelationshipTo(sequence_node, RelTypes.has);
                    sequence_ids[sequence]=curr_node.getId();
                    curr_side=0;
                    position=-1;
                    seq_len=genomes.sequence_length[genome][sequence]-1;
                    initial_kmers();
                tx.success();}
                curr_index=index.find(k_mer.bytes);
                index.get_pointer(byte_pointer,pointer,curr_index);
                finish=false;
                while(!finish)
                {
                    try(Transaction tx = graphDb.beginTx()){
                        for(i=0;i<trsc_limit && !finish;++i)
                        {
                            if (pointer.node_id == -1L) // kmer is new
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
        int i,m,loc,node_len,s_side,d_side,count;
        long seq_len;
        long[] anchor_nodes;
        int[] anchor_positions;
        Node node,neighbor,seq_node;
        StringBuilder nds=new StringBuilder();
        StringBuilder pos=new StringBuilder();
        StringBuilder sds=new StringBuilder();
        String[] ids_list,posis_list;
        String lable, prp;
        boolean found;
        for(genome=1;genome<=genomes.num_genomes;++genome) 
        {
            for(sequence=1;sequence<=genomes.num_sequences[genome];++sequence) 
            {
                //System.out.println("Indexing sequence "+genome+"_"+sequence+"...");
                try(Transaction tx = graphDb.beginTx()){
                seq_len=genomes.sequence_length[genome][sequence]-1;
                seq_node=node=graphDb.findNode(sequence_label, "number", genome+"_"+sequence);
                prp=genome+"_"+sequence;
                node_len=K-1;
                found=true;
                count=0;
                for(loc=0;loc+node_len-K+1<=seq_len && found;)
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
                    if(loc%(seq_len/100+1)==0) 
                        System.out.print((long)loc*100/seq_len+1+"%    \r");
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
        //System.out.println("\n(Running "+cmd+")");
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
        int sq_num1=0,sq_num2=0,ed_num1,ed_num2,K1,K2,ng1,ng2,len1,len2; 
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
                K1=(int)db_node1.getProperty("k_mer_size");
                K2=(int)db_node2.getProperty("k_mer_size");
                sq_num1=(int)db_node1.getProperty("num_nodes");
                sq_num2=(int)db_node2.getProperty("num_nodes");
                ed_num1=(int)db_node1.getProperty("num_edges");
                ed_num2=(int)db_node2.getProperty("num_edges");
                ng1=(int)db_node1.getProperty("num_genomes");
                ng2=(int)db_node2.getProperty("num_genomes");
                knum1=(long)db_node1.getProperty("num_k_mers");
                knum2=(long)db_node2.getProperty("num_k_mers");
                if(K1!=K2 || sq_num1!=sq_num2 || sq_num1!=sq_num2 || ed_num1!=ed_num2 || ng1!=ng2 || knum1!=knum2)
                    return false;
                try{
                    index1=new Index(path1,null,0);
                    index2=new Index(path2,null,0);
                }
                catch(IOException ioe){
                    System.out.println("Failed to open a file!");  
                    System.exit(1);
                }  
                System.out.println("comparing pangenomes...");    
                K=K1;
                dim=(int)Math.ceil(K/4.0);
                Node n1,n2;
                byte[] s1=new byte[21],seq1;
                byte[] s2=new byte[21],seq2;
                kmer k_mer1=new kmer(dim);
                kmer k_mer2=new kmer(dim);
                pointer p1=new pointer();
                pointer p2=new pointer();
                for(k=0;k<knum1;++k)
                {
                    k_mer1=index1.get_kmer(k);
                    k_mer2=index2.get_kmer(k);
                    if(k_mer1.compare(k_mer2,dim)!=0)
                        return false;
                    index1.get_pointer(s1, p1, k);
                    index2.get_pointer(s2, p2, k);
                    if(p1.next_index==-1L && p2.next_index==-1L)
                    {
                        n1=g1.getNodeById(p1.node_id);
                        len1=(int)n1.getProperty("length");
                        n2=g2.getNodeById(p2.node_id);
                        len2=(int)n2.getProperty("length");
                        if(len1!=len2 || (!genomes.compare((int[])n1.getProperty("address"),(int[])n2.getProperty("address"),len1,true) && 
                                          !genomes.compare((int[])n1.getProperty("address"),(int[])n2.getProperty("address"),len1,false)))
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
