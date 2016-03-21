/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pantools;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.io.UnsupportedEncodingException;
import java.lang.management.MemoryUsage;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean; 
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Map;
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
    /*
    Some shared attributes
    */    

    /**
     *
     */
    public static String PATH;

    /**
     *
     */
    public static String GRAPH_DATABASE_PATH="/graph.db/";

    /**
     *
     */
    public static String INDEX_DATABASE_PATH="/index.db/";

    /**
     *
     */
    public static String GENOME_DATABASE_PATH="/genome.db/";

    /**
     *
     */
    public static GraphDatabaseService graphDb;

    /**
     *
     */
    public static index_database indexDb;

    /**
     *
     */
    public static genome_database genomeDb;

    /**
     *
     */
    public final int trsc_limit=1000;    //   The number of transactions to be committed in batch

    /**
     *
     */
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
    private static int seq_nodes;
    private static int gene_nodes;
    private static int num_edges;
    private static int num_bases;


    private kmer fwd_kmer,rev_kmer,k_mer;
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
    private String[][] fwd_locations;
    private String[][] rev_locations;
    private StringBuilder[][] occurences;
    private long[][] sum_lengths;
    private int[][] num_nodes;    
    private int fwd_code,rev_code;
    //private boolean canonical;//,annotate_canonical;
    private boolean finish=false;
    
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
            print_help_comment();
            System.exit(1);
        } 
        Pantools program=new Pantools();
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
                program.annotate(args[2]);
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
                program.retrieve_genes(args[2]);
                break;
            case "retrieve_regions":
                PATH=args[1];
                program.retrieve_regions(args[2]);
                break;
            case "cnv":
                PATH=args[1];
                program.copy_number_variation(args[2]);
                break;
            default:
                print_help_comment();
                System.exit(1);            
        }
        System.out.println("-----------------------------------------------------------------------");        
    }
    /*
    This function builds a pan-genome database out of a set of input genomes.
    fasta_names is a text file containing paths to FASTA files; each in a new line
    */  
    private void build(String file) throws IOException
    {
        int i,j;
        File theDir = new File(PATH);
        if (theDir.exists()) 
            FileUtils.deleteRecursively( new File( PATH+GRAPH_DATABASE_PATH ) ); // deletes old files
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
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+GRAPH_DATABASE_PATH)
                .setConfig( "keep_logical_logs","100M size").newGraphDatabase();
        registerShutdownHook( graphDb );            
        startTime = System.currentTimeMillis();
        seq_nodes=0;
        num_edges=0;
        gene_nodes=0;
        try(Transaction tx = graphDb.beginTx()){
            db_node=graphDb.createNode(pangenome_label);
            db_node.setProperty("k_mer_size", K);
        tx.success();}
        genomeDb=new genome_database(PATH+GENOME_DATABASE_PATH,file,false); // adding = false
        indexDb=new index_database(PATH+INDEX_DATABASE_PATH,file,K,genomeDb);
        fwd_locations=new String[genomeDb.num_genomes+1][];
        rev_locations=new String[genomeDb.num_genomes+1][];
        for(i=1;i<=genomeDb.num_genomes;++i)
        {
            fwd_locations[i]=new String[genomeDb.num_sequences[i]+1];
            rev_locations[i]=new String[genomeDb.num_sequences[i]+1];
            for(j=1;j<=genomeDb.num_sequences[i];++j)
            {
                fwd_locations[i][j]="F"+i+"_"+j;
                rev_locations[i][j]="R"+i+"_"+j;
            }
        }
        construct_pangenome();
        System.out.println("Number of kmers:   "+indexDb.length());
        System.out.println("Number of nodes:   "+seq_nodes);
        System.out.println("Number of edges:   "+num_edges);
        System.out.println("Number of bases:   "+num_bases);
            try(Transaction tx = graphDb.beginTx())
            {
                db_node.setProperty("k_mer_size", K);
                db_node.setProperty("num_k_mers", indexDb.length());
                db_node.setProperty("num_nodes",seq_nodes);
                db_node.setProperty("num_edges",num_edges);
                db_node.setProperty("num_genomes",genomeDb.num_genomes);
                db_node.setProperty("num_genes",gene_nodes);
                db_node.setProperty("num_bases",num_bases);
                tx.success();
            }
        graphDb.shutdown();  
        System.out.println("Total time : "+(System.currentTimeMillis()-startTime)/1000+"."+(System.currentTimeMillis()-startTime)%1000+" seconds");
        print_peak_memory();
        for(i=0;Files.exists(Paths.get(PATH+GRAPH_DATABASE_PATH+"/neostore.transaction.db."+i));++i)
            Files.delete(Paths.get(PATH+GRAPH_DATABASE_PATH+"/neostore.transaction.db."+i));
        System.out.println("graph.db size: "+getFolderSize(new File(PATH+GRAPH_DATABASE_PATH))+" MB");
        System.out.println("index.db size: "+getFolderSize(new File(PATH+INDEX_DATABASE_PATH))+" MB");
        System.out.println("genome.db size: "+getFolderSize(new File(PATH+GENOME_DATABASE_PATH))+" MB");
    }
    /*
    This function adds new genomes to an available pan-genome database.
    fasta_names is a text file containing paths to FASTA files; each in a new line
    */      
    private void add(String file) throws IOException
    {
        int i,j;
        if (new File(PATH+GRAPH_DATABASE_PATH).exists()) 
        {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+GRAPH_DATABASE_PATH)
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
                seq_nodes=(int)db_node.getProperty("num_nodes");
                num_edges=(int)db_node.getProperty("num_edges");
                gene_nodes=(int)db_node.getProperty("num_genes");
                genomeDb=new genome_database(PATH+GENOME_DATABASE_PATH,file,true); // adding = true
                indexDb=new index_database(PATH+INDEX_DATABASE_PATH,file, genomeDb,graphDb);
            tx.success();} 
            fwd_locations=new String[genomeDb.num_genomes+1][];
            rev_locations=new String[genomeDb.num_genomes+1][];
            for(i=1;i<=genomeDb.num_genomes;++i)
            {
                fwd_locations[i]=new String[genomeDb.num_sequences[i]+1];
                rev_locations[i]=new String[genomeDb.num_sequences[i]+1];
                for(j=1;j<=genomeDb.num_sequences[i];++j)
                {
                    fwd_locations[i][j]="F"+i+"_"+j;
                    rev_locations[i][j]="R"+i+"_"+j;
                }
            }
            construct_pangenome();
            System.out.println("Number of kmers:   "+indexDb.length());
            System.out.println("Number of nodes:   "+seq_nodes);
            System.out.println("Number of edges:   "+num_edges);
            System.out.println("Number of bases:   "+num_bases);
            try(Transaction tx = graphDb.beginTx())
            {
                db_node.setProperty("k_mer_size", K);
                db_node.setProperty("num_k_mers", indexDb.length());
                db_node.setProperty("num_nodes",seq_nodes);
                db_node.setProperty("num_edges",num_edges);
                db_node.setProperty("num_genomes",genomeDb.num_genomes);
                db_node.setProperty("num_genes",gene_nodes);
                db_node.setProperty("num_bases",num_bases);
                tx.success();
            }
            graphDb.shutdown(); 
            System.out.println("Total time : "+(System.currentTimeMillis()-startTime)/1000+"."+(System.currentTimeMillis()-startTime)%1000+" seconds");
            print_peak_memory();
        for(i=0;Files.exists(Paths.get(PATH+GRAPH_DATABASE_PATH+"/neostore.transaction.db."+i));++i)
            Files.delete(Paths.get(PATH+GRAPH_DATABASE_PATH+"/neostore.transaction.db."+i));
        System.out.println("graph.db size: "+getFolderSize(new File(PATH+GRAPH_DATABASE_PATH))+" MB");
        System.out.println("index.db size: "+getFolderSize(new File(PATH+INDEX_DATABASE_PATH))+" MB");
        System.out.println("genome.db size: "+getFolderSize(new File(PATH+GENOME_DATABASE_PATH))+" MB");
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
    private void annotate(String gff_names) throws IOException
    {
        int i;
        int begin,end,s,line_len,g_number;
        String seq;
        Node gene_node,genome_node,start_node,stop_node;
        pointer start_ptr,stop_ptr;
        Relationship rel;
        String[] fields;
        String strand,gff_name,line;
        if(new File(PATH+GRAPH_DATABASE_PATH).exists()) 
        {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+GRAPH_DATABASE_PATH)
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
                seq_nodes=(int)db_node.getProperty("num_nodes");
                num_edges=(int)db_node.getProperty("num_edges");
                gene_nodes=(int)db_node.getProperty("num_genes");
                startTime = System.currentTimeMillis();
                genomeDb=new genome_database(PATH+GENOME_DATABASE_PATH);
                //indexDb=new index_database(PATH+INDEX_DATABASE_PATH,null,0,genomeDb);
                try(BufferedReader gff_files = new BufferedReader(new FileReader(gff_names))){
                 for(g_number=1;gff_files.ready();++g_number)
                    {
                        gff_name=gff_files.readLine();
                        if(gff_name.equals(""))
                            continue;
                        genome_node=graphDb.findNode(genome_label, "number", g_number);
                        if(genome_node.hasProperty("annotated"))
                        {
                            System.out.println("genome number "+g_number+" already annotated.");
                            continue;
                        }
                        //g=(int)genome_node.getProperty("number");
                        BufferedReader in = new BufferedReader(new FileReader(gff_name));
                        System.out.println("Annotating genome "+g_number);
                        while(in.ready())
                        {
                            line=in.readLine();
                            if(line.equals(""))
                                continue;
                            fields=line.split("\\t");
                            line_len=fields.length;
                            for(i=0;i<line_len && !fields[i].equals("gene");++i)
                                {}
                            if(i<line_len && fields[i].equals("gene"))
                            {
                                seq=fields[0];
                                s=find_sequence(seq, g_number);
                                if(s>0) // if sequence found
                                {
                                    begin=Integer.parseInt(fields[i+1])-1;
                                    end=Integer.parseInt(fields[i+2])-1;
                                    strand=fields[i+4];
                                    start_ptr=locate(g_number,s,begin);
                                    stop_ptr =locate(g_number,s,end);
                                    //if(start_ptr!=null && stop_ptr!=null)
                                        if( end-begin+1>=K) // genes should not be shorter than K
                                        {
                                            ++gene_nodes;
                                            gene_node=graphDb.createNode(gene_label);
                                            gene_node.setProperty("genome_number", g_number);
                                            gene_node.setProperty("sequence_number", s);
                                            gene_node.setProperty("begin", begin+1);
                                            gene_node.setProperty("end", end+1);
                                            gene_node.setProperty("length", end-begin+1);
                                            gene_node.setProperty("strand", strand);
                                            gene_node.setProperty("title",line);
                                            start_node=graphDb.getNodeById(start_ptr.node_id);
                                            start_node.setProperty("gene_starts", true);
                                            rel=gene_node.createRelationshipTo(start_node, RelTypes.begin);
                                            rel.setProperty("forward", start_ptr.canonical);//start_canonical ?start_point.canonical:1-start_point.canonical);
                                            rel.setProperty("pos", start_ptr.position);
                                            stop_node=graphDb.getNodeById(stop_ptr.node_id);
                                            stop_node.setProperty("gene_stops", true);
                                            rel=gene_node.createRelationshipTo(stop_node, RelTypes.end);
                                            rel.setProperty("forward", stop_ptr.canonical);//start_canonical ?start_point.canonical:1-start_point.canonical);
                                            rel.setProperty("pos", stop_ptr.position);
                                        }
                                        else
                                            System.out.println("feature is less than "+K+" nt long!");
                                }
                                else 
                                {
                                    //System.out.println(seq +" missed in genome "+g); // usually organal genes
                                }
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
            for(i=0;Files.exists(Paths.get(PATH+"/database.db/neostore.transaction.db."+i));++i)
                Files.delete(Paths.get(PATH+"/database.db/neostore.transaction.db."+i));

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
        if (new File(PATH+GRAPH_DATABASE_PATH).exists()) 
        {
            Node group_node=null,gene;
            String line;
            String name=null;
            ResourceIterator<Node> genes;
            int g=0;
            System.out.println("Grouping genes...");
            try(BufferedReader in = new BufferedReader(new FileReader(group_file))){
                graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+GRAPH_DATABASE_PATH)
                        .setConfig( "keep_logical_logs","100M size").newGraphDatabase();  
                registerShutdownHook( graphDb );
                startTime = System.currentTimeMillis();
                try(Transaction tx = graphDb.beginTx())
                {
                    while(in.ready())
                    {
                        line=in.readLine();
                        if(line.equals(""))
                            continue;
                        if(line.startsWith(">"))
                        {
                            name=line.substring(1);
                            if(group_node!=null)
                            {
                                group_node.setProperty("num_genes", g);
                                System.out.println(g+" genes grouped in "+name);
                            }
                            g=0;
                            group_node=graphDb.createNode(group_lable);
                            group_node.setProperty("group_name", name);
                        }
                        else
                        {
                            genes=graphDb.findNodes(gene_label,"title",line);
                            while(genes.hasNext())
                            {
                                gene=genes.next();
                                group_node.createRelationshipTo(gene,RelTypes.has);
                                ++g;
                            }
                            genes.close();
                        }
                    }
                    group_node.setProperty("num_genes", g);
                    System.out.println(g+" genes grouped in "+name);
                    tx.success();
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
            System.out.println("No database found in "+PATH+GRAPH_DATABASE_PATH); 
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
    private void retrieve_genes(String file) throws IOException
    {
        if (new File(PATH+GRAPH_DATABASE_PATH).exists()) 
        {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+GRAPH_DATABASE_PATH)
                    .setConfig( "keep_logical_logs","100M size").newGraphDatabase(); 
            registerShutdownHook( graphDb );
            startTime = System.currentTimeMillis();
            try(Transaction tx = graphDb.beginTx()){
                K=(int)graphDb.findNodes(pangenome_label).next().getProperty("k_mer_size");
            tx.success();}
            BufferedReader in;
            ResourceIterator<Node> gene_nodes;
            Relationship rstart,rstop;
            Node start,stop,gene=null;
            String prp,record;
            String line;
            boolean strand;
            int i,j,begin,end,num_genes=0,genome,sequence;
            StringBuilder gene_seq;
            try{
                in = new BufferedReader(new FileReader(file));
                while (in.ready()) 
                {
                    line=in.readLine();
                    if(line.equals(""))
                        continue;
                    num_genes++;
                }
                in.close();
            }catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            } 
            String[] records=new String[num_genes];
            try(Transaction tx = graphDb.beginTx()){
            try{ 
                in = new BufferedReader(new FileReader(file));
                BufferedWriter out = new BufferedWriter(new FileWriter(file+".fasta"));
                for(i=0;in.ready();)
                {
                    line=in.readLine();
                    if(line.equals(""))
                        continue;                   
                    records[i]=line;
                    ++i;
                }
                Arrays.sort(records);
                genomeDb=new genome_database(PATH+GENOME_DATABASE_PATH);
                for(i=1,j=0,gene_nodes=graphDb.findNodes(gene_label);gene_nodes.hasNext();++i)
                {
                    gene=gene_nodes.next();
                    record=(String)gene.getProperty("title");
                    if(record!=null && Arrays.binarySearch(records, record) >=0)
                    {
                        rstart=gene.getSingleRelationship(RelTypes.begin, Direction.OUTGOING);
                        //rstop=gene.getSingleRelationship(RelTypes.end, Direction.OUTGOING);
                        start=rstart.getEndNode();
                        //stop=rstop.getEndNode();
                        genome=(int)gene.getProperty("genome_number");
                        sequence=(int)gene.getProperty("sequence_number");
                        prp=genome+"_"+sequence;
                        begin=(int)gene.getProperty("begin");
                        end=(int)gene.getProperty("end");
                        strand=gene.getProperty("strand").toString().equals("+");
                        gene_seq=extract_sequence(start,(boolean)rstart.getProperty("forward"),(int)rstart.getProperty("pos"),prp, begin-1, end-1);//
                        if(gene_seq.length() == end-begin+1)
                        {
                            ++j;
                            out.write(">"+record+"\n");
                            if(strand)
                                write_fasta(out,gene_seq.toString(),70);
                            else
                                write_fasta(out,reverse_complement(gene_seq.toString()),70);
                        }
                        else
                        {
                            System.out.println(record);//"Failed to assemble "+i+"\'th gene.");//+(strand?"+":"-"));
                            System.out.println(gene_seq.length()+" != "+(end-begin+1));
                        }
                        gene_seq.setLength(0);
                    }
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
    This function retrieves genes' sequence from the graph
    Input:
        file: A text file containing ID of genes whose sequence should be retrieved, One ID per line
    Output:
        A FASTA file of gene sequences.
    */ 
    private void copy_number_variation(String file) throws IOException
    {
        int i,j;
        if (new File(PATH+GRAPH_DATABASE_PATH).exists()) 
        {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+GRAPH_DATABASE_PATH)
                    .setConfig( "keep_logical_logs","100M size").newGraphDatabase(); 
            registerShutdownHook( graphDb );
            startTime = System.currentTimeMillis();
            try(Transaction tx = graphDb.beginTx()){
                K=(int)graphDb.findNodes(pangenome_label).next().getProperty("k_mer_size");
            tx.success();}
            genomeDb=new genome_database(PATH+GENOME_DATABASE_PATH);
            fwd_locations=new String[genomeDb.num_genomes+1][];
            rev_locations=new String[genomeDb.num_genomes+1][];
            occurences=new StringBuilder[genomeDb.num_genomes+1][];
            sum_lengths=new long[genomeDb.num_genomes+1][];
            num_nodes=new int[genomeDb.num_genomes+1][]; 
            for(i=1;i<=genomeDb.num_genomes;++i)
            {
                sum_lengths[i]=new long[genomeDb.num_sequences[i]+1];
                num_nodes[i]=new int[genomeDb.num_sequences[i]+1];
                fwd_locations[i]=new String[genomeDb.num_sequences[i]+1];
                rev_locations[i]=new String[genomeDb.num_sequences[i]+1];
                occurences[i]=new StringBuilder[genomeDb.num_sequences[i]+1];
                for(j=1;j<=genomeDb.num_sequences[i];++j)
                {
                    fwd_locations[i][j]="F"+i+"_"+j;
                    rev_locations[i][j]="R"+i+"_"+j;
                    occurences[i][j]=new StringBuilder();
                }
            }
            BufferedReader in;
            ResourceIterator<Node> gene_nodes;
            Relationship rstart,rstop;
            Node start,stop,gene=null;
            String prp,record;
            String line;
            int begin,end,num_genes=0,genome,sequence,g,s;
            try{
                in = new BufferedReader(new FileReader(file));
                while (in.ready()) 
                {
                    line=in.readLine();
                    if(line.equals(""))
                        continue;
                    num_genes++;
                }
                in.close();
            }catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            } 
            String[] records=new String[num_genes];
            try(Transaction tx = graphDb.beginTx()){
            try{ 
                in = new BufferedReader(new FileReader(file));
                BufferedWriter out = new BufferedWriter(new FileWriter(file+".fasta"));
                for(i=0;in.ready();)
                {
                    line=in.readLine();
                    if(line.equals(""))
                        continue;                   
                    records[i]=line;
                    ++i;
                }
                Arrays.sort(records);
                genomeDb=new genome_database(PATH+GENOME_DATABASE_PATH);
                for(i=1,j=0,gene_nodes=graphDb.findNodes(gene_label);gene_nodes.hasNext();++i)
                {
                    gene=gene_nodes.next();
                    record=(String)gene.getProperty("title");
                    if(record!=null && Arrays.binarySearch(records, record) >=0)
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
                        out.write(">"+record+"length: "+(end-begin+1)+"\n");
                        out.write(extract_cnv(start,stop,(int)rstart.getProperty("pos"),prp, begin-1, end-1).toString());
                    }
                    else
                    {}//System.out.println(gene_record+" not annotated!");
                    if(i%(num_genes/100+1)==0) 
                        System.out.print((long)i*100/num_genes+1+"%\r");
                }//for i
               // System.out.println(j+" genes retrieved successfully.");
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
    private void retrieve_regions(String file) throws IOException
    {
        if (new File(PATH+GRAPH_DATABASE_PATH).exists()) 
        {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+GRAPH_DATABASE_PATH)
                    .setConfig( "keep_logical_logs","100M size").newGraphDatabase(); 
            registerShutdownHook( graphDb );
            try(Transaction tx = graphDb.beginTx()){
                K=(int)graphDb.findNodes(pangenome_label).next().getProperty("k_mer_size");
            tx.success();}
            String[] fields;
            String line;
            StringBuilder seq;
            Node start_node;
            pointer start_ptr;
            int c,g,s,from, to, num_regions=0;
            try{
                BufferedReader in = new BufferedReader(new FileReader(file));
                while (in.ready()) 
                {
                    line=in.readLine();
                    if(line.equals(""))
                        continue;
                    num_regions++;
                }
                in.close();
            }catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            } 
            //num_regions=Integer.parseInt(Pantools.executeCommand("wc -l "+file).trim().split("\\s")[0]);
            startTime = System.currentTimeMillis();
            try(Transaction tx = graphDb.beginTx()){
            try (BufferedReader in = new BufferedReader(new FileReader(file))) {
                BufferedWriter out = new BufferedWriter(new FileWriter(file+".fasta"));
                //genomeDb=new genome_database(PATH+GENOME_DATABASE_PATH);
                for(c=0;in.ready();)
                {
                    line=in.readLine();
                    if(line.equals(""))
                        continue;
                    fields=line.split("\\s");
                    g= Integer.parseInt(fields[0]);
                    s= Integer.parseInt(fields[1]);
                    from= Integer.parseInt(fields[2]);
                    to= Integer.parseInt(fields[3]);
                    start_ptr=locate(g,s,from);
                    //stop_ptr=locate(g,s,to);
                    start_node=graphDb.getNodeById(start_ptr.node_id);
                    //stop_node =graphDb.getNodeById(stop_ptr.node_id);
                    seq=extract_sequence(start_node, start_ptr.canonical,start_ptr.position,g+"_"+s, from, to);
                    /*if(start_ptr.canonical)
                        seq=extract_sequence(start_node,stop_node, true,from-start_ptr.position,g+"_"+s, from, to);
                    else
                        seq=extract_sequence(start_node,stop_node, true,(int)start_node.getProperty("length")-(from-start_ptr.position)-K,g+"_"+s,from, to);
                   */
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
    private StringBuilder extract_sequence(Node node, boolean forward,int position, String coordinates, int begin, int end) throws UnsupportedEncodingException
    {
        boolean found;
        Node neighbor;
        String lable,name;
        String[] prp={"F"+coordinates,"R"+coordinates};
        int loc,node_len,neighbor_len,s_side,d_side,seq_len=end-begin+1;
        StringBuilder seq=new StringBuilder();
        loc=begin;
        node_len=(int)node.getProperty("length");
        if(forward) 
        {
            if(position+seq_len-1<=node_len-1)
                loc+=append_fwd(seq,(String)node.getProperty("sequence"),position,position+seq_len-1);
            else
                loc+=append_fwd(seq,(String)node.getProperty("sequence"),position,node_len-1);
        }
        else
        {
            if(position-(seq_len-1)>=0)
                loc+=append_rev(seq,(String)node.getProperty("sequence"),position-(seq_len-1),position);
            else
                loc+=append_rev(seq,(String)node.getProperty("sequence"),0,position);
        }
        found=true;
        while(seq.length()<seq_len && found) //loc<end)
            {
                found=false;
                for(Relationship r:node.getRelationships(Direction.OUTGOING))
                {
                    neighbor=r.getEndNode();
                    name=r.getType().name();
                    d_side=(name.charAt(2)-48)%2;
                    lable=(d_side==0?prp[0]:prp[1]);
                    if(neighbor.hasProperty(lable) )
                    {
                        neighbor_len=(int)neighbor.getProperty("length");
                        if(b_search((int[])neighbor.getProperty(lable), loc-K+1)>=0)
                        {
                            found=true;
                            if(d_side==0)
                                if(seq.length()+neighbor_len-K+1>seq_len)
                                    loc+=append_fwd(seq,(String)neighbor.getProperty("sequence"),K-1,seq_len-seq.length()+K-2);
                                else
                                    loc+=append_fwd(seq,(String)neighbor.getProperty("sequence"),K-1,neighbor_len-1);
                            else
                                if(seq.length()+neighbor_len-K+1>seq_len)
                                    loc+=append_rev(seq,(String)neighbor.getProperty("sequence"),neighbor_len-K-(seq_len-seq.length())+1,neighbor_len-K);
                                else
                                    loc+=append_rev(seq,(String)neighbor.getProperty("sequence"),0,neighbor_len-K);
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
                    lable=(s_side==0?prp[1]:prp[0]);
                    if(neighbor.hasProperty(lable) )
                    {
                        neighbor_len=(int)neighbor.getProperty("length");
                        if(b_search((int[])neighbor.getProperty(lable), loc-K+1)>=0)
                        {
                            found=true;
                            if(s_side==1)
                                if(seq.length()+neighbor_len-K+1>seq_len)
                                    loc+=append_fwd(seq,(String)neighbor.getProperty("sequence"),K-1,seq_len-seq.length()+K-2);
                                else
                                    loc+=append_fwd(seq,(String)neighbor.getProperty("sequence"),K-1,neighbor_len-1);
                            else
                                if(seq.length()+neighbor_len-K+1>seq_len)
                                    loc+=append_rev(seq,(String)neighbor.getProperty("sequence"),neighbor_len-K-(seq_len-seq.length())+1,neighbor_len-K);
                                else
                                    loc+=append_rev(seq,(String)neighbor.getProperty("sequence"),0,neighbor_len-K);
                            node=neighbor;
                            node_len=(int)node.getProperty("length");
                            break;
                        }
                    }
                }
            }//for loc
        return seq;
    }
    /*
    This function returns the sequence starting at (start_node, start_side, start_pos) and stopping at stop_node which belongs to 
    genoeme and sequence specified by "coordinates" and is of length end-begin+1
    */
    private StringBuilder extract_cnv(Node start_node, Node stop_node, int start_pos, String coordinates, int begin, int end) throws UnsupportedEncodingException
    {
        boolean found;
        Node node=start_node,neighbor;
        String lable,name;
        int i,j,p,len=0,loc,node_len,neighbor_len,s_side,d_side,seq_len=end-begin+1;
        for(i=1;i<=genomeDb.num_genomes;++i)
            for(j=1;j<=genomeDb.num_sequences[i];++j)
            {
                sum_lengths[i][j]=0;
                num_nodes[i][j]=0;
                occurences[i][j].setLength(0);
            }
        StringBuilder result=new StringBuilder();
        if(start_node.equals(stop_node) && start_pos+end-begin<(int)start_node.getProperty("length")) // if gene is inside the node not a circle which starts and ends in the node
        {
            p=end-begin+1;
            append_statistics(start_node,p);
            len+=p;
        }
        else
        {
            node_len=(int)node.getProperty("length");
            if(start_pos+end-begin<(int)start_node.getProperty("length"))
            {
                p=end-begin+1;
                append_statistics(start_node,p);
                len+=p;
            }
            else
            {
                p=(int)node.getProperty("length")-start_pos;
                append_statistics(start_node,p);
                len+=p;
            }
            loc=begin-start_pos;
            found=true;
            while(len<seq_len && found) 
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
                            if(len+neighbor_len-K+1>seq_len)
                            {
                                append_statistics(neighbor,seq_len-len);
                                len=seq_len;
                            }
                            else
                            {
                                append_statistics(neighbor,neighbor_len-K+1);
                                len+=neighbor_len-K+1;
                            }
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
                            if(len+neighbor_len-K+1>seq_len)
                            {
                                append_statistics(neighbor,seq_len-len);
                                len=seq_len;
                            }
                            else
                            {
                                append_statistics(neighbor,neighbor_len-K+1);
                                len+=neighbor_len-K+1;
                            }
                            node=neighbor;
                            node_len=(int)node.getProperty("length");
                            break;
                        }
                    }
                }
            }//for loc
        }
        for(i=1;i<=genomeDb.num_genomes;++i)
        {
            for(j=1;j<=genomeDb.num_sequences[i];++j)
                if(num_nodes[i][j]>0)
                    result.append("g"+i+"s"+j+" len: "+sum_lengths[i][j]+" num_nodes: "+num_nodes[i][j]+"\n["+occurences[i][j]+ "]\n");///(double)(end-begin+1)
        }        
        return result;
    }
    /*
    This function appends s[from..to] to "seq". 
    */
    private int append_fwd(StringBuilder seq, String s, int from, int to)
    {
        seq.append(s.substring(from,to+1));//genomeDb.get_sequence(addrs[0],addrs[1],addrs[2]+from, to-from+1));
        return to-from+1;
    }
    /*
    This function appends reverse_complement of s[from..to] to "seq". 
    */
    private int append_rev(StringBuilder seq, String s, int from, int to)
    {
        seq.append(reverse_complement(s.substring(from,to+1)));//genomeDb.get_sequence(addrs[0],addrs[1],addrs[2]+from, to-from+1)));
        return to-from+1;
    }
    /*
    This function appends reverse_complement of s[from..to] to "seq". 
    */
    private void append_statistics(Node node, int len)
    {
        int i,j,k,n;
        int[] occ;
        for(i=1;i<=genomeDb.num_genomes;++i)
            for(j=1;j<=genomeDb.num_sequences[i];++j)
            {
                if(node.hasProperty(fwd_locations[i][j]))
                {
                    occ=(int[])node.getProperty(fwd_locations[i][j]);
                    n=occ.length;
                    sum_lengths[i][j]+=len*n;
                    num_nodes[i][j]+=n;
                    for(k=0;k<n;++k)
                        occurences[i][j].append(occ[k]+" ");
                }
                if(node.hasProperty(rev_locations[i][j]))
                {
                    occ=(int[])node.getProperty(rev_locations[i][j]);
                    n=occ.length;
                    sum_lengths[i][j]+=len*n;
                    num_nodes[i][j]+=n;
                    for(k=0;k<n;++k)
                        occurences[i][j].append(occ[k]+" ");
                }
            }
    }
    /*
    This function returns the "location" (node, side, position) where a region located at "position" in "genoeme" and "sequence"
    specified by "coordinates" occurs. It starts the search from "start" node.
    */
    private pointer locate(int g, int s, int position)
    {
        int i,loc,low,high,mid,node_len;
        boolean forward;
        Node node,neighbor;
        String lable;
        boolean found;
        long[] anchor_nodes;
        int[] anchor_positions;
        String anchor_sides;
        String[] prp={"F"+g+"_"+s,"R"+g+"_"+s};
        Node seq_node=graphDb.findNode(sequence_label, "number", g+"_"+s);
        anchor_nodes=(long [])seq_node.getProperty("anchor_nodes");
        anchor_positions=(int [])seq_node.getProperty("anchor_positions");
        anchor_sides=(String)seq_node.getProperty("anchor_sides");
        for(low=0,high=anchor_sides.length()-1,mid=(low+high)/2;low<high;mid=(low+high)/2)
        {
            if(position<anchor_positions[mid])
                high=mid-1;
            else if(position>anchor_positions[mid])
                low=mid+1;
            else
                break;
        }
        if(position<anchor_positions[mid])
            --mid;
        try(Transaction tx = graphDb.beginTx()){
        node=graphDb.getNodeById(anchor_nodes[mid]);
        loc=anchor_positions[mid];
        forward=anchor_sides.charAt(mid)=='0';
        node_len=(int)node.getProperty("length");
        found=true;
        while(loc+node_len-K+1<=position && found) // while have not reached to target
        {
            found=false;
            for(Relationship r:node.getRelationships(Direction.OUTGOING)) // check all the outgoing edges
            {
                neighbor=r.getEndNode();
                forward=(r.getType().name().charAt(2)-48)%2==0;
                lable=forward?prp[0]:prp[1];
                if(neighbor.hasProperty(lable))
                {
                    if(b_search((int[])neighbor.getProperty(lable), loc+node_len-K+1)>=0)
                    {
                        found=true;
                        loc+=node_len-K+1;
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
                forward=(r.getType().name().charAt(2)-48)/2==0;
                lable=forward?prp[1]:prp[0];
                if(neighbor.hasProperty(lable) )
                {
                    if(b_search((int[])neighbor.getProperty(lable), loc+node_len-K+1)>=0)
                    {
                        found=true;
                        loc+=node_len-K+1;
                        node=neighbor;
                        node_len=(int)node.getProperty("length");
                        break;
                    }
                }
            }
        }
        tx.success();}
        return new pointer(node.getId(),forward,forward?position-loc:node_len-1-(position-loc),-1L);
    }
    /*
    This function returns the number of the sequence from genome "g" whose name contains "prefix".
    */    
    int find_sequence(String prefix, int g)
    {
        boolean found=false;
        int s;
        for(found=false,s=1;!found && s<=genomeDb.num_sequences[g];++s)
            if(genomeDb.sequence_titles[g][s].contains(prefix))
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
                case 'A': 
                    rv.append('T');
                    break;
                case 'C': 
                    rv.append('G');
                    break;
                case 'G': 
                    rv.append('C');
                    break;
                case 'T': 
                    rv.append('A');
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
            s_base=s_side==0?genomeDb.get(s_gen,s_seq,s_loc+s_len-K):(3-genomeDb.get(s_gen,s_seq,s_loc+K-1));
            if(s_base<0 || s_base>3)
                s_base=4;
        }
        d_len= (int)des.getProperty("length") ;
        d_add=(int[])des.getProperty("address");
        d_gen=d_add[0];
        d_seq=d_add[1];
        d_loc=d_add[2];
        d_base=d_side==0?genomeDb.get(d_gen,d_seq,d_loc+K-1):(3-genomeDb.get(d_gen,d_seq,d_loc+d_len-K));
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
            r_rt=RelTypes.values()[(3-genomeDb.binary[rt_name.charAt(1)])*20+(3-genomeDb.binary[rt_name.charAt(0)])*4+(rt_name.charAt(2)=='0'?3:(rt_name.charAt(2)-48)%3)];
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
        num_bases+=(K-1);
        int l=(int)node.getProperty("length");
        int i,j,k=0,s_id,s_side=0,d_side=0,gen,seq,loc;
        long inx,split_first_kmer,node_last_kmer=0;
        int[] address;
        Node neighbor;
        RelationshipType rel_type;
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
        if(node.hasProperty("gene_starts"))
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
        if(node.hasProperty("gene_stops"))
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
        int[] node_positions;
        int len;
        int[] split_positions;
        for(i=1;i<=genome;++i) // initial occurence properties for the split_node
            for(j=1;j<=genomeDb.num_sequences[i];++j)
            {
                prp=fwd_locations[i][j];
                if(node.hasProperty(prp) )
                {
                    node_positions=(int[])node.getProperty(prp);
                    len=node_positions.length;
                    split_positions=new int[len];
                    for(k=0;k<len;++k)
                        split_positions[k]=node_positions[k]+pos;
                    split_node.setProperty(prp,split_positions);
                }
                prp=rev_locations[i][j];
                if(node.hasProperty(prp) )
                {
                    node_positions=(int[])node.getProperty(prp);
                    len=node_positions.length;
                    split_positions=new int[len];
                    for(k=0;k<len;++k)
                    {
                        split_positions[k]=node_positions[k];
                        node_positions[k]=node_positions[k]+l-pos-K+1;
                    }
                    split_node.setProperty(prp,split_positions);
                    node.setProperty(prp,node_positions);
                }
            }
        node_last_kmer=indexDb.find(make_kmer(gen,seq,loc+pos-1));
        split_first_kmer=indexDb.get_next_index(node_last_kmer);
        indexDb.put_next_index(-1L, node_last_kmer); 
        split_node.setProperty("first_kmer",split_first_kmer);
        split_node.setProperty("last_kmer",node.getProperty("last_kmer"));
        s_id=(int)split_node.getId();
        for(i=0,inx=split_first_kmer;inx!=-1L;inx=indexDb.get_next_index(inx),++i) // update kmer coordinates
        {
            indexDb.put_node_id(s_id, inx);
            indexDb.put_position(i, inx);
        }  
        for(Relationship r:node.getRelationships(Direction.OUTGOING))
        {
            rel_type=r.getType();
            s_side=(rel_type.name().charAt(2)-48)/2;
            if(s_side==0)
            {
                d_side=(rel_type.name().charAt(2)-48)%2;
                neighbor=r.getEndNode();
                if(neighbor.equals(node))
                    neighbor=d_side==0?node:split_node;
                split_node.createRelationshipTo(neighbor,rel_type);
                r.delete();
            }                
        }
        for(Relationship r:node.getRelationships(Direction.INCOMING))
        {
            rel_type=r.getType();
            d_side=(rel_type.name().charAt(2)-48)%2;
            if(d_side==1)
            {
                s_side=(rel_type.name().charAt(2)-48)/2;
                neighbor=r.getStartNode();
                if(neighbor.equals(node))
                    neighbor=s_side==1?node:split_node;
                neighbor.createRelationshipTo(split_node,rel_type);
                r.delete();
            }                
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
        int begin,len=(int)node.getProperty("length");
        long id=node.getId(),last_kmer=(long)node.getProperty("last_kmer");
        boolean broke=false;
        while(position<seq_len)
        {
            //System.out.println("extend "+position);
            if(genomeDb.get(genome,sequence,position+1)>3)
            {
                ++position;
                begin=position-K+1;
                jump();
                create_degenerate(begin);
                curr_node=degenerate_node;
                //curr_side=0;
                break;
            }
            next_kmer();
            curr_index=indexDb.find(k_mer);
            indexDb.get_pointer(byte_pointer,pointer,curr_index);
            if(pointer.node_id == -1L)
            {
                indexDb.put_next_index(curr_index,last_kmer);
                ++len;
                //node.setProperty("length",(int)node.getProperty("length")+1);
                ++num_bases;
                pointer.node_id=id;
                pointer.canonical=fwd_kmer.canonical;
                pointer.position=len-K;//(int)node.getProperty("length")-K;
                pointer.next_index=-1L;
                indexDb.put_pointer(pointer, curr_index);
                last_kmer=curr_index;
                //node.setProperty("last_kmer",curr_index);
                if(position%(seq_len/100+1)==0) 
                    System.out.print((long)position*100/seq_len+1+"%\r");
            }
            else
            {
                broke=true;
                break;
            }
        }
        node.setProperty("length",len);
        node.setProperty("last_kmer",last_kmer);
        if(!broke && position==seq_len )
            finish=true;
    }
    /*
    This function creates a new node
    */
    void create()
    {
        //System.out.println("create "+position);
        num_bases+=K;
        ++seq_nodes;
        new_node=graphDb.createNode(node_label);
        new_node.setProperty("F"+genome+"_"+sequence,new int[]{position-K+1});
        new_node.setProperty("address",new int[] {genome,sequence,position-K+1});
        new_node.setProperty("length",K);
        new_node.setProperty("last_kmer",curr_index);
        new_node.setProperty("first_kmer",curr_index);
        pointer.node_id=new_node.getId();
        pointer.canonical=fwd_kmer.canonical;
        pointer.position=0;
        pointer.next_index=-1L;
        indexDb.put_pointer(pointer, curr_index);
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
            if(curr_node.equals(node))// && curr_side==0) // is always 0
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
        for(pos=0;pos<=l && position<=seq_len && genomeDb.get(g, s, loc+pos+K-1)==genomeDb.get(genome,sequence,position); ++pos)
        {
            ++position;
            if(position<=seq_len && genomeDb.get(genome,sequence,position)>3)
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
        curr_index=indexDb.find(k_mer);
        //if(curr_index!=-1L)
            indexDb.get_pointer(byte_pointer,pointer,curr_index);
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
            if(curr_node.equals(node))// && curr_side==0)// always 0
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
        for(pos=(int)node.getProperty("length")-K;pos>=0 && position<=seq_len && genomeDb.get(g,s,loc+pos)==genomeDb.get_complement(genome,sequence,position); --pos) 
        {
            ++position;
            if(position<=seq_len && genomeDb.get(genome,sequence,position)>3)
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
        curr_index=indexDb.find(k_mer);
        //if(curr_index!=-1L)
            indexDb.get_pointer(byte_pointer,pointer,curr_index);
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
        fwd_code=genomeDb.get(genome,sequence,position);
        rev_code=3-fwd_code;
        do{
            while(fwd_code>3 && position<seq_len)
            {
                ++position;
                fwd_code=genomeDb.get(genome,sequence,position);
                rev_code=3-fwd_code;
            }
            fwd_kmer=new kmer(K,indexDb.get_pre_len(),indexDb.get_suf_len());
            rev_kmer=new kmer(K,indexDb.get_pre_len(),indexDb.get_suf_len());
            fwd_kmer.next_up_kmer(fwd_code);
            rev_kmer.next_down_kmer(rev_code);
            for(j=0; j<K-1 && position<seq_len;++j)
            {
                ++position;
                fwd_code=genomeDb.get(genome,sequence,position);
                rev_code=3-fwd_code;
                if(fwd_code>3 )
                    break;
                fwd_kmer.next_up_kmer(fwd_code);
                rev_kmer.next_down_kmer(rev_code);
            }
        }while(fwd_code>3 && position<seq_len);
        fwd_kmer.canonical=fwd_kmer.compare(rev_kmer)==-1;
        k_mer=fwd_kmer.canonical?fwd_kmer:rev_kmer;
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
        degenerate_node.setProperty("genome",genome);
        degenerate_node.setProperty("F"+genome+"_"+sequence,new int[]{begin});
        degenerate_node.setProperty("address",new int[] {genome,sequence,begin});
        degenerate_node.setProperty("length",position-begin);
        connect(curr_node,degenerate_node,curr_side,0,0);
        curr_index=indexDb.find(k_mer);//canonical?fwd_kmer:rev_kmer);
        num_bases+=(position-begin);
        curr_side=0;
        indexDb.get_pointer(byte_pointer,pointer,curr_index);
    }
    /*
    This function produce the next "fwd_kmer", "rev_kmer" and the canonical "kmer" from the current ones.
    */ 
    void next_kmer()
    {
        ++position;
        fwd_code=genomeDb.get(genome,sequence,position);
        rev_code=3-fwd_code;
        fwd_kmer.next_up_kmer(fwd_code);
        rev_kmer.next_down_kmer(rev_code);
        if(position%(seq_len/100+1)==0) 
            System.out.print((long)position*100/seq_len+1+"%\r");
        fwd_kmer.canonical=fwd_kmer.compare(rev_kmer)==-1;
        k_mer=fwd_kmer.canonical?fwd_kmer:rev_kmer;
    }
    /*
    This function initialize the first K-mer of the "genome", "sequence" at "position". 
    It might jump over the degenerate regions creating a degenerate_node 
    */ 

    /**
     *
     * @throws IOException
     */
     
    public void initial_kmers() throws IOException
    {
        int i;
        fwd_kmer=new kmer(K,indexDb.get_pre_len(),indexDb.get_suf_len());
        rev_kmer=new kmer(K,indexDb.get_pre_len(),indexDb.get_suf_len());
        for(i=0;i<K && position<seq_len;++i)
        {
            if(genomeDb.get(genome,sequence,position+1)>3 )
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
        fwd_kmer.canonical=fwd_kmer.compare(rev_kmer)==-1;
        k_mer=fwd_kmer.canonical?fwd_kmer:rev_kmer;
    }
    /*
    This function gives the K-mer in "genome", "sequence" at "position". 
    fwd_kmer: Forward sequence of the K-mer
    rev_kmer: Reverse sequence of the K-mer
    kmer  : canonical k-mer
    canonical: a boolean determining if the K-mer is canonical  
    */ 
    kmer make_kmer(int genome, int sequence, int position)
    {
        int j,fwd_code,rev_code;
        kmer fwd_kmer=new kmer(K,indexDb.get_pre_len(),indexDb.get_suf_len());
        kmer rev_kmer=new kmer(K,indexDb.get_pre_len(),indexDb.get_suf_len());
        for(j=0;j<K;++j)
        {
            fwd_code=genomeDb.get(genome,sequence,position+j);
            rev_code=3-fwd_code;
            fwd_kmer.next_up_kmer(fwd_code);
            rev_kmer.next_down_kmer(rev_code);
        }   
        fwd_kmer.canonical=fwd_kmer.compare(rev_kmer)==-1;
        return fwd_kmer.canonical?fwd_kmer:rev_kmer;        
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
        for(genome=genomeDb.previous_num_genomes+1;genome<=genomeDb.num_genomes;++genome) 
        {
            try(Transaction tx = graphDb.beginTx()){
                genome_node=graphDb.createNode(genome_label);
                genome_node.setProperty("number", genome);
                genome_node.setProperty("num_sequences",genomeDb.num_sequences[genome]);
                db_node.createRelationshipTo(genome_node, RelTypes.has);
                System.out.println("Processing genome "+genome+" :             ");
            tx.success();}
            sequence_ids=new long[genomeDb.num_sequences[genome]+1];
            for(sequence=1;sequence<=genomeDb.num_sequences[genome];++sequence) 
            {
                System.out.println("sequence "+sequence+"/"+genomeDb.num_sequences[genome]+" of genome "+genome+"\tlength="+genomeDb.sequence_length[genome][sequence]);
                try(Transaction tx = graphDb.beginTx()){
                    sequence_node=curr_node=graphDb.createNode(sequence_label);
                    sequence_node.setProperty("number",genome+"_"+sequence);
                    sequence_node.setProperty("sequence_name",genomeDb.sequence_titles[genome][sequence]);
                    sequence_node.setProperty("sequence_length",genomeDb.sequence_length[genome][sequence]);
                    genome_node.createRelationshipTo(sequence_node, RelTypes.has);
                    sequence_ids[sequence]=curr_node.getId();
                    curr_side=0;
                    position=-1;
                    seq_len=genomeDb.sequence_length[genome][sequence]-1;
                    initial_kmers();
                tx.success();}
                curr_index=indexDb.find(k_mer);
                indexDb.get_pointer(byte_pointer,pointer,curr_index);
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
                            else if( fwd_kmer.canonical ^ pointer.canonical )// if sides don't agree
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
        index_genomes();
        include_sequences();        
    }
    /*
    This function adds list of "anchor_nodes", "anchor_sides" and "anchor_positions" to each sequence_node.
    These properties facilitate extracting genomic regions. 
    */ 
    void index_genomes() 
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
        for(genome=1;genome<=genomeDb.num_genomes;++genome) 
        {
            for(sequence=1;sequence<=genomeDb.num_sequences[genome];++sequence) 
            {
                //System.out.println("Indexing sequence "+genome+"_"+sequence+"...");
                try(Transaction tx = graphDb.beginTx()){
                seq_len=genomeDb.sequence_length[genome][sequence]-1;
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
    This function adds "sequence" property to the nodes and writes the K-mer index on disk
    */ 
    void include_sequences() throws IOException
    {
        int i,len;
        int[] address;
        byte[] sequence;
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
                    node.setProperty("sequence",genomeDb.get_sequence(address[0], address[1], address[2], len));
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
                    node.setProperty("sequence",genomeDb.get_sequence(address[0], address[1], address[2], len));
                }
            tx.success();}
        nodes.close();
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
    This function compares two graph databases located in paths "path1" and "path2" and return the boolean result
    */    
    private boolean compare_pan(String path1, String path2)
    {
        int sq_num1=0,sq_num2=0,ed_num1,ed_num2,K1,K2,ng1,ng2,len1,len2; 
        long k,knum1=0,knum2;
        index_database indexDb1=null,indexDb2=null;
        genome_database genomeDb1,genomeDb2;
        Node db_node1,db_node2;
        if (new File(path1+GRAPH_DATABASE_PATH).exists() && new File(path2+GRAPH_DATABASE_PATH).exists()) 
        {
            GraphDatabaseService g1 = new GraphDatabaseFactory().newEmbeddedDatabase(path1+GRAPH_DATABASE_PATH );
            GraphDatabaseService g2 = new GraphDatabaseFactory().newEmbeddedDatabase(path2+GRAPH_DATABASE_PATH );
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
                genomeDb1=new genome_database(path1+GENOME_DATABASE_PATH);
                genomeDb2=new genome_database(path2+GENOME_DATABASE_PATH);
                indexDb1=new index_database(path1+INDEX_DATABASE_PATH,null,0,genomeDb1);
                indexDb2=new index_database(path2+INDEX_DATABASE_PATH,null,0,genomeDb2);
                System.out.println("comparing pangenomes...");    
                K=K1;
                Node n1,n2;
                byte[] s1=new byte[21];
                byte[] s2=new byte[21];
                kmer k_mer1=new kmer(K,indexDb1.get_pre_len(),indexDb1.get_suf_len());
                kmer k_mer2=new kmer(K,indexDb2.get_pre_len(),indexDb2.get_suf_len());
                pointer p1=new pointer();
                pointer p2=new pointer();
                for(k=0;k<knum1;++k)
                {
                    k_mer1=indexDb1.get_kmer(k);
                    k_mer2=indexDb2.get_kmer(k);
                    if(k_mer1.compare(k_mer2)!=0)
                        return false;
                    indexDb1.get_pointer(s1, p1, k);
                    indexDb2.get_pointer(s2, p2, k);
                    if(p1.next_index==-1L && p2.next_index==-1L)
                    {
                        n1=g1.getNodeById(p1.node_id);
                        len1=(int)n1.getProperty("length");
                        n2=g2.getNodeById(p2.node_id);
                        len2=(int)n2.getProperty("length");
                        if(len1!=len2 || (!genomeDb1.compare(genomeDb2,(int[])n1.getProperty("address"),(int[])n2.getProperty("address"),len1,true) && 
                                          !genomeDb1.compare(genomeDb2,(int[])n1.getProperty("address"),(int[])n2.getProperty("address"),len1,false)))
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
     private static void print_help_comment()
    {
        System.out.println("Requirements:\n" +
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
        "          java -jar pantools.jar annotate database_path gff_names_file\n" +
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
        "          java -jar pantools.jar group database_path groups_file\n" +
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
        "gff_names_file  : A text file containing paths to GFF files corresponding to the stored genomes,\n"
         + " in the same order. Missing annotations are shown by an empty line.\n" +
        "\n" +
        "annotation_records : A text file containing annotation titles as they appear in gff file.\n" +
        "\n" +
        "regions_file    : A text file containing genome_number, sequence_number, begin and\n"
         + " end of a region in each line seperated by one space \n" +
        "\n" +
        "groups_file : A FASTA file with titles being names given to the groups "
        + "followed by lines containing annotation titles as they appear in gff files.");
    }

    /**
     *
     * @param dir
     * @return
     */
    public static long getFolderSize(File dir) {
    long size = 0;
    for (File file : dir.listFiles()) {
        if (file.isFile()) {
            // System.out.println(file.getName() + " " + file.length());
            size += file.length();
        } else
            size += getFolderSize(file);
    }
    return size/1048576+1;
}     
}
