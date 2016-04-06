/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pantools;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.Console;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.lang.management.MemoryUsage;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean; 
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Map;
import java.util.Map.Entry;
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
import org.neo4j.graphdb.QueryExecutionException;
import org.neo4j.graphdb.Result;


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
        
        begin,end,
        has
    } 

    public static String PATH;
    public static String GRAPH_DATABASE_PATH="/graph.db/";
    public static String INDEX_DATABASE_PATH="/index.db/";
    public static String GENOME_DATABASE_PATH="/genome.db/";
    public static GraphDatabaseService graphDb;
    public static index_database indexDb;
    public static genome_database genomeDb;
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
    private static int seq_nodes;
    private static int degenerate_nodes;
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
    private int fwd_code,rev_code;
    private boolean finish=false;
    private int[] address;
    
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
    public static void main(String[] args) {
        if(args.length < 2 || args[1].equals("--help") || args[1].equals("-h")){
            print_help_comment();
            System.exit(1);
        } 
        Pantools program=new Pantools();
        System.out.println("------------------------------- PanTools -------------------------------");        
        switch (args[0]) {
            case "reconstruct":
                PATH=args[2];
                program.reconstruct_genomes(args[1]);
                break;
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
                if(program.compare_pangenomes(args[1],args[2]))
                    System.out.println("Databases are equal.");
                else
                    System.out.println("Databases are different");
                break;
            case "retrieve":
                PATH=args[2];
                if(args[1].equals("genes"))
                    program.retrieve_genes(args[3]);
                else if(args[1].equals("regions"))
                    program.retrieve_regions(args[3]);
                else
                {
                    print_help_comment();
                    System.exit(1);       
                }
                break;
            case "query":
                PATH=args[1];
                program.run_query();
                break;
            default:
                print_help_comment();
                System.exit(1);            
        }
        System.out.println("-----------------------------------------------------------------------");        
    }
    
    /*
    To reconstructe genomes out of the pan-genome.
    genome_records : a text file containing number and a given name to the genes
    */      
    private void reconstruct_genomes(String genome_records)
    {
        if (new File(PATH+GENOME_DATABASE_PATH).exists()) 
        {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+GRAPH_DATABASE_PATH)
                    .setConfig( "keep_logical_logs","100M size").newGraphDatabase(); 
            registerShutdownHook( graphDb );
            startTime = System.currentTimeMillis();
            genomeDb=new genome_database(PATH+GENOME_DATABASE_PATH, graphDb);
            BufferedReader in;
            BufferedWriter out;
            Node seq_node,start;
            String[] fields;
            String line;
            int g,s;
            StringBuilder seq=new StringBuilder();
            try(Transaction tx = graphDb.beginTx()){
                db_node=graphDb.findNodes(pangenome_label).next();
                if(db_node==null)
                {
                    System.out.println("Can not locate database node!");
                    System.exit(1);
                }
                K = (int)db_node.getProperty("k_mer_size");
                if(genome_records.equals("all"))
                {
                    for(g=1;g<=genomeDb.num_genomes;++g)
                    {
                        System.out.println("Reconstructing genome "+g+"...");
                        try{
                            out = new BufferedWriter(new FileWriter(PATH+"/genome_"+g+".fasta"));  
                            for(s=1;s<=genomeDb.num_sequences[g];++s) 
                            {
                                out.write(">"+genomeDb.sequence_titles[g][s]+"\n");
                                seq_node=graphDb.findNode(DynamicLabel.label( "sequence" ), "number", g+"_"+s);
                                start=seq_node.getRelationships(Direction.OUTGOING).iterator().next().getEndNode();
                                extract_sequence(seq,start, true, 0, g+"_"+s, 0,  (int)genomeDb.sequence_length[g][s]-1);
                                write_fasta(out,seq.toString(),80);
                                seq.setLength(0);
                            }
                            out.close();
                        }catch (IOException e) {
                            System.out.println(e.getMessage());
                            System.exit(1);
                        }                    
                    }
                }
                else
                {
                    try{ 
                        in = new BufferedReader(new FileReader(genome_records));
                        while(in.ready())
                        {
                            line=in.readLine();
                            if(line.equals(""))
                                continue;    
                            fields=line.split(" ");
                            g=Integer.parseInt(fields[0]);
                            System.out.println("Reconstructing genome "+fields[1]+"...");
                            try{
                                out = new BufferedWriter(new FileWriter(PATH+(fields.length>1?"/"+fields[1]:"/genome_"+g)+".fasta"));  
                                for(s=1;s<=genomeDb.num_sequences[g];++s) 
                                {
                                    out.write(">"+genomeDb.sequence_titles[g][s]+"\n");
                                    seq_node=graphDb.findNode(DynamicLabel.label( "sequence" ), "number", g+"_"+s);
                                    start=seq_node.getRelationships(Direction.OUTGOING).iterator().next().getEndNode();
                                    extract_sequence(seq,start, true, 0, g+"_"+s, 0,  (int)genomeDb.sequence_length[g][s]-1);
                                    write_fasta(out,seq.toString(),80);
                                    seq.setLength(0);
                                }
                                out.close();
                            }catch (IOException e) {
                                System.out.println(e.getMessage());
                                System.exit(1);
                            } 
                        }
                        in.close();
                    }catch(IOException ioe){
                       System.out.println("Failed to read file names!");  
                       System.exit(1);
                    }
                }
            tx.success();}
            System.out.println("Genomes are ready in "+PATH);
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
    To builds a pan-genome database out of a set of input genomes.
    file : a text file containing paths to FASTA files
    */  
    private void build(String genome_paths)
    {
        int i,j;
        File theDir = new File(PATH);
        if (theDir.exists()) 
            try{
                FileUtils.deleteRecursively( new File( PATH ) ); // deletes old files
            }catch(IOException ioe)
            { 
            }
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
        genomeDb=new genome_database(PATH+GENOME_DATABASE_PATH,genome_paths);
        indexDb=new index_database(PATH+INDEX_DATABASE_PATH,genome_paths,K,genomeDb);
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
        construct_pangenome(0);
        System.out.println("Number of kmers:   "+indexDb.length());
        System.out.println("Number of nodes:   "+seq_nodes);
        System.out.println("Number of edges:   "+num_edges);
        System.out.println("Number of bases:   "+num_bases);
        System.out.println("Number of degenerate nodes:   "+degenerate_nodes);
        try(Transaction tx = graphDb.beginTx())
        {
            db_node.setProperty("k_mer_size", K);
            db_node.setProperty("num_k_mers", indexDb.length());
            db_node.setProperty("num_nodes",seq_nodes);
            db_node.setProperty("num_degenerate_nodes",degenerate_nodes);
            db_node.setProperty("num_edges",num_edges);
            db_node.setProperty("num_genomes",genomeDb.num_genomes);
            db_node.setProperty("num_genes",gene_nodes);
            db_node.setProperty("num_bases",num_bases);
            tx.success();
        }
        graphDb.shutdown();  
        System.out.println("Total time : "+(System.currentTimeMillis()-startTime)/1000+"."+(System.currentTimeMillis()-startTime)%1000+" seconds");
        print_peak_memory();
        try{
            for(i=0;Files.exists(Paths.get(PATH+GRAPH_DATABASE_PATH+"/neostore.transaction.db."+i));++i)
                Files.delete(Paths.get(PATH+GRAPH_DATABASE_PATH+"/neostore.transaction.db."+i));
        }catch(IOException ioe)
        { 
        }
        System.out.println("graph.db size: "+getFolderSize(new File(PATH+GRAPH_DATABASE_PATH))+" MB");
        System.out.println("index.db size: "+getFolderSize(new File(PATH+INDEX_DATABASE_PATH))+" MB");
        System.out.println("genome.db size: "+getFolderSize(new File(PATH+GENOME_DATABASE_PATH))+" MB");
    }
    /*
    To adds new genomes to an available pan-genome database.
    genome_paths : a text file containing paths to FASTA files
    */      
    private void add(String genome_paths)
    {
        int i,j,g,s,len,previous_num_genomes;
        long byte_number=0;
        Node start,seq_node;
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
                degenerate_nodes=(int)db_node.getProperty("num_degenerate_nodes");
                previous_num_genomes=(int)db_node.getProperty("num_genomes");
                if(!Files.exists(Paths.get(PATH+GENOME_DATABASE_PATH)))
                {
                    genomeDb=new genome_database(PATH+GENOME_DATABASE_PATH, graphDb);
                    StringBuilder seq=new StringBuilder();
                    for(g=1;g<=genomeDb.num_genomes;++g)
                    {
                        for(s=1;s<=genomeDb.num_sequences[g];++s)
                        {
                            seq_node=graphDb.findNode(DynamicLabel.label( "sequence" ), "number", g+"_"+s);
                            start=seq_node.getRelationships(Direction.OUTGOING).iterator().next().getEndNode();
                            extract_sequence(seq,start, true, 0, g+"_"+s, 0,  (int)genomeDb.sequence_length[g][s]-1);
                            len=seq.length();
                            if(len%2==1)
                                --len;
                            for(j=0;j<len;j+=2,++byte_number)
                                genomeDb.genomes_buff[(int)(byte_number/genomeDb.parts_size[0])].put((byte)((genomeDb.binary[seq.charAt(j)]<<4) | genomeDb.binary[seq.charAt(j+1)]));
                            if(len==seq.length()-1)
                            {
                                genomeDb.genomes_buff[(int)(byte_number/genomeDb.parts_size[0])].put((byte)(genomeDb.binary[seq.charAt(len)]<<4 ));
                                ++byte_number;                            
                            }                                  
                        }
                    }
                    //genomeDb.decode_genomes(PATH+GENOME_DATABASE_PATH);
                }
                else
                    genomeDb=new genome_database(PATH+GENOME_DATABASE_PATH);
                genomeDb.add_genomes(PATH+GENOME_DATABASE_PATH, genome_paths);
                indexDb=new index_database(PATH+INDEX_DATABASE_PATH,genome_paths, genomeDb,graphDb,previous_num_genomes);
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
                        drop_property(fwd_locations[i][j]);
                        rev_locations[i][j]="R"+i+"_"+j;
                        drop_property(fwd_locations[i][j]);
                    }
                }
                drop_property("sequence");
                construct_pangenome(previous_num_genomes);
                System.out.println("Number of kmers:   "+indexDb.length());
                System.out.println("Number of nodes:   "+seq_nodes);
                System.out.println("Number of edges:   "+num_edges);
                System.out.println("Number of bases:   "+num_bases);
                System.out.println("Number of degenerate nodes:   "+degenerate_nodes);
                try(Transaction tx = graphDb.beginTx())
                {
                    db_node.setProperty("k_mer_size", K);
                    db_node.setProperty("num_k_mers", indexDb.length());
                    db_node.setProperty("num_nodes",seq_nodes);
                    db_node.setProperty("num_degenerate_nodes",degenerate_nodes);
                    db_node.setProperty("num_edges",num_edges);
                    db_node.setProperty("num_genomes",genomeDb.num_genomes);
                    db_node.setProperty("num_genes",gene_nodes);
                    db_node.setProperty("num_bases",num_bases);
                    tx.success();
                }
                graphDb.shutdown(); 
                System.out.println("Total time : "+(System.currentTimeMillis()-startTime)/1000+"."+(System.currentTimeMillis()-startTime)%1000+" seconds");
                print_peak_memory();
                try{
                    for(i=0;Files.exists(Paths.get(PATH+GRAPH_DATABASE_PATH+"/neostore.transaction.db."+i));++i)
                        Files.delete(Paths.get(PATH+GRAPH_DATABASE_PATH+"/neostore.transaction.db."+i));
                }catch(IOException ioe)
                { 
                }
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
    To add gene_nodes to the graph. 
    gff_paths : a text file containing paths to annotation files 
    */      
    private void annotate(String gff_paths)
    {
        int i;
        int begin,end,s,line_len,g_number;
        String seq;
        Node gene_node,genome_node,start_node,stop_node;
        pointer start_ptr,stop_ptr;
        Relationship rel;
        StringBuilder log=new StringBuilder();
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
                try(BufferedReader gff_files = new BufferedReader(new FileReader(gff_paths))){
                 for(g_number=1;gff_files.ready();++g_number)
                    {
                        gff_name=gff_files.readLine();
                        genome_node=graphDb.findNode(genome_label, "number", g_number);
                        if(gff_name.equals("") || genome_node==null)
                            continue;
                        if(genome_node.hasProperty("annotated"))
                        {
                            log.append("genome number "+g_number+" already annotated.\n");
                            continue;
                        }
                        BufferedReader in = new BufferedReader(new FileReader(gff_name));
                        System.out.println("Annotating genome "+g_number+"...");
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
                                        rel.setProperty("forward", start_ptr.canonical);
                                        rel.setProperty("pos", start_ptr.position);
                                        stop_node=graphDb.getNodeById(stop_ptr.node_id);
                                        stop_node.setProperty("gene_stops", true);
                                        rel=gene_node.createRelationshipTo(stop_node, RelTypes.end);
                                        rel.setProperty("forward", stop_ptr.canonical);
                                        rel.setProperty("pos", stop_ptr.position);
                                    }
                                    else
                                        log.append("This feature is less than "+K+" nt long:\n"+line+"\n");
                                }
                                else 
                                {
                                     log.append(seq+" missed in genome "+g_number+"\n"); // usually organal genes
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
            try{ 
                BufferedWriter out = new BufferedWriter(new FileWriter(PATH+"/annotation.log"));
                out.write(log.toString());
                out.close();
            }
            catch(IOException ioe)
            { 
            }
            graphDb.shutdown(); 
            System.out.println(gene_nodes+" genes annotated.");
            System.out.println("Total time : "+(System.currentTimeMillis()-startTime)/1000+"."+(System.currentTimeMillis()-startTime)%1000+" seconds");
            print_peak_memory();
            try{
                for(i=0;Files.exists(Paths.get(PATH+GRAPH_DATABASE_PATH+"/neostore.transaction.db."+i));++i)
                    Files.delete(Paths.get(PATH+GRAPH_DATABASE_PATH+"/neostore.transaction.db."+i));
            }catch(IOException ioe)
            { 
            }
        }
        else
        {
            System.out.println("No database found in "+PATH); 
            System.exit(1);
        }  
    }
    /*
    To group ortholog/homolog/etc genes. 
    group_names : a text file with lines starting with a group name followed by a colon, followed by genes IDs seperated by one space.
    */     
    private void group(String group_names)
    {
        if (new File(PATH+GRAPH_DATABASE_PATH).exists()) 
        {
            Node group_node=null,gene;
            String line;
            String name=null;
            ResourceIterator<Node> genes;
            int g=0;
            System.out.println("Grouping genes...");
            try(BufferedReader in = new BufferedReader(new FileReader(group_names))){
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
               System.out.println("Failed to open "+group_names);  
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
    To retrieve genes' sequence from the graph
    annotation_records : a text file containing annotation records of the genes to be retrieved
    */ 
    private void retrieve_genes(String annotation_records)
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
            Relationship rstart;
            Node start,gene=null;
            String prp,record;
            String line;
            boolean strand;
            int i,j,begin,end,num_genes=0,genome,sequence;
            StringBuilder gene_seq=new StringBuilder();
            try{
                in = new BufferedReader(new FileReader(annotation_records));
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
                    in = new BufferedReader(new FileReader(annotation_records));
                    for(i=0;in.ready();)
                    {
                        line=in.readLine();
                        if(line.equals(""))
                            continue;                   
                        records[i]=line;
                        ++i;
                    }
                    in.close();
                    Arrays.sort(records);
                }catch(IOException ioe){
                   System.out.println("Failed to read file names!");  
                   System.exit(1);
                }                
                //genomeDb=new genome_database(PATH+GENOME_DATABASE_PATH);
                try{ 
                    BufferedWriter out = new BufferedWriter(new FileWriter(annotation_records+".fasta"));
                    for(i=j=0,gene_nodes=graphDb.findNodes(gene_label);gene_nodes.hasNext();++i)
                    {
                        gene=gene_nodes.next();
                        record=(String)gene.getProperty("title");
                        if(record!=null && Arrays.binarySearch(records, record) >=0)
                        {
                            //System.out.println(record);
                            rstart=gene.getSingleRelationship(RelTypes.begin, Direction.OUTGOING);
                            start=rstart.getEndNode();
                            genome=(int)gene.getProperty("genome_number");
                            sequence=(int)gene.getProperty("sequence_number");
                            prp=genome+"_"+sequence;
                            begin=(int)gene.getProperty("begin");
                            end=(int)gene.getProperty("end");
                            strand=gene.getProperty("strand").toString().equals("+");
                            extract_sequence(gene_seq,start,(boolean)rstart.getProperty("forward"),(int)rstart.getProperty("pos"),prp, begin-1, end-1);//
                            //if(gene_seq.toString().equals(genomeDb.get_sequence(genome, sequence, begin-1, end-begin+1, strand))
                            //|| gene_seq.toString().equals(genomeDb.get_sequence(genome, sequence, begin-1, end-begin+1, !strand)) )//gene_seq.length() == end-begin+1)//
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
                                System.out.println("Failed to assemble:\n"+record);
                            }
                            gene_seq.setLength(0);
                        }
                        if(i%(num_genes/100+1)==0) 
                            System.out.print((long)i*100/num_genes+1+"%\r");
                    }//for i
                    System.out.println(j+" out of "+i+" genes retrieved successfully.");
                    out.close();
                }catch(IOException ioe){
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
    To retrieves genes' sequence from the graph
    Input:
        file: A text file containing ID of genes whose sequence should be retrieved, One ID per line
    Output:
        A FASTA file of gene sequences.
    */ 
    private void run_query()
    {
        int i;
        String query=null;
        BufferedWriter out;
        if (new File(PATH+GRAPH_DATABASE_PATH).exists()) 
        {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(PATH+GRAPH_DATABASE_PATH)
                    .setConfig( "keep_logical_logs","100M size").newGraphDatabase(); 
            registerShutdownHook( graphDb );
            startTime = System.currentTimeMillis();
            for(i=1;;++i)
            {
                Console console = System.console();
                query = console.readLine("Enter a query or type exit : "); 
                if(query.equals("exit"))
                    break;
                try( Transaction ignored = graphDb.beginTx(); Result result = graphDb.execute(query) )
                {
                    try{
                        out = new BufferedWriter(new FileWriter(PATH+"/query"+i+".result"));
                        while ( result.hasNext() )
                        {
                            Map<String,Object> row = result.next();
                            for ( Entry<String,Object> column : row.entrySet() )
                            {
                               out.write(column.getKey() + ": " + column.getValue() + "\n");
                               System.out.println(column.getKey() + ": " + column.getValue());
                            }
                        }
                        out.close();
                    }catch(IOException ioe){
                       System.out.println("Failed to write the result of query!");  
                       System.exit(1);
                    }
                }catch(QueryExecutionException qee){
                       System.out.println(qee.getMessage());  
                }
            }
            graphDb.shutdown(); 
        }
    }    
    /*
    To retrieve a genomic regions from the graph
    region_records : a text file with lines containing genome number, sequence number, start and stop positions 
    */ 
    private void retrieve_regions(String region_records)
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
            StringBuilder seq=new StringBuilder();
            Node start_node;
            pointer start_ptr;
            int c,g,s,from, to, num_regions=0;
            try{
                BufferedReader in = new BufferedReader(new FileReader(region_records));
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
            startTime = System.currentTimeMillis();
            try(Transaction tx = graphDb.beginTx()){
            try (BufferedReader in = new BufferedReader(new FileReader(region_records))) {
                BufferedWriter out = new BufferedWriter(new FileWriter(region_records+".fasta"));
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
                    start_node=graphDb.getNodeById(start_ptr.node_id);
                    extract_sequence(seq,start_node, start_ptr.canonical,start_ptr.position,g+"_"+s, from, to);
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
    To extract and return a sequence belonging to sequence_number starting at coordinate (node, forward, position) of length end-begin+1
    */
    private void extract_sequence(StringBuilder seq,Node node, boolean forward,int position, String sequence_number, int begin, int end)
    {
        boolean found;
        Node neighbor;
        String lable,name;
        String[] prp={"F"+sequence_number,"R"+sequence_number};
        int loc,node_len,neighbor_len,s_side,d_side,seq_len=end-begin+1;
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
        while(seq.length()<seq_len && found)
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
    }
    /*
    To append s[from..to] to "seq". 
    */
    private int append_fwd(StringBuilder seq, String s, int from, int to)
    {
        seq.append(s.substring(from,to+1));//genomeDb.get_sequence(addrs[0],addrs[1],addrs[2]+from, to-from+1));
        return to-from+1;
    }
    /*
    To append reverse_complement of s[from..to] to "seq". 
    */
    private int append_rev(StringBuilder seq, String s, int from, int to)
    {
        seq.append(reverse_complement(s.substring(from,to+1)));//genomeDb.get_sequence(addrs[0],addrs[1],addrs[2]+from, to-from+1)));
        return to-from+1;
    }
    /*
    To return the coordinate pointing to the genomic position (g,s,position)
    */
    private pointer locate(int g, int s, int position)
    {
        int i,loc,low,high,mid=0,node_len;
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
        for(low=0,high=anchor_sides.length()-1,mid=(low+high)/2;low<=high;mid=(low+high)/2)
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
                forward=(r.getType().name().charAt(2)-48)/2==1;
                lable=forward?prp[0]:prp[1];
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
    To return the number of the sequence from genome "g" whose name contains "prefix".
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
    To return the reverse complement of "s".
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
    To create an edge between "src" and "des" nodes. 
    Type of the edge is determined by "s_side", "d_side" and transition bases.
    */     
    void connect(Node src, Node des,int s_side, int d_side)
    {
        RelationshipType f_rt,r_rt;
        int s_base,d_base=0;
        String rt_name;
        boolean has_fwd, has_rev;
        s_base=     get_source_symbol(src,s_side);
        d_base=get_destination_symbol(des,d_side);
        f_rt=RelTypes.values()[s_base*20+d_base*4+s_side*2+d_side];
        if(src.hasLabel(degenerate_label) || des.hasLabel(degenerate_label) 
        || src.hasLabel(sequence_label)   || des.hasLabel(sequence_label) )
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
    */
    int get_source_symbol(Node src, int s_side)
    {
        int s_base,s_len;
        int[] s_add;
        if(src.hasLabel(sequence_label))
            s_base=4;
        else
        {
            s_len= (int)src.getProperty("length");
            s_add=(int[])src.getProperty("address");
            s_base=s_side==0?genomeDb.get_code(s_add[0],s_add[1],s_add[2]+s_len-K):genomeDb.get_complement_code(s_add[0],s_add[1],s_add[2]+K-1);
            if(s_base>3)
                s_base=4;
        }
        return s_base;
    }
    /*
    */
    int get_destination_symbol(Node des,int d_side)
    {
        int d_base,d_len;
        int[] d_add;
        if(des.hasLabel(sequence_label))
            d_base=4;
        else
        {
            d_len= (int)des.getProperty("length") ;
            d_add=(int[])des.getProperty("address");
            d_base=d_side==0?genomeDb.get_code(d_add[0],d_add[1],d_add[2]+K-1):genomeDb.get_complement_code(d_add[0],d_add[1],d_add[2]+d_len-K);
            if(d_base>3)
                d_base=4;
        }
        return d_base;
    }
    /*
    To split 'node' at position 'pos' by creating a 'split_node' as a part has been seperated from 'node'   
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
        address[0]=gen;
        address[1]=seq;
        address[2]=loc+pos;
        split_node.setProperty("address",address);
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
        connect(node,split_node,0,0);
    }
    /*
    To extend a new node till reach to a previously visited K-mer
    */
    void extend(Node node)
    {
        //System.out.println("extend "+position);
        int begin,len=(int)node.getProperty("length");
        long id=node.getId(),last_kmer=(long)node.getProperty("last_kmer");
        boolean broke=false;
        while(position<seq_len-1)
        {
            //System.out.println("extend "+position);
            if(genomeDb.get_code(genome,sequence,position+1)>3)
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
                    System.out.print((long)position*100/seq_len+"%\r");
            }
            else
            {
                broke=true;
                break;
            }
        }
        node.setProperty("length",len);
        node.setProperty("last_kmer",last_kmer);
        if(!broke && position==seq_len-1 )
            finish=true;
    }
    /*
    To create a new node
    */
    void create()
    {
        //System.out.println("create "+position);
        num_bases+=K;
        ++seq_nodes;
        new_node=graphDb.createNode(node_label);
        address[0]=genome;
        address[1]=sequence;
        address[2]=position-K+1;
        new_node.setProperty("address",address);
        new_node.setProperty("length",K);
        new_node.setProperty("last_kmer",curr_index);
        new_node.setProperty("first_kmer",curr_index);
        pointer.node_id=new_node.getId();
        pointer.canonical=fwd_kmer.canonical;
        pointer.position=0;
        pointer.next_index=-1L;
        indexDb.put_pointer(pointer, curr_index);
        connect(curr_node,new_node,curr_side,0);
        curr_node=new_node;
        curr_side=0;
    }
    /*
    To enter a node and follow it in reverse direction till reaching to a split position or end of the node
    */
    void follow_forward()
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
        connect(curr_node,node,curr_side,0);                
        curr_node=node;
        curr_side=0;
        l=(int)curr_node.getProperty("length")-K;
        address=(int[])curr_node.getProperty("address");
        g=address[0];
        s=address[1];
        loc=address[2];
        p=position-K+1;
        for(pos=0;pos<=l && position<=seq_len-1 && genomeDb.get_code(g, s, loc+pos+K-1)==genomeDb.get_code(genome,sequence,position); ++pos)
        {
            ++position;
            if(position<=seq_len-1 && genomeDb.get_code(genome,sequence,position)>3)
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
        if(position==seq_len )
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
        if(degenerated)
        {
            curr_node=degenerate_node;
            curr_side=0;
        }
        else
        {
            curr_index=indexDb.find(k_mer);
            indexDb.get_pointer(byte_pointer,pointer,curr_index);            
        }
    }
    /*
    To enter a node and follow it in forward direction till reaching to a split position or end of the node
    */
    void follow_reverse()
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
        connect(curr_node,node,curr_side,1);                
        curr_node=node;
        curr_side=1;
        p=position-K+1;
        address=(int[])curr_node.getProperty("address");
        g=address[0];
        s=address[1];
        loc=address[2];
        for(pos=(int)node.getProperty("length")-K;pos>=0 && position<=seq_len-1 && genomeDb.get_code(g,s,loc+pos)==genomeDb.get_complement_code(genome,sequence,position); --pos) 
        {
            ++position;
            if(position<=seq_len-1 && genomeDb.get_code(genome,sequence,position)>3)
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
        if(position==seq_len )
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
        if(degenerated)
        {
            curr_node=degenerate_node;
            curr_side=0;
        }
        else
        {
            curr_index=indexDb.find(k_mer);
            indexDb.get_pointer(byte_pointer,pointer,curr_index);            
        }
    }
    /*
    To jump over an ambiguous region.
    position points to the first position which degenerate starts, after jumping it points to the last base of the first K-mer after the ambiguous region 
    */
    void jump()
    {
        int j;
        fwd_code=genomeDb.get_code(genome,sequence,position);
        rev_code=3-fwd_code;
        do{
            while(fwd_code>3 && position<seq_len-1)
            {
                ++position;
                fwd_code=genomeDb.get_code(genome,sequence,position);
                rev_code=3-fwd_code;
            }
            fwd_kmer.reset();//=new kmer(K,indexDb.get_pre_len(),indexDb.get_suf_len());
            rev_kmer.reset();//=new kmer(K,indexDb.get_pre_len(),indexDb.get_suf_len());
            fwd_kmer.next_up_kmer(fwd_code);
            rev_kmer.next_down_kmer(rev_code);
            for(j=0; j<K-1 && position<seq_len-1;++j)
            {
                ++position;
                fwd_code=genomeDb.get_code(genome,sequence,position);
                rev_code=3-fwd_code;
                if(fwd_code>3 )
                    break;
                fwd_kmer.next_up_kmer(fwd_code);
                rev_kmer.next_down_kmer(rev_code);
            }
            if(j==K-1)
            {
                fwd_kmer.canonical=fwd_kmer.compare(rev_kmer)==-1;
                k_mer=fwd_kmer.canonical?fwd_kmer:rev_kmer;
            }
            else if(position==seq_len-1)
                finish=true;
        }while(fwd_code>3 && position<seq_len-1);
    }
    /*
    To create a degenerate node starting at "begin" ending at position-1
    and connects it to previous node and makes the first K-mer after the ambiguous region
    */
    void create_degenerate(int begin)
    {
        //System.out.println("create_degenerate "+position);
        Node node=null;
        boolean exist=false;
        address[0]=genome;
        address[1]=sequence;
        address[2]=begin;        
        for(Relationship r:curr_node.getRelationships(Direction.OUTGOING))
        {
            node=r.getEndNode();
            if(node.hasLabel(degenerate_label) && 
            genomeDb.compare(genomeDb, address, (int[])node.getProperty("address"), (int)node.getProperty("length"), true)) 
            {
                exist=true;
                degenerate_node=node;
                break;
            }
        }      
        if(!exist)
        {
            ++degenerate_nodes;
            degenerate_node=graphDb.createNode(degenerate_label);
            //degenerate_node.setProperty("genome",genome);
            degenerate_node.setProperty("address",address);
            degenerate_node.setProperty("length",position-begin);
            connect(curr_node,degenerate_node,curr_side,0);
            num_bases+=(position-begin);
        }
        if(!finish)
        {
            curr_index=indexDb.find(k_mer);
            indexDb.get_pointer(byte_pointer,pointer,curr_index);      
        }
    }
    /*
    To produce the next "fwd_kmer", "rev_kmer" and the canonical "kmer" from the current ones.
    */ 
    void next_kmer()
    {
        ++position;
        fwd_code=genomeDb.get_code(genome,sequence,position);
        rev_code=3-fwd_code;
        fwd_kmer.next_up_kmer(fwd_code);
        rev_kmer.next_down_kmer(rev_code);
        if(position%(seq_len/100+1)==0) 
            System.out.print((long)position*100/seq_len+1+"%\r");
        fwd_kmer.canonical=fwd_kmer.compare(rev_kmer)==-1;
        k_mer=fwd_kmer.canonical?fwd_kmer:rev_kmer;
    }
    /*
    To initialize the first K-mer of the "genome", "sequence" at "position". 
    It might jump over the degenerate regions creating a degenerate_node 
    */ 

    /**
     *
     * @throws IOException
     */
     
    public void initial_kmers()
    {
        int i;
        fwd_kmer.reset();//=new kmer(K,indexDb.get_pre_len(),indexDb.get_suf_len());
        rev_kmer.reset();//=new kmer(K,indexDb.get_pre_len(),indexDb.get_suf_len());
        for(i=0;i<K && position<seq_len-1;++i)
        {
            if(genomeDb.get_code(genome,sequence,position+1)>3 )
            {
                ++position;
                jump();
                create_degenerate(0);
                curr_node=degenerate_node;
                break;
            }
            next_kmer();
        }   
        fwd_kmer.canonical=fwd_kmer.compare(rev_kmer)==-1;
        k_mer=fwd_kmer.canonical?fwd_kmer:rev_kmer;
    }
    /*
    To give the K-mer in genomic position (genome,sequence,position) 
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
            fwd_code=genomeDb.get_code(genome,sequence,position+j);
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
    void construct_pangenome(int previous_num_genomes)
    {
        int i;
        Node genome_node,sequence_node;
        long[] sequence_ids;
        phaseTime = System.currentTimeMillis();
        fwd_kmer=new kmer(K,indexDb.get_pre_len(),indexDb.get_suf_len());
        rev_kmer=new kmer(K,indexDb.get_pre_len(),indexDb.get_suf_len()); 
        address=new int[3];
        for(genome=previous_num_genomes+1;genome<=genomeDb.num_genomes;++genome) 
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
                    sequence_node.setProperty("sequence_title",genomeDb.sequence_titles[genome][sequence]);
                    sequence_node.setProperty("sequence_length",genomeDb.sequence_length[genome][sequence]);
                    genome_node.createRelationshipTo(sequence_node, RelTypes.has);
                    sequence_ids[sequence]=curr_node.getId();
                    curr_side=0;
                    position=-1;
                    seq_len=genomeDb.sequence_length[genome][sequence];
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
                    connect(curr_node,sequence_node,curr_side,0);// to point to the last k-mer of the sequence located in the other strand
                tx.success();}
            }//sequences
            try(Transaction tx = graphDb.beginTx()){
                genome_node.setProperty("sequence_ids", sequence_ids);
            tx.success();}
            System.out.println("Running time : " + (System.currentTimeMillis() - phaseTime) / 1000 + " seconds");
        }//genomes
        include_sequences();        
        index_genomes();
    }
    /*
    To add list of "anchor_nodes", "anchor_sides" and "anchor_positions" to each sequence_node.
    These properties facilitate extracting genomic regions. 
    */ 
    void index_genomes() 
    {
        int i,len,m,neighbor_side,count;
        long seq_len;
        long[] anchor_nodes;
        int[] anchor_positions;
        int[] initial_coordinate=new int[1];
        Node node,neighbor,seq_node;
        StringBuilder nds=new StringBuilder();
        StringBuilder pos=new StringBuilder();
        StringBuilder sds=new StringBuilder();
        String[] ids_list,posis_list;
        String prp;
        int[] positions;
        int[] new_positions;
        boolean found;
        for(address[0]=1;address[0]<=genomeDb.num_genomes;++address[0]) 
        {
            for(address[1]=1;address[1]<=genomeDb.num_sequences[address[0]];++address[1]) 
            {
                System.out.println("Indexing sequence "+address[1]+"/"+genomeDb.num_sequences[address[0]]+" of genome "+address[0]+"...                        ");
                try(Transaction tx = graphDb.beginTx()){
                seq_len=genomeDb.sequence_length[address[0]][address[1]]-1;
                node=seq_node=graphDb.findNode(sequence_label, "number", address[0]+"_"+address[1]);
                found=true;
                count=0;
                for(address[2]=0;address[2]+K-1<seq_len && found;) // K-1 bases of the last node not added
                {
                    found=false;
                    //System.out.println("out");
                    for(Relationship r:node.getRelationships(Direction.OUTGOING))
                    {
                        neighbor=r.getEndNode();
                        neighbor_side=(r.getType().name().charAt(2)-48)%2;
                        if(neighbor.hasProperty("address") && 
                        genomeDb.compare(genomeDb, (int[])neighbor.getProperty("address"), address, (int)neighbor.getProperty("length"), neighbor_side==0) ) 
                        {
                            found=true;
                            prp=neighbor_side==0?fwd_locations[address[0]][address[1]]:rev_locations[address[0]][address[1]];
                            if(neighbor.hasProperty(prp))
                            {
                                positions=(int[])neighbor.getProperty(prp);
                                len = positions.length;
                                new_positions=new int[len+1];
                                for(i=0;i<len;++i)
                                    new_positions[i]=positions[i];
                                new_positions[i]=address[2];
                                neighbor.setProperty(prp, new_positions);
                            }
                            else
                            {
                                initial_coordinate[0]=address[2];
                                neighbor.setProperty(prp,initial_coordinate);
                            }
                            if(count % anchor_distance==0)
                            {
                                nds.append(neighbor.getId()).append(" ");
                                sds.append(neighbor_side);
                                pos.append(address[2]).append(" ");
                            }
                            count++;
                            address[2]=address[2]+(int)neighbor.getProperty("length")-K+1;
                            node=neighbor;
                            break;
                        }
                        
                    }
                    if(!found)
                    {
                    //System.out.println("in ");
                        for(Relationship r:node.getRelationships(Direction.INCOMING))
                        {
                            neighbor=r.getStartNode();
                            neighbor_side=(r.getType().name().charAt(2)-48)/2;
                            if(neighbor.hasProperty("address") && genomeDb.compare(genomeDb, (int[])neighbor.getProperty("address"), address, (int)neighbor.getProperty("length"), neighbor_side==1)  ) 
                            {
                                found=true;
                                prp=neighbor_side==1?fwd_locations[address[0]][address[1]]:rev_locations[address[0]][address[1]];
                                if(neighbor.hasProperty(prp))
                                {
                                    positions=(int[])neighbor.getProperty(prp);
                                    len = positions.length;
                                    new_positions=new int[len+1];
                                    for(i=0;i<len;++i)
                                        new_positions[i]=positions[i];
                                    new_positions[i]=address[2];
                                    neighbor.setProperty(prp, new_positions);
                                }
                                else
                                {
                                    initial_coordinate[0]=address[2];
                                    neighbor.setProperty(prp,initial_coordinate);
                                }
                                if(count % anchor_distance==0)
                                {
                                    nds.append(neighbor.getId()).append(" ");
                                    sds.append(1-neighbor_side);
                                    pos.append(address[2]).append(" ");
                                }
                                count++;
                                address[2]=address[2]+(int)neighbor.getProperty("length")-K+1;
                                node=neighbor;
                                break;
                            }
                        }
                    }
                    if(address[2]%(seq_len/100+1)==0) 
                        System.out.print((long)address[2]*100/seq_len+1+"%    \r");
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
    To add "sequence" property to the nodes
    */ 
    void include_sequences()
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
                    node.setProperty("sequence",genomeDb.get_sequence(address[0], address[1], address[2], len,true));
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
                    node.setProperty("sequence",genomeDb.get_sequence(address[0], address[1], address[2], len,true));
                }
            tx.success();}
        nodes.close();
    }
    /*
    To add "sequence" property to the nodes
    */ 
    void drop_property(String prp)
    {
        int i;
        byte[] sequence;
        ResourceIterator<Node> nodes;
        Node node;
        try(Transaction tx = graphDb.beginTx()){
        nodes = graphDb.findNodes( node_label );
        tx.success();}
        while(nodes.hasNext())
            try(Transaction tx = graphDb.beginTx()){
                for(i=0;i<trsc_limit && nodes.hasNext();++i)
                {
                    node=nodes.next();
                    node.removeProperty(prp);
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
                    node.removeProperty(prp);
                }
            tx.success();}
        nodes.close();
    }    /*
    To register the action to be taken if the program halts unexpectedly
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
    To compare two graph databases located in paths "path1" and "path2"
    */    
    private boolean compare_pangenomes(String path1, String path2)
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
                indexDb1=new index_database(path1+INDEX_DATABASE_PATH,genomeDb1);
                indexDb2=new index_database(path2+INDEX_DATABASE_PATH,genomeDb2);
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
    To write the "seq" in a FASTA file with lines justified to lines "length" long 
    */
    private void write_fasta(BufferedWriter fasta_file, String seq, int length)
    {
        int i;
        try{
            for(i=1;i<=seq.length();++i)
            {
                fasta_file.write(seq.charAt(i-1));
                if(i%length==0)
                    fasta_file.write("\n");
            }
            fasta_file.write("\n");
        }catch( IOException ioe)
        {
            
        }
        
    }
    /*
    To provide a simple binary search
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
    To calculates and prints peak of memory usage of the program in mega bytes.
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
    /*
    To print a simple manual for the users
    */    
     private static void print_help_comment()
    {
        System.out.println("Requirements:\n" +
        "\n" +
        "- KMC: add paths to kmc and kmc_tools executables to the OS path environment variable.\n" +
        "- neo4j-community-2.3.1-unix\n" +
        "- Java Virtual Machine jdk1.7 : add the corresponding path to java executable\n" +
        "\n" +
        "To run the program:\n" +
        "java  [-server] [-XX:+UseConcMarkSweepGC]  [-Xmx(a number followed by g or m)] -jar ./pantools/dist/pantools.jar command arguments\n" +
        "\n" +
        "\n" +
        "List of commands:\n" +
        "\n" +
        "1. reconstruct:\nTo reconstruct all or a set of genomes out of the pan-genome\njava -jar pantools.jar reconstruct all [or genome_names_file] database_path\n" +
        "\n" +
        "2. build:\nTo build a pan-genome out of a set of genomes\njava -jar pantools.jar build K database_path fasta_names_file\n" +
        "\n" +
        "3. annotate:\nTo add annotations to a pan-genome\njava -jar pantools.jar annotate database_path gff_names_file\n" +
        "\n" +
        "4. add:\nTo add new genomes to an available pan-genome\njava -jar pantools.jar add database_path fasta_names_file\n" +
        "\n" +
        "5. retrieve genes:\nTo extract sequence of some annotated genes\njava -jar pantools.jar retrieve_genes database_path annotation_records\n" +
        "\n" +
        "6. retrieve regions:\nTo extract region sequence\njava -jar pantools.jar retrieve_regions database_path regions_file\n" +
        "        \n" +
        "7. group:\nTo group some genes by adding group nodes pointing to them\njava -jar pantools.jar group database_path groups_file\n" +
        "\n" +
        "8. compare:\nTo compare two pan-genomes\njava -jar pantools.jar compare database_1_path database_2_path\n" +
        "\n" +
        "8. query:\nTo compare write Cypher queries and receive the results\njava -jar pantools.jar query database_path\n" +
        "\n\n" +
        "Description of commands' arguments:\n" +
        "\n" +
        "a) K :\nSize of K for de Bruijn graph\n" +
        "\n" +
        "b) genome_names_file:\nA text file containing genome_number and genome_name seperated by a single space in each line\n" +
        "\n" +
        "c) fasta_names_file:\nA text file containing paths to FASTA files; each in one line\n" +
        "\n" +
        "d) gff_names_file:\nA text file containing paths to GFF files corresponding to the stored genomes,\n"
         + "in the same order. Missing annotations are shown by an empty line.\n" +
        "\n" +
        "e) annotation_records:\nA text file containing annotation titles as they appear in gff file.\n" +
        "\n" +
        "f) regions_file:\nA text file containing genome_number, sequence_number, begin and\n"
         + "end of a region in each line seperated by one space \n" +
        "\n" +
        "g) groups_file:\nA FASTA file with titles being names given to the groups "
        + "followed by lines containing annotation titles as they appear in gff files.");
    }
    /*
    To return size of a folder
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
