/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package index;

import sequence.SequenceDatabase;
import static index.kmer.adjust_fwd_kmer;
import static index.kmer.compare_suffix;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.neo4j.io.fs.FileUtils;
//import org.apache.commons.io.FileUtils;
import static pantools.Pantools.MAX_TRANSACTION_SIZE;
import static pantools.Pantools.cores;
import static pantools.Pantools.degenerate_label;
import static pantools.Pantools.labels;
import static pantools.Pantools.executeCommand;
import static pantools.Pantools.nucleotide_label;

/**
 * Implements all the functionality to work with a KMC-based kmer index database. 
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public final class IndexDatabase {

    public int MAX_BYTE_COUNT = 1 << 24;// 16 MB 
    public int header_pos;
    public final String INFORMATION_FILE_EXTENTION = ".info";
    public final String PREFIX_FILE_EXTENTION = ".kmc_pre";
    public final String SUFFIX_FILE_EXTENTION = ".kmc_suf";
    public final String POINTER_FILE_EXTENTION = ".pointer";

    public static int K;
    public int mode;
    public int ctr_size;
    public int pre_len;
    public int min_count;
    public int max_count;
    public long kmers_num;
    public int id_len;
    public int offset_len;

    public int suf_len;
    public  int poniter_length;
    public long[] prefix_ptr;
    public RandomAccessFile suf_file;
    public RandomAccessFile pre_file;
    public MappedByteBuffer[] suf_buff;
    public RandomAccessFile ptr_file;
    public MappedByteBuffer[] ptr_buff;
    public MappedByteBuffer pre_buff[];
    public String KMC = "kmc";


    /**
     * Mounts an available index database to the index database object
     * 
     * @param index_path Path to the index database 
     */
    public IndexDatabase(String index_path, String index_name) {
        int i, j, number_of_pages, record_size, page_size, full_page_size;
        long number_of_prefixes;
        try {
            //System.out.println("Mounting index " + index_path);
            pre_file = new RandomAccessFile(index_path + index_name + PREFIX_FILE_EXTENTION, "r");
            BufferedReader in = new BufferedReader(new FileReader(index_path + index_name + INFORMATION_FILE_EXTENTION));
            K = Integer.parseInt(in.readLine().trim().split(":")[1]);
            mode = Integer.parseInt(in.readLine().trim().split(":")[1]);
            ctr_size = Integer.parseInt(in.readLine().trim().split(":")[1]);
            pre_len = Integer.parseInt(in.readLine().trim().split(":")[1]);
            min_count = Integer.parseInt(in.readLine().trim().split(":")[1]);
            max_count = Integer.parseInt(in.readLine().trim().split(":")[1]);
            kmers_num = Long.valueOf(in.readLine().trim().split(":")[1]);
            suf_len = Integer.parseInt(in.readLine().trim().split(":")[1]);
            id_len = Integer.parseInt(in.readLine().trim().split(":")[1]);
            offset_len = Integer.parseInt(in.readLine().trim().split(":")[1]);
            poniter_length = Integer.parseInt(in.readLine().trim().split(":")[1]);
            in.close();
            //System.out.println("Indexing " + kmers_num + " kmers...                    ");
        // load the prefix file into the memory    
            //System.out.println("Loading existing prefixes in memory...                    ");
            number_of_prefixes = 1 << (2 * pre_len);
            prefix_ptr = new long[(int)number_of_prefixes];
            record_size = 8; // in bytes
            full_page_size = MAX_BYTE_COUNT / record_size * record_size; // in bytes
            number_of_pages = (int) ((number_of_prefixes * record_size) / full_page_size + ((number_of_prefixes * record_size) % full_page_size == 0 ? 0 : 1));
            pre_buff = new  MappedByteBuffer[number_of_pages];
            //System.out.println("pre_len "  + pre_len);
            //System.out.println("number_of_prefixes  "+ number_of_prefixes);
            //System.out.println("full_page_size  "+ full_page_size);
            //System.out.println("number_of_pages  "+ number_of_pages);
            for (i = 0; i < number_of_pages; ++i) {
                page_size = (int) (i == number_of_pages - 1 ? (number_of_prefixes * record_size) % full_page_size : full_page_size); // in bytes
                page_size = page_size == 0 ? full_page_size : page_size; // in bytes
                //System.out.println("page_size  "+ page_size);
                pre_buff[i] = pre_file.getChannel().map(FileChannel.MapMode.READ_ONLY, ((long)full_page_size) * i + 4, page_size);
                pre_buff[i].position(0);
                for (j = 0; j < page_size / 8; ++j) 
                    prefix_ptr[(int)(((long)full_page_size) * i / 8 + j)] = read_prefix(pre_buff[i]);
                pre_buff[i] = null;
            }
            pre_buff = null;
            pre_file.close();
        // mapping suffix file into the memory
            //System.out.println("Mapping existing suffix file in memory...                    ");
            record_size = ctr_size + suf_len / 4; // in bytes
            full_page_size = MAX_BYTE_COUNT / record_size * record_size; // in bytes
            number_of_pages = (int) ((kmers_num * record_size) / full_page_size + ((kmers_num * record_size) % full_page_size == 0 ? 0 : 1));
            suf_file = new RandomAccessFile(index_path + index_name + SUFFIX_FILE_EXTENTION, "r");
            suf_buff = new MappedByteBuffer[number_of_pages];
            for (i = 0; i < number_of_pages; ++i) {
                page_size = (int) (i == number_of_pages - 1 ? (kmers_num * record_size) % full_page_size : full_page_size); // in bytes
                page_size = page_size == 0 ? full_page_size : page_size; // in bytes
                suf_buff[i] = suf_file.getChannel().map(FileChannel.MapMode.READ_ONLY, ((long)full_page_size) * i + 4, page_size);
            }
        // mapping pointers file into the memory
            //System.out.println("Mapping existing pointers file in memory...                    ");
            record_size = poniter_length; // in bytes
            full_page_size = MAX_BYTE_COUNT / record_size * record_size; // in bytes
            number_of_pages = (int) ((kmers_num * record_size) / full_page_size + ((kmers_num * record_size) % full_page_size == 0 ? 0 : 1));
            ptr_file = new RandomAccessFile(index_path + index_name + POINTER_FILE_EXTENTION, "rw");
            ptr_buff = new MappedByteBuffer[number_of_pages];
            for (i = 0; i < number_of_pages; ++i) {
                page_size = (int) (i == number_of_pages - 1 ? (kmers_num * record_size) % full_page_size : full_page_size); // in bytes
                page_size = page_size == 0 ? full_page_size : page_size; // in bytes
                ptr_buff[i] = ptr_file.getChannel().map(FileChannel.MapMode.READ_WRITE, ((long)full_page_size) * i + 4, page_size);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
    }

    /**
     * Creates a new index database
     * 
     * @param index_path Path to the index database 
     * @param genomes_path_file A text file containing path to the genomes
     * @param K_size Size of K
     * @param genomeDb The genome database
     */
    public IndexDatabase(String index_path, String genomes_path_file, SequenceDatabase genomeDb, int k) {
        long p;
        int i, j, number_of_pages, record_size, page_size, full_page_size;
        long number_of_prefixes, longest_scaffold = 0;
        System.out.println("Creating index in " + index_path);
        try {
            Files.createDirectory(Paths.get(index_path));
            find_k(index_path, genomes_path_file, genomeDb, k);
        /*
        All integers in the KMC output files are stored in LSB (least significant byte first) format.
            
        A .kmc pre file contains, in the order:
        [marker] : 4 bytes with the letters: KMCP.
        [data]: Here are k-mer prefix data. More precisely, it is an array of uint64 elements, of size 4^prefix length. Position x
                in the array stores the number of different k-mers (in binary) whose prefix is less
                than x. DNA symbols are encoded as follows: A   0, C   1, G   2, T   3. For example, if the queried
                k-mer is ACGACTGA and lut prefix length = 4, then we cut off the first 4 symbols, i.e., the prefix ACGA,
                which is interpreted in binary as x = 24 (since 0 (* 2^6 + 1 * 2^4 + 2 * 2^2 + 0 * 2^0 = 24). Now we look into
                “data” at locations x and x+1 to read, e.g., 1523 and 1685. This means that in the file .kmc suf in records from
                1523 to 1685􀀀1 there are suffixes of the k-mers with prefix ACGA. Having got this range, we can now binary
                search the suffix CTGA.            
        [header] : 
        [header position] : An integer consisting of the last 4 bytes in the file. Its contains the relative position of the beginning of the
        field [header].
        [marker] : another copy, to signal the file is not truncated.
        
        After opening the file, one should do the following:
        1. Read the first 4 bytes and check if they contain the letters KMCP.
        2. Jump to position 4 bytes back from end of file (excluding the marker) and read the header position x.
        3. Jump to position x + 4 (excluding the marker) bytes back from end of file and read the header.
        4. Read [data].
            
        A .kmc suf file contains, in order:
        [marker] : 4 bytes with the letters: KMCS
        [data] :An array record t records[total kmers].
                total kmers is a value taken from the .kmc pre file.
                record t is a type representing a k-mer. Its first field is the k-mer suffix string, stored on (kmer length 􀀀
                lut prefix length)=4 bytes. The next field is counter size, with the number of bytes used by the counter,
                which is either a 1 : : : 4-byte integer, or a 4-byte float.            
                The k-mers are stored with their leftmost symbol first, packed into bytes. For example, CCACAAAT is
                represented as 0x51 (for CCAC), 0x03 (for AAAT). Integers are stored according to the LSB (little endian)
                convention, floats are stored in the same way as they are stored in the memory.
        [marker] (another copy, to signal the file is not truncated).
        */
            pre_file = new RandomAccessFile(index_path + "sorted" + PREFIX_FILE_EXTENTION, "r");
            pre_file.seek(pre_file.length() - 8);
            header_pos = read_int(pre_file);
            pre_file.seek(pre_file.length() - (8 + header_pos));
        // read the header of the index    
            K = read_int(pre_file);
        //  k-mer length   
            mode = read_int(pre_file); 
        // 0 (occurrence count) or 1 (counting according to Quake quality)
            ctr_size = read_int(pre_file); 
        // counter field size: for mode 0 it is 1, 2, 3, or 4; for mode 1 it is always 4
            pre_len = read_int(pre_file);  
        // prefix length such that suffix length is divisible by 4. 
        // Max prefix length is limited up to 15 in KMC (exact value is calculated in such a way that summary size of kmc_pre and kmc_suf should be minimal)
            min_count = read_int(pre_file);
        // minimum number of k-mer occurrences to be written in the database    
            max_count = read_int(pre_file);
        // maximum number of k-mer occurrences to be written in the database    
            kmers_num = read_long(pre_file);
        // total number of k-mers in the database    
            suf_len = K - pre_len;
            for (i = 1; i <= genomeDb.num_genomes; ++i)
                for (j = 1; j <= genomeDb.num_sequences[i]; ++j)
                    if (genomeDb.sequence_length[i][j] > longest_scaffold)
                        longest_scaffold = genomeDb.sequence_length[i][j];
            id_len = (int)Math.round(Math.ceil( Math.log(kmers_num) / Math.log(2) / 8));
            offset_len = (int)Math.round(Math.ceil( Math.log(longest_scaffold) / Math.log(2) / 8));
            poniter_length = 2 * id_len + offset_len + 1;
            write_info(index_path, "sorted");
            System.out.println("Indexing kmers...                    ");
        // load the prefix file into the memory    
            //System.out.println("Loading prefixes in memory...                    ");
            number_of_prefixes = 1 << (2 * pre_len);
            prefix_ptr = new long[(int)number_of_prefixes];
            record_size = 8; // in bytes
            full_page_size = MAX_BYTE_COUNT / record_size * record_size; // in bytes
            number_of_pages = (int) ((number_of_prefixes * record_size) / full_page_size + ((number_of_prefixes * record_size) % full_page_size == 0 ? 0 : 1));
            pre_buff = new  MappedByteBuffer[number_of_pages];
            for (i = 0; i < number_of_pages; ++i) {
                page_size = (int) (i == number_of_pages - 1 ? (number_of_prefixes * record_size) % full_page_size : full_page_size); // in bytes
                page_size = page_size == 0 ? full_page_size : page_size; // in bytes
                //System.out.println("page_size  "+ page_size);
                pre_buff[i] = pre_file.getChannel().map(FileChannel.MapMode.READ_ONLY, ((long)full_page_size) * i + 4, page_size);
                for (j = 0; j < page_size / 8; ++j) 
                    prefix_ptr[(int)(((long)full_page_size) * i / 8 + j)] = read_prefix(pre_buff[i]);
                pre_buff[i] = null;
            }
            pre_buff = null;
            pre_file.close();
        // mapping suffix file into the memory
            //System.out.println("Mapping suffixes file to memory...                    ");
            record_size = ctr_size + suf_len / 4; // in bytes
            full_page_size = MAX_BYTE_COUNT / record_size * record_size; // in bytes
            number_of_pages = (int) ((kmers_num * record_size) / full_page_size + ((kmers_num * record_size) % full_page_size == 0 ? 0 : 1));
            suf_file = new RandomAccessFile(index_path + "sorted" + SUFFIX_FILE_EXTENTION, "r");
            suf_buff = new MappedByteBuffer[number_of_pages];
            for (i = 0; i < number_of_pages; ++i) {
                page_size = (int) (i == number_of_pages - 1 ? (kmers_num * record_size) % full_page_size : full_page_size); // in bytes
                page_size = page_size == 0 ? full_page_size : page_size; // in bytes
                suf_buff[i] = suf_file.getChannel().map(FileChannel.MapMode.READ_ONLY, ((long)full_page_size) * i + 4, page_size);
            }
        // mapping pointers file into the memory
            //System.out.println("Mapping pointers file to memory...                    ");
            record_size = poniter_length; // in bytes
            full_page_size = MAX_BYTE_COUNT / record_size * record_size; // in bytes
            number_of_pages = (int) ((kmers_num * record_size) / full_page_size + ((kmers_num * record_size) % full_page_size == 0 ? 0 : 1));
            ptr_file = new RandomAccessFile(index_path + "sorted" + POINTER_FILE_EXTENTION, "rw");
            ptr_buff = new MappedByteBuffer[number_of_pages];
            for (i = 0; i < number_of_pages; ++i) {
                page_size = (int) (i == number_of_pages - 1 ? (kmers_num * record_size) % full_page_size : full_page_size); // in bytes
                page_size = page_size == 0 ? full_page_size : page_size; // in bytes
                ptr_buff[i] = ptr_file.getChannel().map(FileChannel.MapMode.READ_WRITE, ((long)full_page_size) * i + 4, page_size);
                for (j = 0; j < page_size; ++j)
                    ptr_buff[i].put((byte)-1);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
    }
       
    /**
     * Updates an available index to include new genomes
     * 
     * @param index_path Path to the index database 
     * @param genomes_path_file A text file containing path to the new genomes
     * @param genomeDb The genome database
     * @param graphDb The graph database
     * @param previous_num_genomes Number of the genomes available in the index
     */
    public IndexDatabase(String index_path, String genomes_path_file, SequenceDatabase genomeDb, GraphDatabaseService graphDb, int previous_num_genomes) {
        int cores = Runtime.getRuntime().availableProcessors() / 2 + 1;
        int i, j, p, number_of_pages, record_size, page_size, full_page_size;
        long number_of_prefixes, longest_scaffold = 0;
        long c_index, p_index, l;
        Node node;
        IndexDatabase old_darabase;
        ResourceIterator<Node> nodes_iterator;
        // move current index files to directory old_darabase
        try{
            Files.move(Paths.get(index_path + "sorted" + PREFIX_FILE_EXTENTION), Paths.get(index_path + "old" + PREFIX_FILE_EXTENTION));
            Files.move(Paths.get(index_path + "sorted" + SUFFIX_FILE_EXTENTION), Paths.get(index_path + "old" + SUFFIX_FILE_EXTENTION));
            Files.move(Paths.get(index_path + "sorted" + POINTER_FILE_EXTENTION), Paths.get(index_path + "old" + POINTER_FILE_EXTENTION));
            Files.move(Paths.get(index_path + "sorted" + INFORMATION_FILE_EXTENTION), Paths.get(index_path + "old" + INFORMATION_FILE_EXTENTION));
        // load old_darabase
            old_darabase = new IndexDatabase(index_path, "old");
        // make new index for new genomes
            System.out.println("Running KMC with K = " + old_darabase.K + " ...                      ");
            executeCommand(KMC + " -cs127 -r -k" + old_darabase.K + " -t" + cores + " -ci1 -fm " + (genomeDb.num_genomes - previous_num_genomes > 1 ? "@" + 
                    genomes_path_file.trim() : genomeDb.genome_names[previous_num_genomes + 1]) + 
                    " " + index_path + "new " + index_path);
        // merge two indeces    
            executeCommand(KMC + "_tools union " + index_path + "old " + index_path + "new " + index_path + "sorted");
            pre_file = new RandomAccessFile(index_path + "sorted" + PREFIX_FILE_EXTENTION, "r");
            pre_file.seek(pre_file.length() - 8);
            header_pos = read_int(pre_file);
            pre_file.seek(pre_file.length() - 8 - header_pos);
        // read the merged index properties    
            K = read_int(pre_file);
            mode = read_int(pre_file);
            ctr_size = read_int(pre_file);
            pre_len = read_int(pre_file);
            min_count = read_int(pre_file);
            max_count = read_int(pre_file);
            kmers_num = read_long(pre_file);
            suf_len = K - pre_len;
            for (i = 1; i <= genomeDb.num_genomes; ++i)
                for (j = 1; j <= genomeDb.num_sequences[i]; ++j)
                    if (genomeDb.sequence_length[i][j] > longest_scaffold)
                        longest_scaffold = genomeDb.sequence_length[i][j];
            id_len = (int)Math.round(Math.ceil( Math.log(kmers_num) / Math.log(2) / 8));
            offset_len = (int)Math.round(Math.ceil( Math.log(longest_scaffold) / Math.log(2) / 8));
            poniter_length = 2 * id_len + offset_len + 1;
            write_info(index_path, "sorted");
            System.out.println("number of available kmers:\t" + old_darabase.kmers_num);
            System.out.println("number of new kmers:\t" + (kmers_num - old_darabase.kmers_num));
            System.out.println("Total number of kmers:\t" + kmers_num);
        // load the prefix file into the memory    
            //System.out.println("loading new prefixes in memory...                    ");
            number_of_prefixes = 1 << (2 * pre_len);
            prefix_ptr = new long[(int)number_of_prefixes];
            record_size = 8; // in bytes
            full_page_size = MAX_BYTE_COUNT / record_size * record_size; // in bytes
            number_of_pages = (int) ((number_of_prefixes * record_size) / full_page_size + ((number_of_prefixes * record_size) % full_page_size == 0 ? 0 : 1));
            pre_buff = new  MappedByteBuffer[number_of_pages];
            for (i = 0; i < number_of_pages; ++i) {
                page_size = (int) (i == number_of_pages - 1 ? (number_of_prefixes * record_size) % full_page_size : full_page_size); // in bytes
                page_size = page_size == 0 ? full_page_size : page_size; // in bytes
                pre_buff[i] = pre_file.getChannel().map(FileChannel.MapMode.READ_ONLY, ((long)full_page_size) * i + 4, page_size);
                for (j = 0; j < page_size / 8; ++j) 
                    prefix_ptr[(int)(((long)full_page_size) * i / 8 + j)] = read_prefix(pre_buff[i]);
                pre_buff[i] = null;
            }
            pre_buff = null;
            pre_file.close();
        // mapping suffix file into the memory
            //System.out.println("Mapping new suffix file to memory...                    ");
            record_size = ctr_size + suf_len / 4;
            full_page_size = MAX_BYTE_COUNT / record_size * record_size; // in bytes
            number_of_pages = (int) ((kmers_num * record_size) / full_page_size + ((kmers_num * record_size) % full_page_size == 0 ? 0 : 1));
            suf_file = new RandomAccessFile(index_path + "sorted" + SUFFIX_FILE_EXTENTION, "r");
            suf_buff = new MappedByteBuffer[number_of_pages];
            for (i = 0; i < number_of_pages; ++i) {
                page_size = (int) (i == number_of_pages - 1 ? (kmers_num * record_size) % full_page_size : full_page_size); // in bytes
                page_size = page_size == 0 ? full_page_size : page_size; // in bytes
                suf_buff[i] = suf_file.getChannel().map(FileChannel.MapMode.READ_ONLY, ((long)full_page_size) * i + 4, page_size);
            }
        // mapping pointers file into the memory
            //System.out.println("Mapping new pointer files to memory...                    ");
            record_size = poniter_length; // in bytes
            full_page_size = MAX_BYTE_COUNT / record_size * record_size; // in bytes
            number_of_pages = (int) ((kmers_num * record_size) / full_page_size + ((kmers_num * record_size) % full_page_size == 0 ? 0 : 1));
            ptr_file = new RandomAccessFile(index_path + "sorted" + POINTER_FILE_EXTENTION, "rw");
            ptr_buff = new MappedByteBuffer[number_of_pages];
            for (i = 0; i < number_of_pages; ++i) {
                page_size = (int) (i == number_of_pages - 1 ? (kmers_num * record_size) % full_page_size : full_page_size); // in bytes
                page_size = page_size == 0 ? full_page_size : page_size; // in bytes
                ptr_buff[i] = ptr_file.getChannel().map(FileChannel.MapMode.READ_WRITE, ((long)full_page_size) * i + 4, page_size);
                for (j = 0; j < page_size; ++j)
                    ptr_buff[i].put((byte)-1);
            }
        // adjusting available pointers
            System.out.println("Updating kmer index...                    ");
            kmer old_kmer = new kmer(K, old_darabase.pre_len);
            kmer new_kmer = new kmer(K, pre_len);
            IndexScanner old_scanner = new IndexScanner(old_darabase);
            IndexScanner new_scanner = new IndexScanner(this);
            IndexPointer ptr = new IndexPointer();
            try(Transaction tx = graphDb.beginTx()){
                nodes_iterator = graphDb.findNodes(nucleotide_label);
                for(;nodes_iterator.hasNext();){
                    node = nodes_iterator.next();
                    if (!node.hasLabel(degenerate_label)){
                        l = (long)node.getProperty("first_kmer");
                        old_scanner.get_kmer(old_kmer, l);
                        adjust_fwd_kmer(new_kmer, old_kmer);
                        p_index = new_scanner.find(new_kmer);
                        old_scanner.get_pointer(ptr,l);
                        new_scanner.put_pointer(ptr,p_index);
                        node.setProperty("first_kmer", p_index);
                        for(l = old_scanner.get_next_index(l); l != -1L; l = old_scanner.get_next_index(l)){
                            old_scanner.get_kmer(old_kmer, l);
                            adjust_fwd_kmer(new_kmer, old_kmer);
                            c_index = new_scanner.find(new_kmer);
                            old_scanner.get_pointer(ptr,l);
                            new_scanner.put_pointer(ptr,c_index);
                            new_scanner.put_next_index(c_index,p_index);
                            p_index = c_index;
                        }
                        new_scanner.put_next_index(-1L,p_index);
                        node.setProperty("last_kmer", p_index);
                    }
                }
                nodes_iterator.close();
                tx.success();
            }
            old_darabase.close();
            Files.delete(Paths.get(index_path + "old" + PREFIX_FILE_EXTENTION));
            Files.delete(Paths.get(index_path + "old" + SUFFIX_FILE_EXTENTION));
            Files.delete(Paths.get(index_path + "old" + INFORMATION_FILE_EXTENTION));
            Files.delete(Paths.get(index_path + "old" + POINTER_FILE_EXTENTION));
            Files.delete(Paths.get(index_path + "new" + PREFIX_FILE_EXTENTION));
            Files.delete(Paths.get(index_path + "new" + SUFFIX_FILE_EXTENTION));
        } catch (IOException e) {
            System.out.println("Failed to make index!\n" + e.getMessage());
            System.exit(1);
        }
    }

    public void find_k(String index_path, String genomes_path_file, SequenceDatabase genomeDb, int k){
        double p = 0;
        String output;
        try{
        Runtime.getRuntime().exec(KMC); // to check if kmc is reachable
        if (k > -1){ // K is not given by the user, then calculate the optimal K
            K = k;
            if (K % 2 == 0) // Even values are problamatic to localization process 
                K += 1;
            executeCommand(KMC + " -cs127 -k" + K + " -t" + cores + " -ci1 -fm " + 
                    (genomeDb.num_genomes > 1 ? "@" + genomes_path_file.trim() : genomeDb.genome_names[1]) + " " + index_path + "/kmers " + index_path);            
            if (new File(index_path + "/kmers.kmc_pre").exists() && new File(index_path + "/kmers.kmc_suf").exists()) {
                output = executeCommand(KMC + "_tools sort " + index_path + "/kmers " + index_path + "/sorted");
            // Small databases are usually sorted already    
                if (output.startsWith("This database contains sorted k-mers already!")) {
                    new File(index_path + "/kmers.kmc_pre").renameTo(new File(index_path + "/sorted.kmc_pre"));
                    new File(index_path + "/kmers.kmc_suf").renameTo(new File(index_path + "/sorted.kmc_suf"));
                } else {
                    new File(index_path + "/kmers.kmc_pre").delete();
                    new File(index_path + "/kmers.kmc_suf").delete();
                }            
            } else {
                System.out.println("No kmc index found in " + index_path);
                System.exit(1);
            }
        } else {
            K = Math.round((float)(Math.log(genomeDb.number_of_bytes * 200)/Math.log(4)));
            K += (K % 2 == 0 ? 1 : 0);
            do{
                K -= 2;
        System.out.println("K = " +  K);
                executeCommand(KMC + " -cs127 -k" + K + " -t" + cores + " -ci1 -fm " +            
                        (genomeDb.num_genomes > 1 ? "@" + genomes_path_file.trim() : genomeDb.genome_names[1]) + " " + index_path + "/kmers " + index_path);
                if (new File(index_path + "/kmers.kmc_pre").exists() && new File(index_path + "/kmers.kmc_suf").exists()) {
                    output = executeCommand(KMC + "_tools sort " + index_path + "/kmers " + index_path + "/sorted");
                // Small databases are usually sorted already    
                    if (output.startsWith("This database contains sorted k-mers already!")) {
                        new File(index_path + "/kmers.kmc_pre").renameTo(new File(index_path + "/sorted.kmc_pre"));
                        new File(index_path + "/kmers.kmc_suf").renameTo(new File(index_path + "/sorted.kmc_suf"));
                    } else {
                        new File(index_path + "/kmers.kmc_pre").delete();
                        new File(index_path + "/kmers.kmc_suf").delete();
                    }            
                    pre_file = new RandomAccessFile(index_path + "/sorted.kmc_pre", "r");
                    pre_file.seek(pre_file.length() - 8);
                    header_pos = read_int(pre_file);
                    pre_file.seek(pre_file.length() - (8 + header_pos));
                    // read the header of the index    
                    K = read_int(pre_file);
                    //  k-mer length
                    mode = read_int(pre_file);
                    // 0 (occurrence count) or 1 (counting according to Quake quality)
                    ctr_size = read_int(pre_file);
                    // counter field size: for mode 0 it is 1, 2, 3, or 4; for mode 1 it is always 4
                    pre_len = read_int(pre_file);
                    // prefix length such that suffix length is divisible by 4.
                    // Max prefix length is limited up to 15 in KMC (exact value is calculated in such a way that summary size of kmc_pre and kmc_suf should be minimal)
                    min_count = read_int(pre_file);
                    // minimum number of k-mer occurrences to be written in the database    
                    max_count = read_int(pre_file);
                    // maximum number of k-mer occurrences to be written in the database
                    kmers_num = read_long(pre_file);  
                    p = (double)kmers_num / (1l  << (2 * K)) * 2;// 2n / 4 ^ K
                    pre_file.close();
                } else {
                    System.out.println("No kmc index found in " + index_path);
                    System.exit(1);
                }
            } while (p < 0.01);
            K += 2;
            executeCommand(KMC + " -cs127 -k" + K + " -t" + cores + " -ci1 -fm " +            
                    (genomeDb.num_genomes > 1 ? "@" + genomes_path_file.trim() : genomeDb.genome_names[1]) + " " + index_path + "/kmers " + index_path);
            if (new File(index_path + "/kmers.kmc_pre").exists() && new File(index_path + "/kmers.kmc_suf").exists()) {
                output = executeCommand(KMC + "_tools sort " + index_path + "/kmers " + index_path + "/sorted");
            // Small databases are usually sorted already    
                if (output.startsWith("This database contains sorted k-mers already!")) {
                    new File(index_path + "/kmers.kmc_pre").renameTo(new File(index_path + "/sorted.kmc_pre"));
                    new File(index_path + "/kmers.kmc_suf").renameTo(new File(index_path + "/sorted.kmc_suf"));
                } else {
                    new File(index_path + "/kmers.kmc_pre").delete();
                    new File(index_path + "/kmers.kmc_suf").delete();
                } 
            } else {
                System.out.println("No kmc index found in " + index_path);
                System.exit(1);
            }
        }
        System.out.println("K = " +  K);
        } catch (IOException ex){
            System.err.println("Failed to find K!");
        }
    }

    public void write_info(String index_path, String index_name) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(index_path + index_name + INFORMATION_FILE_EXTENTION));
            out.write("K:" + K + "\n");
            out.write("mode:" + mode + "\n");
            out.write("ctr_size:" + ctr_size + "\n");
            out.write("pre_len:" + pre_len + "\n");
            out.write("min_count:" + min_count + "\n");
            out.write("max_count:" + max_count + "\n");
            out.write("kmers_num:" + kmers_num + "\n");
            out.write("suf_len:" + suf_len + "\n");
            out.write("id_len:" + id_len + "\n");
            out.write("offset_len:" + offset_len + "\n");
            out.write("poniter_length:" + poniter_length + "\n");
            out.close();
        } catch (IOException ioe) {
            System.out.println(ioe.getMessage());
            System.exit(1);
        }
    }

    /**
     * Gives the length of index. 
     * 
     * @return Number of kmers
     */
    public long length() {
        return kmers_num;
    }

    /**
     * Gives the prefix length of kmers.
     * 
     * @return Prefix length of kmers
     */
    public int get_pre_len() {
        return pre_len;
    }

    /**
     * Gives the suffix length of kmers.
     * 
     * @return suffix length of kmers
     */
    public int get_suf_len() {
        return suf_len;
    }

    /**
     * Closes the index database
     */
    public void close() {
        try {
            ptr_file.close();
            suf_file.close();
            int k;
            for (k = 0; k < ptr_buff.length; ++k) {
                ptr_buff[k] = null;
            }
            for (k = 0; k < suf_buff.length; ++k) {
                suf_buff[k] = null;
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
    }
    
    /**
     * Reads a long integer from the memory mapped buffer.
     * 
     * @param buff The file associated to the byte stream
     * @return The long integer value
     * @throws IOException 
     */
    private long read_prefix(MappedByteBuffer buff){
        int i;
        long number = 0;
        for (i = 0; i < 8; ++i)
            number += ((long)(buff.get() & 0x00FF)) << (8 * i);
        return number;
    }    

    /**
     * Reads an integer from the byte stream.
     * 
     * @param file The file associated to the byte stream
     * @return The integer value
     * @throws IOException 
     */
    private int read_int(RandomAccessFile file){
        int number = 0;
        for (int i = 0; i < 4; ++i) {
            try {
                number += ((int)(file.readByte() & 0x00FF) << (8 * i));
            } catch (IOException ex) {
                Logger.getLogger(IndexDatabase.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return number;
    }

    /**
     * Reads a long integer from the byte stream.
     * 
     * @param file The file associated to the byte stream
     * @return The long integer value
     * @throws IOException 
     */
    private long read_long(RandomAccessFile file) throws IOException {
        long number = 0;
        for (int i = 0; i < 8; ++i) {
            number += ((long)(file.readByte() & 0x00FF) << (8 * i));
        }
        return number;
    }

    /**
     * Converts the byte array to an integer value.
     * @param b The byte array
     * @param offset Index of the first byte in the array
     * @return The integer
     */
     public int b2i(byte[] b, int offset) {
        return (b[offset + 3] & 0xFF)
                | ((b[offset + 2] & 0xFF) << 8)
                | ((b[offset + 1] & 0xFF) << 16)
                | ((b[offset] & 0xFF) << 24);
    }

    /**
     * Converts an integer to a byte array.
     * @param value The integer value
     * @return The byte array
     */ 
    public static void i2b(int value, byte[] data) {
        data[0] = (byte) ((value >> 24) & 0xFF);
        data[1] = (byte) ((value >> 16) & 0xFF);
        data[2] = (byte) ((value >> 8) & 0xFF);
        data[3] = (byte) (value & 0xFF);
    }

    /**
     * Converts the byte array to a long integer value.
     * @param b The byte array
     * @param offset Index of the first byte in the array
     * @return The long integer
     */
    public static long b2l(byte[] b, int offset) {
        return (  (long)  b[offset + 7] & 0x00FF)
                | (long)((b[offset + 6] & 0x00FF) << 8)
                | (long)((b[offset + 5] & 0x00FF) << 16)
                | (long)((b[offset + 4] & 0x00FF) << 24)
                | (long)((b[offset + 3] & 0x00FF) << 32)
                | (long)((b[offset + 2] & 0x00FF) << 40)
                | (long)((b[offset + 1] & 0x00FF) << 48)
                | (long)((b[offset] & 0x00FF) << 56);
    }

    /**
     * Converts a long integer to a byte array.
     * @param value The long value
     * @return The byte array
     */
    public static void l2b(long value, byte[] data) {
        data[0] = (byte) ((value >> 56) & 0x00FF);
        data[1] = (byte) ((value >> 48) & 0x00FF);
        data[2] = (byte) ((value >> 40) & 0x00FF);
        data[3] = (byte) ((value >> 32) & 0x00FF);
        data[4] = (byte) ((value >> 24) & 0x00FF);
        data[5] = (byte) ((value >> 16) & 0x00FF);
        data[6] = (byte) ((value >> 8) & 0x00FF);
        data[7] = (byte) (value & 0x00FF);
    }
    
    public int get_K(){
        return K;
    }
    
    public int get_pointer_length(){
        return poniter_length;
    }
}
