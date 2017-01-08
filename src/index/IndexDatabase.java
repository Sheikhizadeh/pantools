/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package index;

import genome.SequenceDatabase;
import java.io.File;
import java.io.IOException;
import org.neo4j.graphdb.DynamicLabel;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import static pantools.Pantools.executeCommand;

/**
 * Implements all the functionality to work with a KMC-based kmer index database. 
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class IndexDatabase {

    private kmer key;
    private long[] prefix_ptr;
    private int ptr_parts_num;
    private long[] ptr_parts_size;
    private int suf_parts_num;
    private long[] suf_parts_size;
    private int MAX_BYTE_COUNT = 2100000000; // The maximum number of bytes addressable by an integer
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
    public static int POINTER_LENGTH = 21; // The length of a poniter in bytes
    private RandomAccessFile suf_file;
    private RandomAccessFile pre_file;
    private MappedByteBuffer[] suf_buff;
    private RandomAccessFile ptr_file;
    private MappedByteBuffer[] ptr_buff;

    /**
     * Mounts an available index database to the index database object
     * 
     * @param index_path Path to the index database 
     */
    public IndexDatabase(String index_path) {
        int k;
        long p;
        try {
            pre_file = new RandomAccessFile(index_path + "/sorted.kmc_pre", "r");
            pre_file.seek(pre_file.length() - 8);
            header_pos = read_int(pre_file);
            pre_file.seek(pre_file.length() - 8 - header_pos);
        // read the index properties    
            K = read_int(pre_file);
            mode = read_int(pre_file);
            ctr_size = read_int(pre_file);
            pre_len = read_int(pre_file);
            min_count = read_int(pre_file);
            max_count = read_int(pre_file);
            kmers_num = read_long(pre_file);
            suf_len = K - pre_len;
            key = new kmer(K, pre_len, suf_len);
            System.out.println("Indexing " + kmers_num + " kmers...                    ");
        // load the prefix file into the memory    
            pre_file.seek(4);
            int q, len = 1 << (2 * pre_len);
            prefix_ptr = new long[len];
            MappedByteBuffer pre_buff;
            for (q = 0, p = 0; p < 8; ++p) {
                pre_buff = pre_file.getChannel().map(FileChannel.MapMode.READ_ONLY, 4 + p * len, len);
                for (k = 0; k < len / 8; ++k, ++q) {
                    prefix_ptr[q] = read_long(pre_buff);
                }
                pre_buff = null;
            }
            pre_file.close();
        // mapping suffix file into the memory
            suf_rec_size = ctr_size + suf_len / 4;
            MAX_BYTE_COUNT = MAX_BYTE_COUNT / suf_rec_size * suf_rec_size;
            suf_parts_num = (int) ((kmers_num * suf_rec_size) % MAX_BYTE_COUNT == 0 ? (kmers_num * suf_rec_size) / MAX_BYTE_COUNT : (kmers_num * suf_rec_size) / MAX_BYTE_COUNT + 1);
            suf_parts_size = new long[suf_parts_num];
            suf_file = new RandomAccessFile(index_path + "/sorted.kmc_suf", "r");
            suf_buff = new MappedByteBuffer[suf_parts_num];
            for (k = 0; k < suf_parts_num; ++k) {
                suf_parts_size[k] = (int) (k == suf_parts_num - 1 ? (kmers_num * suf_rec_size) % MAX_BYTE_COUNT : MAX_BYTE_COUNT);
                suf_buff[k] = suf_file.getChannel().map(FileChannel.MapMode.READ_ONLY, 4 + k * suf_parts_size[0], suf_parts_size[k]);
            }
        // mapping pointers file into the memory
            MAX_BYTE_COUNT = MAX_BYTE_COUNT / POINTER_LENGTH * POINTER_LENGTH;
            ptr_parts_num = (int) ((kmers_num * POINTER_LENGTH) % MAX_BYTE_COUNT == 0 ? (kmers_num * POINTER_LENGTH) / MAX_BYTE_COUNT : (kmers_num * POINTER_LENGTH) / MAX_BYTE_COUNT + 1);
            ptr_parts_size = new long[ptr_parts_num];
            ptr_file = new RandomAccessFile(index_path + "/pointers.db", "rw");
            ptr_buff = new MappedByteBuffer[ptr_parts_num];
            for (k = 0; k < ptr_parts_num; ++k) {
                ptr_parts_size[k] = (int) (k == ptr_parts_num - 1 ? (kmers_num * POINTER_LENGTH) % MAX_BYTE_COUNT : MAX_BYTE_COUNT);
                ptr_buff[k] = ptr_file.getChannel().map(FileChannel.MapMode.READ_WRITE, k * ptr_parts_size[0], ptr_parts_size[k]);
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
    public IndexDatabase(String index_path, String genomes_path_file, int K_size, SequenceDatabase genomeDb) {
        int cores = Runtime.getRuntime().availableProcessors() / 2 + 1;
        int k;
        long p;
        IndexPointer null_pointer = new IndexPointer();
        try {
            Files.createDirectory(Paths.get(index_path));
            System.out.println("Running KMC2...                      ");
            executeCommand("kmc -r -k" + K_size + " -t" + cores + " -m" + 
                (Runtime.getRuntime().maxMemory() / 1073741824L) + " -ci1 -fm " + 
                (genomeDb.num_genomes > 1 ? "@" + genomes_path_file.trim() : genomeDb.genome_names[1]) + " " + index_path + "/kmers " + index_path);
            System.out.println("Sorting kmers...                      ");
            String output = executeCommand("kmc_tools sort " + index_path + "/kmers " + index_path + "/sorted");
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
            pre_file.seek(pre_file.length() - 8 - header_pos);
            // read the index properties    
            K = read_int(pre_file);
            mode = read_int(pre_file);
            ctr_size = read_int(pre_file);
            pre_len = read_int(pre_file);
            min_count = read_int(pre_file);
            max_count = read_int(pre_file);
            kmers_num = read_long(pre_file);
            suf_len = K - pre_len;
            key = new kmer(K, pre_len, suf_len);
            System.out.println("Indexing " + kmers_num + " kmers...                    ");
            // load the prefix file into the memory    
            pre_file.seek(4);
            int q, len = 1 << (2 * pre_len);
            prefix_ptr = new long[len];
            MappedByteBuffer pre_buff;
            for (q = 0, p = 0; p < 8; ++p) {
                pre_buff = pre_file.getChannel().map(FileChannel.MapMode.READ_ONLY, 4 + p * len, len);
                for (k = 0; k < len / 8; ++k, ++q) {
                    prefix_ptr[q] = read_long(pre_buff);
                }
                pre_buff = null;
            }
            pre_file.close();

            // mapping suffix file into the memory
            suf_rec_size = ctr_size + suf_len / 4;
            MAX_BYTE_COUNT = MAX_BYTE_COUNT / suf_rec_size * suf_rec_size;
            suf_parts_num = (int) ((kmers_num * suf_rec_size) % MAX_BYTE_COUNT == 0 ? (kmers_num * suf_rec_size) / MAX_BYTE_COUNT : (kmers_num * suf_rec_size) / MAX_BYTE_COUNT + 1);
            suf_parts_size = new long[suf_parts_num];
            suf_file = new RandomAccessFile(index_path + "/sorted.kmc_suf", "r");
            suf_buff = new MappedByteBuffer[suf_parts_num];
            for (k = 0; k < suf_parts_num; ++k) {
                suf_parts_size[k] = (int) (k == suf_parts_num - 1 ? (kmers_num * suf_rec_size) % MAX_BYTE_COUNT : MAX_BYTE_COUNT);
                suf_buff[k] = suf_file.getChannel().map(FileChannel.MapMode.READ_ONLY, 4 + k * suf_parts_size[0], suf_parts_size[k]);
            }
            // mapping pointers file into the memory
            MAX_BYTE_COUNT = MAX_BYTE_COUNT / POINTER_LENGTH * POINTER_LENGTH;
            ptr_parts_num = (int) ((kmers_num * POINTER_LENGTH) % MAX_BYTE_COUNT == 0 ? (kmers_num * POINTER_LENGTH) / MAX_BYTE_COUNT : (kmers_num * POINTER_LENGTH) / MAX_BYTE_COUNT + 1);
            ptr_parts_size = new long[ptr_parts_num];
            ptr_file = new RandomAccessFile(index_path + "/pointers.db", "rw");
            ptr_buff = new MappedByteBuffer[ptr_parts_num];
            for (k = 0; k < ptr_parts_num; ++k) {
                ptr_parts_size[k] = (int) (k == ptr_parts_num - 1 ? (kmers_num * POINTER_LENGTH) % MAX_BYTE_COUNT : MAX_BYTE_COUNT);
                ptr_buff[k] = ptr_file.getChannel().map(FileChannel.MapMode.READ_WRITE, k * ptr_parts_size[0], ptr_parts_size[k]);
            }
            for (p = 0; p < kmers_num; ++p) {
                put_pointer(null_pointer, p);
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
        int i, k, p, trsc, seq_nodes;
        long c_index, p_index, l, new_kmers_num;
        IndexPointer null_pointer = new IndexPointer();
        Node node;
        ResourceIterator<Node> nodes;
        Path old_index_folder;
        // move current index files to directory old_index
        try {
            old_index_folder = Files.createDirectory(Paths.get(index_path+"/old_index"));
            Files.move(Paths.get(index_path + "/sorted.kmc_pre"), Paths.get(index_path + "/old_index/sorted.kmc_pre"));
            Files.move(Paths.get(index_path + "/sorted.kmc_suf"), Paths.get(index_path + "/old_index/sorted.kmc_suf"));
            Files.move(Paths.get(index_path + "/pointers.db"), Paths.get(index_path + "/old_index/pointers.db"));
        // load old_index
            IndexDatabase old_index = new IndexDatabase(index_path + "/old_index");
            K = old_index.K;
        // make new index for new genomes
            System.out.println("Running KMC2...                      ");
            executeCommand("kmc -r -k" + K + " -t" + cores + " -m" + (Runtime.getRuntime().maxMemory() / 1073741824L) + 
                    " -ci1 -fm " + (genomeDb.num_genomes - previous_num_genomes > 1 ? "@" + 
                    genomes_path_file.trim() : genomeDb.genome_names[previous_num_genomes + 1]) + 
                    " " + index_path + "/new_kmers " + index_path);
        // merge two indeces    
            executeCommand("kmc_tools union " + index_path + "/old_index/sorted " + index_path + "/new_kmers " + index_path + "/sorted");
            pre_file = new RandomAccessFile(index_path + "/sorted.kmc_pre", "r");
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
            new_kmers_num = read_long(pre_file);
            kmers_num += new_kmers_num;
            System.out.println(new_kmers_num + " new kmers generated.                    ");
        // load the prefix file into the memory    
            pre_file.seek(4);
            int q, len = 1 << (2 * pre_len);
            prefix_ptr = new long[len];
            MappedByteBuffer pre_buff;
            for (q = 0, p = 0; p < 8; ++p) {
                pre_buff = pre_file.getChannel().map(FileChannel.MapMode.READ_ONLY, 4 + p * len, len);
                for (k = 0; k < len / 8; ++k, ++q) {
                    prefix_ptr[q] = read_long(pre_buff);
                    //System.out.println(q+" "+prefix_ptr[q]);
                }
                pre_buff = null;
            }
            pre_file.close();
        // mapping suffix file into the memory    
            suf_len = K - pre_len;
            suf_rec_size = ctr_size + suf_len / 4;
            MAX_BYTE_COUNT = MAX_BYTE_COUNT / suf_rec_size * suf_rec_size;
            suf_parts_num = (int) ((kmers_num * suf_rec_size) % MAX_BYTE_COUNT == 0 ? (kmers_num * suf_rec_size) / MAX_BYTE_COUNT : (kmers_num * suf_rec_size) / MAX_BYTE_COUNT + 1);
            suf_parts_size = new long[suf_parts_num];
            suf_file = new RandomAccessFile(index_path + "/sorted.kmc_suf", "r");
            suf_buff = new MappedByteBuffer[suf_parts_num];
            key = new kmer(K, pre_len, suf_len);
            for (k = 0; k < suf_parts_num; ++k) {
                suf_parts_size[k] = (int) (k == suf_parts_num - 1 ? (kmers_num * suf_rec_size) % MAX_BYTE_COUNT : MAX_BYTE_COUNT);
                suf_buff[k] = suf_file.getChannel().map(FileChannel.MapMode.READ_ONLY, 4 + k * suf_parts_size[0], suf_parts_size[k]);
            }
        // mapping pointers file into the memory    
            MAX_BYTE_COUNT = MAX_BYTE_COUNT / POINTER_LENGTH * POINTER_LENGTH;
            ptr_parts_num = (int) ((kmers_num * POINTER_LENGTH) % MAX_BYTE_COUNT == 0 ? (kmers_num * POINTER_LENGTH) / MAX_BYTE_COUNT : (kmers_num * POINTER_LENGTH) / MAX_BYTE_COUNT + 1);
            ptr_parts_size = new long[ptr_parts_num];
            ptr_file = new RandomAccessFile(index_path + "/pointers.db", "rw");
            ptr_buff = new MappedByteBuffer[ptr_parts_num];
            for (k = 0; k < ptr_parts_num; ++k) {
                ptr_parts_size[k] = (int) (k == ptr_parts_num - 1 ? (kmers_num * POINTER_LENGTH) % MAX_BYTE_COUNT : MAX_BYTE_COUNT);
                ptr_buff[k] = ptr_file.getChannel().map(FileChannel.MapMode.READ_WRITE, k * ptr_parts_size[0], ptr_parts_size[k]);
            }
            for (p = 0; p < kmers_num; ++p) {
                put_pointer(null_pointer, p);
            }
        // adjusting available pointers
            try (Transaction tx = graphDb.beginTx()) {
                nodes = graphDb.findNodes(DynamicLabel.label("node"));
                seq_nodes = (int) graphDb.findNodes(DynamicLabel.label("pangenome")).next().getProperty("num_nodes");
                tx.success();
            }
            IndexPointer ptr = new IndexPointer();
            System.out.println("Updating kmer index...                    ");
            for (i = 0; nodes.hasNext();) {
                try (Transaction tx = graphDb.beginTx()) {
                    for (trsc = 0; trsc < 10000 && nodes.hasNext(); ++trsc, ++i) {
                        node = nodes.next();
                        l = (long) node.getProperty("first_kmer");
                        p_index = find(old_index.get_kmer(l));
                        old_index.get_pointer(ptr, l);
                        put_pointer(ptr, p_index);
                        node.setProperty("first_kmer", p_index);
                        for (l = old_index.get_next_index(l); l != -1L; l = old_index.get_next_index(l)) {
                            c_index = find(old_index.get_kmer(l));
                            old_index.get_pointer(ptr, l);
                            put_pointer(ptr, c_index);
                            put_next_index(c_index, p_index);
                            p_index = c_index;
                        }
                        put_next_index(-1L, p_index);
                        node.setProperty("last_kmer", p_index);
                        if (i % (seq_nodes / 100 + 1) == 0) {
                            System.out.print((long) i * 100 / seq_nodes + 1 + "%    \r");
                        }
                    }
                    tx.success();
                }
            }
            old_index.close();
            Files.delete(Paths.get(index_path + "/old_index/sorted.kmc_suf"));
            Files.delete(Paths.get(index_path + "/new_kmers.kmc_pre"));
            Files.delete(Paths.get(index_path + "/new_kmers.kmc_suf"));
            Files.delete(Paths.get(index_path + "/old_index/pointers.db"));
            Files.delete(Paths.get(index_path + "/old_index/sorted.kmc_pre"));
            //Files.delete(old_index_folder);
        } catch (IOException e) {
            System.out.println(e.getMessage() + "\nFailed to make index!");
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
            for (k = 0; k < ptr_parts_num; ++k) {
                ptr_buff[k] = null;
            }
            for (k = 0; k < suf_parts_num; ++k) {
                suf_buff[k] = null;
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
    }

    /**
     * Reads an integer from the byte stream.
     * 
     * @param file The file associated to the byte stream
     * @return The integer value
     * @throws IOException 
     */
    private int read_int(RandomAccessFile file) throws IOException {
        int number = 0;
        for (int i = 0; i < 4; ++i) {
            number += (file.readByte() & 0xFF) << (8 * i);
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
            number += ((long) (file.readByte() & 0xFF)) << (8 * i);
        }
        return number;
    }

    /**
     * Reads a long integer from the memory mapped buffer.
     * 
     * @param buff The file associated to the byte stream
     * @return The long integer value
     * @throws IOException 
     */
    private long read_long(MappedByteBuffer buff) throws IOException {
        long number = 0;
        for (int i = 0; i < 8; ++i) {
            number += ((long) (buff.get() & 0xFF)) << (8 * i);
        }
        return number;
    }

    /**
     * Reads a kmer from th index database.
     * 
     * @param number The number of kmer in the index
     * @return The kmer object
     */
    public kmer get_kmer(long number) {
        kmer k_mer = new kmer(K, pre_len, suf_len);
        boolean found = false;
        int low = 0, high = (1 << (2 * pre_len)) - 1, mid = 0;
        while (low <= high && !found) {
            mid = (low + high) / 2;
            if (number < prefix_ptr[mid]) {
                high = mid - 1;
            } else if (number > prefix_ptr[mid]) {
                low = mid + 1;
            } else {
                found = true;
            }
        }
        if (!found) {
            mid = high;
        } else {
            for (; mid < prefix_ptr.length - 1 && prefix_ptr[mid + 1] == prefix_ptr[mid]; ++mid);
        }
        k_mer.prefix = mid;
        suf_buff[(int) (number * suf_rec_size / suf_parts_size[0])].position((int) (number * suf_rec_size % suf_parts_size[0]));
        suf_buff[(int) (number * suf_rec_size / suf_parts_size[0])].get(k_mer.suffix, 0, suf_len / 4);
        return k_mer;
    }

    /**
     * Reads a pointer from the kmer index.
     * @param poniter The pointer
     * @param number Number of the pointer.
     */
    public void get_pointer(IndexPointer poniter, long number) {
        poniter.node_id = ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].getLong((int) (number * POINTER_LENGTH % ptr_parts_size[0]));
        poniter.canonical = ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].get((int) (number * POINTER_LENGTH % ptr_parts_size[0] + 8)) == 0;
        poniter.position = ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].getInt((int) (number * POINTER_LENGTH % ptr_parts_size[0] + 9));
        poniter.next_index = ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].getLong((int) (number * POINTER_LENGTH % ptr_parts_size[0] + 13));
    }

    /**
     * Writes a pointer in the index database
     * @param poniter The pointer
     * @param number Number of the pointer.
     */
    public void put_pointer(IndexPointer poniter, long number) {
        ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].putLong((int) (number * POINTER_LENGTH % ptr_parts_size[0]), poniter.node_id);
        ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].put((int) (number * POINTER_LENGTH % ptr_parts_size[0] + 8), (byte) (poniter.canonical ? 0 : 1));
        ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].putInt((int) (number * POINTER_LENGTH % ptr_parts_size[0] + 9), poniter.position);
        ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].putLong((int) (number * POINTER_LENGTH % ptr_parts_size[0] + 13), poniter.next_index);
    }

    /**
     * Reads the node id of a kmer from the index database.
     * @param number Number of the kmer in the index
     * @return The id of node
     */
    public long get_node_id(long number) {
        return ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].getLong((int) (number * POINTER_LENGTH % ptr_parts_size[0]));
    }

    /**
     * Writes the node id of a kmer into the index database.
     * @param node_id The node id to be written
     * @param number Number of the kmer in the index
     */
    public void put_node_id(long node_id, long number) {
        ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].putLong((int) (number * POINTER_LENGTH % ptr_parts_size[0]), node_id);
    }

    /**
     * Looks if a kmer has been canonical at the first visit.
     * @param number Number of the kmer in the index
     * @return The canonical status
     */
    public boolean get_canonical(long number) {
        return ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].get((int) (number * POINTER_LENGTH % ptr_parts_size[0] + 8)) == 0;
    }

    /**
     * Writes the canonical status of a kmer at the first visit.
     * @param canonical The canonical status to be written ( 0 if it was canonical, 1 otherwise )
     * @param number Number of the kmer in the index
     */
    public void put_canonical(boolean canonical, long number) {
        ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].put((int) (number * POINTER_LENGTH % ptr_parts_size[0] + 8), (byte) (canonical ? 0 : 1));
    }

    /**
     * Reads the position of a kmer in the node where it occurs.
     * @param number Number of the kmer in the index
     * @return The position in the node
     */
    public int get_position(long number) {
        return ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].getInt((int) (number * POINTER_LENGTH % ptr_parts_size[0] + 9));
    }

    /**
     * Writes the position of a kmer in the node where it occurs.
     * @param position The node id to be written
     * @param number Number of the kmer in the index
     */
    public void put_position(int position, long number) {
        ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].putInt((int) (number * POINTER_LENGTH % ptr_parts_size[0] + 9), position);
    }

    /**
     * Reads the number of the next in the index database.
     * @param number Number of the kmer in the index
     * @return The number of the next kmer
     */
    public long get_next_index(long number) {
        return ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].getLong((int) (number * POINTER_LENGTH % ptr_parts_size[0] + 13));
    }

    /**
     * Writes the number of the next in the index database.
     * @param next_index The node id to be written
     * @param number Number of the kmer in the index
     */
    public void put_next_index(long next_index, long number) {
        ptr_buff[(int) (number * POINTER_LENGTH / ptr_parts_size[0])].putLong((int) (number * POINTER_LENGTH % ptr_parts_size[0] + 13), next_index);
    }

    /**
     * Finds the number of a kmer in the database.
     * 
     * @param k_mer The kmer object
     * @return The number of kmer
     */
    public long find(kmer k_mer) {
        long low, mid, high;
        int j, comp;
        low = prefix_ptr[k_mer.prefix];
        if (k_mer.prefix == prefix_ptr.length - 1) {
            high = kmers_num - 1;
        } else {
            high = prefix_ptr[k_mer.prefix + 1] - 1;
        }
        while (low <= high) {
            mid = (low + high) / 2;
            j = (int) (mid * suf_rec_size / suf_parts_size[0]);
            suf_buff[j].position((int) (mid * suf_rec_size % suf_parts_size[0]));
            suf_buff[j].get(key.suffix, 0, suf_len / 4);
            comp = k_mer.compare_suffix(key);
            if (comp == -1) {
                high = mid - 1;
            } else if (comp == 1) {
                low = mid + 1;
            } else {
                return mid;
            }
        }
        return -1L; // not found
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
    public byte[] i2b(int value) {
        byte[] data = new byte[4];
        data[0] = (byte) ((value >> 24) & 0xFF);
        data[1] = (byte) ((value >> 16) & 0xFF);
        data[2] = (byte) ((value >> 8) & 0xFF);
        data[3] = (byte) (value & 0xFF);
        return data;
    }

    /**
     * Converts the byte array to a long integer value.
     * @param b The byte array
     * @param offset Index of the first byte in the array
     * @return The long integer
     */
    public long b2l(byte[] b, int offset) {
        return (b[offset + 7] & 0xFF)
                | ((b[offset + 6] & 0xFF) << 8)
                | ((b[offset + 5] & 0xFF) << 16)
                | ((b[offset + 4] & 0xFF) << 24)
                | ((b[offset + 3] & 0xFF) << 32)
                | ((b[offset + 2] & 0xFF) << 40)
                | ((b[offset + 1] & 0xFF) << 48)
                | ((b[offset] & 0xFF) << 56);
    }

    /**
     * Converts a long integer to a byte array.
     * @param value The long value
     * @return The byte array
     */
    public byte[] l2b(long value) {
        byte[] data = new byte[8];
        data[0] = (byte) ((value >> 56) & 0xFF);
        data[1] = (byte) ((value >> 48) & 0xFF);
        data[2] = (byte) ((value >> 40) & 0xFF);
        data[3] = (byte) ((value >> 32) & 0xFF);
        data[4] = (byte) ((value >> 24) & 0xFF);
        data[5] = (byte) ((value >> 16) & 0xFF);
        data[6] = (byte) ((value >> 8) & 0xFF);
        data[7] = (byte) (value & 0xFF);
        return data;
    }
}
