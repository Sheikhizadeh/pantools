/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package genome;

//import java.util.regex.*;
import java.util.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Transaction;
import static pantools.Pantools.genome_label;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.sequence_label;
import static pantools.Pantools.write_fasta;

/**
 * Implements all the functionality to work with a 4-bit compressed sequence database. 
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class SequenceDatabase {

    private final String INFO_FILE = "/genomes.info";
    private final String DB_FILE = "/genomes.db";
    public long num_bytes;
    public int num_genomes;
    public int num_sequences[];        // Number of sequences in each genome    
    public String[][] sequence_titles; // Name of sequences for each genome
    public String[] genome_names;     // Paths to the genome FASTA files
    public long genome_length[];    // Length of sequences for each genome
    public long sequence_length[][];    // Length of sequences for each genome
    public long sequence_start[][];    // Length of sequences for each genome    
    private RandomAccessFile genomes_file;
    public MappedByteBuffer[] genomes_buff;
    private int parts_num;
    public long[] parts_size;
    public int max_byte = 2100000000;
    public char sym[];
    public int[] binary;
    public int[] complement;
    private String db_path;

    /**
     * Initialize sym, binary and complement arrays.
     * sym: All the nucleotide IUPAC symbols.
     * binary: The binary code for IUPAC symbols.
     * complement: The binary code for the complement of every binary code 
     */
    public void initalize() {
        sym = new char[]{'A', 'C', 'G', 'T', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N'};
        complement = new int[]{3, 2, 1, 0, 9, 8, 6, 7, 5, 4, 13, 12, 11, 10, 14};
        binary = new int[256];
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
    }

    /**
     * Mounts the genome database to the database object.
     * @param path Path to the genome database
     */
    public SequenceDatabase(String path) {
        int g, s, k;
        BufferedReader in;
        db_path = path;
        initalize();
        try {
            in = new BufferedReader(new FileReader(path + INFO_FILE));
            num_bytes = Long.valueOf(in.readLine().split(":")[1]);
            num_genomes = Integer.parseInt(in.readLine().split(":")[1]);
            genome_names = new String[num_genomes + 1];
            genome_length = new long[num_genomes + 1];
            sequence_titles = new String[num_genomes + 1][];
            sequence_length = new long[num_genomes + 1][];
            sequence_start = new long[num_genomes + 1][];
            num_sequences = new int[num_genomes + 1];
            for (g = 1; g <= num_genomes; ++g) {
                genome_length[g] = Long.valueOf(in.readLine().split(":")[1]);
                num_sequences[g] = Integer.parseInt(in.readLine().split(":")[1]);
                sequence_titles[g] = new String[num_sequences[g] + 1];
                sequence_length[g] = new long[num_sequences[g] + 1];
                sequence_start[g] = new long[num_sequences[g] + 1];
                for (s = 1; s <= num_sequences[g]; ++s) {
                    sequence_titles[g][s] = in.readLine().split(":")[1];
                    sequence_length[g][s] = Long.valueOf(in.readLine().split(":")[1]);
                    sequence_start[g][s] = Long.valueOf(in.readLine().split(":")[1]);
                }
            }
            in.close();
            if (Files.exists(Paths.get(path + DB_FILE))) {
                genomes_file = new RandomAccessFile(path + DB_FILE, "r");
                parts_num = (int) (num_bytes % max_byte == 0 ? num_bytes / max_byte : num_bytes / max_byte + 1);
                parts_size = new long[parts_num];
                genomes_buff = new MappedByteBuffer[parts_num];
                for (k = 0; k < parts_num; ++k) {
                    parts_size[k] = (int) (k == parts_num - 1 ? num_bytes % max_byte : max_byte);
                    genomes_buff[k] = genomes_file.getChannel().map(FileChannel.MapMode.READ_ONLY, k * parts_size[0], parts_size[k]);
                }
            } else {
                System.out.println("genome database not found!");
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
    }

    /**
     * Creates a new genome database.
     * @param path Path of the database
     * @param genome_paths_file A text file containing path to the genomes
     */
    public SequenceDatabase(String path, String genome_paths_file) {
        int g;
        BufferedReader in;
        String line;
        List<String> genome_list = new LinkedList();
        db_path = path;
        initalize();
        num_genomes = 0;
        new File(path).mkdir();
        num_bytes = 0;
        // count number of genomes    
        try {
            in = new BufferedReader(new FileReader(genome_paths_file));
            while (in.ready()) {
                line = in.readLine();
                if (line.equals("")) {
                    continue;
                }
                genome_list.add(line);
                num_genomes++;
            }
            in.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
        genome_names = new String[num_genomes + 1];
        genome_length = new long[num_genomes + 1];
        sequence_titles = new String[num_genomes + 1][];
        sequence_length = new long[num_genomes + 1][];
        sequence_start = new long[num_genomes + 1][];
        num_sequences = new int[num_genomes + 1];
        Iterator<String> itr = genome_list.iterator();
        // count number of sequences of genomes    
        for (g = 1; itr.hasNext(); ++g) {
            num_sequences[g] = 0;
            genome_names[g] = itr.next();
            try {
                in = new BufferedReader(new FileReader(genome_names[g]));
                while (in.ready()) {
                    line = in.readLine();
                    if (line.equals("")) {
                        continue;
                    }
                    if (line.charAt(0) == '>') {
                        num_sequences[g]++;
                    }
                }
                in.close();
            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            }
            sequence_titles[g] = new String[num_sequences[g] + 1];
            sequence_length[g] = new long[num_sequences[g] + 1];
            sequence_start[g] = new long[num_sequences[g] + 1];
        }
        code_genomes(path, 0);
        write_info();
    }

    /**
     * Reconstructs a genome database from the pangenome.
     * 
     * @param path Path of the genome database
     * @param graphDb The graph database object
     */
    public SequenceDatabase(String path, GraphDatabaseService graphDb) {
        new File(path).mkdir();
        Node db_node, seq_node, gen_node;
        int k, j, g, s;
        db_path = path;
        num_bytes = 0;
        initalize();
        try (Transaction tx = graphDb.beginTx()) {
            db_node = graphDb.findNodes(pangenome_label).next();
            num_genomes = (int) db_node.getProperty("num_genomes");
            genome_length = new long[num_genomes + 1];
            sequence_titles = new String[num_genomes + 1][];
            sequence_length = new long[num_genomes + 1][];
            sequence_start = new long[num_genomes + 1][];
            num_sequences = new int[num_genomes + 1];
            for (g = 1; g <= num_genomes; ++g) {
                gen_node = graphDb.findNode(genome_label, "number", g);
                num_sequences[g] = (int) gen_node.getProperty("num_sequences");
                sequence_titles[g] = new String[num_sequences[g] + 1];
                sequence_length[g] = new long[num_sequences[g] + 1];
                sequence_start[g] = new long[num_sequences[g] + 1];
                for (s = 1; s <= num_sequences[g]; ++s) {
                    seq_node = graphDb.findNode(sequence_label, "number", g + "_" + s);
                    sequence_titles[g][s] = (String) seq_node.getProperty("sequence_title");
                    sequence_length[g][s] = (long) seq_node.getProperty("sequence_length");
                    sequence_start[g][s] = num_bytes;
                    num_bytes += sequence_length[g][s] % 2 == 0 ? sequence_length[g][s] / 2 : sequence_length[g][s] / 2 + 1;
                }
            }
            try {
                genomes_file = new RandomAccessFile(path + DB_FILE, "rw");
            } catch (FileNotFoundException ioe) {
                System.out.println(ioe.getMessage());
                System.exit(1);
            }
            tx.success();
        }
        parts_num = (int) (num_bytes % max_byte == 0 ? num_bytes / max_byte : num_bytes / max_byte + 1);
        parts_size = new long[parts_num];
        genomes_buff = new MappedByteBuffer[parts_num];
        try {
            for (k = 0; k < parts_num; ++k) {
                parts_size[k] = (int) (k == parts_num - 1 ? num_bytes % max_byte : max_byte);
                genomes_buff[k] = genomes_file.getChannel().map(FileChannel.MapMode.READ_WRITE, k * parts_size[0], parts_size[k]);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
        write_info();
    }

    /**
     * Compresses genomes in a binary database, each nucleotide in 4 bits.
     * 
     * @param path Path of the genome database
     * @param previous_num_genomes The number of the genomes were already in the genome database
     */
    public void code_genomes(String path, int previous_num_genomes) {
        String line;
        char carry;
        boolean havecarry;
        long size = 0, byte_number;
        int j, g, s, len, k;
        BufferedReader in;
        byte_number = previous_num_genomes == 0 ? 0 : num_bytes;
        initalize();
        System.out.println("Reading " + (num_genomes - previous_num_genomes) + " genome(s)...");
        for (g = previous_num_genomes + 1; g <= num_genomes; ++g) {
            try {
                in = new BufferedReader(new FileReader(genome_names[g]));
                s = 0;
                while (in.ready()) {
                    line = in.readLine();
                    if (line.equals("")) {
                        continue;
                    }
                    if (line.charAt(0) == '>') {
                        ++s;
                        sequence_titles[g][s] = line.substring(1);
                        if (size % 2 == 1) {
                            ++size;
                        }
                        sequence_start[g][s] = num_bytes + size / 2;
                    } else {
                        len = line.length();
                        sequence_length[g][s] += len;
                        genome_length[g] += len;
                        size += len;
                    }
                }
                if (size % 2 == 1) {
                    ++size;
                }
                in.close();
            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            }
        }
        num_bytes += size / 2;
        try {
            genomes_file = new RandomAccessFile(path + DB_FILE, "rw");
        } catch (FileNotFoundException ioe) {
            System.out.println(ioe.getMessage());
            System.exit(1);
        }
        parts_num = (int) (num_bytes % max_byte == 0 ? num_bytes / max_byte : num_bytes / max_byte + 1);
        parts_size = new long[parts_num];
        genomes_buff = new MappedByteBuffer[parts_num];
        try {
            for (k = 0; k < parts_num; ++k) {
                parts_size[k] = (int) (k == parts_num - 1 ? num_bytes % max_byte : max_byte);
                genomes_buff[k] = genomes_file.getChannel().map(FileChannel.MapMode.READ_WRITE, k * parts_size[0], parts_size[k]);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
        genomes_buff[(int) (byte_number / parts_size[0])].position((int) (byte_number % parts_size[0]));
        for (g = previous_num_genomes + 1; g <= num_genomes; ++g) {
            try {
                in = new BufferedReader(new FileReader(genome_names[g]));
                carry = ' ';
                havecarry = false;
                s = 0;
                while (in.ready()) {
                    line = in.readLine().toUpperCase();
                    if (line.equals("")) {
                        continue;
                    }
                    if (line.charAt(0) != '>' && havecarry) {
                        line = carry + line;
                    }
                    if (line.charAt(0) == '>') {
                        if (havecarry) {
                            genomes_buff[(int) (byte_number / parts_size[0])].put((byte) (binary[carry] << 4));
                            ++byte_number;
                        }
                        havecarry = false;
                        ++s;
                        System.out.print("Sequence " + s + "/" + num_sequences[g] + " of genome " + g + " : " + genome_names[g]+"\t\t\t\t\t\r");
                    } else {
                        len = line.length();
                        havecarry = (len % 2 == 1);
                        if (havecarry) {
                            carry = line.charAt(len - 1);
                            --len;
                        }
                        for (j = 0; j < len; j += 2, ++byte_number) {
                            genomes_buff[(int) (byte_number / parts_size[0])].put((byte)((binary[line.charAt(j)] << 4) | binary[line.charAt(j + 1)]));
                        }
                    }
                }
                if (havecarry) {
                    genomes_buff[(int) (byte_number / parts_size[0])].put((byte) (binary[carry] << 4));
                    ++byte_number;
                }
                havecarry = false;
                in.close();
            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            }
            System.out.println();
        }
    }

    /**
     * Make FASTA file of a genome from the genome database.
     * 
     * @param path Path of the genome database
     * @param genome_number Number of the genome in the database
     * @param genome_name Name of the resulting FASTA file
     */
    public void decode_genome(String path, int genome_number, String genome_name) {
        int s;
        BufferedWriter out;
        initalize();
        System.out.println("Decoding " + num_genomes + " genome(s) into FASTA files...");
        try {
            out = new BufferedWriter(new FileWriter(db_path + genome_name + ".fasta"));
            for (s = 1; s <= num_sequences[genome_number]; ++s) {
                out.write(">" + sequence_titles[genome_number][s] + "\n");
                write_fasta(out, get_sequence(genome_number, s, 0, (int) sequence_length[genome_number][s], true), 80);
            }
            out.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
    }

    /**
     * Adds new genomes to the genome database.
     * 
     * @param path Path of the genome database
     * @param genome_paths_file A text file containing path to the genomes
     */
    public void add_genomes(String path, String genome_paths_file) {
        int g, s, i, previous_num_genomes = 0;
        BufferedReader in;
        String line;
        List<String> genome_list = new LinkedList();
        initalize();
        try {
            // count number of new genomes
            in = new BufferedReader(new FileReader(genome_paths_file));
            while (in.ready()) {
                line = in.readLine();
                if (line.equals("")) {
                    continue;
                }
                genome_list.add(line);
                num_genomes++;
            }
            in.close();
            in = new BufferedReader(new FileReader(path + INFO_FILE));
            num_bytes = Long.valueOf(in.readLine().split(":")[1]);
            previous_num_genomes = Integer.parseInt(in.readLine().split(":")[1]);
            genome_names = new String[num_genomes + 1];
            genome_length = new long[num_genomes + 1];
            sequence_titles = new String[num_genomes + 1][];
            sequence_length = new long[num_genomes + 1][];
            sequence_start = new long[num_genomes + 1][];
            num_sequences = new int[num_genomes + 1];
            for (g = 1; g <= previous_num_genomes; ++g) {
                genome_length[g] = Long.valueOf(in.readLine().split(":")[1]);
                num_sequences[g] = Integer.parseInt(in.readLine().split(":")[1]);
                sequence_titles[g] = new String[num_sequences[g] + 1];
                sequence_length[g] = new long[num_sequences[g] + 1];
                sequence_start[g] = new long[num_sequences[g] + 1];
                for (s = 1; s <= num_sequences[g]; ++s) {
                    sequence_titles[g][s] = in.readLine().split(":")[1];
                    sequence_length[g][s] = Long.valueOf(in.readLine().split(":")[1]);
                    sequence_start[g][s] = Long.valueOf(in.readLine().split(":")[1]);
                }
            }
            in.close();
            Iterator<String> itr = genome_list.iterator();
            for (g = previous_num_genomes + 1; itr.hasNext(); ++g) {
                num_sequences[g] = 0;
                genome_names[g] = itr.next();
                // count number of sequences
                in = new BufferedReader(new FileReader(genome_names[g]));
                while (in.ready()) {
                    line = in.readLine();
                    if (line.equals("")) {
                        continue;
                    }
                    if (line.charAt(0) == '>') {
                        num_sequences[g]++;
                    }
                }
                in.close();
                sequence_titles[g] = new String[num_sequences[g] + 1];
                sequence_length[g] = new long[num_sequences[g] + 1];
                sequence_start[g] = new long[num_sequences[g] + 1];
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
        code_genomes(path, previous_num_genomes);
        write_info();
    }

    /**
     * Closes the database.
     */
    public void close() {
        try {
            genomes_file.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
        for (int k = 0; k < parts_num; ++k) {
            genomes_buff[k] = null;
        }
    }

    /**
     * Writes the genomes.info file to disk.
     */
    public void write_info() {
        int g, s;
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(db_path + INFO_FILE));
            out.write("number_of_bytes:" + num_bytes + "\n");
            out.write("number_of_genomes:" + num_genomes + "\n");
            for (g = 1; g <= num_genomes; ++g) {
                out.write("genome length:" + genome_length[g] + "\n");
                out.write("number_of_sequences:" + num_sequences[g] + "\n");
                for (s = 1; s <= num_sequences[g]; ++s) {
                    out.write("sequence title:" + sequence_titles[g][s] + "\n");
                    out.write("sequence length:" + sequence_length[g][s] + "\n");
                    out.write("sequence start:" + sequence_start[g][s] + "\n");
                }
            }
            out.close();
        } catch (IOException ioe) {
            System.out.println(ioe.getMessage());
            System.exit(1);
        }
    }

    /**
     * Returns the nucleotide at a specified genomic position.
     * @param g Genome number 
     * @param s Sequence number
     * @param p Base position
     * @return The base 
     */
    public char get_symbol(int g, int s, int p) {
        if (p < sequence_length[g][s]) {
            byte b;
            long position = sequence_start[g][s] + p / 2;
            b = genomes_buff[(int) (position / parts_size[0])].get((int) (position % parts_size[0]));
            if (p % 2 == 0) {
                return sym[(b >> 4) & 0x0f];
            } else {
                return sym[b & 0x0f];
            }
        } else {
            return 0;
        }
    }

    /**
     * Returns the complement of a nucleotide at a specified genomic position.
     * @param g Genome number 
     * @param s Sequence number
     * @param p Base position
     * @return The base 
     */
    public char get_complement_symbol(int g, int s, int p) {
        if (p < sequence_length[g][s]) {
            byte b;
            long position = sequence_start[g][s] + p / 2;
            b = genomes_buff[(int) (position / parts_size[0])].get((int) (position % parts_size[0]));
            if (p % 2 == 0) {
                return sym[complement[(b >> 4) & 0x0f]];
            } else {
                return sym[complement[(b & 0x0f)]];
            }
        } else {
            return 0;
        }
    }

    /**
     * Returns the binary code of anucleotide at a specified genomic position.
     * @param g Genome number 
     * @param s Sequence number
     * @param p Base position
     * @return The base 
     */
    public int get_code(int g, int s, int p) {
        if (p < sequence_length[g][s] && p > -1) {
            byte b;
            long position = sequence_start[g][s] + p / 2;
            b = genomes_buff[(int) (position / parts_size[0])].get((int) (position % parts_size[0]));
            if (p % 2 == 0) {
                return (b >> 4) & 0x0f;
            } else {
                return (b & 0x0f);
            }
        } else {
            System.out.println("Wrong genomic position: " + p);
            return -1;
        }
    }

    /**
     * Returns the binary code of complement of a nucleotide at a specified genomic position.
     * @param g Genome number 
     * @param s Sequence number
     * @param p Base position
     * @return The base 
     */    
    public int get_complement_code(int g, int s, int p) {
        if (p < sequence_length[g][s]) {
            byte b;
            long position = sequence_start[g][s] + p / 2;
            b = genomes_buff[(int) (position / parts_size[0])].get((int) (position % parts_size[0]));
            if (p % 2 == 0) {
                return complement[(b >> 4) & 0x0f];
            } else {
                return complement[(b & 0x0f)];
            }
        } else {
            return -1;
        }
    }

    /**
     * Retrieves a genomic region from the genome database.
     * 
     * @param g Genome number
     * @param s Sequence number
     * @param p Start Position of the region
     * @param l Length of the region
     * @param direction specifies the direction, True for forward and False for reverse
     * @return The genomic region
     */
    public StringBuilder get_sequence(int g, int s, int p, int l, boolean direction) {
        int i;
        StringBuilder seq = new StringBuilder();
        if (p < 0) // take the part is available at the start of the sequence
        {
            l += p;
            p = 0;
        }
        if (p + l > sequence_length[g][s]) // take the part is available at the end of the sequence
        {
            l = (int) sequence_length[g][s] - p;
        }
        if (direction) {
            for (i = 0; i < l; ++i) {
                seq.append(get_symbol(g, s, p + i));
            }
        } else {
            for (i = l - 1; i >= 0 && p + i < sequence_length[g][s]; --i) {
                seq.append(get_complement_symbol(g, s, p + i));
            }
        }
        return seq;
    }

    /**
     * Determines the identity of two genomic regions.
     * 
     * @param database The second database we are making a comparison with. 
     * @param a1 Genomic address in the first database
     * @param a2 Genomic address in the second database
     * @param offset1 Distance to the start position of the first sequence
     * @param offset2 Distance to the start position of the second sequence
     * @param len The length of sequences
     * @param direction Direction of the second sequence True for forward and False for reverse
     * @return The result of the comparison
     */
    public boolean compare(SequenceDatabase database, int[] a1, int[] a2, int offset1, int offset2, int len, boolean direction) {
        if (a1[2] + offset1 + len - 1 >= sequence_length[a1[0]][a1[1]] || a2[2] + offset2 + len - 1 >= database.sequence_length[a2[0]][a2[1]]) {
            return false;
        }
        int i;
        boolean equal;
        if (direction) {
            for (equal = true, i = 0; i < len && equal; ++i) {
                if (get_code(a1[0], a1[1], a1[2] + offset1 + i) != database.get_code(a2[0], a2[1], a2[2] + offset2 + i)) {
                    equal = false;
                }
            }
        } else {
            for (equal = true, i = 0; i < len && equal; ++i) {
                if (get_code(a1[0], a1[1], a1[2] + offset1 + i) != database.get_complement_code(a2[0], a2[1], a2[2] + offset2 + len - i - 1)) {
                    equal = false;
                }
            }
        }
        return equal;
    }
}
