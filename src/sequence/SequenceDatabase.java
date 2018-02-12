/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sequence;

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
    public String[][] sequence_qualities; // FASTQ quality strings
    public String[] genome_names;     // Paths to the genome FASTA files
    public long genome_length[];    // Length of sequences for each genome
    public long sequence_length[][];    // Length of sequences for each genome
    public long sequence_offset[][];    // Cummulative length of previous sequences
    public long sequence_start[][];    // Length of sequences for each genome    
    private RandomAccessFile genomes_file;
    public MappedByteBuffer[] genomes_buff;
    private int parts_num;
    public long[] parts_size;
    public int max_byte = 100000000;
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
            sequence_qualities = new String[num_genomes + 1][];
            sequence_length = new long[num_genomes + 1][];
            sequence_offset = new long[num_genomes + 1][];
            sequence_start = new long[num_genomes + 1][];
            num_sequences = new int[num_genomes + 1];
            for (g = 1; g <= num_genomes; ++g) {
                genome_names[g] = in.readLine().split(":")[1];
                genome_length[g] = Long.valueOf(in.readLine().split(":")[1]);
                num_sequences[g] = Integer.parseInt(in.readLine().split(":")[1]);
                sequence_titles[g] = new String[num_sequences[g] + 1];
                sequence_qualities[g] = new String[num_sequences[g] + 1];
                sequence_length[g] = new long[num_sequences[g] + 1];
                sequence_offset[g] = new long[num_sequences[g] + 1];
                sequence_start[g] = new long[num_sequences[g] + 1];
                for (s = 1; s <= num_sequences[g]; ++s) {
                    sequence_titles[g][s] = in.readLine().split(":")[1];
                    sequence_qualities[g][s] = in.readLine().split(":")[1];
                    sequence_length[g][s] = Long.valueOf(in.readLine().split(":")[1]);
                    sequence_offset[g][s] = Long.valueOf(in.readLine().split(":")[1]);
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
        String line, file_type;
        String[] fields;
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
                line = in.readLine().trim();
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
        sequence_qualities = new String[num_genomes + 1][];
        sequence_length = new long[num_genomes + 1][];
        sequence_offset = new long[num_genomes + 1][];
        sequence_start = new long[num_genomes + 1][];
        num_sequences = new int[num_genomes + 1];
        Iterator<String> itr = genome_list.iterator();
        // count number of sequences of genomes    
        for (g = 1; itr.hasNext(); ++g) {
            num_sequences[g] = 0;
            genome_names[g] = itr.next();
            try {
                in = new BufferedReader(new FileReader(genome_names[g]));
                fields = genome_names[g].split("\\.");
                file_type = fields[fields.length - 1].toLowerCase();
                if (file_type.equals("fasta") || file_type.equals("fa") || file_type.equals("fna") || file_type.equals("fn")){
                    while (in.ready()) {
                        line = in.readLine().trim();
                        if (line.equals("")) 
                            continue;
                        if (line.charAt(0) == '>' || line.charAt(0) == '+') 
                            num_sequences[g]++;
                    }
                }else if (file_type.equals("fastq") || file_type.equals("fq") || file_type.equals("fnq") || file_type.equals("q")){
                    while (in.ready()) {
                        line = in.readLine().trim();
                        if (line.equals(""))
                            continue;
                        num_sequences[g]++;
                    }
                    num_sequences[g] /= 4;
                } else {
                    System.out.println(genome_names[g] + " does not have a valid extention (fasta, fa, fna, fn, fastq, fq, fnq, q)");
                    System.exit(1);
                }
                in.close();
            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            }
            sequence_titles[g] = new String[num_sequences[g] + 1];
            sequence_qualities[g] = new String[num_sequences[g] + 1];
            sequence_length[g] = new long[num_sequences[g] + 1];
            sequence_offset[g] = new long[num_sequences[g] + 1];
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
            genome_names = new String[num_genomes + 1];
            genome_length = new long[num_genomes + 1];
            sequence_titles = new String[num_genomes + 1][];
            sequence_qualities = new String[num_genomes + 1][];
            sequence_length = new long[num_genomes + 1][];
            sequence_offset = new long[num_genomes + 1][];
            sequence_start = new long[num_genomes + 1][];
            num_sequences = new int[num_genomes + 1];
            for (g = 1; g <= num_genomes; ++g) {
                gen_node = graphDb.findNode(genome_label, "number", g);
                genome_names[g] = (String) gen_node.getProperty("path");
                num_sequences[g] = (int) gen_node.getProperty("num_sequences");
                sequence_titles[g] = new String[num_sequences[g] + 1];
                sequence_qualities[g] = new String[num_sequences[g] + 1];
                sequence_length[g] = new long[num_sequences[g] + 1];
                sequence_offset[g] = new long[num_sequences[g] + 1];
                sequence_start[g] = new long[num_sequences[g] + 1];
                for (s = 1; s <= num_sequences[g]; ++s) {
                    seq_node = graphDb.findNode(sequence_label, "identifier", g + "_" + s);
                    sequence_titles[g][s] = (String) seq_node.getProperty("title");
                    sequence_qualities[g][s] = (String) seq_node.getProperty("quality","*");
                    sequence_length[g][s] = (long) seq_node.getProperty("length");
                    sequence_offset[g][s] = (long) seq_node.getProperty("offset");
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
        String line, file_type;
        String[] fields;
        char carry;
        boolean havecarry;
        long size = 0, byte_number;
        int j, g, s, len, k;
        BufferedReader in;
        byte_number = previous_num_genomes == 0 ? 0 : num_bytes;
        initalize();
        System.out.println("Reading " + (num_genomes - previous_num_genomes) + " genome(s)...");
        try {
            for (g = previous_num_genomes + 1; g <= num_genomes; ++g) {
                in = new BufferedReader(new FileReader(genome_names[g]));
                fields = genome_names[g].split("\\.");
                file_type = fields[fields.length - 1].toLowerCase();
                if (file_type.equals("fasta") || file_type.equals("fa") || file_type.equals("fna") || file_type.equals("fn")){
                    s = 0;
                    sequence_offset[g][0] = 0;
                    sequence_length[g][0] = 0;
                    while (in.ready()) {
                        line = in.readLine().trim();
                        if (line.equals("")) {
                            continue;
                        }
                        if (line.charAt(0) == '>') {
                            ++s;
                            sequence_titles[g][s] = line.substring(1);
                            sequence_qualities[g][s] = "*";
                            sequence_offset[g][s] = sequence_offset[g][s - 1] + sequence_length[g][s - 1];
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
                }else if (file_type.equals("fastq") || file_type.equals("fq") || file_type.equals("fnq") || file_type.equals("q")){
                    
                    sequence_offset[g][0] = 0;
                    sequence_length[g][0] = 0;
                    for (s = 1; in.ready(); ++s) {
                    // read title    
                        line = in.readLine().trim();
                        sequence_titles[g][s] = line.substring(1);
                        sequence_offset[g][s] = sequence_offset[g][s - 1] + sequence_length[g][s - 1];
                        if (size % 2 == 1) {
                            ++size;
                        }
                        sequence_start[g][s] = num_bytes + size / 2;
                    // read sequence
                        line = in.readLine().trim();
                        sequence_length[g][s] += line.length();
                        genome_length[g] += line.length();
                        size += line.length();
                    // read +    
                        in.readLine();
                    // read quality    
                        sequence_qualities[g][s] = in.readLine().trim();
                    }
                    if (size % 2 == 1)
                        ++size;
                } else {
                    System.out.println(genome_names[g] + " does not have a valid extention (fasta, fa, fna, fn, fastq, fq, fnq, q)");
                    System.exit(1);
                }
                in.close();
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
            for (k = 0; k < parts_num; ++k) {
                parts_size[k] = (int) (k == parts_num - 1 ? num_bytes % max_byte : max_byte);
                genomes_buff[k] = genomes_file.getChannel().map(FileChannel.MapMode.READ_WRITE, k * parts_size[0], parts_size[k]);
            }
            genomes_buff[(int) (byte_number / parts_size[0])].position((int) (byte_number % parts_size[0]));
            for (g = previous_num_genomes + 1; g <= num_genomes; ++g) {
                in = new BufferedReader(new FileReader(genome_names[g]));
                fields = genome_names[g].split("\\.");
                file_type = fields[fields.length - 1].toLowerCase();
                if (file_type.equals("fasta") || file_type.equals("fa") || file_type.equals("fna") || file_type.equals("fn")){
                    carry = ' ';
                    havecarry = false;
                    s = 0;
                    while (in.ready()) {
                        line = in.readLine().trim().toUpperCase();
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
                            //System.out.println("Reading sequence " + s + "/" + num_sequences[g] + " of genome " + g + " : " + genome_names[g]);
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
                }else if (file_type.equals("fastq") || file_type.equals("fq") || file_type.equals("fnq") || file_type.equals("q")){
                    carry = ' ';
                    havecarry = false;
                    while (in.ready()) {
                    // read title
                        in.readLine();
                        if (havecarry) {
                            genomes_buff[(int) (byte_number / parts_size[0])].put((byte) (binary[carry] << 4));
                            ++byte_number;
                        }
                        havecarry = false;
                    //read sequence    
                        line = in.readLine().trim().toUpperCase();
                        len = line.length();
                        havecarry = (len % 2 == 1);
                        if (havecarry) {
                            carry = line.charAt(len - 1);
                            --len;
                        }
                        for (j = 0; j < len; j += 2, ++byte_number) {
                            genomes_buff[(int) (byte_number / parts_size[0])].put((byte)((binary[line.charAt(j)] << 4) | binary[line.charAt(j + 1)]));
                        }
                    // read +
                        in.readLine();
                    // read quality
                        in.readLine();
                    }
                    if (havecarry) {
                        genomes_buff[(int) (byte_number / parts_size[0])].put((byte) (binary[carry] << 4));
                        ++byte_number;
                    }
                    havecarry = false;
                }else{
                    System.err.println(genome_names[g] + " should have one of these extensions: fasta, fa, fna, fn for FATSA or"
                            + "fastq, fq, fnq, q for FASTQ.");
                    continue;
                }
                in.close();
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
    }

    /**
     * Make FASTA file of a genome from the genome database.
     * 
     * @param path Path of the genome database
     * @param genome_number Number of the genome in the database
     * @param genome_name Name of the resulting FASTA file
     */
    /*public void decode_genome(String path, int genome_number, String genome_name) {
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
    }*/

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
                line = in.readLine().trim();
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
            sequence_qualities = new String[num_genomes + 1][];
            sequence_length = new long[num_genomes + 1][];
            sequence_offset = new long[num_genomes + 1][];
            sequence_start = new long[num_genomes + 1][];
            num_sequences = new int[num_genomes + 1];
            for (g = 1; g <= previous_num_genomes; ++g) {
                genome_names[g] = in.readLine();
                genome_length[g] = Long.valueOf(in.readLine().split(":")[1]);
                num_sequences[g] = Integer.parseInt(in.readLine().split(":")[1]);
                sequence_titles[g] = new String[num_sequences[g] + 1];
                sequence_qualities[g] = new String[num_sequences[g] + 1];
                sequence_length[g] = new long[num_sequences[g] + 1];
                sequence_offset[g] = new long[num_sequences[g] + 1];
                sequence_start[g] = new long[num_sequences[g] + 1];
                for (s = 1; s <= num_sequences[g]; ++s) {
                    sequence_titles[g][s] = in.readLine().split(":")[1];
                    sequence_qualities[g][s] = in.readLine().split(":")[1];
                    sequence_length[g][s] = Long.valueOf(in.readLine().split(":")[1]);
                    sequence_offset[g][s] = Long.valueOf(in.readLine().split(":")[1]);
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
                    line = in.readLine().trim();
                    if (line.equals("")) {
                        continue;
                    }
                    if (line.charAt(0) == '>') {
                        num_sequences[g]++;
                    }
                }
                in.close();
                sequence_titles[g] = new String[num_sequences[g] + 1];
                sequence_qualities[g] = new String[num_sequences[g] + 1];
                sequence_length[g] = new long[num_sequences[g] + 1];
                sequence_offset[g] = new long[num_sequences[g] + 1];
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
                out.write("genome name:" + genome_names[g] + "\n");
                out.write("genome length:" + genome_length[g] + "\n");
                out.write("number_of_sequences:" + num_sequences[g] + "\n");
                for (s = 1; s <= num_sequences[g]; ++s) {
                    out.write("sequence title:" + sequence_titles[g][s] + "\n");
                    out.write("sequence quality:" + sequence_qualities[g][s] + "\n");
                    out.write("sequence length:" + sequence_length[g][s] + "\n");
                    out.write("sequence offset:" + sequence_offset[g][s] + "\n");
                    out.write("sequence start:" + sequence_start[g][s] + "\n");
                }
            }
            out.close();
        } catch (IOException ioe) {
            System.out.println(ioe.getMessage());
            System.exit(1);
        }
    }
}
