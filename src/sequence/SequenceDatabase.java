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
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.zip.GZIPInputStream;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Transaction;
import static pantools.Pantools.genome_label;
import static pantools.Pantools.labels;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.sequence_label;

/**
 * Implements all the functionality to work with a 4-bit compressed sequence database. 
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class SequenceDatabase {

    private final String INFO_FILE = "/sequences.info";
    private final String DB_FILE = "/sequences.db";
    public long number_of_bytes;
    public int num_genomes;
    public int num_sequences[];        // Number of sequences in each genome    
    public String[][] sequence_titles; // Name of sequences for each genome
    //public String[][] sequence_qualities; // FASTQ quality strings
    public String[] genome_names;     // Paths to the genome FASTA files
    public long genome_length[];    // Length of sequences for each genome
    public long sequence_length[][];    // Length of sequences for each genome
    public long sequence_offset[][];    // Cummulative length of previous sequences
    public long sequence_start[][];    // Length of sequences for each genome    
    private RandomAccessFile genomes_file;
    public MappedByteBuffer[] genomes_buff;
    public int MAX_BYTE_COUNT = 1 << 24;// 16 MB 
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
        int g, s, i, number_of_pages, page_size;
        String[] fields;
        BufferedReader in;
        db_path = path;
        initalize();
        try {
            in = new BufferedReader(new FileReader(path + INFO_FILE));
            number_of_bytes = Long.valueOf(in.readLine().split(":")[1]);
            num_genomes = Integer.parseInt(in.readLine().split(":")[1]);
            genome_names = new String[num_genomes + 1];
            genome_length = new long[num_genomes + 1];
            sequence_titles = new String[num_genomes + 1][];
            //sequence_qualities = new String[num_genomes + 1][];
            sequence_length = new long[num_genomes + 1][];
            sequence_offset = new long[num_genomes + 1][];
            sequence_start = new long[num_genomes + 1][];
            num_sequences = new int[num_genomes + 1];
            for (g = 1; g <= num_genomes; ++g) {
                genome_names[g] = in.readLine().split(":")[1];
                genome_length[g] = Long.valueOf(in.readLine().split(":")[1]);
                num_sequences[g] = Integer.parseInt(in.readLine().split(":")[1]);
                sequence_titles[g] = new String[num_sequences[g] + 1];
                //sequence_qualities[g] = new String[num_sequences[g] + 1];
                sequence_length[g] = new long[num_sequences[g] + 1];
                sequence_offset[g] = new long[num_sequences[g] + 1];
                sequence_start[g] = new long[num_sequences[g] + 1];
                in.readLine(); // skip the top row
                for (s = 1; s <= num_sequences[g]; ++s) {
                    fields = in.readLine().split("\t");
                    sequence_titles[g][s] = fields[0];
                    sequence_length[g][s] = Long.valueOf(fields[1]);
                    sequence_offset[g][s] = Long.valueOf(fields[2]);
                    sequence_start[g][s] = Long.valueOf(fields[3]);
                    //sequence_qualities[g][s] = fields[4];
                }
            }
            in.close();
            if (Files.exists(Paths.get(path + DB_FILE))) {
                genomes_file = new RandomAccessFile(path + DB_FILE, "r");
                number_of_pages = (int) (number_of_bytes % MAX_BYTE_COUNT == 0 ? number_of_bytes / MAX_BYTE_COUNT : number_of_bytes / MAX_BYTE_COUNT + 1);
                genomes_buff = new MappedByteBuffer[number_of_pages];
                for (i = 0; i < number_of_pages; ++i) {
                    page_size = (int) (i == number_of_pages - 1 ? number_of_bytes % MAX_BYTE_COUNT : MAX_BYTE_COUNT);
                    page_size = page_size == 0 ? MAX_BYTE_COUNT : page_size;
                    genomes_buff[i] = genomes_file.getChannel().map(FileChannel.MapMode.READ_ONLY, ((long)i) * MAX_BYTE_COUNT, page_size);
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
        number_of_bytes = 0;
        // count number of genomes    
        try {
            in = open_file(genome_paths_file);
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
        //sequence_qualities = new String[num_genomes + 1][];
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
                if (is_fasta(genome_names[g])){
                    in = open_file(genome_names[g]);
                    while ((line = in.readLine()) != null){
                        if (line.equals("")) 
                            continue;
                        if (line.charAt(0) == '>')
                            num_sequences[g]++;
                    }
                    in.close();
                }else if (is_fastq(genome_names[g])){
                    in = open_file(genome_names[g]);
                    while ((line = in.readLine()) != null){
                        if (line.equals(""))
                            continue;
                        num_sequences[g]++;
                    }
                    num_sequences[g] /= 4;
                    in.close();
                } else {
                    System.out.println(genome_names[g] + " is not a proper FASTA of FAASTQ file");
                    System.exit(1);
                }
            } catch (IOException e) {
                System.out.println(e.getMessage());
                System.exit(1);
            }
            sequence_titles[g] = new String[num_sequences[g] + 1];
            //sequence_qualities[g] = new String[num_sequences[g] + 1];
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
        int i, number_of_pages, g, s, page_size;
        db_path = path;
        number_of_bytes = 0;
        initalize();
        try (Transaction tx = graphDb.beginTx()) {
            db_node = graphDb.findNodes(pangenome_label).next();
            num_genomes = (int) db_node.getProperty("num_genomes");
            genome_names = new String[num_genomes + 1];
            genome_length = new long[num_genomes + 1];
            sequence_titles = new String[num_genomes + 1][];
            //sequence_qualities = new String[num_genomes + 1][];
            sequence_length = new long[num_genomes + 1][];
            sequence_offset = new long[num_genomes + 1][];
            sequence_start = new long[num_genomes + 1][];
            num_sequences = new int[num_genomes + 1];
            for (g = 1; g <= num_genomes; ++g) {
                gen_node = graphDb.findNode(genome_label, "number", g);
                genome_names[g] = (String) gen_node.getProperty("path");
                num_sequences[g] = (int) gen_node.getProperty("num_sequences");
                sequence_titles[g] = new String[num_sequences[g] + 1];
                //sequence_qualities[g] = new String[num_sequences[g] + 1];
                sequence_length[g] = new long[num_sequences[g] + 1];
                sequence_offset[g] = new long[num_sequences[g] + 1];
                sequence_start[g] = new long[num_sequences[g] + 1];
                for (s = 1; s <= num_sequences[g]; ++s) {
                    seq_node = graphDb.findNode(sequence_label, "identifier", g + "_" + s);
                    sequence_titles[g][s] = (String) seq_node.getProperty("title");
                    //sequence_qualities[g][s] = (String) seq_node.getProperty("quality","*");
                    sequence_length[g][s] = (long) seq_node.getProperty("length");
                    sequence_offset[g][s] = (long) seq_node.getProperty("offset");
                    sequence_start[g][s] = number_of_bytes;
                    number_of_bytes += sequence_length[g][s] % 2 == 0 ? sequence_length[g][s] / 2 : sequence_length[g][s] / 2 + 1;
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
        number_of_pages = (int) (number_of_bytes % MAX_BYTE_COUNT == 0 ? number_of_bytes / MAX_BYTE_COUNT : number_of_bytes / MAX_BYTE_COUNT + 1);
        genomes_buff = new MappedByteBuffer[number_of_pages];
        try {
            for (i = 0; i < number_of_pages; ++i) {
                page_size = (int) (i == number_of_pages - 1 ? number_of_bytes % MAX_BYTE_COUNT : MAX_BYTE_COUNT);
                page_size = page_size == 0 ? MAX_BYTE_COUNT : page_size;
                genomes_buff[i] = genomes_file.getChannel().map(FileChannel.MapMode.READ_WRITE, ((long)i) * MAX_BYTE_COUNT, page_size);
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
        int i, j, g, s, len, number_of_pages, page_size;
        BufferedReader in;
        byte_number = previous_num_genomes == 0 ? 0 : number_of_bytes;
        initalize();
        try {
            for (g = previous_num_genomes + 1; g <= num_genomes; ++g) {
                System.out.println("Reading " + genome_names[g] + " ...");
                if (is_fasta(genome_names[g])){
                    in = open_file(genome_names[g]);
                    s = 0;
                    sequence_offset[g][0] = 0;
                    sequence_length[g][0] = 0;
                    while ((line = in.readLine()) != null){
                        if (line.equals("")) 
                            continue;
                        if (line.charAt(0) == '>') {
                            ++s;
                            sequence_titles[g][s] = line.substring(1);
                            //sequence_qualities[g][s] = "*";
                            sequence_offset[g][s] = sequence_offset[g][s - 1] + sequence_length[g][s - 1];
                            if (size % 2 == 1) {
                                ++size;
                            }
                            sequence_start[g][s] = number_of_bytes + size / 2;
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
                }else if (is_fastq(genome_names[g])){
                    in = open_file(genome_names[g]);
                    sequence_offset[g][0] = 0;
                    sequence_length[g][0] = 0;
                    for (s = 1;(line = in.readLine()) != null; ++s){
                    // read title    
                        if (line.equals(""))
                            continue;
                        sequence_titles[g][s] = line.substring(1);
                        sequence_offset[g][s] = sequence_offset[g][s - 1] + sequence_length[g][s - 1];
                        if (size % 2 == 1) {
                            ++size;
                        }
                        sequence_start[g][s] = number_of_bytes + size / 2;
                    // read sequence
                        line = in.readLine();
                        sequence_length[g][s] += line.length();
                        genome_length[g] += line.length();
                        size += line.length();
                    // read +    
                        in.readLine();
                    // read quality    
                        in.readLine();
                        //sequence_qualities[g][s] = "*";
                    }
                    if (size % 2 == 1)
                        ++size;
                    in.close();
                } else {
                    System.out.println(genome_names[g] + " is not a proper FASTA of FAASTQ file");
                    System.exit(1);
                }
            }
            number_of_bytes += size / 2;
            try {
                genomes_file = new RandomAccessFile(path + DB_FILE, "rw");
            } catch (FileNotFoundException ioe) {
                System.out.println(ioe.getMessage());
                System.exit(1);
            }
            number_of_pages = (int) (number_of_bytes / MAX_BYTE_COUNT + (number_of_bytes % MAX_BYTE_COUNT == 0 ? 0 : 1));
            genomes_buff = new MappedByteBuffer[number_of_pages];
            for (i = 0; i < number_of_pages; ++i) {
                page_size = (int) (i == number_of_pages - 1 ? number_of_bytes % MAX_BYTE_COUNT : MAX_BYTE_COUNT);
                page_size = page_size == 0 ? MAX_BYTE_COUNT : page_size;
                genomes_buff[i] = genomes_file.getChannel().map(FileChannel.MapMode.READ_WRITE, ((long)i) * MAX_BYTE_COUNT, page_size);
            }
            genomes_buff[(int) (byte_number / MAX_BYTE_COUNT)].position((int) (byte_number % MAX_BYTE_COUNT));
            for (g = previous_num_genomes + 1; g <= num_genomes; ++g) {
                if (is_fasta(genome_names[g])){
                    in = open_file(genome_names[g]);
                    carry = ' ';
                    havecarry = false;
                    s = 0;
                    while ((line = in.readLine()) != null){
                        line = line.toUpperCase();
                        if (line.equals("")) {
                            continue;
                        }
                        if (line.charAt(0) != '>' && havecarry) {
                            line = carry + line;
                        }
                        if (line.charAt(0) == '>') {
                            if (havecarry) {
                                genomes_buff[(int) (byte_number / MAX_BYTE_COUNT)].put((byte) (binary[carry] << 4));
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
                                genomes_buff[(int) (byte_number / MAX_BYTE_COUNT)].put((byte)((binary[line.charAt(j)] << 4) | binary[line.charAt(j + 1)]));
                            }
                        }
                    }
                    if (havecarry) {
                        genomes_buff[(int) (byte_number / MAX_BYTE_COUNT)].put((byte) (binary[carry] << 4));
                        ++byte_number;
                    }
                    havecarry = false;
                    in.close();
                }else if (is_fastq(genome_names[g])){
                    in = open_file(genome_names[g]);
                    carry = ' ';
                    havecarry = false;
                    while ((line = in.readLine()) != null){
                    // read title
                        if (havecarry) {
                            genomes_buff[(int) (byte_number / MAX_BYTE_COUNT)].put((byte) (binary[carry] << 4));
                            ++byte_number;
                        }
                        havecarry = false;
                    //read sequence    
                        line = in.readLine().toUpperCase();
                        len = line.length();
                        havecarry = (len % 2 == 1);
                        if (havecarry) {
                            carry = line.charAt(len - 1);
                            --len;
                        }
                        for (j = 0; j < len; j += 2, ++byte_number) {
                            genomes_buff[(int) (byte_number / MAX_BYTE_COUNT)].put((byte)((binary[line.charAt(j)] << 4) | binary[line.charAt(j + 1)]));
                        }
                    // read +
                        in.readLine();
                    // read quality
                        in.readLine();
                    }
                    if (havecarry) {
                        genomes_buff[(int) (byte_number / MAX_BYTE_COUNT)].put((byte) (binary[carry] << 4));
                        ++byte_number;
                    }
                    havecarry = false;
                    in.close();
                }else{
                    System.out.println(genome_names[g] + " is not a proper FASTA of FAASTQ file");
                    continue;
                }
                in.close();
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
    }

    public boolean is_fasta(String file_name){
        try{
            BufferedReader in = open_file(file_name);
            String line;
            while ((line = in.readLine()) != null){
                if (line.equals("")) 
                    continue;
                else {
                    in.close();
                    return line.charAt(0) == '>';
                }            
            }
        } catch (IOException ex){
            System.err.println(ex.getMessage());
        }
        return false;
    }

    public boolean is_fastq(String file_name){
        try{
            BufferedReader in = open_file(file_name);
            String line;
            while ((line = in.readLine()) != null){
                if (line.equals("")) 
                    continue;
                else {
                    in.close();
                    return line.charAt(0) == '@';
                }
            }
        } catch (IOException ex){
            System.err.println(ex.getMessage());
        }
        return false;
    }
    
    private BufferedReader open_file(String filename){
        try{        
            String[] fields = filename.split("\\.");
            String file_type = fields[fields.length - 1].toLowerCase();
            if (file_type.equals("gz") || file_type.equals("gzip"))
                return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filename)), "UTF-8"));                    
            else 
                return new BufferedReader(new BufferedReader(new FileReader(filename)));                    
        } catch (IOException ex){
            System.out.println(ex.getMessage());
            return null;
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
        String[] fields;
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
            number_of_bytes = Long.valueOf(in.readLine().split(":")[1]);
            previous_num_genomes = Integer.parseInt(in.readLine().split(":")[1]);
            genome_names = new String[num_genomes + 1];
            genome_length = new long[num_genomes + 1];
            sequence_titles = new String[num_genomes + 1][];
            //sequence_qualities = new String[num_genomes + 1][];
            sequence_length = new long[num_genomes + 1][];
            sequence_offset = new long[num_genomes + 1][];
            sequence_start = new long[num_genomes + 1][];
            num_sequences = new int[num_genomes + 1];
            for (g = 1; g <= previous_num_genomes; ++g) {
                genome_names[g] = in.readLine();
                genome_length[g] = Long.valueOf(in.readLine().split(":")[1]);
                num_sequences[g] = Integer.parseInt(in.readLine().split(":")[1]);
                sequence_titles[g] = new String[num_sequences[g] + 1];
                //sequence_qualities[g] = new String[num_sequences[g] + 1];
                sequence_length[g] = new long[num_sequences[g] + 1];
                sequence_offset[g] = new long[num_sequences[g] + 1];
                sequence_start[g] = new long[num_sequences[g] + 1];
                in.readLine(); // skip the top row
                for (s = 1; s <= num_sequences[g]; ++s) {
                    fields = in.readLine().split("\t");
                    sequence_titles[g][s] = fields[0];
                    sequence_length[g][s] = Long.valueOf(fields[1]);
                    sequence_offset[g][s] = Long.valueOf(fields[2]);
                    sequence_start[g][s] = Long.valueOf(fields[3]);
                    //sequence_qualities[g][s] = fields[4];
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
                //sequence_qualities[g] = new String[num_sequences[g] + 1];
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
        for (int k = 0; k < genomes_buff.length; ++k) {
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
            out.write("number_of_bytes:" + number_of_bytes + "\n");
            out.write("number_of_genomes:" + num_genomes + "\n");
            for (g = 1; g <= num_genomes; ++g) {
                out.write("genome name:" + genome_names[g] + "\n");
                out.write("genome length:" + genome_length[g] + "\n");
                out.write("number_of_sequences:" + num_sequences[g] + "\n");
                out.write("sequence_title\tsequence_quality\tsequence_length\tsequence_offset\tsequence_start\n" );
                for (s = 1; s <= num_sequences[g]; ++s) 
                    out.write(sequence_titles[g][s] + "\t" + sequence_length[g][s] + "\t" + sequence_offset[g][s] + "\t" + sequence_start[g][s] + "\n");
            }
            out.close();
        } catch (IOException ioe) {
            System.out.println(ioe.getMessage());
            System.exit(1);
        }
    }
}
