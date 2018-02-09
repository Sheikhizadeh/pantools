/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sequence;

import index.kmer;
import static pantools.Pantools.DEBUG;

/**
 *
 * @author siavash
 */
public class SequenceScanner {
    private int to_genome;
    private int to_sequence;
    private int genome;
    private int sequence;
    private int position;
    private int K;
    private kmer curr_kmer;
    private long curr_index;
    SequenceDatabase database;
    
    public SequenceScanner(SequenceDatabase db, int fg, int tg, int fs, int ts, int k, int pre_len){
        database = db;
        K = k;
        curr_kmer = new kmer(K,pre_len);
        genome = fg;
        to_genome = tg;
        sequence = fs;
        to_sequence = ts;
    }
    
    public int get_genome(){
        return genome;
    }
    public int get_sequence(){
        return sequence;
    }
    public int get_position(){
        return position;
    }
    public kmer get_curr_kmer(){
        return curr_kmer;
    }
    public long get_curr_index(){
     return curr_index;   
    }   
    public int[] get_address(){
        return new int[]{genome, sequence, position};
    }
    
    public long get_sequence_length(){
        return database.sequence_length[genome][sequence];
    }
    
    public long get_sequence_length(int g, int s){
        return database.sequence_length[g][s];
    }

    public int get_code(int offset) {
        if (position + offset < database.sequence_length[genome][sequence] && position + offset > -1) {
            byte b;
            long pos = database.sequence_start[genome][sequence] + (position + offset) / 2;
            b = database.genomes_buff[(int) (pos / database.parts_size[0])].get((int) (pos % database.parts_size[0]));
            if ((position + offset) % 2 == 0) {
                return (b >> 4) & 0x0f;
            } else {
                return (b & 0x0f);
            }
        } else {
            System.out.println("Wrong genomic position: " + (position + offset));
            return -1;
        }
    }

    public void set_genome(int g){
        genome = g;
    }
    public void set_sequence(int s){
        sequence = s;
    }
    public void set_position(int p){
        position = p;
    } 
    public void set_curr_index(long ci){
     curr_index = ci;   
    }

    public boolean end_of_scan(){
        return genome > to_genome;
    }      
    public boolean end_of_genome(){
        return sequence > to_sequence;
    }    
    public boolean end_of_sequence(){
        return position >= get_sequence_length() - 1;
    }
    public boolean start_of_sequence(){
        return position <= 0;
    }
    
    public void next_genome(){
        ++genome;
        sequence = 1;
        if(!end_of_scan())
            to_sequence = database.num_sequences[genome];
    }
    public void next_sequence(){
        ++sequence;
    }
    public void next_position(){
        ++position;
    }
    public void previous_position(){
        --position;
    }

    
    public boolean initialize_left_kmer(int start, int stop) {
        if (DEBUG) System.out.println("initialize_left_kmer");
        int i;
        curr_kmer.reset();
        position = start - 1;
        for (i = 0; i < K && position < stop; ++i) {
            if (get_code(1) > 3) {
                next_position();
                if (DEBUG) System.out.println("jump_forward");
                jump_forward();
                return false;
            }
            next_position();
            curr_kmer.next_kmer(get_code(0));
        }
        if (DEBUG) System.out.println(curr_kmer.toString());
        return true;
    }

    public boolean initialize_right_kmer(int start, int stop) {
        if (DEBUG) System.out.println("initialize_right_kmer min = " + start);
        int i;
        curr_kmer.reset();
        position = stop + 1;
        for (i = 0; i < K && position > start; ++i) {
            if (get_code(-1) > 3) {
                previous_position();
                if (DEBUG) System.out.println("jump_backward");
                jump_backward();
                return false;
            }
            previous_position();
            curr_kmer.prev_kmer(get_code(0));
        }
        if (DEBUG) System.out.println(curr_kmer.toString());
        return true;
    }    


    /**
     * Jump over an ambiguous region; at first, position points to the first position which degenerate starts, 
     * after jumping it points to the last base of the first K-mer after the ambiguous region. 
     */
    public void jump_forward() {
        int j;
        int base_code = get_code(0);
        curr_kmer.reset();
        do {
            while (base_code > 3 && !end_of_sequence()) {
                next_position();
                base_code = get_code(0);
            }
            curr_kmer.next_kmer(base_code);
            for (j = 0; j < K - 1 && !end_of_sequence(); ++j) {
                next_position();
                base_code = get_code(0);
                if (base_code > 3) {
                    break;
                }
                curr_kmer.next_kmer(base_code);
            }
        } while (base_code > 3 && !end_of_sequence());
        if (DEBUG) System.out.println(curr_kmer.toString());
    }

    /**
     * Jump over an ambiguous region; at first, position points to the first position which degenerate starts, 
     * after jumping it points to the last base of the first K-mer after the ambiguous region. 
     */
    public void jump_backward() {
        int j;
        int base_code = get_code(0);
        curr_kmer.reset();
        do {
            while (base_code > 3 && !start_of_sequence()) {
                previous_position();
                base_code = get_code(0);
            }
            curr_kmer.prev_kmer(base_code);
            for (j = 0; j < K - 1 && !start_of_sequence(); ++j) {
                previous_position();
                base_code = get_code(0);
                if (base_code > 3) {
                    break;
                }
                curr_kmer.prev_kmer(base_code);
            }
        } while (base_code > 3 && !start_of_sequence());
        if (DEBUG) System.out.println(curr_kmer.toString());
    }
    
    /**
     * Returns the nucleotide at a specified genomic position.
     * @param g Genome number 
     * @param s Sequence number
     * @param p Base position
     * @return The base 
     */
    public char get_symbol(int g, int s, int p) {
        if (p < database.sequence_length[g][s]) {
            byte b;
            long pos = database.sequence_start[g][s] + p / 2;
            b = database.genomes_buff[(int) (pos / database.parts_size[0])].get((int) (pos % database.parts_size[0]));
            if (p % 2 == 0) {
                return database.sym[(b >> 4) & 0x0f];
            } else {
                return database.sym[b & 0x0f];
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
        if (p < database.sequence_length[g][s]) {
            byte b;
            long pos = database.sequence_start[g][s] + p / 2;
            b = database.genomes_buff[(int) (pos / database.parts_size[0])].get((int) (pos % database.parts_size[0]));
            if (p % 2 == 0) {
                return database.sym[database.complement[(b >> 4) & 0x0f]];
            } else {
                return database.sym[database.complement[(b & 0x0f)]];
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
        if (p < database.sequence_length[g][s] && p > -1) {
            byte b;
            long pos = database.sequence_start[g][s] + p / 2;
            b = database.genomes_buff[(int) (pos / database.parts_size[0])].get((int) (pos % database.parts_size[0]));
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
        if (p < database.sequence_length[g][s]) {
            byte b;
            long pos = database.sequence_start[g][s] + p / 2;
            b = database.genomes_buff[(int) (pos / database.parts_size[0])].get((int) (pos % database.parts_size[0]));
            if (p % 2 == 0) {
                return database.complement[(b >> 4) & 0x0f];
            } else {
                return database.complement[(b & 0x0f)];
            }
        } else {
            return -1;
        }
    }

    public int get_complement_current_code(int offset) {
        if (position + offset < database.sequence_length[genome][sequence]) {
            byte b;
            long pos = database.sequence_start[genome][sequence] + (position + offset) / 2;
            b = database.genomes_buff[(int) (pos / database.parts_size[0])].get((int) (pos % database.parts_size[0]));
            if ((position + offset) % 2 == 0) {
                return database.complement[(b >> 4) & 0x0f];
            } else {
                return database.complement[(b & 0x0f)];
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
    public void get_sequence_string(StringBuilder seq, int g, int s, int p, int l, boolean direction) {
        int i;
        seq.setLength(0);
        if (p >= 0 && p + l <= database.sequence_length[g][s]){
            if (direction) {
                for (i = 0; i < l; ++i) {
                    seq.append(get_symbol(g, s, p + i));
                }
            } else {
                for (i = l - 1; i >= 0; --i) {
                    seq.append(get_complement_symbol(g, s, p + i));
                }
            }
        }
    } 
    public void get_sequence_string(StringBuilder seq, int[] adderess, boolean direction) {
        int i, g, s, p, l;
        g = adderess[0];
        s = adderess[1];
        p =adderess[2];
        l = adderess[3] - adderess[2] + 1;
        seq.setLength(0);
        if (p >= 0 && p + l <= database.sequence_length[g][s]){
            if (direction) {
                for (i = 0; i < l; ++i) {
                    seq.append(get_symbol(g, s, p + i));
                }
            } else {
                for (i = l - 1; i >= 0; --i) {
                    seq.append(get_complement_symbol(g, s, p + i));
                }
            }
        }
    } 

    public void get_sequence_quality(StringBuilder quality, int g, int s) {
        quality.setLength(0);
        quality.append(database.sequence_qualities[g][s]);
    } 

    public void get_sequence_title(StringBuilder title, int g, int s) {
        title.setLength(0);
        title.append(database.sequence_titles[g][s]);
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
    public boolean compare(int[] a1, int[] a2, int offset1, int offset2, int len, boolean direction) {
        if (a1[2] + offset1 + len - 1 >= database.sequence_length[a1[0]][a1[1]] || 
                a2[2] + offset2 + len - 1 >= database.sequence_length[a2[0]][a2[1]]) {
            return false;
        }
        int i;
        boolean equal;
        if (direction) {
            for (equal = true, i = 0; i < len && equal; ++i) {
                if (get_code(a1[0], a1[1], a1[2] + offset1 + i) != get_code(a2[0], a2[1], a2[2] + offset2 + i)) {
                    equal = false;
                }
            }
        } else {
            for (equal = true, i = 0; i < len && equal; ++i) {
                if (get_code(a1[0], a1[1], a1[2] + offset1 + i) != get_complement_code(a2[0], a2[1], a2[2] + offset2 + len - i - 1)) {
                    equal = false;
                }
            }
        }
        return equal;
    }  
    
    
    /**
     * Makes a kmer located at a specific genomic position.
     * 
     * @param genome The genome number
     * @param sequence The sequence number
     * @param position The position of the kmer in the sequence
     * @return kmer The canonical form of the kmer 
     */
    public kmer make_kmer(int genome, int sequence, int position) {
        int j,fwd_code;
        kmer tmp_kmer=new kmer(curr_kmer);
        for(j = 0; j < K; ++j)
        {
            fwd_code=get_code(genome,sequence,position+j);
            tmp_kmer.next_kmer(fwd_code);
        }   
        return tmp_kmer;             
    }    
}
