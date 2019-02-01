/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sequence;

import index.IndexDatabase;
import index.IndexScanner;
import index.kmer;
import static pantools.Pantools.DEBUG;
import static pantools.Pantools.complement;
import static pantools.Pantools.sym;

/**
 *
 * @author siavash
 */
public class SequenceScanner {
    private int genome;
    private int sequence;
    private int position;
    private int K;
    private kmer curr_kmer;
    private long curr_index;
    SequenceDatabase database;
    
    public SequenceScanner(SequenceDatabase db, int g, int s, int k, int pre_len){
        genome = g;
        sequence = s;
        database = db;
        position = 0;
        K = k;
        curr_kmer = new kmer(K, pre_len);
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
    public String get_title(){
        return database.sequence_titles[genome][sequence];
    }    
    public long get_offset(){
        return database.sequence_offset[genome][sequence];
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
            b = database.genomes_buff[(int) (pos / database.MAX_BYTE_COUNT)].get((int) (pos % database.MAX_BYTE_COUNT));
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
        position = 0;
    }
    public void set_position(int p){
        position = p;
    } 
    public void set_curr_index(long ci){
     curr_index = ci;   
    }

    public boolean end_of_scan(){
        return genome > database.num_genomes;
    }      
    public boolean end_of_genome(){
        return sequence > database.num_sequences[genome];
    }    
    
    public void next_genome(){
        ++genome;
        sequence = 1;
    }
    public void next_sequence(){
        ++sequence;
        position = 0;
    }
    public void next_position(){
        ++position;
    }
    public void previous_position(){
        --position;
    }

    /*
    * get a kmer from the sequence starting at "start".
    * @param start The start position of the kmer
    * @retutn The position of the first base after the kmer or
    * MININT if the is no kmer at that position or
    * The negative value of current position if it passess a degenerate region.
    */
    public int initialize_kmer(int start) {
        if (DEBUG) System.out.println("initialize_kmer at " + start);
        int i;
        curr_kmer.reset();
        position = start - 1;
        for (i = 0; i < K && position < get_sequence_length() - 1; ++i) {
            next_position();
            if (get_code(0) > 3) {
                if (DEBUG) System.out.println("jump_forward");
                if (jump_forward())
                    return - (position - K + 1); // success but passed a degenerate region
                else
                    return Integer.MIN_VALUE; // failed
            } else
                curr_kmer.next_kmer(get_code(0));
        }
        if (DEBUG) System.out.println(curr_kmer.toString());
        if (i == K)
            return position - K + 1; // success
        else 
            return Integer.MIN_VALUE; // failed
    }

    public int next_kmer() {
        if (position < get_sequence_length() - 1){
            next_position();
            if (get_code(0) > 3){
                if (jump_forward())
                    return - (position - K + 1); // success but passed a degenerate region
                else
                    return Integer.MIN_VALUE; // failed
            } else {
                curr_kmer.next_kmer(get_code(0));
                return position - K + 1;
            }
        } else
            return Integer.MIN_VALUE; // failed
    }

    /**
     * Jump over an ambiguous region; at first, position points to the first position which degenerate starts, 
     * after jumping it points to the last base of the first K-mer after the ambiguous region. 
     */
    public boolean jump_forward() {
        int j;
        int base_code = get_code(0);
        do {
            curr_kmer.reset();
            while (base_code > 3 && position < get_sequence_length() - 1) {
                next_position();
                base_code = get_code(0);
            }
            curr_kmer.next_kmer(base_code);
            for (j = 0; j < K - 1 && position < get_sequence_length() - 1; ++j) {
                next_position();
                base_code = get_code(0);
                if (base_code > 3) {
                    break;
                }
                curr_kmer.next_kmer(base_code);
            }
        } while (base_code > 3 && position < get_sequence_length() - 1);
        if (DEBUG) System.out.println(curr_kmer.toString());
        return j == K - 1; // got valid kmer
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
            b = database.genomes_buff[(int) (pos / database.MAX_BYTE_COUNT)].get((int) (pos % database.MAX_BYTE_COUNT));
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
        if (p < database.sequence_length[g][s]) {
            byte b;
            long pos = database.sequence_start[g][s] + p / 2;
            b = database.genomes_buff[(int) (pos / database.MAX_BYTE_COUNT)].get((int) (pos % database.MAX_BYTE_COUNT));
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
        if (p < database.sequence_length[g][s] && p > -1) {
            byte b;
            long pos = database.sequence_start[g][s] + p / 2;
            b = database.genomes_buff[(int) (pos / database.MAX_BYTE_COUNT)].get((int) (pos % database.MAX_BYTE_COUNT));
            if (p % 2 == 0) {
                return (b >> 4) & 0x0f;
            } else {
                return (b & 0x0f);
            }
        } else {
            System.out.println(p + " is not in range 0.." + database.sequence_length[g][s]);
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
            b = database.genomes_buff[(int) (pos / database.MAX_BYTE_COUNT)].get((int) (pos % database.MAX_BYTE_COUNT));
            if (p % 2 == 0) {
                return complement[(b >> 4) & 0x0f];
            } else {
                return complement[(b & 0x0f)];
            }
        } else {
            return -1;
        }
    }

    public int get_complement_current_code(int offset) {
        if (position + offset < database.sequence_length[genome][sequence]) {
            byte b;
            long pos = database.sequence_start[genome][sequence] + (position + offset) / 2;
            b = database.genomes_buff[(int) (pos / database.MAX_BYTE_COUNT)].get((int) (pos % database.MAX_BYTE_COUNT));
            if ((position + offset) % 2 == 0) {
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
    public void get_sub_sequence(StringBuilder seq, int g, int s, int p, int l, boolean direction) {
        int i;
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
        } else
            System.err.println("Reading out of range!");
    } 

    public void get_sub_sequence(StringBuilder seq, int[] adderess, boolean direction) {
        int i, g, s, p, l;
        g = adderess[0];
        s = adderess[1];
        p =adderess[2];
        l = adderess[3] - adderess[2] + 1;
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
        } else
            System.err.println("Reading out of range!");
    } 

    public void get_complete_sequence(StringBuilder seq, int g, int s, boolean direction) {
        int i, l;
        l = (int)database.sequence_length[g][s];
        if (direction) {
            for (i = 0; i < l; ++i) {
                seq.append(get_symbol(g, s, i));
            }
        } else {
            for (i = l - 1; i >= 0; --i) {
                seq.append(get_complement_symbol(g, s, i));
            }
        }
    } 
    
    public void get_complete_sequence(StringBuilder seq, int[] adderess, boolean direction) {
        int i, g, s, p, l;
        g = adderess[0];
        s = adderess[1];
        l = (int)database.sequence_length[g][s];
        if (direction) {
            for (i = 0; i < l; ++i) {
                seq.append(get_symbol(g, s, i));
            }
        } else {
            for (i = l - 1; i >= 0; --i) {
                seq.append(get_complement_symbol(g, s, i));
            }
        }
    } 
    
    public void get_sequence_quality(StringBuilder quality, int g, int s) {
        //quality.append(database.sequence_qualities[g][s]);
    } 

    public void get_sequence_title(StringBuilder title, int g, int s) {
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
        tmp_kmer.reset();
        for(j = 0; j < K; ++j) {
            fwd_code=get_code(genome, sequence, position+j);
            tmp_kmer.next_kmer(fwd_code);
        }   
        return tmp_kmer;             
    }   
    
    public long find_curr_kmer(IndexScanner inx){
        curr_index = inx.find(curr_kmer);
        return curr_index;
    }
        
}
