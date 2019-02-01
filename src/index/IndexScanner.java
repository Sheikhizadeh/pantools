/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package index;

import static index.kmer.compare_suffix;
import java.io.IOException;
import java.nio.MappedByteBuffer;

/**
 * Implements all the functionality to work with a KMC-based kmer index database. 
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public final class IndexScanner {

    public kmer key;
    public IndexDatabase database;
    public int suffix_record_size;
    public int full_suffix_page_size;
    public int full_pointer_page_size;

    public IndexScanner(IndexDatabase db){
        database = db;   
        suffix_record_size = db.ctr_size + db.suf_len / 4; // in bytes
        full_suffix_page_size = db.MAX_BYTE_COUNT / suffix_record_size * suffix_record_size ; // in bytes
        full_pointer_page_size = db.MAX_BYTE_COUNT / db.poniter_length * db.poniter_length; // in bytes
        key=new kmer(db.K, db.pre_len);
    }
    
    /**
     * Gives the length of index. 
     * 
     * @return Number of kmers
     */
    public long length() {
        return database.kmers_num;
    }

    /**
     * Gives the prefix length of kmers.
     * 
     * @return Prefix length of kmers
     */
    public int get_pre_len() {
        return database.pre_len;
    }

    /**
     * Gives the suffix length of kmers.
     * 
     * @return suffix length of kmers
     */
    public int get_suf_len() {
        return database.suf_len;
    }
    
    /**
     * Reads a long integer from the memory mapped buffer.
     * 
     * @param buff The file associated to the byte stream
     * @return The long integer value
     * @throws IOException 
     */
    private int read_offset(MappedByteBuffer buff, int pos) {
        int i, offset = 0;
            for (i = 0; i < database.offset_len; ++i, ++pos) {
                offset = offset << 8;
                offset = offset | (int)(buff.get(pos) & 0x00FF);
            }
        return offset == (1l << 8 * database.offset_len) - 1 ? -1 : offset;
    }
    
    /**
     * Reads a long integer from the memory mapped buffer.
     * 
     * @param buff The file associated to the byte stream
     * @return The long integer value
     * @throws IOException 
     */
    private long read_node_id(MappedByteBuffer buff, int pos) {
        long i, id = 0;
            for (i = 0; i < database.id_len; ++i, ++pos) {
                id = id << 8;
                id = id | (long)(buff.get(pos) & 0x00FF);
            }
        return id == (1l << 8 * database.id_len) - 1 ? -1l : id;
    }

    /**
     * Reads a long integer from the memory mapped buffer.
     * 
     * @param buff The file associated to the byte stream
     * @return The long integer value
     * @throws IOException 
     */
    private void write_offset(MappedByteBuffer buff, int pos, int offset) {
        for (int i = database.offset_len - 1; i >= 0; --i, ++pos) {
            buff.put(pos, (byte)((offset >> (8 * i)) & 0x00FF));
        }
    }
    
    /**
     * Reads a long integer from the memory mapped buffer.
     * 
     * @param buff The file associated to the byte stream
     * @return The long integer value
     * @throws IOException 
     */
    private void write_node_id(MappedByteBuffer buff, int pos, long id) {
        for (int i = database.id_len - 1; i >= 0; --i, ++pos) {
            buff.put(pos, (byte)((id >> (8 * i)) & 0x00FF));
        }
    }    

    /**
     * Reads a kmer from th index database.
     * 
     * @param number The number of kmer in the index
     * @return The kmer object
     */
    public void get_kmer(kmer key, long number) {
        boolean found = false;
        int i, low = 0, high = database.prefix_ptr.length - 1, mid = 0, pos;
        while (low <= high && !found) {
            mid = (low + high) / 2;
            if (number < database.prefix_ptr[mid]) {
                high = mid - 1;
            } else if (number > database.prefix_ptr[mid]) {
                low = mid + 1;
            } else {
                found = true;
            }
        }
        if (!found) {
            mid = high;
        } else {
            for (; mid < database.prefix_ptr.length - 1 && database.prefix_ptr[mid + 1] == database.prefix_ptr[mid]; ++mid);
        }
        key.set_fwd_prefix(mid);
        pos = (int) (number * suffix_record_size % full_suffix_page_size);
        get_suffix(database.suf_buff[(int) (number * suffix_record_size / full_suffix_page_size)], pos, key.get_fwd_suffix());
    }
    
    public void get_suffix(MappedByteBuffer buff, int pos, byte[] des){
        for (int i = 0; i < database.suf_len / 4; ++i, ++pos)
            des[i] = buff.get(pos);
    }

    public int get_kmer_count(long number) {
        database.suf_buff[(int) (number * suffix_record_size / full_suffix_page_size)].position((int) (number * suffix_record_size % full_suffix_page_size) + database.suf_len / 4);
        return database.suf_buff[(int) (number * suffix_record_size / full_suffix_page_size)].get() & 0x00FF;
    }

    /**
     * Reads a pointer from the kmer index.
     * @param poniter The pointer
     * @param number Number of the pointer.
     */
    public void get_pointer(IndexPointer poniter, long number){
        MappedByteBuffer buff = database.ptr_buff[(int) (number * database.poniter_length / full_pointer_page_size)];
        int position = (int)(number * database.poniter_length % full_pointer_page_size);
        poniter.node_id = read_node_id(buff, position);
        position += database.id_len;
        poniter.offset = read_offset(buff, position);
        position += database.offset_len;
        poniter.canonical = buff.get(position) == 0;
        position += 1;
        poniter.next_index = read_node_id(buff, position);
    }
    
    /**
     * Writes a pointer in the index database
     * @param poniter The pointer
     * @param number Number of the pointer.
     */
    public void put_pointer(IndexPointer poniter, long number){
        MappedByteBuffer buff = database.ptr_buff[(int) (number * database.poniter_length / full_pointer_page_size)];
        int position = (int)(number * database.poniter_length % full_pointer_page_size);
        write_node_id(buff, position, poniter.node_id);
        position += database.id_len;
        write_offset(buff, position, poniter.offset);
        position += database.offset_len;
        buff.put(position, (byte)(poniter.canonical ? 0 : 1));
        position += 1;
        write_node_id(buff, position, poniter.next_index);
    }

    public void put_pointer(long node_id, int offset, boolean canonical,long next_index, long number){
        MappedByteBuffer buff = database.ptr_buff[(int) (number * database.poniter_length / full_pointer_page_size)];
        int position = (int)(number * database.poniter_length % full_pointer_page_size);
        write_node_id(buff, position, node_id);
        position += database.id_len;
        write_offset(buff, position, offset);
        position += database.offset_len;
        buff.put(position, (byte)(canonical ? 0 : 1));
        position += 1;
        write_node_id(buff, position, next_index);
    }
    
    /**
     * Reads the node id of a kmer from the index database.
     * @param number Number of the kmer in the index
     * @return The id of node
     */
    public long get_node_id(long number) {
        MappedByteBuffer buff = database.ptr_buff[(int) (number * database.poniter_length / full_pointer_page_size)];
        return read_node_id(buff, (int)(number * database.poniter_length % full_pointer_page_size));
    }

    /**
     * Writes the node id of a kmer into the index database.
     * @param node_id The node id to be written
     * @param number Number of the kmer in the index
     */
    public void put_node_id(long node_id, long number) {
        MappedByteBuffer buff = database.ptr_buff[(int) (number * database.poniter_length / full_pointer_page_size)];
        write_node_id(buff, (int)(number * database.poniter_length % full_pointer_page_size), node_id);
    }

    /**
     * Reads the position of a kmer in the node where it occurs.
     * @param number Number of the kmer in the index
     * @return The position in the node
     */
    public int get_position(long number){
        MappedByteBuffer buff = database.ptr_buff[(int) (number * database.poniter_length / full_pointer_page_size)];
        return read_offset(buff, (int)(number * database.poniter_length % full_pointer_page_size) + database.id_len);
    }

    /**
     * Writes the position of a kmer in the node where it occurs.
     * @param position The node id to be written
     * @param number Number of the kmer in the index
     */
    public void put_position(int position, long number) {
        MappedByteBuffer buff = database.ptr_buff[(int) (number * database.poniter_length / full_pointer_page_size)];
        write_offset(buff, (int) (number * database.poniter_length % full_pointer_page_size) + database.id_len, position);
    }

    /**
     * Looks if a kmer has been canonical at the first visit.
     * @param number Number of the kmer in the index
     * @return The canonical status
     */
    public boolean get_canonical(long number) {
        MappedByteBuffer buff = database.ptr_buff[(int) (number * database.poniter_length / full_pointer_page_size)];
        return buff.get((int)(number * database.poniter_length % full_pointer_page_size) + database.id_len + database.offset_len) == 0;
    }

    /**
     * Writes the canonical status of a kmer at the first visit.
     * @param canonical The canonical status to be written ( 0 if it was canonical, 1 otherwise )
     * @param number Number of the kmer in the index
     */
    public void put_canonical(boolean canonical, long number) {
        MappedByteBuffer buff = database.ptr_buff[(int) (number * database.poniter_length / full_pointer_page_size)];
        buff.put((int)(number * database.poniter_length % full_pointer_page_size) + database.id_len + database.offset_len, (byte) (canonical ? 0 : 1));
    }
    
    public long get_next_index(long number) {
        MappedByteBuffer buff = database.ptr_buff[(int) (number * database.poniter_length / full_pointer_page_size)];
        return read_node_id(buff, (int)(number * database.poniter_length % full_pointer_page_size) + database.id_len + database.offset_len + 1);
    }

    public void put_next_index(long next_index, long number) {
        MappedByteBuffer buff = database.ptr_buff[(int) (number * database.poniter_length / full_pointer_page_size)];
        write_node_id(buff, (int)(number * database.poniter_length % full_pointer_page_size) + database.id_len + database.offset_len + 1, next_index);
    }  
    
    /**
     * Finds the number of a kmer in the database.
     * 
     * @param k_mer The kmer object
     * @return The number of kmer
     */
    public long find(kmer k_mer) {
        long low, mid, high;
        int j, comp, prefix = k_mer.get_canonical_prefix(), pos;
        byte[] suffix = k_mer.get_canonical_suffix();
        low = database.prefix_ptr[prefix];
        if (prefix == database.prefix_ptr.length - 1) {
            high = database.kmers_num - 1;
        } else {
            high = database.prefix_ptr[prefix + 1] - 1;
        }
        while (low <= high) {
            mid = (low + high) / 2;
            j = (int) (mid * suffix_record_size / full_suffix_page_size);
            pos = (int) (mid * suffix_record_size % full_suffix_page_size);
            get_suffix(database.suf_buff[(int) (mid * suffix_record_size / full_suffix_page_size)], pos, key.get_canonical_suffix());
            comp = compare_suffix(suffix, key.get_canonical_suffix());
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
    
    public int get_K(){
        return database.K;
    }
    
    public int get_pointer_length(){
        return database.poniter_length;
    }
}
