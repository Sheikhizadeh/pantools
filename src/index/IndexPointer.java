/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package index;

/**
 *
 * @author siavashs
 */
public class IndexPointer {

    public long node_id;
    public int position;
    public boolean canonical;
    public long next_index;

    public IndexPointer() {
        node_id = -1L;
        canonical = false;
        position = -1;
        next_index = -1L;
    }

    public IndexPointer(long id, boolean f, int p, long n) {
        node_id = id;
        canonical = f;
        position = p;
        next_index = n;
    }

    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

    void reset(int d) {
        node_id = -1L;
        canonical = false;
        position = -1;
        next_index = -1L;
    }

}
