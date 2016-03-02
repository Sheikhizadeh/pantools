/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pantools;

import pantools.*;
import java.nio.MappedByteBuffer;

/**
 *
 * @author siavashs
 */
public class pointer {
    public long node_id;
    public int position;
    public long next_index;
    public byte format;
    public pointer()
    {
        node_id=-1L;
        format=-1;
        position=-1;
        next_index=-1L;
    }
    public pointer(long id,byte f,int p,long n)
    {
        node_id=id;
        format=f;
        position=p;
        next_index=n;
    }
    //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    void reset(int d) 
    {
        node_id=-1L;
        format=-1;
        position=-1;
        next_index=-1L;
    }
    
}
