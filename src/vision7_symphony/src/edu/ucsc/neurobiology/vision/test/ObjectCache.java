package edu.ucsc.neurobiology.vision.test;

import java.io.*;
import java.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ObjectCache {

    public static interface ObjectSource {
        public Object getObject(int objectID) throws IOException;
    }


    private static class CachedObject {
        public final Object object;
        public final long time;


        public CachedObject(Object object) {
            this.object = object;
            this.time = System.currentTimeMillis();
        }
    }


    private ObjectSource objectSource;
    private HashMap<Integer, CachedObject> cache;
    private int maxCacheSize;


    public ObjectCache(ObjectSource objectSource, int maxCacheSize) {
        this.objectSource = objectSource;
        cache = new HashMap(maxCacheSize);
        this.maxCacheSize = maxCacheSize;
    }


    public Object getObject(int objectID) throws NoSuchElementException {
        Integer id = new Integer(objectID);

        if (cache.containsKey(id)) {
            return ( (CachedObject) cache.get(id)).object;
        } else {
            Object o;
            try {
                o = objectSource.getObject(objectID);
            } catch (IOException e) {
                throw new NoSuchElementException("ID " + objectID);
            }
            if (o == null) {
                throw new NoSuchElementException("ID " + objectID);
            }

            putObject(id, o);

            return o;
        }
    }


    private void putObject(Integer newId, Object object) {
        if (cache.size() >= maxCacheSize) {
            long oldestTime = Long.MAX_VALUE;
            Integer oldestID = null;

            for (Integer id : cache.keySet()) {
                long time = cache.get(id).time;
                if (time < oldestTime) {
                    oldestTime = time;
                    oldestID = id;
                }
            }

            cache.remove(oldestID);
        }

        cache.put(newId, new CachedObject(object));
    }
}
