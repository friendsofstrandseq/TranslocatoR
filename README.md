# TranslocatoR

TranslocatoR finds translocations in MosaiCatcher-processed data. It can be used for both reciprocal and non-reciprocal translocations.

## Basics
In order to run TranslocatoR you must provide it with a path to the MosaiCatcher output folder that contains your sample(s) of choice. 
TranslocatoR assumes that the input directory has the exact structure of the MosaiCatcher output, so if you want to test something or supply different files, you'll have to mimic the
following structure:  
```
|__<your folder>  
    |__strand_states  
    |   |__<sample ID>  
    |       |__<input file>
    |__segmentation
    |   |__<sample ID>
    |       |__<input file>
    |__counts
        |__<sample ID>
            |__<input file>
```
## Options
You can run TranslocatoR with different options. ```pq``` takes the strand state for
the end of each arm and compares them all to each other. ```majority``` does the same
but for the majority state of each chromosome. ```segments``` automatically identifies
all recurring segments in a library and compares them to each other. This is useful
for very complex events.

## Providing your own files
If you suspect a non-reciprocal translocation, you can pass TranslocatoR a file
```trfile```that contains the manually-identified strand state of the translocation for each
cell. This file must contain a 'sample' column, a 'cell' column, and a 'state' 
column. An example:

| sample | cell | state |
|--------|------|-------|
|RPE-BM510|BM510_20306| C|
|RPE-BM510|BM510_20310| W|
|RPE-BM510|BM510_20315| W|

Each file cannot contain more than one sample, but you can provide different files
for different samples.

If you want to examine one or more specific regions, you can provide a ```regions``` file
containing four columns: 'sample', 'cell', 'start' and 'end'.
This file can contains different regions for different samples, just make sure your sample names match those used by MosaiCatcher.
