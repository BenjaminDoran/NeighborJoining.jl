#!/usr/bin/env bash

# generate trees
for i in {2..10}
do
    # simple generate tree and distance matrix representation
    gotree generate yuletree --seed 42 -l $((2**i)) | gotree matrix -o dist_mtx_$((2**i)).txt;

    # # if instead I want to generate alignments as intermediaries
    # gotree generate yuletree --seed 42 -l $((2**i)) | -o tree_$((2**i)); 
    # seq-gen -mGTR -z 42 -l $((i*10)) < tmp_tree.nw > tmp_align.phy;
    # goalign compute distance -p -t10 go-m f81 -i tmp_align.phy -o dist_mtx_$((2**i)).phy 
done
# clean up tmp files
rm -f tmp_*