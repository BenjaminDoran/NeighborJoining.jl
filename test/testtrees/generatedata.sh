#!/usr/bin/env bash

# simple generate tree and distance matrix representation
# gotree generate yuletree --seed 42 -l $((2**i)) | gotree matrix -o dist_mtx_$((2**i)).txt;
rm -f tree_varieties.nw
for i in $(seq 3 5 100)
do
    gotree generate yuletree -l$((i)) >> tree_varieties.nw
done
for i in $(seq 3 5 100)
do
    gotree generate uniformtree -l$((i)) >> tree_varieties.nw
done
for i in $(seq 3 5 100)
do
    gotree generate caterpillartree -l$((i)) >> tree_varieties.nw
done
for i in $(seq 2 7)
do
    gotree generate balancedtree -d$((i)) >> tree_varieties.nw
done
for i in $(seq 3 5 100)
do
    gotree generate startree -l$((i)) >> tree_varieties.nw
done



# # if instead I want to generate alignments as intermediaries
# gotree generate yuletree --seed 42 -l $((2**i)) | -o tree_$((2**i)); 
# seq-gen -mGTR -z 42 -l $((i*10)) < tmp_tree.nw > tmp_align.phy;
# goalign compute distance -p -t10 go-m f81 -i tmp_align.phy -o dist_mtx_$((2**i)).phy 