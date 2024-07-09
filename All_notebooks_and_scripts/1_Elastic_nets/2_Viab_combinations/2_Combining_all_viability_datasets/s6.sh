#!/bin/bash
for k0 in {0..6}; do
    for k1 in {0..6}; do
        for k2 in {6..6}; do
            sbatch runner6.sh "$k0" "$k1" "$k2"
        done
    done
done

# 6 6 4

##!/bin/bash
#for k0 in {0..6}; do
#    for k1 in {0..6}; do
#        for k2 in {0..5}; do
#            sbatch runner6.sh "$k0" "$k1" "$k2"
#        done
#    done
#done