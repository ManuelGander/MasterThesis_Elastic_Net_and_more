#!/bin/bash
for k0 in {0..6}; do
    for k1 in {0..6}; do
        for k2 in {0..5}; do
            sbatch runner6.sh "$k0" "$k1" "$k2"
        done
    done
done

##!/bin/bash
#for k0 in {0..6}; do
#    for k1 in {0..6}; do
#        for k2 in {0..5}; do
#            sbatch runner6.sh "$k0" "$k1" "$k2"
#        done
#    done
#done