#!/bin/bash
for k0 in {0..6}; do
    for k1 in {0..6}; do
        for k2 in {0..7}; do
            for k3 in {0..5}; do
                sbatch runner0.sh "$k0" "$k1" "$k2" "$k3"
            done
        done
    done
done

#for k0 in {0..6}; do
#    for k1 in {0..6}; do
#        for k2 in {0..7}; do
#            for k3 in {0..5}; do
#                sbatch runner0.sh "$k0" "$k1" "$k2" "$k3"
#            done
 #       done
#    done
#done