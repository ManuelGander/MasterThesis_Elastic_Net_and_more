#!/bin/bash
for k0 in {0..0}; do
    for k1 in {0..0}; do
        for k2 in {0..0}; do
            for k3 in {0..0}; do
            sbatch runner8.sh "$k0" "$k1" "$k2" "$k3"
            done
        done
    done
done