folder=/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.02.28_CJ_batch87
gene=CDKN2A

# getting the tab character in grep from a shell script is a bit tricky: $(printf '\t')
grep "Modified\|$gene;\|$gene$(printf '\t')"  $folder/evidence.txt > $folder/evidence_$gene.txt

# Manual steps:
# 1. Replace empty cells in "Intensity" and "Reporter intensity corrected" with zeroes
# 2. Duplicate "Reporter intensity corrected" columns as "Reporter intensity" and "Reporter intensity count" columns
# 3. Rename "Batch" column to "Experiment"

if true; then
    python -m picked_group_fdr --mq_evidence $folder/evidence_$gene.txt \
        --protein_groups_out $folder/proteinGroups_$gene.txt \
        --fasta /media/kusterlab/internal_projects/active/TOPAS/Databases/uniprot_proteome_up000005640_03112020.fasta \
        --enzyme trypsinp \
        --min-length 6 \
        --cleavages 3 \
        --lfq_min_peptide_ratios 1 \
        --suppress_missing_peptide_warning \
        --methods picked_protein_group_mq_input \
        --do_quant \
        | tee $folder/proteinGroups_$gene.log
fi
