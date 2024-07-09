import sys

import pandas as pd
from picked_group_fdr import digest

"""
python create_ibaq_table.py --gene_level --fasta /media/kusterlab/internal_projects/active/TOPAS/Databases/uniprot_proteome_up000005640_03112020.fasta --enzyme trypsinp --output_file /media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/MT/2022.01.10_iBAQ_tables/ibaq_table.tsv
"""

def main(argv):
    args = parseArgs(argv)
    
    parseId = digest.parseUntilFirstSpace
    if args.gene_level:
        parseId = digest.parseGeneNameFunc
    elif args.fasta_use_uniprot_id:
        parseId = digest.parseUniProtId

    minLenIbaq = max([6, args.min_length])
    maxLenIbaq = min([30, args.max_length])

    pre, not_post = digest.getCleavageSites(args.enzyme)

    peptideToProteinMapIbaq = digest.getPeptideToProteinMap(
            args.fasta, db = 'target', digestion = args.digestion, 
            min_len = minLenIbaq, max_len = maxLenIbaq, 
            pre = pre, not_post = not_post, miscleavages = 0,
            methionineCleavage = False, specialAAs = list(args.special_aas),
            parseId = parseId, useHashKey = (args.digestion == "none"))
        
    numIbaqPeptidesPerProtein = digest.getNumPeptidesPerProtein(peptideToProteinMapIbaq)
        
    df = pd.DataFrame(numIbaqPeptidesPerProtein.items(), columns=['Gene name', 'num iBAQ peptides'])
    
    df.to_csv(args.output_file, sep='\t', index=False)


def parseArgs(argv):
    import argparse
    apars = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    apars.add_argument('--fasta', required=True, default=None, metavar = "F",
                         help='''Fasta file to create mapping from peptides to proteins.''')

    apars.add_argument('--output_file', required=True, default=None, metavar = "F",
                         help='''Output tsv file.''')
    
    apars.add_argument('--fasta_use_uniprot_id',
                         help='''Parse protein identifiers in the fasta file as UniProt IDs, 
                                 i.e. Q9UM47 for the protein identifier sp|Q9UM47|NOTC3_HUMAN''',
                         action='store_true')
    
    apars.add_argument('--gene_level',
                         help='''Do quantification on gene-level instead of on protein group level''',
                         action='store_true')                         
    
    digest.addArguments(apars)
                                
    # ------------------------------------------------
    args = apars.parse_args(argv)
    
    return args


if __name__ == "__main__":
    main(sys.argv[1:])
