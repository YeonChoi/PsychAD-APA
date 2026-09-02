import pandas as pd
import numpy as np

def make_pas_weights(df_pas_gene, genes_col='mapped_genes', gene_splitter=';'):
    """Compute per-gene PAS order and uniform weight (0-1) along the 3' UTR.
    
    df_pas_gene: PAS info table in bed12 style (cols '0'-'12') + 'mapped_genes' column.
    """
    df = df_pas_gene.copy()
    df['list_mapped_genes'] = df[genes_col].str.split(gene_splitter)
    df_pas_exp = df.explode('list_mapped_genes')
    
    all_gene_names = set(df_pas_exp['list_mapped_genes'])
    print(f'Exploded from {len(df_pas_gene)} lines to {len(df_pas_exp)}, '
          f'for {len(all_gene_names)} genes')
    
    # Strand check
    if (df_pas_exp.groupby('list_mapped_genes')['5'].nunique() != 1).any():
        strand_counts = df_pas_exp.groupby('list_mapped_genes')['5'].nunique()
        mixed_genes = strand_counts[strand_counts != 1].index
        print(df_pas_exp[df_pas_exp['list_mapped_genes'].isin(mixed_genes)])
        raise ValueError(f'{len(mixed_genes)} genes on mixed strands')
    
    # + strand: sort by end ascending (proximal first)
    df_pos = df_pas_exp[df_pas_exp['5']=='+'].sort_values(['list_mapped_genes', '2', '1'])
    # - strand: sort by start descending (proximal = larger coord)
    df_neg = df_pas_exp[df_pas_exp['5']=='-'].sort_values(['list_mapped_genes', '1', '2'], ascending=[True, False, False])
    df_cat = pd.concat([df_pos, df_neg], ignore_index=True)
    
    # Assign order (1-indexed) and weight (0-1 uniform) per gene
    grp = df_cat.groupby('list_mapped_genes', sort=False)
    df_cat['pas_order'] = grp.cumcount() + 1
    n_pas = grp['pas_order'].transform('size')
    df_cat['pas_weight'] = np.where(n_pas == 1, 1.0, (df_cat['pas_order'] - 1) / (n_pas - 1).clip(lower=1))

    df_cat['gene_name'] = df_cat['list_mapped_genes']
    del df_cat['list_mapped_genes']
    
    return df_cat