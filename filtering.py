import numpy as np


def filter_by_relative_threshold(allseqs, threshold=0.01):
    thresholded_seqs = {}
    for file, df in allseqs.items():
        thresholded_seqs[file] = df.drop(df.loc[df.value < threshold].index)
        if df[df.value >= threshold].shape[0] < 1:
            thresholded_seqs[file] = df.sort_values('value', ascending = False).head(1)		

    return thresholded_seqs



def filter_by_coverage(allseqs, threshold=0.5):
    thresholded_seqs = {}
    for file, df in allseqs.items():
        df_sorted = df.sort_values('value', ascending = False)
        cumsum = np.cumsum(df_sorted['value'])
        cov_ix = np.searchsorted(cumsum, threshold, side='right')+1

        thresholded_seqs[file] = df_sorted.head(cov_ix)		

    return thresholded_seqs