import numpy as np


def toKmers(df, kmer_size=5):
    counter = {}
    for _, row in df.iterrows():
        seq = row['seq']
        value = row['value']
        if kmer_size:
            seq_d = seq*2
            kmers = [seq_d[i:i+kmer_size] for i in range(len(seq))]
        else:
            kmers = [seq]
        
        for kmer in kmers:
            counter[kmer] = counter.get(kmer, 0.0) + value

    return counter


def weighted_overlap_kmer(df1, df2, kmer_size=5, weighted='linear'):
    # options are weighted = 'linear', 'none', 'log'   

    kmers1 = toKmers(df1, kmer_size=kmer_size)
    kmers2 = toKmers(df2, kmer_size=kmer_size)

    common_kmers = kmers1.keys() & kmers2.keys()
    kmers_min = [min(kmers1[kmer], kmers2[kmer]) for kmer in common_kmers]
    
    all_kmers = kmers1.keys() | kmers2.keys()
    kmers_max = [max(kmers1.get(kmer, 0), kmers2.get(kmer, 0)) for kmer in all_kmers]

    if weighted == 'linear':
        # standard weighting
        sum_overlap = sum(kmers_min)
        sum_all = sum(kmers_max)
    
    elif weighted == 'log':
        # logarithmic weighting
        sum_overlap = sum(np.log(kmers_min))
        sum_all = sum(np.log(kmers_max))
    
    elif weighted == 'none':
        # binary weighting
        sum_overlap = len(kmers_min)
        sum_all = len(kmers_max)

    if sum_all == 0:
        return 0
    else:
        return sum_overlap / sum_all
