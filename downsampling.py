import numpy as np
import pandas as pd
import multiprocessing as mp
import collections
rng = np.random.default_rng()

# local imports
import fileio
import filtering
import kmers
import analysis

# logging related
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


NS = [100000, 50000, 25000, 10000, 5000, 2500, 1000]
ITERS = 30

settings = dict(kmer_size=8, weighting='linear', threshold=0.7)



def make_global(allSeqs):
    global allSeqs_global
    allSeqs_global = allSeqs




def downsampler(input_seqs, count=10000):
    newSeqs = {}
    for filename, df in input_seqs.items():

        if not count > df['count'].sum():
            counts = collections.Counter(
                rng.choice(
                    df.seq,
                    size=count,
                    p=df.value,
                    replace=True,
                )
            )
            df_new = pd.DataFrame.from_dict(counts, orient="index", columns=["count"])
            df_new.reset_index(inplace=True)
            df_new.rename(columns={"index": "seq"}, inplace=True)
        
        else:
            df_new = df.copy()

        df_new.drop(df_new.loc[df_new['count'] == 0].index, inplace=True)
        df_new['value'] = df_new['count']/df_new['count'].sum()
        newSeqs[filename] = df_new
        
    return newSeqs


def run(input_seqs, kmer_size=5, weighting="linear", threshold=0.1, count=10000):

    ds = []
    logger.info(f"Running for N={count} with {ITERS} iterations ...")
        
    for i in range(ITERS):

        downsampled_seqs = downsampler(input_seqs, count=count)

        if threshold < 0.1:
            # use only sequences with a value of more than threshold 
            allSeqs = filtering.filter_by_relative_threshold(downsampled_seqs, threshold=threshold)
        else:
            # use most frequent sequences up to a certain coverage
            allSeqs = filtering.filter_by_coverage(downsampled_seqs, threshold=threshold)

        # perform overlap analysis
        overlaps = np.zeros((len(allSeqs), len(allSeqs)))
        for i1, (file1, seqs1) in enumerate(allSeqs.items()):
            for i2, (file2, seqs2) in enumerate(allSeqs.items()):
                overlaps[i1, i2] = kmers.weighted_overlap_kmer(seqs1, seqs2, kmer_size=kmer_size, weighted=weighting)

        full_positive_matrix, full_negative_matrix = analysis.generate_truthtables()
        posdata = np.where(full_positive_matrix, overlaps, np.nan)
        negdata = np.where(full_negative_matrix, overlaps, np.nan)

        threshold_width = np.nanmin(posdata) - np.nanmax(negdata)
        ds.append(threshold_width)

    return ds


def run_parallel(N):
    result = run(allSeqs_global, **settings, count=N)
    return N, result



if __name__ == '__main__':

    # read all data
    allSeqs = fileio.readData("./data/")

    with mp.Pool(4, initializer=make_global, initargs=(allSeqs,)) as pool:
        data = pool.map(run_parallel, NS)

        pd.DataFrame.from_dict(dict(data)).to_csv(f"./results_downsampling/N{ITERS}.csv", sep=",", index=False)
