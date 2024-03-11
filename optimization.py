import numpy as np
import multiprocessing as mp
import itertools

# local imports
import fileio
import filtering
import kmers

# logging related
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)



DATA_DIR = "./data/"

kmer_sizes = [4, 6, 8, 10, 12]
weightings = ['log', 'linear', 'none']
thresholds = [0.01, 0.005, 0.001, 0.5, 0.7, 0.9]



def make_global(allSeqs):
    global allSeqs_global
    allSeqs_global = allSeqs



def run(input_seqs, kmer_size, weighting, threshold):

    if threshold < 0.1:
        # use only sequences with a value of more than threshold 
        allSeqs = filtering.filter_by_relative_threshold(input_seqs, threshold=threshold)
    else:
        # use most frequent sequences up to a certain coverage
        allSeqs = filtering.filter_by_coverage(input_seqs, threshold=threshold)

    # perform overlap analysis
    overlaps = np.zeros((len(allSeqs), len(allSeqs)))
    for i1, (file1, seqs1) in enumerate(allSeqs.items()):
        for i2, (file2, seqs2) in enumerate(allSeqs.items()):
            overlaps[i1, i2] = kmers.weighted_overlap_kmer(seqs1, seqs2, kmer_size=kmer_size, weighted=weighting)


    filename = f"./results/KMR{kmer_size}_WGT{weighting}_TSH{str(threshold).replace('.', '-')}.csv"

    np.savetxt(filename, overlaps, delimiter=",")


def run_parallel(data):
    kmer, weighting, threshold = data
    logger.info(f"Running {kmer}, {weighting}, {threshold}")
    run(allSeqs_global, kmer, weighting, threshold)



if __name__ == '__main__':

    # read all data
    allSeqs = fileio.readData(DATA_DIR)

    # parameter grid
    params_grid = list(itertools.product(kmer_sizes, weightings, thresholds))

    with mp.Pool(24, initializer=make_global, initargs=(allSeqs,)) as pool:
        pool.map(run_parallel, params_grid)
