import gzip
import collections
import pandas as pd

# logging related
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def fastq_to_list(fastq_filename):
    seqs = []
    with gzip.open(fastq_filename, "rt") as f:
        for i,line in enumerate(f):
            if (i + 3) % 4 == 0:
                seqs += [line[0:-1]]
    return seqs



def readData(rootDir):


    files = [
        "Sp1-Infw4-Inrv5-1_S5_L001_R1_001",
        "Sp1-Infw4-Inrv5-2_S6_L001_R1_001",
        "Infw4-Inrv5-1_S1_L001_R1_001",
        "Infw4-Inrv5-2_S2_L001_R1_001",
        "Sp1-amp-Infw4-Inrv5-1_S11_L001_R1_001",
        "Sp1-amp-Infw4-Inrv5-2_S12_L001_R1_001",
        "0_S1_L001_R1_001",
        "1_S2_L001_R1_001",
        "2_S3_L001_R1_001",
        "3_S4_L001_R1_001",
        "Sp1-Infw4-1-Inrv5-1_S9_L001_R1_001",
        "Sp1-Infw4-1-Inrv5-2_S10_L001_R1_001",
        "Sp1-Infw4-Inrv5-2-1_S11_L001_R1_001",
        "Sp1-Infw4-Inrv5-2-2_S12_L001_R1_001",
        "Sp1-Infw4-3-Inrv5-3-1_S13_L001_R1_001",
        "Sp1-Infw4-3-Inrv5-3-2_S14_L001_R1_001",
        "Lib4-Sp1-G2-infw4-4-inrv5-4-1_S13_L001_R1_001",
        "Lib4-Sp1-G2-infw4-4-inrv5-4-2_S14_L001_R1_001",
        "Sp1-Infw4-Inrv8-1_S9_L001_R1_001",
        "Sp1-Infw4-Inrv8-2_S10_L001_R1_001",
        "8_S9_L001_R1_001",
        "9_S10_L001_R1_001",
        "Infw1-Inrv1-2_S3_L001_R1_001",
        "Sp1-Infw1-Inrv1-1_S7_L001_R1_001",
        "Sp1-Infw1-Inrv1-2_S8_L001_R1_001",
        "Sp2-Infw4-Inrv5-1_S1_L001_R1_001",
        "Sp2-Infw4-Inrv5-2_S2_L001_R1_001",
        "Sp2-Infw1-Inrv1-1_S3_L001_R1_001",
        "Sp2-Infw1-Inrv1-2_S4_L001_R1_001",
        "Sp3-Infw4-Inrv5-1_S3_L001_R1_001",
        "Sp3-Infw4-Inrv5-2_S4_L001_R1_001",
        "Sp-e9-Inrv5-Inrv6-1_S5_L001_R1_001",
        "Sp-e9-Inrv5-Inrv6-2_S6_L001_R1_001",
        "Sp-e9-Inrv5-Inrv5-1_S7_L001_R1_001",
        "Sp-e9-Inrv5-Inrv5-2_S8_L001_R1_001",
        "Sp-e10-Infw6-Inrv7-1_S15_L001_R1_001",
        "Sp-e10-Infw6-Inrv7-2_S16_L001_R1_001",
        "Sp-e10-Infw6-Inrv6-1_S13_L001_R1_001",
        "Sp-e10-Infw6-Inrv6-2_S14_L001_R1_001",
    ]

    files = [rootDir + f + ".filtered.fastq.gz" for f in files]

    allseqs = {}
    for i, filename in enumerate(files):
        counter = collections.Counter(fastq_to_list(filename))
        df = pd.DataFrame.from_dict(counter, orient='index').reset_index().rename(columns={'index': 'seq', 0: 'count'})
        df['value'] = df['count']/df['count'].sum()
        allseqs[filename] = df.sort_values(by='value', ascending=False)
        logger.debug(f"Read {i+1}/{len(files)} {filename}.")

    logger.info(f"Read total of {len(allseqs)} files ...")
    return allseqs