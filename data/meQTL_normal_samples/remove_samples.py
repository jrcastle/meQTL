#!/opt/anaconda/anaconda2/bin/python
import pandas as pd

OUTDIR = './'

infiles = [
    '../cov_all.txt',
    '../meth_all_fixed.txt',
    '../ImputationResults_all_loci_maf_geno_hwe_EURandAFR_1000G_finalAD.tab.txt'
]

outfiles = [
    OUTDIR + 'cov.txt',
    OUTDIR + 'meth.txt',
    OUTDIR + 'ImputationResults_all_loci_maf_geno_hwe_EURandAFR_1000G_finalAD.tab.txt'
]

samples_already_cut = [
    '1551-06G_GC502', '1552-12A_GC518', '1551-12C_GC290', '1546-04D_GC21',
    '1551-03G_GB501', '1549-04G_GC195', '1546-08G_GC42',  '1551-05G_GB502',
    '1546-02D_GC20'
]

  
normal_file = open('../normal_samples.txt', 'r')
normal_samples = normal_file.readlines()

samples_to_keep = normal_samples
for i in range(0, len(samples_to_keep)):
    samples_to_keep[i] = samples_to_keep[i].replace('\n','')
#ENDFOR

samples_to_keep = [x for x in samples_to_keep if x not in samples_already_cut]


i = 0
for f in infiles:

    print "Loading " + f
    first = []
    if f == '../ImputationResults_all_loci_maf_geno_hwe_EURandAFR_1000G_finalAD.tab.txt':
        first = ['IID']
    if f == '../meth_all_fixed.txt':
        first = ['pos']
    if f == '../cov_all.txt':
        first = ['id']

    df = pd.read_table(
        f,
        sep = '\t',
        header = 0,
        usecols = first + samples_to_keep,
        na_filter = False,
        dtype = object
    )

    print "Saving "+ outfiles[i]
    df.to_csv(
        outfiles[i],
        index = False,
        sep = '\t'
    )
    i = i+1

#ENDFOR
