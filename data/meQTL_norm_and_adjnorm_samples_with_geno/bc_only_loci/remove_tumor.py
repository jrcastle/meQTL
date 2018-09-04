#!/opt/anaconda/anaconda2/bin/python
import pandas as pd

OUTDIR = './'

infiles = [
    '../../cov_all.txt',
    '../../meth_all_fixed.txt',
    '../../ImputationResults_bc_loci_EURandAFR_1000G_finalAD_all.tab.txt'
]

outfiles = [
    OUTDIR + 'cov.txt',
    OUTDIR + 'meth.txt',
    OUTDIR + 'ImputationResults_bc_loci_EURandAFR_1000G_finalAD.tab.txt'
]

samples_already_cut = [
    '1551-06G_GC502', '1546-07H_GB48', '1552-12A_GC518', '1551-12C_GC290', '1546-04D_GC21',
    '1551-03G_GB501', '1549-04G_GC195', '1546-08G_GC42', '1551-05G_GB502', '1546-02D_GC20'
]

first = [
    ['id'],
    ['pos'],
    ['IID']
]
  
normal_file = open('../../normal_samples.txt', 'r')
adjnormal_file = open('../../adjnormal_samples.txt', 'r')
normal_samples = normal_file.readlines()
adjnormal_samples = adjnormal_file.readlines()

samples_to_keep = normal_samples + adjnormal_samples
for i in range(0, len(samples_to_keep)):
    samples_to_keep[i] = samples_to_keep[i].replace('\n','')
#ENDFOR

samples_to_keep = [x for x in samples_to_keep if x not in samples_already_cut]


i = 0
for f in infiles:

    if i < 2:
        i = i + 1 
        continue
    #ENDIF
    
    print "Loading " + f
    naFlag = True
    if f == infiles[2]:
        naFlag = False
    #ENDIF
    
    df = pd.read_table(
        f,
        sep = '\t',
        header = 0,
        usecols = first[i] + samples_to_keep,
        na_filter = naFlag,
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
