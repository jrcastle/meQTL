#!/opt/anaconda/anaconda2/bin/python
import pandas as pd

OUTDIR = 'rm_adjnormal/'

infiles = [
#    '../cov_all.txt',
#    '../meth_all_fixed.txt',
    '../ImputationResults_EURandAFR_1000G_finalAD_all.tab.txt',
]

outfiles = [
#    OUTDIR + 'cov.txt',
#    OUTDIR + 'meth.txt',
    OUTDIR + 'ImputationResults_EURandAFR_1000G_finalAD.tab.txt'
]

samples_already_cut = [
    '1548-09C_GB120', '1551-06G_GC502', '1547-11F_GB89',  '1552-12A_GC518', '1547-09G_GB94',
    '1548-07A_GB105', '1551-12C_GC290', '1546-09D_GB24',  '1548-05D_GB124', '1547-03H_GB97',
    '1546-04D_GC21',  '1547-11E_GB81',  '1551-03G_GB501', '1547-09F_GB87',  '1549-04G_GC195',
    '1546-08G_GC42',  '1548-11C_GB121', '1548-07C_GB118', '1547-01D_GB70',  '1552-12C_GB711',
    '1551-05G_GB502', '1547-01G_GB90',  '1547-09D_GB74',  '1547-01E_GB76',  '1546-02D_GC20'
]



  
normal_file = open('normal_samples.txt', 'r')
tumor_file = open('tumor_samples.txt', 'r')
normal_samples = normal_file.readlines()
tumor_samples = tumor_file.readlines()

samples_to_keep = normal_samples + tumor_samples
for i in range(0, len(samples_to_keep)):
    samples_to_keep[i] = samples_to_keep[i].replace('\n','')
#ENDFOR

samples_to_keep = [x for x in samples_to_keep if x not in samples_already_cut]
print len(samples_to_keep)

i = 0
for f in infiles:

    print "Loading " + f
    first = []
    if f == 'ImputationResults_EURandAFR_1000G_finalAD_all.tab.txt':
        first = ['IID']
    if f == 'meth_all_fixed.txt':
        first = ['pos']
    if f == 'cov_all.txt':
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
