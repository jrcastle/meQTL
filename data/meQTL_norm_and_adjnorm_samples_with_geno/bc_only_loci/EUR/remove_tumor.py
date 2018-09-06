#!/opt/anaconda/anaconda2/bin/python
import pandas as pd

OUTDIR = './'

infiles = [
    '../../../cov_all.txt',
    '../../../meth_all_fixed.txt',
    '../../../ImputationResults_bc_loci_EUR_1000G_finalAD_all.tab.txt'
]

outfiles = [
    OUTDIR + 'cov.txt',
    OUTDIR + 'meth.txt',
    OUTDIR + 'ImputationResults_bc_loci_EUR_1000G_finalAD.tab.txt'
]

samples_already_cut = [
    '1548-10C_GC120', '1551-06G_GC502', '1548-12B_GC113', '1546-07H_GB48',  '1552-12A_GC518',
    '1551-06B_GC277', '1548-06H_GC151', '1546-05C_GB16',  '1551-12C_GC290', '1546-11B_GB12',
    '1546-04D_GC21',  '1546-01B_GB7',   '1551-03G_GB501', '1550-02A_GC208', '1551-06F_GC347',
    '1551-05F_GB347', '1551-10F_GC349', '1547-02H_GC96',  '1547-02D_GC70',  '1546-11G_GB44',
    '1549-04G_GC195', '1550-12H_GC265', '1546-08G_GC42',  '1548-02F_GC135', '1548-02G_GC142',
    '1546-12B_GC12',  '1550-12G_GC256', '1546-06B_GC9',   '1551-05G_GB502', '1550-04H_GC259',
    '1549-04C_GC168', '1547-12H_GC101', '1551-09F_GB349', '1546-02D_GC20'
]

normal_file = open('../../../normal_samples.txt', 'r')
adjnormal_file = open('../../../adjnormal_samples.txt', 'r')
normal_samples = normal_file.readlines()
adjnormal_samples = adjnormal_file.readlines()

samples_to_keep = normal_samples + adjnormal_samples
for i in range(0, len(samples_to_keep)):
    samples_to_keep[i] = samples_to_keep[i].replace('\n','')
#ENDFOR

samples_to_keep = [x for x in samples_to_keep if x not in samples_already_cut]


i = 0
for f in infiles:

    first = []
    first.append( open(f, 'r').readline().split('\t')[0] )
    
    naFlag = True
    if "Imputation" in f:
        naFlag = False
    #ENDIF

    print "Loading " + f
    df = pd.read_table(
        f,
        sep = '\t',
        header = 0,
        usecols = first + samples_to_keep,
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
