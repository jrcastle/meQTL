#!/opt/anaconda/anaconda2/bin/python
import pandas as pd

infiles = [
    'cov.txt',
    'meth.txt',
    'ImputationResults_all_loci_EURandAFR_1000G_finalAD.tab.txt'
]

outfiles = [
    'cov',
    'meth',
    'ImputationResults_all_loci_EURandAFR_1000G_finalAD.tab'
]

# GET LIST OF EUR SAMPLES
EUR_file = open('EUR_samples.txt', 'r')
EUR_samples = EUR_file.readlines()
EUR_samples_removed = [
    '1551-07F_GB348', '1548-05A_GB104', '1547-05F_GB85', '1547-09E_GB80', '1547-03G_GB91',
    '1548-01A_GB102', '1547-05E_GB78', '1547-05C_GB67', '1548-07A_GB105', '1547-07E_GB79',
    '1551-12C_GC290', '1552-12C_GB711', '1547-11C_GB69', '1546-09D_GB24', '1548-09B_GB112',
    '1547-03D_GB71', '1547-03H_GB97', '1546-04D_GC21', '1547-11E_GB81', '1551-03G_GB501',
    '1546-03B_GB8', '1547-09D_GB74', '1548-01B_GB108', '1547-03F_GB84', '1548-07D_GB125',
    '1547-05G_GB92', '1549-04G_GC195', '1548-09A_GB106', '1547-11F_GB89', '1547-09G_GB94',
    '1546-08G_GC42', '1548-03B_GB109', '1548-03A_GB103', '1547-09F_GB87', '1548-07C_GB118',
    '1547-03E_GB77', '1546-01D_GB20', '1547-09C_GB68', '1548-03D_GB123', '1548-01C_GB115',
    '1547-11D_GB75', '1551-05G_GB502', '1548-03C_GB116', '1548-05D_GB124', '1548-07B_GB111',
    '1548-05B_GB110', '1552-01C_GB700', '1547-01E_GB76', '1546-02D_GC20', '1551-06G_GC502'
]

for i in range(0, len(EUR_samples)):
    EUR_samples[i] = EUR_samples[i].replace('\n','')
#ENDFOR

EUR_samples_to_keep = [x for x in EUR_samples if x not in EUR_samples_removed]

# GET LIST OF AFR SAMPLES
AFR_file = open('AFR_samples.txt', 'r')
AFR_samples = AFR_file.readlines()
AFR_samples_removed = [
    '1548-11B_GB113', '1546-07H_GB48', '1547-01G_GB90', '1548-11A_GB107', '1548-11C_GB121',
    '1546-01H_GB45', '1547-01H_GB96', '1547-09B_GB62', '1552-12A_GC518', '1547-01D_GB70',
    '1548-09C_GB120', '1547-03C_GB65', '1547-11B_GB63', '1547-11H_GB101'
]

for i in range(0, len(AFR_samples)):
    AFR_samples[i] = AFR_samples[i].replace('\n','')
#ENDFOR

AFR_samples_to_keep = [x for x in AFR_samples if x not in AFR_samples_removed]

i = 0
for f in infiles:
    print "Loading " + f
    first  = []
    naflag = True
    if f == 'ImputationResults_EURandAFR_1000G_finalAD.tab.txt':
        first = ['IID']
        naflag = False
    if f == 'meth.txt':
        first = ['pos']
    if f == 'cov.txt':
        first = ['id']

    df_EUR = pd.read_table(
        f,
        header = 0,
        sep = '\t',
        usecols = first + EUR_samples_to_keep,
        na_filter = naflag
    )
    df_EUR = df_EUR[ first + EUR_samples_to_keep ]
    
    df_AFR = pd.read_table(
        f,
        header = 0,
        sep = '\t',
        usecols	= first	+ AFR_samples_to_keep,
        na_filter = naflag
    )
    df_AFR = df_AFR[ first + AFR_samples_to_keep ]

    df_EUR.to_csv(
        outfiles[i] + '_EUR.txt',
        sep = '\t',
        index = False
    )

    df_AFR.to_csv(
	outfiles[i] + '_AFR.txt',
	sep = '\t',
	index =	False
    )

    i = i+1
#END FOR
