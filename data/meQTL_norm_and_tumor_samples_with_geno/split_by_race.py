#!/opt/anaconda/anaconda2/bin/python
import pandas as pd

infiles = [
    'cov.txt',
    'meth.txt',                                                                                                                                                                                                                     
    'ImputationResults_EURandAFR_1000G_finalAD.tab.txt'
]

outfiles = [
    'cov',
    'meth',
    'ImputationResults_1000G_finalAD.tab'
]

# GET LIST OF EUR SAMPLES
EUR_file = open('EUR_samples.txt', 'r')
EUR_samples = EUR_file.readlines()
EUR_samples_removed = [
    '1546-01G_GB39', '1546-03C_GB15', '1546-03E_GB27', '1546-07G_GB42', '1546-05H_GB47',
    '1548-07A_GB105', '1551-12C_GC290', '1552-12C_GB711', '1546-09A_GB5', '1546-11F_GB38',
    '1551-01F_GB344', '1546-09D_GB24', '1546-07C_GB17', '1547-03H_GB97', '1546-05A_GB3',
    '1546-04D_GC21', '1547-11E_GB81', '1551-03G_GB501', '1547-09D_GB74', '1546-09C_GB18',
    '1546-07B_GB10', '1546-01C_GB14', '1546-01F_GB32', '1546-09H_GB49', '1546-01A_GB1',
    '1546-09E_GB30', '1549-04G_GC195', '1546-07A_GB4', '1547-11F_GB89', '1547-09G_GB94',
    '1547-01B_GB58', '1546-05E_GB28', '1546-08G_GC42', '1547-09F_GB87', '1546-01E_GB26',
    '1548-07C_GB118', '1547-03B_GB59', '1546-11E_GB31', '1546-05F_GB34', '1546-03H_GB46',
    '1546-07F_GB35', '1546-05G_GB41', '1547-09A_GB56', '1551-05G_GB502', '1548-05D_GB124',
    '1546-03G_GB40', '1547-11A_GB57', '1546-03D_GB21', '1546-09F_GB37', '1547-01E_GB76',
    '1546-05D_GB22', '1547-05B_GB60', '1546-02D_GC20', '1551-06G_GC502'
]
for i in range(0, len(EUR_samples)):
    EUR_samples[i] = EUR_samples[i].replace('\n','')
#ENDFOR

EUR_samples_to_keep = [x for x in EUR_samples if x not in EUR_samples_removed]

# GET LIST OF AFR SAMPLES
AFR_file = open('AFR_samples.txt', 'r')
AFR_samples = AFR_file.readlines()
AFR_samples_removed = [
    '1546-05C_GB16', '1546-07H_GB48', '1551-05F_GB347', '1547-01G_GB90', '1546-01B_GB7',
    '1551-09F_GB349', '1548-11C_GB121', '1546-11G_GB44', '1552-12A_GC518', '1547-01D_GB70',
    '1548-09C_GB120', '1546-11B_GB12'
]
for i in range(0, len(AFR_samples)):
    AFR_samples[i] = AFR_samples[i].replace('\n','')
#ENDFOR

AFR_samples_to_keep = [x for x in AFR_samples if x not in AFR_samples_removed]

i = 0
for f in infiles:
    print "Loading " + f
    first = []
    if f == 'ImputationResults_EURandAFR_1000G_finalAD.tab.txt':
        first = ['IID']
    if f == 'meth.txt':
        first = ['pos']
    if f == 'cov.txt':
        first = ['id']

    df_EUR = pd.read_table(
        f,
        header = 0,
        sep = '\t',
        usecols = first + EUR_samples_to_keep
    )
    df_EUR = df_EUR[ first + EUR_samples_to_keep ]
    
    df_AFR = pd.read_table(
        f,
        header = 0,
        sep = '\t',
        usecols	= first	+ AFR_samples_to_keep
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
