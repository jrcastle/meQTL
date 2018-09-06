#!/opt/anaconda/anaconda2/bin/python
import pandas as pd

OUTDIR = './'

infiles = [
    '../../../cov_all.txt',
    '../../../meth_all_fixed.txt',
    '../../../ImputationResults_bc_loci_AFR_1000G_finalAD_all.tab.txt'
]

outfiles = [
    OUTDIR + 'cov.txt',
    OUTDIR + 'meth.txt',
    OUTDIR + 'ImputationResults_bc_loci_AFR_1000G_finalAD.tab.txt'
]

samples_already_cut = [
    '1551-06G_GC502', '1546-07H_GB48',  '1552-12A_GC518', '1551-12C_GC290',
    '1546-04D_GC21',  '1551-03G_GB501', '1549-04G_GC195', '1546-08G_GC42',
    '1551-05G_GB502', '1546-02D_GC20'
]

samples_already_cut = samples_already_cut + [
    '1546-08B_GC10',  '1546-01G_GB39',  '1550-12F_GC249', '1546-01A_GB1',   '1546-03C_GB15',  '1550-02E_GC238', '1547-02B_GC58',  '1550-08G_GC253',
    '1550-04B_GC218', '1546-10E_GC30',  '1546-03E_GB27',  '1548-08D_GC125', '1552-05A_GB515', '1546-07G_GB42',  '1546-02C_GC14',  '1551-08H_GC509',
    '1547-10A_GC56',  '1551-06C_GC285', '1547-05B_GB60',  '1550-10B_GC222', '1551-11H_GB511', '1546-05H_GB47',  '1549-08B_GC164', '1551-04D_GC293',
    '1548-12D_GC127', '1551-02D_GC291', '1549-04A_GC156', '1550-08C_GC229', '1546-12F_GC38',  '1547-12G_GC95',  '1547-12E_GC81',  '1551-02A_GC266',
    '1552-04A_GC514', '1546-06H_GC47',  '1547-12C_GC69',  '1546-04E_GC27',  '1546-09A_GB5',   '1547-06D_GC72',  '1552-01A_GB512', '1546-11F_GB38',
    '1551-02C_GC283', '1551-07G_GB503', '1546-09E_GB30',  '1549-04H_GC203', '1551-01F_GB344', '1548-06B_GC110', '1546-10C_GC18',  '1550-02C_GC225',
    '1548-04G_GC143', '1549-04F_GC189', '1549-10D_GC179', '1551-02F_GC344', '1549-12D_GC181', '1550-08D_GC235', '1546-07C_GB17',  '1548-04C_GC116',
    '1551-08F_GC348', '1547-04H_GC97',  '1552-08A_GC516', '1550-06B_GC220', '1547-10F_GC87',  '1549-12F_GC193', '1548-06A_GC104', '1549-12G_GC200',
    '1549-06F_GC190', '1548-06G_GC144', '1548-02B_GC108', '1550-04D_GC233', '1546-05F_GB34',  '1549-08C_GC170', '1546-01C_GB14',  '1549-12E_GC187',
    '1549-02D_GC173', '1549-04D_GC174', '1546-09H_GB49',  '1549-10G_GC199', '1547-11A_GB57',  '1549-10B_GC165', '1550-06A_GC210', '1546-09C_GB18',
    '1550-12D_GC237', '1550-02G_GC250', '1550-10A_GC213', '1546-02A_GC1',   '1546-03G_GB40',  '1546-10F_GC37',  '1550-02H_GC257', '1551-02B_GC274',
    '1546-07A_GB4',   '1549-08H_GC205', '1550-10C_GC230', '1549-08F_GC191', '1552-06A_GC515', '1550-06H_GC260', '1547-01B_GB58',  '1546-05E_GB28',
    '1548-02A_GC102', '1546-06D_GC22',  '1546-04G_GC40',  '1551-10G_GC504', '1550-10E_GC242', '1546-05A_GB3',   '1551-02G_GC500', '1546-01E_GB26',
    '1551-06H_GC508', '1548-12G_GC147', '1548-06D_GC124', '1547-03B_GB59',  '1546-11E_GB31',  '1546-01F_GB32',  '1547-10E_GC80',  '1547-06H_GC98',
    '1546-07B_GB10',  '1547-02F_GC82',  '1546-03H_GB46',  '1547-06G_GC92',  '1546-07F_GB35',  '1549-02G_GC194', '1551-11G_GB505', '1546-05G_GB41',
    '1551-07H_GB509', '1551-03H_GB507', '1551-10A_GC272', '1550-06E_GC240', '1546-04A_GC2',   '1551-05H_GB508', '1549-02F_GC188', '1548-04E_GC129',
    '1546-08D_GC23',  '1549-08D_GC178', '1547-08G_GC93',  '1547-02E_GC76',  '1551-01G_GB500', '1547-12D_GC75',  '1547-08E_GC79',  '1550-06C_GC228',
    '1546-12D_GC25',  '1546-03D_GB21',  '1552-10A_GC517', '1550-10D_GC236', '1546-08E_GC29',  '1548-02H_GC148', '1551-02E_GC337', '1550-04C_GC226',
    '1548-08E_GC131', '1550-12C_GC231', '1552-08B_GC522', '1546-08F_GC35',  '1547-04G_GC91',  '1547-09A_GB56',  '1551-08G_GC503', '1548-04F_GC136',
    '1546-09F_GB37',  '1546-04F_GC33',  '1549-10C_GC171', '1548-02C_GC115', '1547-10C_GC68',  '1549-12H_GC207', '1546-05D_GB22',  '1549-10A_GC159'
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
