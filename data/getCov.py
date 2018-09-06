#!/opt/anaconda/anaconda2/bin/python

infile1 = open('matchedID.txt','r')
infile2 = open('ImputationResults_all_loci_EURandAFR_1000G_finalAD.tab.txt','r')

line = infile2.readline()
line = line.rstrip().split('\t')[1:]

SnpNameToAge = {}
SnpNameToRace = {}
SnpNameToN = {}
SnpNameToT = {}
SnpNameTo032 = {}
SnpNameTo049 = {}
SnpNameTo052 = {}
SnpNameTo056 = {}
SnpNameTo058 = {}
SnpNameTo061 = {}
SnpNameTo066 = {}
SnpNameTo067 = {}
SnpNameTo070 = {}
SnpNameTo076 = {}
SnpNameTo087 = {}
SnpNameTo089 = {}
SnpNameTo091 = {}

while True:
	line1 = infile1.readline()
	if line1 == '':
		break
	line1 = line1.rstrip().split('\t')
	SnpName = line1[6]
	age = line1[3]
	race = line1[4]
	num = float(line1[0])

	SnpNameToAge[SnpName] = age
	SnpNameToRace[SnpName] = race

	Type = line1[2]

	if 'N' in Type:
		SnpNameToN[SnpName] = '1'
	else:
		SnpNameToN[SnpName] = '0'
	if 'T' in Type:
		SnpNameToT[SnpName] = '1'
	else:
		SnpNameToT[SnpName] = '0'

        if num <= 48 :
                HS032 = '1'
        else:
                HS032 = '0'
                
        if 96 >= num > 48:
                HS049 = '1'
        else:
                HS049 = '0'

        if 144 >= num > 96:
                HS056 = '1'
        else:
                HS056 = '0'

        if  192 >= num > 144:
                HS052 = '1'
        else:
                HS052 = '0'

        if 240 >= num > 192:
                HS058 = '1'
        else:
                HS058 = '0'

        if  288 >= num > 240:
                HS061 = '1'
        else:
                HS061 = '0'

        if  336 >= num > 288:
                HS066 = '1'
        else:
                HS066 = '0'

        if  384 >= num > 336:
                HS067 = '1'
        else:
                HS067 = '0'

        if  432 >= num > 384:
		HS070 = '1'
	else:
             	HS070 = '0'

        if  480 >= num > 432:
                HS076 = '1'
        else:
                HS076 = '0'

        if  528 >= num > 480:
                HS087 = '1'
        else:
                HS087 = '0'

        if  576 >= num > 528:
		HS089 = '1'
	else:
                HS089 = '0'                

        if  624 >= num > 576:
                HS091 = '1'
        else:
                HS091 = '0'
                
	SnpNameTo032[SnpName] = HS032
        SnpNameTo049[SnpName] = HS049
        SnpNameTo052[SnpName] = HS052
        SnpNameTo056[SnpName] = HS056
        SnpNameTo058[SnpName] = HS058
        SnpNameTo061[SnpName] = HS061
        SnpNameTo066[SnpName] = HS066
        SnpNameTo067[SnpName] = HS067
        SnpNameTo070[SnpName] = HS070
        SnpNameTo076[SnpName] = HS076
        SnpNameTo087[SnpName] = HS087
        SnpNameTo089[SnpName] = HS089
        SnpNameTo091[SnpName] = HS091
        
outfile = open('cov_all.txt','w')

line1 = 'id' + '\t' + '\t'.join(line) + '\n'
age = []
race = []
T = []
N = []
HS032 = []
HS049 = []
HS052 = []
HS056 = []
HS058 = []
HS061 = []
HS066 = []
HS067 = []
HS070 = []
HS076 = []
HS087 = []
HS089 = []
HS091 = []

for i in line:
	age_temp = SnpNameToAge[i]
	age.append(age_temp)
	race_temp = SnpNameToRace[i]
	race.append(race_temp)
	T_temp = SnpNameToT[i]
	T.append(T_temp)
	N_temp = SnpNameToN[i]
	N.append(N_temp)
	batch_temp = SnpNameTo032[i]
	HS032.append(batch_temp)
	batch_temp = SnpNameTo049[i]
	HS049.append(batch_temp)
        batch_temp = SnpNameTo052[i]
        HS052.append(batch_temp)
	batch_temp = SnpNameTo056[i]
	HS056.append(batch_temp)
	batch_temp = SnpNameTo058[i]
	HS058.append(batch_temp)
	batch_temp = SnpNameTo061[i]
	HS061.append(batch_temp)
        batch_temp = SnpNameTo066[i]
        HS066.append(batch_temp)
        batch_temp = SnpNameTo067[i]
        HS067.append(batch_temp)
        batch_temp = SnpNameTo070[i]
        HS070.append(batch_temp)
        batch_temp = SnpNameTo076[i]
        HS076.append(batch_temp)
        batch_temp = SnpNameTo087[i]
        HS087.append(batch_temp)
        batch_temp = SnpNameTo089[i]
        HS089.append(batch_temp)
        batch_temp = SnpNameTo091[i]
        HS091.append(batch_temp)
        
line2 = 'age' + '\t' + '\t'.join(age) + '\n'
line3 = 'race' + '\t' + '\t'.join(race) + '\n'
line4 = 'N' + '\t' + '\t'.join(N) + '\n'
line5 = 'T' + '\t' + '\t'.join(T) + '\n'
line6  = 'HS032' + '\t' + '\t'.join(HS032) + '\n'
line7  = 'HS049' + '\t' + '\t'.join(HS049) + '\n'
line8  = 'HS052' + '\t' + '\t'.join(HS052) + '\n'
line9  = 'HS056' + '\t' + '\t'.join(HS056) + '\n'
line10 = 'HS058' + '\t' + '\t'.join(HS058) + '\n'
line11 = 'HS061' + '\t' + '\t'.join(HS061) + '\n'
line12 = 'HS066' + '\t' + '\t'.join(HS066) + '\n'
line13 = 'HS067' + '\t' + '\t'.join(HS067) + '\n'
line14 = 'HS070' + '\t' + '\t'.join(HS070) + '\n'
line15 = 'HS076' + '\t' + '\t'.join(HS076) + '\n'
line16 = 'HS087' + '\t' + '\t'.join(HS087) + '\n'
line17 = 'HS089' + '\t' + '\t'.join(HS089) + '\n'
line18 = 'HS091' + '\t' + '\t'.join(HS091) + '\n'

outfile.write(line1)
outfile.write(line2)
outfile.write(line3)
outfile.write(line4)
outfile.write(line5)
outfile.write(line6)
outfile.write(line7)
outfile.write(line8)
outfile.write(line9)
outfile.write(line10)
outfile.write(line11)
outfile.write(line12)
outfile.write(line13)
outfile.write(line14)
outfile.write(line15)
outfile.write(line16)
outfile.write(line17)
outfile.write(line18)
