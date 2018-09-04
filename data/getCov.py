infile1 = open('matchedID.txt','r')
infile2 = open('ImputationResults_all_loci_EURandAFR_1000G_finalAD.tab.txt','r')

line = infile2.readline()
line = line.rstrip().split('\t')[1:]

SnpNameToAge = {}
SnpNameToRace = {}
SnpNameToN = {}
SnpNameToT = {}
SnpNameTo061 = {}
SnpNameTo058 = {}
SnpNameTo052 = {}
SnpNameTo056 = {}
SnpNameTo049 = {}
SnpNameTo032 = {}

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
	
	if num > 240:
		HS061 = '1'
	else:
		HS061 = '0'

	if 240 >= num >192:
		HS058 = '1'
	else:
		HS058 = '0'

	if  192 >= num > 144:
		HS052 = '1'
	else:
		HS052 = '0'
	if 144 >= num > 96:
		HS056 = '1'
	else:
		HS056 = '0'		
	if 96 >= num > 48:
		HS049 = '1'
	else:
		HS049 = '0'
	if num <= 48 :
		HS032 = '1'
	else:
		HS032 = '0'		
	SnpNameTo061[SnpName] = HS061
	SnpNameTo058[SnpName] = HS058
	SnpNameTo052[SnpName] = HS052
	SnpNameTo056[SnpName] = HS056
	SnpNameTo049[SnpName] = HS049
	SnpNameTo032[SnpName] = HS032
print SnpNameTo032
outfile = open('cov_all.txt','w')

line1 = 'id' + '\t' + '\t'.join(line) + '\n'
age = []
race = []
T = []
N = []
HS032 = []
HS049 = []
HS056 = []
HS052 = []
HS058 = []
HS061 = []

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
	batch_temp = SnpNameTo056[i]
	HS056.append(batch_temp)
	batch_temp = SnpNameTo052[i]
	HS052.append(batch_temp)
	batch_temp = SnpNameTo058[i]
	HS058.append(batch_temp)
	batch_temp = SnpNameTo061[i]
	HS061.append(batch_temp)

line2 = 'age' + '\t' + '\t'.join(age) + '\n'
line3 = 'race' + '\t' + '\t'.join(race) + '\n'
line4 = 'HS032' + '\t' + '\t'.join(HS032) + '\n'
line5 = 'HS049' + '\t' + '\t'.join(HS049) + '\n'
line6 = 'HS056' + '\t' + '\t'.join(HS056) + '\n'
line7 = 'HS052' + '\t' + '\t'.join(HS052) + '\n'
line8 = 'HS058' + '\t' + '\t'.join(HS058) + '\n'
line9 = 'HS061' + '\t' + '\t'.join(HS061) + '\n'
line10 = 'N' + '\t' + '\t'.join(N) + '\n'
line11 = 'T' + '\t' + '\t'.join(T) + '\n'

outfile.write(line1)
outfile.write(line2)
outfile.write(line3)
outfile.write(line10)
#outfile.write(line11)
#outfile.write(line4)
#outfile.write(line5)
#outfile.write(line6)
#outfile.write(line7)
#outfile.write(line8)
#outfile.write(line9)
