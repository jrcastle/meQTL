#!/opt/anaconda/anaconda2/bin/python
import os
fSNPname = "SNPs.txt"
fSNPlocName = "SNPloc.txt"

fSNPs = open(fSNPname, "r")
if os.path.isfile( fSNPlocName ):
    command = "rm " + fSNPlocName
    os.system( command )
#ENDIF

fOUT  = open( fSNPlocName, "w+" )

fOUT.write("snp\tchr\tpos\n")

for line in fSNPs.readlines():
    line = line.replace("\n", "")
    chr = line.split(":")[0]
    pos = line.split(":")[1].split("_")[0]
    fOUT.write(str(line) + "\t" + 'chr' + str(chr) + "\t" + str(pos) + "\n")

