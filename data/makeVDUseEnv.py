#!/opt/anaconda/anaconda2/bin/python
VDfile = open('DNAMethylation_VDUse.csv', 'r')
outfile = open('VDUse_all.txt', 'w+')

finalSampleList = []
finalVDUseList  = []
    
for l in VDfile.readlines():
    l = l.replace('\r\n', '')
    sample = l.split(',')[1].replace(" ", "")
    VD_use = l.split(',')[16].replace(" ", "")
    if VD_use.lower() == "yes" or VD_use.lower() == "no":
        finalSampleList.append( sample )
        use = 0
        if VD_use.lower() == "yes":
            use = 1
        #ENDIF
        finalVDUseList.append(use)
    #ENDIF
#ENDFOR


# Write samples
outfile.write('IID')
for i in finalSampleList:
    outfile.write('\t' + i)


# Write VD Use
outfile.write('\n')
outfile.write('VDUse')
for i in finalVDUseList:
    outfile.write('\t' + str(i))
    


