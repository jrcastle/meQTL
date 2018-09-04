#!/opt/anaconda/anaconda2/bin/python
imp_file = 'test1.txt'
cov_file = 'test2.txt'
meth_file = 'test3.txt'

f_imp = open(imp_file, 'r')
f_cov = open(cov_file, 'r')
f_meth = open(meth_file, 'r')

imp_samples = f_imp.readlines()
cov_samples = f_cov.readlines()
meth_samples = f_meth.readlines()

imp_samples = imp_samples[1:]
cov_samples = cov_samples[1:]
meth_samples = meth_samples[1:]

for i in range(0, len(imp_samples)):
    imp_samples[i] = imp_samples[i].replace('\n','')
    cov_samples[i] = cov_samples[i].replace('\n','')
    meth_samples[i] = meth_samples[i].replace('\n','')

    #print str(imp_samples[i]) + '\t' + str(cov_samples[i]) + '\t' + str(meth_samples[i])
    #continue
    
    if (imp_samples[i] != cov_samples[i]) or (imp_samples[i] != meth_samples[i]):
        print "ALIGNMENT ERROR!"
    #ENDIF
    
#ENDFOR

print "DONE!"
        

               
