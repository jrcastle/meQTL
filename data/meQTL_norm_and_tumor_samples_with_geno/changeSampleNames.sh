#!/bin/sh
echo 'sed ... meth_all.txt > tmp.txt'
sed '1 s/X//g' ../meth_all.txt > ./tmp.txt
echo 'sed ... tmp.txt > meth_all_fixed.txt'
sed '1 s/\./-/g' ./tmp.txt > ./meth_all_fixed.txt
echo 'rm tmp.txt'
rm tmp.txt
echo 'python -u remove_adjnormal.py'
python -u remove_adjnormal.py
echo "DONE!"
