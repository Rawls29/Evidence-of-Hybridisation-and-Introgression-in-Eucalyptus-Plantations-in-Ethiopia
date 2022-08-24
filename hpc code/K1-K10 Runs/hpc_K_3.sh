#!/bin/bash
#PBS -l walltime=11:00:00
#PBS -l select=1:ncpus=1:mem=1gb
module load anaconda3/personal
echo "Hope this works!"
cp $HOME/pops/mainparams .
cp $HOME/pops/extraparams .
cp $HOME/pops/Eucalyptus_STRUCTURE_input_MAF_0.05.txt .
structure -K 3 -o MAF_K_3_1
structure -K 3 -o MAF_K_3_2
structure -K 3 -o MAF_K_3_3
mv MAF_K_3_1_f $HOME/pops
mv MAF_K_3_2_f $HOME/pops
mv MAF_K_3_3_f $HOME/pops
echo "Successful?"