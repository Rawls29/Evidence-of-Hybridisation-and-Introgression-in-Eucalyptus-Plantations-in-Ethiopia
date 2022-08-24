#!/bin/bash
#PBS -l walltime=11:00:00
#PBS -l select=1:ncpus=1:mem=1gb
module load anaconda3/personal
echo "Hope this works!"
cp $HOME/pops/mainparams .
cp $HOME/pops/extraparams .
cp $HOME/pops/Eucalyptus_STRUCTURE_input_MAF_0.05.txt .
structure -K 2 -o MAF_K_2_1
structure -K 2 -o MAF_K_2_2
structure -K 2 -o MAF_K_2_3
mv MAF_K_2_1_f $HOME/pops
mv MAF_K_2_2_f $HOME/pops
mv MAF_K_2_3_f $HOME/pops
echo "Successful?"