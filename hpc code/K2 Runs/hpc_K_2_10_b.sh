#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=1:mem=1gb
module load anaconda3/personal
echo "Hope this works!"
cp $HOME/K_2/mainparams .
cp $HOME/K_2/extraparams .
cp $HOME/K_2/Eucalyptus_STRUCTURE_input_MAF_0.05.txt .
structure -K 2 -o MAF_K_2_1_b
structure -K 2 -o MAF_K_2_2_b
structure -K 2 -o MAF_K_2_3_b
structure -K 2 -o MAF_K_2_4_b
structure -K 2 -o MAF_K_2_5_b
structure -K 2 -o MAF_K_2_6_b
structure -K 2 -o MAF_K_2_7_b
structure -K 2 -o MAF_K_2_8_b
structure -K 2 -o MAF_K_2_9_b
structure -K 2 -o MAF_K_2_10_b
mv MAF_K_2_1_b_f $HOME/K_2
mv MAF_K_2_2_b_f $HOME/K_2
mv MAF_K_2_3_b_f $HOME/K_2
mv MAF_K_2_4_b_f $HOME/K_2
mv MAF_K_2_5_b_f $HOME/K_2
mv MAF_K_2_6_b_f $HOME/K_2
mv MAF_K_2_7_b_f $HOME/K_2
mv MAF_K_2_8_b_f $HOME/K_2
mv MAF_K_2_9_b_f $HOME/K_2
mv MAF_K_2_10_b_f $HOME/K_2
echo "Successful?"