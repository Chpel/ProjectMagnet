#!/bin/bash
for ((i=41 ; i <90 ; i+=4))
do
  for len in 100 300 400 500 1000
  do
  sbatch --time=1-0:0 --wrap="srun ../a.out $len $i 0"
  done
done
