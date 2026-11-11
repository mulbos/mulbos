#!/bin/bash
#$ -S /bin/bash # Use the Bash shell
#$ -cwd # Use the current directory as working directory
#$ -N Mulbos # Name of the job
#$ -V # Export all variables to compute nodes
#$ -R y # Reserve slots
# from here are the programs/commands to run
# the output is written into mulbos.oNNN (output)
# and mulbos.eNNN (error if there is)
#echo "This is run XX"
./mulbos
