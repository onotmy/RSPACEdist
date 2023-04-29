#!/bin/csh
#$ -cwd -pe mpi 16
setenv OMP_NUM_THREADS 2
mv $USER.$PE.$JOB_ID.1 hosts
set TEST = ""
foreach aaaa (`cat hosts`)
@ I++
set TEST = (${TEST} `echo $aaaa`)
end
rm hosts
set i=0
while ( $i < $NSLOTS )
@ i = $i + $OMP_NUM_THREADS
echo $TEST[$i] >> hosts
end
@ i = $NSLOTS / $OMP_NUM_THREADS

mpirun -np $i -machinefile hosts ./kukan43 > output.txt

