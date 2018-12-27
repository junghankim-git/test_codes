#!/bin/bash

function make_run_sh {

rm -rf run.sh
nproc=$1
maxcpusnode=$2
EXPDIR=$3
SYSTEM=$4

nnode=$(($nproc/$maxcpusnode))
extra=$(($nproc%$maxcpusnode))
if [ $extra -gt 0 ]; then
  nnode=$(($nnode+1))
fi
echo $nproc

if [ $SYSTEM = "gaon2" ]; then

echo "system: KIAPS"

cat <<EOF >./run.sh
#!/bin/bash
#PBS -l select=$nnode:ncpus=$maxcpusnode:mpiprocs=$maxcpusnode:ompthreads=1
#PBS -N log_$nproc
#PBS -q normal
#PBS -W block=true
#PBS -V
#PBS -j oe

cd $EXPDIR

mpirun -np $nproc -f \$PBS_NODEFILE ./run.exe | tee exp.txt
EOF

else

echo "system: KMA"

cat <<EOF >./run.sh
#!/bin/bash
#PBS -l select=$nnode:ncpus=$maxcpusnode:mpiprocs=$maxcpusnode:ompthreads=1
#PBS -N log_$nproc
#PBS -q normal@uri
#PBS -W block=true
#PBS -V
#PBS -j oe

cd $EXPDIR

aprun -n $nproc ./run.exe | tee exp.txt
EOF

fi

}

###############################################
# Main
###############################################

nprocs=(1 11 20 21 50 51 75 101)
#nprocs=(1 11)

expdir=$1
outdir=$2
system="gaon2"

maxcpus=-1
if [ $system = "gaon2" ]; then
maxcpus=16
fi
if [ $system = "nuri" ]; then
maxcpus=24
fi

rm main_log.txt
touch main_log.txt
echo "PROGRAM start" >> main_log.txt
for ip in ${nprocs[*]}
do
    
    #rm -rf $outdir/*.nc
    rm -rf /scratch/jhkim/data/multi_io_cp/*.nc
    make_run_sh $ip $maxcpus $expdir $system
    qsub run.sh
  
done
echo "PROGRAM  end " >> main_log.txt
