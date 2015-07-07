#!/bin/bash

# settings from input

wtime=${1:-72:00:00}

# settings for files

binary=tm1dNNN.GF

inpfile=tm1dNNN.inp
orgfile=$inpfile"-org"

# settings for directories

currdir=`pwd`

binarydir=$HOME/Projects/tm1DNNN/EXE

submitdir=$currdir
[ -d ${submitdir} ] || mkdir ${submitdir}

jobdir=${submitdir}

[ -d ${jobdir} ] || mkdir ${jobdir}

# compute jobs for gate sweep
# this will use individual jobs along each x-line

for imodel in 0 1 2 3
do

for energy in 0.0 1.0
do

for irange in 1 2 3 4 5 6 7 8 9 10 20 50 100 200 500 1000
do

jobname="NN_M$imodel_E$energy_R$irange"
echo $jobname

jobfile=`printf "$jobname.sh"`

targetdir=${submitdir}/../${jobname}
[ -d ${targetdir} ] || mkdir ${targetdir}

tmpdir=/storage/disqs/phsht/RUNS/${jobname}

echo "binarydir=" $binarydir " submitdir=" $submitdir 
echo "tmpdir=" $tmpdir " targetdir=" $targetdir

# settings for parallel submission

cat > ${jobdir}/${jobfile} << EOD
#!/bin/bash --login
#PBS -l pvmem=400mb
##PBS -M r.roemer@warwick.ac.uk
#PBS -m a
#PBS -r y
#PBS -V
##PBS -k oe
#PBS -j oe

#       The jobname
#PBS -N ${jobname}

#       The total number of parallel tasks for your job.
#       This is the sum of the number of parallel tasks required by each
#       of the aprun commands you are using. In this example we have
#       ${mppwidth} tasks
##PBS -l mppwidth=${mppwidth}
#PBS -l nodes=1:ppn=1

#       Specify how many processes per node.
##PBS -l mppnppn=32

#       Specify the wall clock time required for your job.
#PBS -l walltime=${wtime}

#       Specify which budget account that your job will be charged to.
##PBS -A midpluswarwick
  
# Make sure any symbolic links are resolved to absolute path
export PBS_O_WORKDIR=\$(readlink -f \$PBS_O_WORKDIR)

# The base directory is the directory that the job was submitted from.
basedir=\$PBS_O_WORKDIR
echo "basedir=" \${basedir}

# The binary directory is the directory where the scripts are
echo "binarydir=" ${binarydir}

# do the run in a TMPdir for safekeeping
#tmpdir=$tmpbasedir/`basename ${jobname} .sh`

[ -d ${tmpdir} ] || mkdir ${tmpdir}

cp ${binarydir}/${binary} ${tmpdir}/

# construct the input file

cd ${tmpdir}
echo "hostname=" $HOSTNAME
pwd

rm -rf $inpfile

touch $inpfile

echo "ISeed        =            1277">> $inpfile   
echo "NOfIter      =       100000000">> $inpfile   
echo "NOfOrtho     =               2">> $inpfile   
echo "NOfPrint     =          100000">> $inpfile        
echo "NOfGamma     =               1">> $inpfile   
echo "IModelFlag   =               $imodel">> $inpfile   
echo "IBCFlag      =               0">> $inpfile   
echo "IRNGFlag     =               0">> $inpfile   
echo "IKeepFlag    =               0">> $inpfile   
echo "IWriteFlag   =               1">> $inpfile   
echo "ISortFlag    =               0">> $inpfile   
echo "IFluxFlag    =               0">> $inpfile   
echo "Range0       =         $irange">> $inpfile   
echo "Range1       =         $irange">> $inpfile   
echo "dRange       =               1">> $inpfile   
echo "DiagDis0     =        0.250000">> $inpfile   
echo "DiagDis1     =        10.10000">> $inpfile   
echo "dDiagDis     =        0.250000">> $inpfile   
echo "Energy0      =        $energy">> $inpfile   
echo "Energy1      =        $energy">> $inpfile   
echo "dEnergy      =        0.500000">> $inpfile   
echo "IKappaFlag   =               1">> $inpfile   
echo "Kappa        =        10.00000">> $inpfile   
echo "Epsilon      =           0.001">> $inpfile   
echo "IMidDump     =               0">> $inpfile   
echo "IWalltime    =              14">> $inpfile   
echo "IRestart     =               0">> $inpfile   

cat $inpfile
ls -al \${tmpdir}

# Make sure ranks are numbered appropriately
export MPICH_RANK_REORDER_METHOD=1 
export MPICH_PTL_MATCH_OFF=1
export MPICH_FAST_MEMCPY=1
export IPATH_NO_CPUAFFINITY=1
  
echo "--- starting the SIM run"
time ./${binary}

# copy the result of the run in the detination for safekeeping
#targetdir=${submitdir}/`basename ${job_file} .sh`

[ -d ${targetdir} ] || mkdir ${targetdir}

cp -vr *.raw ${targetdir}
cp -vr *.inp ${targetdir}

wait
#exit 0
EOD

chmod 755 ${jobdir}/${jobfile}
(cd ${jobdir} ; qsub -q serial ./${jobfile})
#(cd ${jobdir} ; qsub -q parallel ./${jobfile})
#(cd ${jobdir} ; qsub -q devel ./${jobfile})
#(cd ${jobdir} ; qsub -q taskfarm ${jobfile})
#(cd ${jobdir} ; qsub ./${jobfile})

done
done
done

