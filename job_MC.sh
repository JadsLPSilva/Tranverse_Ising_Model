#!/bin/bash
#PBS -N IT_16 			#Job Name
#PBS -j oe			#
#PBS -l nodes=1:ppn=1 		#Number of nodes and processor per nodes(ppn)
#PBS -l walltime=120:00:00	
#=========================================================================================
#
# Job.sh - Job Model for TORQUE manager
# 
# Jadson Lucas
# 
# Universidade Federal do Piauí(2019)
#
#=========================================================================================


####################################
######### Initialization  ##########
####################################

EXECFILE="isingt2d"                #Executable file
PARAMLIST="readisingt.in"              #Parameter file
#MODFILE="systemvariables.mod"          #variable file
##SRCFILE="isingt.F90"   #Source file
IFILELIST="$EXECFILE $PARAMLIST $MODFILE"	#Files required to run the code
OFILELIST="*.dat"                       #Output files

####################################
######### Bookeeping ###############
####################################
echo    "Executando no  : $HOSTNAME"
echo -n "Data              : "
date
echo    "Job ID            : $PBS_JOBID"
echo -n "Diretorio         : "
pwd

####################################
######### Preparing to run #########
####################################

# creating directory
#JRM ++++++++++++++++++++++++++
WRKDIR=$SCRATCH/$PBS_JOBID
mkdir -p $WRKDIR

## Transfer the necessary entries and files 
## to the /scratch directory
cd $PBS_O_WORKDIR
# JRM ++++++++++++++++++++++++

#mkdir $PBS_O_WORKDIR/$PBS_JOBID
#cd $PBS_O_WORKDIR/$PBS_JOBID
##mkdir /tmp/$PBS_JOBID
##cd /tmp/$PBS_JOBID

# Copying files to directory
echo "Copying the files to the directory provided"
for i in $IFILELIST
do
#        cp $PBS_O_WORKDIR/$i .
#JRM +++++++++++++++++++++++++++++
cp -r cp $PBS_O_WORKDIR/$i   $WRKDIR/

cd $WRKDIR
#JRM +++++++++++++++++++++++++++++

done

####################################
######### Running the program ######
####################################

#source /opt/intel/composerxe/bin/compilervars.sh intel64
#ifort  -o $EXECFILE $SRCFILE -O2

time ./$EXECFILE $PARAMLIST -p 8

####################################
######### Finishing ################
####################################

#copying output files to home directory
##for i in $OFILELIST
##do
##        cp $i $PBS_O_WORKDIR
##done

#JRm+++++++++

cp $WRKDIR/$OFILELIST $PBS_O_WORKDIR/


APAGA_SCRATCH=Y
## Delete the directory that the job ran
if [ "$APAGA_SCRATCH" = "Y" ]; then

    rm -rf $WRKDIR

else

    echo -e "\nO directory \e[00;31m$WRKDIR\e[00m must be removed manually
    to avoid problems for other jobs and/or users. \n"

fi

#JRM ++++++++++++++++++



echo -n "Job finished in : "
date
