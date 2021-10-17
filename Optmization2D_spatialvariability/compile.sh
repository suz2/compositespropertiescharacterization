#!/bin/bash

# To compile umat with Calculix, make sure all necessary .f90 files are 
# included within the /source folder and all input files are included
# in the parent folder


# Step1: cd to the calculix source folder and remove the /source folder
#        which contains the old source codes for umat. Remove the file
#        UMATSRC.inc as well to remove the name of the old source files

#		 The calculix source folder could vary based on installation,
#		 the user should make necessary adjustments.
SOLVER_NAME='CalculiX'
SOLVER_DIR1='/home/zsu/calculix/'
SOLVER_DIR2='/ccx_2.17/src'
SOLVER_DIR=${SOLVER_DIR1}/${SOLVER_NAME}${SOLVER_DIR2}
##cd /home/zimu/MUMSresearch/calculix/${SOLVER_NAME}/ccx_2.17/src
cd ${SOLVER_DIR}
rm -r source
rm UMATSRC.inc


# Step2: cd back to parent folder, create UMATSRC.inc, which contains 
#        all source file names
cd -
> UMATSRC.inc # clear file content
echo "UMATSRC = \\" >> UMATSRC.inc
for name in ./source/*
do
	filename=`basename $name`
	echo "$filename \\" >> UMATSRC.inc
done


# Step3: cp the source folder and UMATSRC.inc to calculix source folder
#cp -r source /home/zimu/MUMSresearch/calculix/${SOLVER_NAME}/ccx_2.17/src
#cp UMATSRC.inc /home/zimu/MUMSresearch/calculix/${SOLVER_NAME}/ccx_2.17/src
cp -r source ${SOLVER_DIR}
cp UMATSRC.inc ${SOLVER_DIR}

# Step4: run make.sh, which excutes the Makefile and makes necessary
#        changes
#cd /home/zimu/MUMSresearch/calculix/${SOLVER_NAME}/ccx_2.17/src
cd ${SOLVER_DIR}
./make.sh
#make

# Now the Calculix is compiled with our own umat
echo 'compiling for the Calculix with umat is finished'

