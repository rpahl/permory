#!/bin/bash
#
# Copyright (c) 2010 Roman Pahl
#               2011 Volker Steiß
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)
#
# Example MPI submission script for use with Marbuger Rechen Cluster (MaRC).
#
# Example submition to parallel queue starting 10 processes using MPI:
# $ qsub -cwd -l h_rt=1:0:0 -pe orte 10 -q parallel@@nodes_ng -R y mpi_marc.sh
#

#############################    Configuration    #############################
##
## Edit the configuration to your needs.
##

# random seed
SEED=123456789

# permutations
PERMS=50000

## Path to PERMORY executable. The binary will be copied to $EXEC_DIR
PERMORY_BIN=permory

## Data source path which contains the data files to use.
DATA_SRC=/home/${USER}/data/

## List of data files.
FILES="chr1.slide.gz chr2.slide.gz chr3.slide.gz chr4.slide.gz chr5.slide.gz"

## Temporary directory for job execution. This place needs to be accassable to
## all nodes. This directory will be created at start, filled with the data
## files and deleted after execution.
EXEC_DIR=/glusterfs/distrib/${USER}/${JOB_ID}

## processes to start
PROCESSES=${NSLOTS}
## Workarround for missing -npernode option in mpirun before version 1.3.
## Uncomment this line to start two processes per node.
# PROCESSES=$((${NSLOTS}*2))

## Prefix for output filenames.
OUT_PREFIX=mpi_${PERMS}_${PROCESSES}p

## Memorize current working directory. The results will be copied here.
OLD_DIR=`pwd`


## build option string to pass to PERMORY
OPTIONS="--seed ${SEED} --nperm ${PERMS} -o ${OUT_PREFIX}"

#############################      Execution      #############################

## Create temporary dir
mkdir -p "${EXEC_DIR}"
cd "${EXEC_DIR}"

## Copy binary
cp "${PERMORY_BIN}" ./permory

## copy data
for FILE in ${FILES}
do
    cp "${DATA_SRC}"/${FILE} ./
done

## run
mpirun -np ${PROCESSES} -machinefile ${TMPDIR}/machines ./permory ${OPTIONS} ${FILES}

cp -f "${OUT_PREFIX}".* "${OLD_DIR}"

# Clean up temporary dir
rm -rf "${EXEC_DIR}"
