#!/bin/bash

# This script runs a Gaussian job.
#
# Usage: run_gaussian.sh /scratch/dae23_gaussian_0000020002
#
# This script assumes the job will be called gaussian.gjf and the output should go into gaussian.out.
# You can modify this script to suit your system's requirements.

# read command line arguments
jobDirectory=$1

# ensure directory is empty
cd $jobDirectory
#export GAUSS_SCRDIR=$jobDirectory

# run Gaussian
g16 gaussian.gjf gaussian.out 

# the job directory will be cleaned up by Jprogdyn, so we don't need to do anything for that
