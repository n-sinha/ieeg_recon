#!/bin/bash

# make a loop that runs though each directory in /project/davis_group_1/nishants/ieeg_recon/data/BIDS
for dir in /project/davis_group_1/nishants/ieeg_recon/data/BIDS/*; do

    # get the rid from the directory name
    rid=$(basename $dir)
    bsub -q bsc_normal "sh /project/davis_group_1/nishants/ieeg_recon/jobs/job.sh $rid"

done