#!/bin/bash

bsub -W 100 -R rhel60 -o test_out.txt python skim_one.py $1 $2
