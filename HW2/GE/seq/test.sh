#!/bin/bash

input_mat_dir="./mat_in"

for entry in "$input_mat_dir"/*
do 
    command="./gauss $entry"
    echo $entry 
    for i in `seq 1 10`;do
        eval $command
    done
done
