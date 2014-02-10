#!/bin/bash

input_mat_dir="./mat_in"
num_thread_ls=(1 2 4 8 16)

for entry in "$input_mat_dir"/*
do 
    for num_thread in ${num_thread_ls[*]}; do  
        command="./ge_pthread $entry -p $num_thread"
        for i in `seq 1 10`;do
            eval $command
        done
    done
done
