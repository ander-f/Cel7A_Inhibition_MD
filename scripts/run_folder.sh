#!/bin/bash

# Check if the correct number of arguments are passed
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <cosolvent> <concentration> <replica> <protein>"
    exit 1
fi

# Assign variables from command-line arguments
cosolvent=$1
concentration=$2
replica=$3
protein=$4

# Return to the scripts directory and run Julia
julia packmolinputcreator_${protein}.jl "$cosolvent" "$concentration" "$replica" "$protein"

# Copy required files to the directories
workdir=./..
\cp -f "cm_water_${protein}.jl"            "$workdir/$cosolvent/$concentration/$replica/$protein/"  # copy the .jl of the water to dir
\cp -f "cm_${cosolvent}_${protein}.jl"     "$workdir/$cosolvent/$concentration/$replica/$protein/"  # copy the .jl of the cosolvent to dir
