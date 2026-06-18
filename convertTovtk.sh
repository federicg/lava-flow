#!/bin/bash

# Check that NPROC was passed as argument
if [ $# -ne 2 ]; then
    echo "Usage: $0 <NPROC> <folder_results>"
    echo "Example: $0 4 results"
    exit 1
fi

NPROC=$1
CURRENT_WD=$2

ORIGINAL_WD=$(pwd)

cd "$CURRENT_WD" || { echo "Failed to cd into results"; exit 1; }

"${ORIGINAL_WD}/_deps/octave/install/bin/octave" --no-gui <<EOF
addpath(canonicalize_file_name('${ORIGINAL_WD}/../scripts/m/'))
export_tmesh_data('swe_h_%4.4d_%4.4d', ...
    {'swe_h_%4.4d_%4.4d', 'swe_Ux_%4.4d_%4.4d', 'swe_Uy_%4.4d_%4.4d', 'swe_Th_%4.4d_%4.4d', 'swe_Z_%4.4d_%4.4d'}, ...
    {'h', 'Ux', 'Uy', 'Th', 'Z'}, ...
    {}, {}, 'sweref', 0:1000, $NPROC)
EOF
