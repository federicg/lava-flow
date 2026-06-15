#!/bin/bash

# Check that NPROC was passed as argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 <NPROC>"
    echo "Example: $0 4"
    exit 1
fi

NPROC=$1

cd results || { echo "Failed to cd into results"; exit 1; }

../_deps/octave/install/bin/octave --no-gui <<EOF
addpath(canonicalize_file_name('../../scripts/m/'))
export_tmesh_data('swe_h_%4.4d_%4.4d', ...
    {'swe_h_%4.4d_%4.4d', 'swe_Ux_%4.4d_%4.4d', 'swe_Uy_%4.4d_%4.4d', 'swe_Th_%4.4d_%4.4d', 'swe_Z_%4.4d_%4.4d'}, ...
    {'h', 'Ux', 'Uy', 'Th', 'Z'}, ...
    {}, {}, 'sweref', 0:1000, $NPROC)
EOF
