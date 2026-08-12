#!/bin/bash

mkdir -p pdbqt
mkdir -p logs
mkdir -p failed

# Number of CPU cores
NCPU=$(sysctl -n hw.ncpu)

echo "========================================"
echo "Preparing ligands"
echo "CPU cores: $NCPU"
echo "========================================"

# Reset failed list
> failed/failed_ligands.txt

find ligands -name "*.sdf" | \
parallel --bar -j "$NCPU" '

sdf="{}"
base=$(basename "$sdf" .sdf)

output="pdbqt/${base}.pdbqt"
log="logs/${base}.log"

# Skip already prepared ligands
if [ -f "$output" ]; then
    exit 0
fi

# Convert existing 3D SDF → PDBQT
# IMPORTANT: no --gen3d
obabel "$sdf" \
    -O "$output" \
    -h \
    -p 7.4 \
    > "$log" 2>&1

status=$?

if [ $status -ne 0 ] || [ ! -s "$output" ]; then

    echo "$base" >> failed/failed_ligands.txt

    rm -f "$output"

fi
'

echo
echo "========================================"
echo "Ligand preparation completed"
echo "========================================"

TOTAL=$(find ligands -name "*.sdf" | wc -l)
SUCCESS=$(find pdbqt -name "*.pdbqt" | wc -l)
FAILED=$(wc -l < failed/failed_ligands.txt 2>/dev/null || echo 0)

echo "Input SDF files : $TOTAL"
echo "PDBQT generated : $SUCCESS"
echo "Failed          : $FAILED"
echo "========================================"