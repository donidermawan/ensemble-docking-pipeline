#!/bin/bash

RECEPTORS=("1VYW" "7RWF" "5A14" "5IF1" "1W98")

for pdb in "${RECEPTORS[@]}"
do
    echo "========================================"
    echo "Preparing receptor: ${pdb}"
    echo "========================================"

    mk_prepare_receptor.py \
        --read_pdb "${pdb}.pdb" \
        --allow_bad_res \
        -o "${pdb}" \
        -p \
        > "${pdb}_prepare.log" 2>&1

    # Check output
    if [ -f "${pdb}.pdbqt" ]; then
        echo "✓ Successfully prepared ${pdb}"
    else
        echo "✗ Failed to prepare ${pdb}"
        echo "Check log: ${pdb}_prepare.log"
    fi

    echo ""
done

echo "========================================"
echo "Receptor preparation completed."
echo "========================================"