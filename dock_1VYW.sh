#!/bin/bash

set -u

# ============================================================
# SETTINGS
# ============================================================

RECEPTOR="1VYW.pdbqt"
CONFIG="vina_1VYW.txt"

LIGAND_DIR="pdbqt"

OUTPUT="vina_1VYW_summary.csv"
RANKED="vina_1VYW_ranked.csv"
JOBLOG="vina_1VYW_parallel.log"

# Use all CPU cores except one
NCPU=$(sysctl -n hw.ncpu)
JOBS=$((NCPU - 1))

if [ "$JOBS" -lt 1 ]; then
    JOBS=1
fi

# ============================================================
# CHECK FILES
# ============================================================

if [ ! -f "$RECEPTOR" ]; then
    echo "ERROR: Receptor not found: $RECEPTOR"
    exit 1
fi

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Configuration file not found: $CONFIG"
    exit 1
fi

if [ ! -d "$LIGAND_DIR" ]; then
    echo "ERROR: Ligand directory not found: $LIGAND_DIR"
    exit 1
fi

# ============================================================
# START
# ============================================================

echo "=============================================="
echo "AutoDock Vina ultra-fast screening"
echo "=============================================="
echo "Receptor      : $RECEPTOR"
echo "Config        : $CONFIG"
echo "Ligand folder : $LIGAND_DIR"
echo "CPU cores     : $NCPU"
echo "Parallel jobs : $JOBS"
echo "Vina CPU/job  : 1"
echo "=============================================="

date
echo

TOTAL=$(find "$LIGAND_DIR" -name "*.pdbqt" | wc -l | tr -d ' ')

echo "Ligands found: $TOTAL"
echo

if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: No PDBQT ligands found."
    exit 1
fi

# ============================================================
# REMOVE OLD RESULTS
# ============================================================

rm -f "$OUTPUT"
rm -f "$RANKED"
rm -f "$JOBLOG"

# ============================================================
# RUN DOCKING
#
# IMPORTANT:
# The entire GNU Parallel output is redirected to the CSV.
# No individual CSV files are generated.
# ============================================================

{
    echo "Ligand,Affinity"

    find "$LIGAND_DIR" -name "*.pdbqt" -print0 | \
    parallel \
        --null \
        --bar \
        --joblog "$JOBLOG" \
        -j "$JOBS" '

        ligand="{}"
        name=$(basename "$ligand" .pdbqt)

        affinity=$(
            vina \
                --receptor "'"$RECEPTOR"'" \
                --ligand "$ligand" \
                --config "'"$CONFIG"'" \
                --cpu 1 \
                2>/dev/null |
            awk '"'"'
                /^[[:space:]]*1[[:space:]]+/ {
                    print $2
                    exit
                }
            '"'"'
        )

        if [ -z "$affinity" ]; then
            affinity="FAILED"
        fi

        printf "%s,%s\n" "$name" "$affinity"

    '

} > "$OUTPUT"

# ============================================================
# SORT RESULTS
# ============================================================

echo
echo "Sorting results..."

{
    head -n 1 "$OUTPUT"

    tail -n +2 "$OUTPUT" |
    awk -F',' '
        $2 ~ /^-?[0-9]+([.][0-9]+)?$/ {
            print
        }
    ' |
    sort -t',' -k2,2n

} > "$RANKED"

# ============================================================
# STATISTICS
# ============================================================

SUCCESS=$(tail -n +2 "$OUTPUT" |
    awk -F',' '
        $2 ~ /^-?[0-9]+([.][0-9]+)?$/ {
            count++
        }
        END {
            print count+0
        }
    ')

FAILED=$(tail -n +2 "$OUTPUT" |
    awk -F',' '
        $2 == "FAILED" {
            count++
        }
        END {
            print count+0
        }
    ')

# ============================================================
# FINISH
# ============================================================

echo
echo "=============================================="
echo "SCREENING COMPLETED"
echo "=============================================="
echo "Total ligands : $TOTAL"
echo "Successful    : $SUCCESS"
echo "Failed        : $FAILED"
echo
echo "Raw results   : $OUTPUT"
echo "Ranked results: $RANKED"
echo "Parallel log  : $JOBLOG"
echo "=============================================="

echo
echo "Top 10 ligands:"
head -n 11 "$RANKED"

echo
date