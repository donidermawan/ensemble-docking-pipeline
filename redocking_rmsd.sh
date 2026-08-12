#!/bin/bash

# ============================================================
# CDK2 REFERENCE-INHIBITOR REDOCKING + RMSD VALIDATION
# + RECEPTOR-LIGAND COMPLEX GENERATION
# + MEAN RMSD CALCULATION
# ============================================================
#
# Systems:
#
#   1. TW8672      – 7RWF
#   2. PNU-292137  – 1VYW
#   3. K03861      – 5A14
#
# Docking:
#   AutoDock Vina 1.2.7
#
# Validation:
#   Heavy-atom RMSD against experimental crystal ligand pose
#   using Open Babel obrms
#
# RMSD:
#   All generated Vina poses are evaluated.
#   Mean RMSD is calculated across all reported poses.
#
# Additional output:
#   Receptor + best Vina ligand pose (MODEL 1)
#   combined into a complex PDB.
#
# ============================================================

set -u


# ============================================================
# SOFTWARE
# ============================================================

VINA="vina"
OBABEL="obabel"
OBRMS="obrms"


# ============================================================
# CHECK SOFTWARE
# ============================================================

echo "Checking software..."

if ! command -v "$VINA" >/dev/null 2>&1; then
    echo "ERROR: Vina not found."
    exit 1
fi

if ! command -v "$OBABEL" >/dev/null 2>&1; then
    echo "ERROR: Open Babel not found."
    exit 1
fi

if ! command -v "$OBRMS" >/dev/null 2>&1; then
    echo "ERROR: obrms not found."
    exit 1
fi

echo "Vina:"
"$VINA" --version

echo
echo "Open Babel:"
"$OBABEL" -V

echo
echo "obrms:"
"$OBRMS" --version 2>/dev/null || true

echo


# ============================================================
# OUTPUT DIRECTORIES
# ============================================================

mkdir -p results

mkdir -p results/7RWF_TW8672
mkdir -p results/1VYW_PNU292137
mkdir -p results/5A14_K03861


# ============================================================
# FUNCTION: CHECK INPUT
# ============================================================

check_file() {

    FILE="$1"

    if [ ! -f "$FILE" ]; then

        echo
        echo "ERROR: Required file not found:"
        echo "  $FILE"

        exit 1

    fi

}


# ============================================================
# FUNCTION: CALCULATE MEAN RMSD
# ============================================================
#
# Input:
#   RMSD output file from obrms
#
# Example:
#
# RMSD crystal.pdb:ligand.pdb 0.46079
# RMSD crystal.pdb:ligand.pdb 1.78574
# RMSD crystal.pdb:ligand.pdb 0.621234
#
# Output:
#   Mean RMSD in Angstrom
#
# ============================================================

calculate_mean_rmsd() {

    RMSD_FILE="$1"

    awk '
    /^RMSD/ {
        sum += $NF
        n++
    }
    END {
        if (n > 0)
            printf "%.6f", sum/n
        else
            print "NA"
    }
    ' "$RMSD_FILE"

}


# ============================================================
# FUNCTION:
# CREATE RECEPTOR-LIGAND COMPLEX
# ============================================================

create_complex() {

    RECEPTOR_PDBQT="$1"
    VINA_OUT="$2"
    OUTPUT_DIR="$3"
    PREFIX="$4"


    echo
    echo "Creating receptor-ligand complex..."
    echo


    # --------------------------------------------------------
    # Convert receptor PDBQT -> PDB
    # --------------------------------------------------------

    RECEPTOR_PDB="${OUTPUT_DIR}/${PREFIX}_receptor.pdb"

    echo "Converting receptor..."
    echo "  $RECEPTOR_PDBQT"
    echo "->"
    echo "  $RECEPTOR_PDB"


    "$OBABEL" \
        "$RECEPTOR_PDBQT" \
        -O "$RECEPTOR_PDB" \
        >/dev/null 2>&1


    if [ $? -ne 0 ] || [ ! -s "$RECEPTOR_PDB" ]; then

        echo "ERROR: Failed to convert receptor to PDB."
        exit 1

    fi


    # --------------------------------------------------------
    # Extract MODEL 1 from Vina output
    # --------------------------------------------------------

    BEST_POSE_PDBQT="${OUTPUT_DIR}/${PREFIX}_best_pose.pdbqt"


    echo
    echo "Extracting best Vina pose (MODEL 1)..."


    awk '
        /^MODEL[[:space:]]+1/ {
            capture=1
        }

        capture {
            print
        }

        /^ENDMDL/ && capture {
            exit
        }
    ' "$VINA_OUT" > "$BEST_POSE_PDBQT"


    if [ ! -s "$BEST_POSE_PDBQT" ]; then

        echo "ERROR: Could not extract MODEL 1 from Vina output."
        exit 1

    fi


    # --------------------------------------------------------
    # Convert best ligand pose PDBQT -> PDB
    # --------------------------------------------------------

    BEST_POSE_PDB="${OUTPUT_DIR}/${PREFIX}_best_pose.pdb"


    echo
    echo "Converting best ligand pose to PDB..."


    "$OBABEL" \
        "$BEST_POSE_PDBQT" \
        -O "$BEST_POSE_PDB" \
        >/dev/null 2>&1


    if [ $? -ne 0 ] || [ ! -s "$BEST_POSE_PDB" ]; then

        echo "ERROR: Failed to convert ligand pose to PDB."
        exit 1

    fi


    # --------------------------------------------------------
    # Create complex
    # --------------------------------------------------------

    COMPLEX_PDB="${OUTPUT_DIR}/${PREFIX}_complex.pdb"


    echo
    echo "Combining receptor + ligand..."


    {
        cat "$RECEPTOR_PDB"

        echo "TER"

        awk '
            /^ATOM/ ||
            /^HETATM/ ||
            /^CONECT/ {
                print
            }
        ' "$BEST_POSE_PDB"

        echo "END"

    } > "$COMPLEX_PDB"


    if [ ! -s "$COMPLEX_PDB" ]; then

        echo "ERROR: Complex PDB was not created."
        exit 1

    fi


    echo
    echo "Complex created:"
    echo "  $COMPLEX_PDB"


    rm -f "$BEST_POSE_PDBQT"

}


# ============================================================
# VARIABLES FOR SUMMARY
# ============================================================

TW_AFFINITY="NA"
PNU_AFFINITY="NA"
K_AFFINITY="NA"

TW_MEAN_RMSD="NA"
PNU_MEAN_RMSD="NA"
K_MEAN_RMSD="NA"


# ============================================================
# 1. TW8672 – 7RWF
# ============================================================

echo
echo "============================================================"
echo "1/3  TW8672 – 7RWF"
echo "============================================================"
echo


check_file "7RWF/7RWF.pdbqt"
check_file "7RWF/TW8672.pdbqt"
check_file "7RWF/TW8672_crystal.pdb"


echo "Running AutoDock Vina..."
echo


"$VINA" \
    --receptor "7RWF/7RWF.pdbqt" \
    --ligand "7RWF/TW8672.pdbqt" \
    --center_x -22.195 \
    --center_y 2.569 \
    --center_z 17.819 \
    --size_x 24 \
    --size_y 24 \
    --size_z 24 \
    --exhaustiveness 8 \
    --num_modes 3 \
    --energy_range 3 \
    --out "results/7RWF_TW8672/TW8672_7RWF_out.pdbqt" \
    > "results/7RWF_TW8672/TW8672_7RWF.log" 2>&1


if [ $? -ne 0 ]; then

    echo "ERROR: Vina docking failed for TW8672–7RWF."
    cat "results/7RWF_TW8672/TW8672_7RWF.log"
    exit 1

fi


echo "Docking completed."


# ------------------------------------------------------------
# Extract best affinity
# ------------------------------------------------------------

TW_AFFINITY=$(awk '
/^[[:space:]]*1[[:space:]]+/ {
    print $2
    exit
}
' "results/7RWF_TW8672/TW8672_7RWF.log")


if [ -z "$TW_AFFINITY" ]; then
    TW_AFFINITY="NA"
fi


echo "Best Vina affinity: ${TW_AFFINITY} kcal/mol"


# ------------------------------------------------------------
# Convert docked poses to SDF
# ------------------------------------------------------------

"$OBABEL" \
    "results/7RWF_TW8672/TW8672_7RWF_out.pdbqt" \
    -O "results/7RWF_TW8672/TW8672_7RWF_out.sdf" \
    >/dev/null 2>&1


# ------------------------------------------------------------
# Convert crystal ligand to SDF
# ------------------------------------------------------------

"$OBABEL" \
    "7RWF/TW8672_crystal.pdb" \
    -O "results/7RWF_TW8672/TW8672_crystal.sdf" \
    >/dev/null 2>&1


# ------------------------------------------------------------
# RMSD
# ------------------------------------------------------------

echo
echo "Calculating TW8672–7RWF RMSD..."


"$OBRMS" \
    -f \
    -m \
    "results/7RWF_TW8672/TW8672_crystal.sdf" \
    "results/7RWF_TW8672/TW8672_7RWF_out.sdf" \
    > "results/7RWF_TW8672/TW8672_RMSD.txt" 2>&1


echo
echo "TW8672–7RWF individual RMSDs:"
cat "results/7RWF_TW8672/TW8672_RMSD.txt"


# ------------------------------------------------------------
# Calculate mean RMSD
# ------------------------------------------------------------

TW_MEAN_RMSD=$(calculate_mean_rmsd \
    "results/7RWF_TW8672/TW8672_RMSD.txt")


echo
echo "TW8672–7RWF mean RMSD:"
echo "  ${TW_MEAN_RMSD} Å"


# ------------------------------------------------------------
# CREATE COMPLEX
# ------------------------------------------------------------

create_complex \
    "7RWF/7RWF.pdbqt" \
    "results/7RWF_TW8672/TW8672_7RWF_out.pdbqt" \
    "results/7RWF_TW8672" \
    "TW8672_7RWF"


# ============================================================
# 2. PNU-292137 – 1VYW
# ============================================================

echo
echo "============================================================"
echo "2/3  PNU-292137 – 1VYW"
echo "============================================================"
echo


check_file "1VYW/1VYW.pdbqt"
check_file "1VYW/PNU-292137.pdbqt"
check_file "1VYW/PNU-292137_crystal.pdb"


echo "Running AutoDock Vina..."
echo


"$VINA" \
    --receptor "1VYW/1VYW.pdbqt" \
    --ligand "1VYW/PNU-292137.pdbqt" \
    --center_x -8.929 \
    --center_y 205.749 \
    --center_z 111.606 \
    --size_x 24 \
    --size_y 24 \
    --size_z 24 \
    --exhaustiveness 8 \
    --num_modes 3 \
    --energy_range 3 \
    --out "results/1VYW_PNU292137/PNU292137_1VYW_out.pdbqt" \
    > "results/1VYW_PNU292137/PNU292137_1VYW.log" 2>&1


if [ $? -ne 0 ]; then

    echo "ERROR: Vina docking failed for PNU-292137–1VYW."
    cat "results/1VYW_PNU292137/PNU292137_1VYW.log"
    exit 1

fi


echo "Docking completed."


# ------------------------------------------------------------
# Extract affinity
# ------------------------------------------------------------

PNU_AFFINITY=$(awk '
/^[[:space:]]*1[[:space:]]+/ {
    print $2
    exit
}
' "results/1VYW_PNU292137/PNU292137_1VYW.log")


if [ -z "$PNU_AFFINITY" ]; then
    PNU_AFFINITY="NA"
fi


echo "Best Vina affinity: ${PNU_AFFINITY} kcal/mol"


# ------------------------------------------------------------
# Convert docked poses
# ------------------------------------------------------------

"$OBABEL" \
    "results/1VYW_PNU292137/PNU292137_1VYW_out.pdbqt" \
    -O "results/1VYW_PNU292137/PNU292137_1VYW_out.sdf" \
    >/dev/null 2>&1


# ------------------------------------------------------------
# Convert crystal ligand
# ------------------------------------------------------------

"$OBABEL" \
    "1VYW/PNU-292137_crystal.pdb" \
    -O "results/1VYW_PNU292137/PNU292137_crystal.sdf" \
    >/dev/null 2>&1


# ------------------------------------------------------------
# RMSD
# ------------------------------------------------------------

echo
echo "Calculating PNU-292137–1VYW RMSD..."


"$OBRMS" \
    -f \
    -m \
    "results/1VYW_PNU292137/PNU292137_crystal.sdf" \
    "results/1VYW_PNU292137/PNU292137_1VYW_out.sdf" \
    > "results/1VYW_PNU292137/PNU292137_RMSD.txt" 2>&1


echo
echo "PNU-292137–1VYW individual RMSDs:"
cat "results/1VYW_PNU292137/PNU292137_RMSD.txt"


# ------------------------------------------------------------
# Mean RMSD
# ------------------------------------------------------------

PNU_MEAN_RMSD=$(calculate_mean_rmsd \
    "results/1VYW_PNU292137/PNU292137_RMSD.txt")


echo
echo "PNU-292137–1VYW mean RMSD:"
echo "  ${PNU_MEAN_RMSD} Å"


# ------------------------------------------------------------
# CREATE COMPLEX
# ------------------------------------------------------------

create_complex \
    "1VYW/1VYW.pdbqt" \
    "results/1VYW_PNU292137/PNU292137_1VYW_out.pdbqt" \
    "results/1VYW_PNU292137" \
    "PNU292137_1VYW"


# ============================================================
# 3. K03861 – 5A14
# ============================================================

echo
echo "============================================================"
echo "3/3  K03861 – 5A14"
echo "============================================================"
echo


check_file "5A14/5A14.pdbqt"
check_file "5A14/K03861.pdbqt"
check_file "5A14/K03861_crystal.pdb"


echo "Running AutoDock Vina..."
echo


"$VINA" \
    --receptor "5A14/5A14.pdbqt" \
    --ligand "5A14/K03861.pdbqt" \
    --center_x 22.287 \
    --center_y 8.979 \
    --center_z 10.665 \
    --size_x 24 \
    --size_y 24 \
    --size_z 24 \
    --exhaustiveness 8 \
    --num_modes 3 \
    --energy_range 3 \
    --out "results/5A14_K03861/K03861_5A14_out.pdbqt" \
    > "results/5A14_K03861/K03861_5A14.log" 2>&1


if [ $? -ne 0 ]; then

    echo "ERROR: Vina docking failed for K03861–5A14."
    cat "results/5A14_K03861/K03861_5A14.log"
    exit 1

fi


echo "Docking completed."


# ------------------------------------------------------------
# Extract affinity
# ------------------------------------------------------------

K_AFFINITY=$(awk '
/^[[:space:]]*1[[:space:]]+/ {
    print $2
    exit
}
' "results/5A14_K03861/K03861_5A14.log")


if [ -z "$K_AFFINITY" ]; then
    K_AFFINITY="NA"
fi


echo "Best Vina affinity: ${K_AFFINITY} kcal/mol"


# ------------------------------------------------------------
# Convert docked poses
# ------------------------------------------------------------

"$OBABEL" \
    "results/5A14_K03861/K03861_5A14_out.pdbqt" \
    -O "results/5A14_K03861/K03861_5A14_out.sdf" \
    >/dev/null 2>&1


# ------------------------------------------------------------
# Convert crystal ligand
# ------------------------------------------------------------

"$OBABEL" \
    "5A14/K03861_crystal.pdb" \
    -O "results/5A14_K03861/K03861_crystal.sdf" \
    >/dev/null 2>&1


# ------------------------------------------------------------
# RMSD
# ------------------------------------------------------------

echo
echo "Calculating K03861–5A14 RMSD..."


"$OBRMS" \
    -f \
    -m \
    "results/5A14_K03861/K03861_crystal.sdf" \
    "results/5A14_K03861/K03861_5A14_out.sdf" \
    > "results/5A14_K03861/K03861_RMSD.txt" 2>&1


echo
echo "K03861–5A14 individual RMSDs:"
cat "results/5A14_K03861/K03861_RMSD.txt"


# ------------------------------------------------------------
# Mean RMSD
# ------------------------------------------------------------

K_MEAN_RMSD=$(calculate_mean_rmsd \
    "results/5A14_K03861/K03861_RMSD.txt")


echo
echo "K03861–5A14 mean RMSD:"
echo "  ${K_MEAN_RMSD} Å"


# ------------------------------------------------------------
# CREATE COMPLEX
# ------------------------------------------------------------

create_complex \
    "5A14/5A14.pdbqt" \
    "results/5A14_K03861/K03861_5A14_out.pdbqt" \
    "results/5A14_K03861" \
    "K03861_5A14"


# ============================================================
# CREATE SUMMARY CSV
# ============================================================

SUMMARY="results/redocking_summary.csv"


echo
echo "Creating summary table..."


cat > "$SUMMARY" << EOF
System,Receptor,Ligand,Vina_Affinity_kcal_mol,Mean_RMSD_A
TW8672-7RWF,7RWF,TW8672,$TW_AFFINITY,$TW_MEAN_RMSD
PNU-292137-1VYW,1VYW,PNU-292137,$PNU_AFFINITY,$PNU_MEAN_RMSD
K03861-5A14,5A14,K03861,$K_AFFINITY,$K_MEAN_RMSD
EOF


# ============================================================
# FINAL
# ============================================================

echo
echo "============================================================"
echo "REDOCKING + RMSD + COMPLEX GENERATION COMPLETED"
echo "============================================================"


echo
echo "Summary:"
echo "  $SUMMARY"


echo
echo "Mean RMSD:"
echo "  TW8672–7RWF      : $TW_MEAN_RMSD Å"
echo "  PNU-292137–1VYW  : $PNU_MEAN_RMSD Å"
echo "  K03861–5A14      : $K_MEAN_RMSD Å"


echo
echo "Complex PDB files:"
echo "  results/7RWF_TW8672/TW8672_7RWF_complex.pdb"
echo "  results/1VYW_PNU292137/PNU292137_1VYW_complex.pdb"
echo "  results/5A14_K03861/K03861_5A14_complex.pdb"


echo
echo "RMSD files:"
echo "  results/7RWF_TW8672/TW8672_RMSD.txt"
echo "  results/1VYW_PNU292137/PNU292137_RMSD.txt"
echo "  results/5A14_K03861/K03861_RMSD.txt"


echo
echo "============================================================"

date

echo "============================================================"