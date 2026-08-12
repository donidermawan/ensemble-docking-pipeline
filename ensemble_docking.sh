#!/bin/bash

set -u

# ============================================================
# CDK2 ENSEMBLE DOCKING + CONSENSUS RANKING
# ============================================================
#
# Initial selection:
#
#       Affinity <= -10.00 kcal/mol
#
# Ensemble receptors:
#
#       1VYW
#       7RWF
#       5A14
#       5IF1
#       1W98
#
# Parallelization:
#
#       12 simultaneous Vina processes per receptor
#       1 CPU per Vina process
#
# Output:
#
#       - Selected ligand list
#       - Individual receptor affinity CSV
#       - Individual receptor ranked CSV
#       - Docked PDBQT poses
#       - Receptor-ligand PDB complexes
#       - Complete ensemble summary
#       - Consensus-ranked ligands
#       - All-ligand statistical ranking
#
# Consensus ranking:
#
#       Primary:
#           Mean affinity
#
#       Secondary:
#           Median affinity
#
#       Tertiary:
#           Standard deviation
#
# More negative affinity = more favorable predicted binding.
#
# ============================================================


# ============================================================
# SETTINGS
# ============================================================

LIGAND_DIR="pdbqt"

SCREENING_FILE="vina_results_1VYW_1-7629.csv"

THRESHOLD="-10.00"

JOBS=12


# ============================================================
# FIVE CDK2 RECEPTORS
# ============================================================

RECEPTORS=(
    "1VYW"
    "7RWF"
    "5A14"
    "5IF1"
    "1W98"
)

N_RECEPTORS=${#RECEPTORS[@]}


# ============================================================
# OUTPUT DIRECTORIES
# ============================================================

BASE_OUT="ensemble_results"

SELECTED_DIR="${BASE_OUT}/selected_ligands"
DOCKING_DIR="${BASE_OUT}/docking"
COMPLEX_DIR="${BASE_OUT}/complexes"
LOG_DIR="${BASE_OUT}/logs"
SUMMARY_DIR="${BASE_OUT}/summaries"
RECEPTOR_PDB_DIR="${BASE_OUT}/receptor_pdb"


mkdir -p "$SELECTED_DIR"
mkdir -p "$DOCKING_DIR"
mkdir -p "$COMPLEX_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$SUMMARY_DIR"
mkdir -p "$RECEPTOR_PDB_DIR"


# ============================================================
# OUTPUT FILES
# ============================================================

SELECTED_CSV="${SELECTED_DIR}/selected_ligands.csv"
SELECTED_LIST="${SELECTED_DIR}/selected_ligands.txt"

ENSEMBLE="${SUMMARY_DIR}/ensemble_docking_summary.csv"

RANKED_OUTPUT="${SUMMARY_DIR}/ensemble_docking_ranked.csv"

ALL_RANKED_OUTPUT="${SUMMARY_DIR}/ensemble_docking_all_ligands_ranked.csv"


# ============================================================
# START
# ============================================================

echo "=============================================="
echo "CDK2 ENSEMBLE DOCKING"
echo "=============================================="
echo

date

echo

echo "Ensemble receptors:"
printf "  %s\n" "${RECEPTORS[@]}"

echo
echo "Number of receptors: ${N_RECEPTORS}"

echo


# ============================================================
# CHECK INPUTS
# ============================================================

if [ ! -f "$SCREENING_FILE" ]; then

    echo "ERROR: Screening file not found:"
    echo "  $SCREENING_FILE"

    exit 1

fi


if [ ! -d "$LIGAND_DIR" ]; then

    echo "ERROR: Ligand directory not found:"
    echo "  $LIGAND_DIR"

    exit 1

fi


# ============================================================
# CHECK RECEPTORS AND CONFIGURATION FILES
# ============================================================

echo "Checking receptor files..."
echo


for receptor in "${RECEPTORS[@]}"
do

    if [ ! -f "${receptor}.pdbqt" ]; then

        echo "ERROR: Missing receptor:"
        echo "  ${receptor}.pdbqt"

        exit 1

    fi


    if [ ! -f "vina_${receptor}.txt" ]; then

        echo "ERROR: Missing configuration:"
        echo "  vina_${receptor}.txt"

        exit 1

    fi

done


echo "All receptor/configuration files found."

echo


# ============================================================
# SELECT LIGANDS
# ============================================================

echo "=============================================="
echo "INITIAL LIGAND SELECTION"
echo "=============================================="

echo

echo "Screening file:"
echo "  $SCREENING_FILE"

echo

echo "Selection threshold:"
echo "  Affinity <= ${THRESHOLD} kcal/mol"

echo


rm -f "$SELECTED_CSV"
rm -f "$SELECTED_LIST"


awk -F',' -v threshold="$THRESHOLD" '

NR > 1 {

    ligand=$1
    affinity=$2

    gsub(/[[:space:]]/, "", ligand)
    gsub(/[[:space:]]/, "", affinity)

    if (affinity ~ /^-?[0-9]+([.][0-9]+)?$/) {

        if ((affinity + 0) <= (threshold + 0)) {

            print ligand "," affinity

        }

    }

}

' "$SCREENING_FILE" > "$SELECTED_CSV"


cut -d',' -f1 "$SELECTED_CSV" > "$SELECTED_LIST"


SELECTED_COUNT=$(wc -l < "$SELECTED_LIST" | tr -d ' ')


echo "Selected ligands: $SELECTED_COUNT"

echo


if [ "$SELECTED_COUNT" -eq 0 ]; then

    echo "ERROR: No ligands passed the threshold."

    exit 1

fi


echo "First 10 selected ligands:"
head -n 10 "$SELECTED_CSV"

echo


# ============================================================
# PREPARE RECEPTOR PDB FILES
# ============================================================

echo "=============================================="
echo "PREPARING RECEPTOR PDB FILES"
echo "=============================================="

echo


for receptor in "${RECEPTORS[@]}"
do

    receptor_pdb="${RECEPTOR_PDB_DIR}/${receptor}.pdb"


    if [ ! -f "$receptor_pdb" ]; then

        echo "Converting ${receptor}.pdbqt -> ${receptor}.pdb"


        obabel \
            "${receptor}.pdbqt" \
            -O "$receptor_pdb" \
            >/dev/null 2>&1


        if [ $? -ne 0 ]; then

            echo "ERROR converting ${receptor}.pdbqt"

            exit 1

        fi

    else

        echo "Already exists: $receptor_pdb"

    fi

done


echo

echo "Receptor PDB files ready."

echo


# ============================================================
# DOCKING
# ============================================================

for receptor in "${RECEPTORS[@]}"
do

    echo
    echo "=============================================="
    echo "DOCKING AGAINST: ${receptor}"
    echo "=============================================="
    echo

    date

    echo


    RECEPTOR_FILE="${receptor}.pdbqt"

    CONFIG_FILE="vina_${receptor}.txt"

    OUT_DIR="${DOCKING_DIR}/${receptor}"

    PDB_DIR="${COMPLEX_DIR}/${receptor}"

    RECEPTOR_LOG_DIR="${LOG_DIR}/${receptor}"

    SUMMARY="${SUMMARY_DIR}/vina_results_${receptor}.csv"

    SORTED_SUMMARY="${SUMMARY_DIR}/vina_results_${receptor}_ranked.csv"

    PARALLEL_LOG="${RECEPTOR_LOG_DIR}/parallel.log"


    mkdir -p "$OUT_DIR"
    mkdir -p "$PDB_DIR"
    mkdir -p "$RECEPTOR_LOG_DIR"


    rm -f "$SUMMARY"
    rm -f "$SORTED_SUMMARY"


    echo "Ligand,Affinity" > "$SUMMARY"


    # ========================================================
    # EXPORT VARIABLES FOR GNU PARALLEL
    # ========================================================

    export LIGAND_DIR
    export RECEPTOR_FILE
    export CONFIG_FILE
    export OUT_DIR
    export RECEPTOR_LOG_DIR


    # ========================================================
    # RUN VINA
    # ========================================================

    echo "Running ${JOBS} Vina processes in parallel..."

    echo


    parallel \
        --bar \
        --line-buffer \
        -j "$JOBS" \
        --joblog "$PARALLEL_LOG" \
        --env LIGAND_DIR \
        --env RECEPTOR_FILE \
        --env CONFIG_FILE \
        --env OUT_DIR \
        --env RECEPTOR_LOG_DIR \
        '
        ligand_name="{}"

        ligand_file="${LIGAND_DIR}/${ligand_name}.pdbqt"


        if [ ! -f "$ligand_file" ]; then

            printf "%s,FAILED\n" "$ligand_name"

            exit 0

        fi


        vina_output="${OUT_DIR}/${ligand_name}_out.pdbqt"


        vina \
            --receptor "$RECEPTOR_FILE" \
            --ligand "$ligand_file" \
            --config "$CONFIG_FILE" \
            --cpu 1 \
            --out "$vina_output" \
            > "${RECEPTOR_LOG_DIR}/${ligand_name}.log" 2>&1


        affinity=$(
            awk '"'"'
                /^[[:space:]]*1[[:space:]]+/ {
                    print $2
                    exit
                }
            '"'"' "${RECEPTOR_LOG_DIR}/${ligand_name}.log"
        )


        if [ -z "$affinity" ]; then

            affinity="FAILED"

        fi


        printf "%s,%s\n" "$ligand_name" "$affinity"

        ' \
        :::: "$SELECTED_LIST" \
        >> "$SUMMARY"


    # ========================================================
    # SORT INDIVIDUAL RECEPTOR RESULTS
    # ========================================================

    {

        head -n 1 "$SUMMARY"


        tail -n +2 "$SUMMARY" |

        awk -F',' '
            $2 ~ /^-?[0-9]+([.][0-9]+)?$/ {
                print
            }
        ' |

        sort -t',' -k2,2n


    } > "$SORTED_SUMMARY"


    # ========================================================
    # CREATE PDB COMPLEXES
    # ========================================================

    echo

    echo "Creating PDB complexes for ${receptor}..."

    echo


    export RECEPTOR_PDB_DIR

    export PDB_DIR


    find "$OUT_DIR" \
        -name "*_out.pdbqt" \
        -print0 |

    parallel \
        --null \
        --bar \
        -j "$JOBS" \
        '
        pdbqt="{}"

        filename=$(basename "$pdbqt" _out.pdbqt)

        receptor_pdb="${RECEPTOR_PDB_DIR}/'"${receptor}"'.pdb"

        complex_dir="$PDB_DIR"


        # ----------------------------------------------------
        # Extract MODEL 1
        # ----------------------------------------------------

        pose_pdbqt="${complex_dir}/${filename}_pose1.pdbqt"


        awk '"'"'
            /^MODEL[[:space:]]+1/ {
                capture=1
            }

            capture {
                print
            }

            /^ENDMDL/ && capture {
                exit
            }

        '"'"' "$pdbqt" > "$pose_pdbqt"


        # ----------------------------------------------------
        # Convert ligand pose
        # ----------------------------------------------------

        pose_pdb="${complex_dir}/${filename}_ligand.pdb"


        obabel \
            "$pose_pdbqt" \
            -O "$pose_pdb" \
            >/dev/null 2>&1


        # ----------------------------------------------------
        # Combine receptor + ligand
        # ----------------------------------------------------

        complex_pdb="${complex_dir}/${filename}_complex.pdb"


        {

            cat "$receptor_pdb"

            echo "TER"


            awk '"'"'
                /^ATOM/ || /^HETATM/ || /^CONECT/ {
                    print
                }
            '"'"' "$pose_pdb"


            echo "END"


        } > "$complex_pdb"


        # ----------------------------------------------------
        # Remove temporary files
        # ----------------------------------------------------

        rm -f "$pose_pdbqt"

        rm -f "$pose_pdb"

        '


    # ========================================================
    # RECEPTOR STATISTICS
    # ========================================================

    SUCCESS=$(awk -F',' '

        NR > 1 &&
        $2 ~ /^-?[0-9]+([.][0-9]+)?$/ {

            count++

        }

        END {

            print count+0

        }

    ' "$SUMMARY")


    FAILED=$(awk -F',' '

        NR > 1 &&
        $2 == "FAILED" {

            count++

        }

        END {

            print count+0

        }

    ' "$SUMMARY")


    PDBQT_COUNT=$(find "$OUT_DIR" \
        -name "*_out.pdbqt" |
        wc -l |
        tr -d ' ')


    PDB_COUNT=$(find "$PDB_DIR" \
        -name "*_complex.pdb" |
        wc -l |
        tr -d ' ')


    # ========================================================
    # DISPLAY RECEPTOR RESULTS
    # ========================================================

    echo

    echo "----------------------------------------------"

    echo "Results for ${receptor}"

    echo "----------------------------------------------"

    echo "Selected ligands : $SELECTED_COUNT"

    echo "Successful       : $SUCCESS"

    echo "Failed            : $FAILED"

    echo "PDBQT poses       : $PDBQT_COUNT"

    echo "PDB complexes     : $PDB_COUNT"

    echo

    echo "CSV:"

    echo "  $SUMMARY"

    echo

    echo "Ranked CSV:"

    echo "  $SORTED_SUMMARY"

    echo

    echo "Complex folder:"

    echo "  $PDB_DIR"

    echo

    echo "----------------------------------------------"


done


# ============================================================
# GENERATE ENSEMBLE SUMMARY + CONSENSUS RANKING
# ============================================================

echo

echo "=============================================="

echo "GENERATING ENSEMBLE SUMMARY + CONSENSUS RANKING"

echo "=============================================="

echo


python3 - \
    "$SELECTED_CSV" \
    "$SUMMARY_DIR" \
    "$ENSEMBLE" \
    "$RANKED_OUTPUT" \
    "$ALL_RANKED_OUTPUT" \
    << 'PY'


import csv
import statistics
import sys
from pathlib import Path


# ============================================================
# COMMAND-LINE ARGUMENTS
# ============================================================

selected_file = Path(sys.argv[1])

summary_dir = Path(sys.argv[2])

output_file = Path(sys.argv[3])

ranked_output = Path(sys.argv[4])

all_ranked_output = Path(sys.argv[5])


# ============================================================
# FIVE ENSEMBLE RECEPTORS
# ============================================================

receptors = [
    "1VYW",
    "7RWF",
    "5A14",
    "5IF1",
    "1W98"
]


N_RECEPTORS = len(receptors)


# ============================================================
# READ SELECTED LIGANDS
# ============================================================

initial = {}


with open(selected_file, newline="") as f:

    reader = csv.reader(f)


    for row in reader:

        if len(row) >= 2:

            ligand = row[0].strip()

            affinity = row[1].strip()


            if not ligand:
                continue


            try:

                initial[ligand] = float(affinity)

            except ValueError:

                pass


print(
    f"Selected ligands loaded: {len(initial)}"
)


# ============================================================
# READ INDIVIDUAL RECEPTOR RESULTS
# ============================================================

results = {}


for receptor in receptors:

    file = summary_dir / f"vina_results_{receptor}.csv"


    results[receptor] = {}


    if not file.exists():

        print(
            f"WARNING: Missing receptor file: {file}"
        )

        continue


    with open(file, newline="") as f:

        reader = csv.DictReader(f)


        for row in reader:

            ligand = row.get(
                "Ligand",
                ""
            ).strip()


            affinity = row.get(
                "Affinity",
                ""
            ).strip()


            if not ligand:
                continue


            try:

                results[receptor][ligand] = float(
                    affinity
                )

            except (ValueError, TypeError):

                results[receptor][ligand] = None


    print(
        f"{receptor}: "
        f"{len(results[receptor])} docking results loaded"
    )


# ============================================================
# PREPARE ENSEMBLE DATA
# ============================================================

ensemble_data = []


for ligand in initial.keys():

    receptor_values = {}

    successful = 0


    for receptor in receptors:

        affinity = results[receptor].get(
            ligand
        )


        if affinity is None:

            receptor_values[receptor] = None

        else:

            receptor_values[receptor] = affinity

            successful += 1


    # --------------------------------------------------------
    # Numerical values
    # --------------------------------------------------------

    values = [
        value
        for value in receptor_values.values()
        if value is not None
    ]


    # --------------------------------------------------------
    # Statistics
    # --------------------------------------------------------

    if values:

        mean_affinity = statistics.mean(values)

        median_affinity = statistics.median(values)

        best_affinity = min(values)

        worst_affinity = max(values)


        if len(values) >= 2:

            sd_affinity = statistics.stdev(values)

        else:

            sd_affinity = 0.0


    else:

        mean_affinity = None

        median_affinity = None

        best_affinity = None

        worst_affinity = None

        sd_affinity = None


    # --------------------------------------------------------
    # Number of receptors <= -10 kcal/mol
    # --------------------------------------------------------

    strong_hits = sum(
        1
        for value in values
        if value <= -10.0
    )


    # --------------------------------------------------------
    # Consensus classification
    # --------------------------------------------------------

    if successful == 0:

        consensus = "FAILED"

    elif (
        successful == N_RECEPTORS
        and
        strong_hits == N_RECEPTORS
    ):

        consensus = "Strong_consensus"

    elif strong_hits >= 4:

        consensus = "High_consensus"

    elif strong_hits >= 3:

        consensus = "Moderate_consensus"

    elif strong_hits >= 2:

        consensus = "Weak_consensus"

    elif strong_hits == 1:

        consensus = "Single_receptor_hit"

    else:

        consensus = "No_consensus"


    ensemble_data.append({

        "Ligand": ligand,

        "Initial_1VYW_Affinity":
            initial[ligand],

        "values":
            receptor_values,

        "Successful_Receptors":
            successful,

        "Mean_Affinity":
            mean_affinity,

        "Median_Affinity":
            median_affinity,

        "Best_Affinity":
            best_affinity,

        "Worst_Affinity":
            worst_affinity,

        "SD_Affinity":
            sd_affinity,

        "Strong_Hits":
            strong_hits,

        "Consensus":
            consensus

    })


# ============================================================
# SORT FUNCTION FOR LIGAND NAMES
# ============================================================

def ligand_sort_key(item):

    ligand = item["Ligand"]

    last_part = ligand.split("_")[-1]


    if last_part.isdigit():

        return (
            0,
            int(last_part)
        )


    return (
        1,
        ligand
    )


# ============================================================
# WRITE COMPLETE ENSEMBLE SUMMARY
# ============================================================

header = [

    "Ligand",

    "Initial_1VYW_Affinity",

    "1VYW",

    "7RWF",

    "5A14",

    "5IF1",

    "1W98",

    "Mean_Affinity",

    "Median_Affinity",

    "Best_Affinity",

    "Worst_Affinity",

    "SD_Affinity",

    "Successful_Receptors",

    "Receptors_<=_-10",

    "Consensus"

]


with open(
    output_file,
    "w",
    newline=""
) as f:

    writer = csv.writer(f)

    writer.writerow(header)


    for item in sorted(
        ensemble_data,
        key=ligand_sort_key
    ):

        row = [

            item["Ligand"],

            f'{item["Initial_1VYW_Affinity"]:.3f}'

        ]


        # ----------------------------------------------------
        # Receptor affinities
        # ----------------------------------------------------

        for receptor in receptors:

            value = item["values"][receptor]


            if value is None:

                row.append("FAILED")

            else:

                row.append(
                    f"{value:.3f}"
                )


        # ----------------------------------------------------
        # Ensemble statistics
        # ----------------------------------------------------

        mean_value = item["Mean_Affinity"]

        median_value = item["Median_Affinity"]

        best_value = item["Best_Affinity"]

        worst_value = item["Worst_Affinity"]

        sd_value = item["SD_Affinity"]


        row.extend([

            f"{mean_value:.3f}"
            if mean_value is not None
            else "NA",

            f"{median_value:.3f}"
            if median_value is not None
            else "NA",

            f"{best_value:.3f}"
            if best_value is not None
            else "NA",

            f"{worst_value:.3f}"
            if worst_value is not None
            else "NA",

            f"{sd_value:.3f}"
            if sd_value is not None
            else "NA",

            item["Successful_Receptors"],

            item["Strong_Hits"],

            item["Consensus"]

        ])


        writer.writerow(row)


print()

print(
    "Created complete ensemble summary:"
)

print(
    f"  {output_file}"
)


# ============================================================
# CONSENSUS RANKING
# ============================================================
#
# ONLY ligands successfully docked against ALL FIVE
# receptors are included.
#
# Ranking:
#
#   1. Mean affinity
#   2. Median affinity
#   3. SD
#
# More negative mean = better.
#
# ============================================================

complete = [

    item

    for item in ensemble_data

    if item["Successful_Receptors"] == N_RECEPTORS

]


print()

print(
    f"Ligands successfully docked against "
    f"all {N_RECEPTORS} receptors: "
    f"{len(complete)}"
)


# ============================================================
# SORT CONSENSUS HITS
# ============================================================

complete.sort(

    key=lambda x: (

        x["Mean_Affinity"],

        x["Median_Affinity"],

        x["SD_Affinity"]

    )

)


# ============================================================
# WRITE CONSENSUS RANKING
# ============================================================

ranked_header = [

    "Rank",

    "Ligand",

    "Initial_1VYW_Affinity",

    "Mean_Affinity",

    "Median_Affinity",

    "Best_Affinity",

    "Worst_Affinity",

    "SD_Affinity",

    "Successful_Receptors",

    "Receptors_<=_-10",

    "Consensus",

    "1VYW",

    "7RWF",

    "5A14",

    "5IF1",

    "1W98"

]


with open(
    ranked_output,
    "w",
    newline=""
) as f:

    writer = csv.writer(f)

    writer.writerow(ranked_header)


    for rank, item in enumerate(
        complete,
        start=1
    ):

        row = [

            rank,

            item["Ligand"],

            f'{item["Initial_1VYW_Affinity"]:.3f}',

            f'{item["Mean_Affinity"]:.3f}',

            f'{item["Median_Affinity"]:.3f}',

            f'{item["Best_Affinity"]:.3f}',

            f'{item["Worst_Affinity"]:.3f}',

            f'{item["SD_Affinity"]:.3f}',

            item["Successful_Receptors"],

            item["Strong_Hits"],

            item["Consensus"]

        ]


        for receptor in receptors:

            value = item["values"][receptor]

            row.append(
                f"{value:.3f}"
            )


        writer.writerow(row)


print()

print(
    "Created consensus ranking:"
)

print(
    f"  {ranked_output}"
)


# ============================================================
# ALL-LIGAND STATISTICAL RANKING
# ============================================================
#
# Includes ligands with at least one successful docking.
#
# This is useful for investigating partial docking results.
#
# ============================================================

all_valid = [

    item

    for item in ensemble_data

    if item["Successful_Receptors"] > 0

]


all_valid.sort(

    key=lambda x: (

        x["Mean_Affinity"],

        x["Median_Affinity"]
        if x["Median_Affinity"] is not None
        else 999999,

        x["SD_Affinity"]
        if x["SD_Affinity"] is not None
        else 999999

    )

)


with open(
    all_ranked_output,
    "w",
    newline=""
) as f:

    writer = csv.writer(f)


    writer.writerow([

        "Rank",

        "Ligand",

        "Initial_1VYW_Affinity",

        "Mean_Affinity",

        "Median_Affinity",

        "Best_Affinity",

        "Worst_Affinity",

        "SD_Affinity",

        "Successful_Receptors",

        "Receptors_<=_-10",

        "Consensus"

    ])


    for rank, item in enumerate(
        all_valid,
        start=1
    ):

        writer.writerow([

            rank,

            item["Ligand"],

            f'{item["Initial_1VYW_Affinity"]:.3f}',

            f'{item["Mean_Affinity"]:.3f}',

            f'{item["Median_Affinity"]:.3f}',

            f'{item["Best_Affinity"]:.3f}',

            f'{item["Worst_Affinity"]:.3f}',

            f'{item["SD_Affinity"]:.3f}',

            item["Successful_Receptors"],

            item["Strong_Hits"],

            item["Consensus"]

        ])


print()

print(
    "Created all-ligand statistical ranking:"
)

print(
    f"  {all_ranked_output}"
)


# ============================================================
# FINAL PYTHON STATISTICS
# ============================================================

strong_consensus = [

    item

    for item in ensemble_data

    if item["Consensus"] == "Strong_consensus"

]


high_consensus = [

    item

    for item in ensemble_data

    if item["Consensus"] == "High_consensus"

]


print()

print("==============================================")

print("CONSENSUS STATISTICS")

print("==============================================")

print()

print(
    f"Total selected ligands       : {len(initial)}"
)

print(
    f"Complete 5/5 docking         : {len(complete)}"
)

print(
    f"Strong consensus (5/5 <=-10) : "
    f"{len(strong_consensus)}"
)

print(
    f"High consensus (>=4 <=-10)   : "
    f"{len(high_consensus)}"
)

print()


PY


# ============================================================
# VERIFY ENSEMBLE OUTPUT
# ============================================================

echo

echo "=============================================="

echo "VERIFYING ENSEMBLE OUTPUT"

echo "=============================================="

echo


if [ ! -f "$ENSEMBLE" ]; then

    echo "ERROR: Ensemble summary was not created:"
    echo "  $ENSEMBLE"

    exit 1

fi


if [ ! -f "$RANKED_OUTPUT" ]; then

    echo "ERROR: Consensus ranking was not created:"
    echo "  $RANKED_OUTPUT"

    exit 1

fi


if [ ! -f "$ALL_RANKED_OUTPUT" ]; then

    echo "ERROR: All-ligand ranking was not created:"
    echo "  $ALL_RANKED_OUTPUT"

    exit 1

fi


ROWS=$(tail -n +2 "$ENSEMBLE" | wc -l | tr -d ' ')

RANKED_ROWS=$(tail -n +2 "$RANKED_OUTPUT" | wc -l | tr -d ' ')

ALL_ROWS=$(tail -n +2 "$ALL_RANKED_OUTPUT" | wc -l | tr -d ' ')


echo "Ensemble summary:"
echo "  $ENSEMBLE"
echo "  Ligands: $ROWS"

echo

echo "Consensus ranking:"
echo "  $RANKED_OUTPUT"
echo "  5/5 ligands: $RANKED_ROWS"

echo

echo "All-ligand ranking:"
echo "  $ALL_RANKED_OUTPUT"
echo "  Ligands with >=1 successful receptor: $ALL_ROWS"


# ============================================================
# DISPLAY TOP 20 CONSENSUS LIGANDS
# ============================================================

echo

echo "=============================================="

echo "TOP 20 CONSENSUS LIGANDS"

echo "=============================================="

echo


head -n 21 "$RANKED_OUTPUT"


# ============================================================
# FINAL SUMMARY
# ============================================================

echo

echo "=============================================="

echo "CDK2 ENSEMBLE DOCKING COMPLETED"

echo "=============================================="

echo


echo "Initial screening:"

echo "  Screening file   : $SCREENING_FILE"

echo "  Selected ligands : $SELECTED_COUNT"

echo "  Threshold        : <= ${THRESHOLD} kcal/mol"


echo

echo "Ensemble receptors:"

printf "  %s\n" "${RECEPTORS[@]}"


echo

echo "Output files:"

echo


echo "1. Selected ligands:"

echo "   $SELECTED_CSV"


echo

echo "2. Individual receptor results:"

echo "   ${SUMMARY_DIR}/vina_results_*.csv"


echo

echo "3. Individual receptor ranked results:"

echo "   ${SUMMARY_DIR}/vina_results_*_ranked.csv"


echo

echo "4. Complete ensemble summary:"

echo "   $ENSEMBLE"


echo

echo "5. Consensus ranking (5/5 receptors):"

echo "   $RANKED_OUTPUT"


echo

echo "6. All-ligand statistical ranking:"

echo "   $ALL_RANKED_OUTPUT"


echo

echo "7. Docking poses:"

echo "   $DOCKING_DIR/"


echo

echo "8. PDB complexes:"

echo "   $COMPLEX_DIR/"


echo

echo "9. Receptor PDB files:"

echo "   $RECEPTOR_PDB_DIR/"


echo

echo "10. Logs:"

echo "    $LOG_DIR/"


echo

echo "=============================================="

date