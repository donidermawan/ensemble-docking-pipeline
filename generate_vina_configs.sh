#!/bin/bash

# ============================================================
# Generate AutoDock Vina configuration files
# for ensemble docking against five CDK2 receptor structures
# ============================================================

SIZE_X=24
SIZE_Y=24
SIZE_Z=24

EXHAUSTIVENESS=8
NUM_MODES=3
ENERGY_RANGE=3

# ============================================================
# 1VYW
# ============================================================

cat > vina_1VYW.txt << EOF
receptor = 1VYW.pdbqt

center_x = -8.929
center_y = 205.749
center_z = 111.606

size_x = ${SIZE_X}
size_y = ${SIZE_Y}
size_z = ${SIZE_Z}

exhaustiveness = ${EXHAUSTIVENESS}
num_modes = ${NUM_MODES}
energy_range = ${ENERGY_RANGE}
EOF

# ============================================================
# 7RWF
# ============================================================

cat > vina_7RWF.txt << EOF
receptor = 7RWF.pdbqt

center_x = -22.195
center_y = 2.569
center_z = 17.819

size_x = ${SIZE_X}
size_y = ${SIZE_Y}
size_z = ${SIZE_Z}

exhaustiveness = ${EXHAUSTIVENESS}
num_modes = ${NUM_MODES}
energy_range = ${ENERGY_RANGE}
EOF

# ============================================================
# 5A14
# ============================================================

cat > vina_5A14.txt << EOF
receptor = 5A14.pdbqt

center_x = 22.287
center_y = 8.979
center_z = 10.665

size_x = ${SIZE_X}
size_y = ${SIZE_Y}
size_z = ${SIZE_Z}

exhaustiveness = ${EXHAUSTIVENESS}
num_modes = ${NUM_MODES}
energy_range = ${ENERGY_RANGE}
EOF

# ============================================================
# 5IF1
# ============================================================

cat > vina_5IF1.txt << EOF
receptor = 5IF1.pdbqt

center_x = -73.670
center_y = -55.376
center_z = 19.600

size_x = ${SIZE_X}
size_y = ${SIZE_Y}
size_z = ${SIZE_Z}

exhaustiveness = ${EXHAUSTIVENESS}
num_modes = ${NUM_MODES}
energy_range = ${ENERGY_RANGE}
EOF

# ============================================================
# 1W98
# ============================================================

cat > vina_1W98.txt << EOF
receptor = 1W98.pdbqt

center_x = 24.006
center_y = 20.821
center_z = -14.838

size_x = ${SIZE_X}
size_y = ${SIZE_Y}
size_z = ${SIZE_Z}

exhaustiveness = ${EXHAUSTIVENESS}
num_modes = ${NUM_MODES}
energy_range = ${ENERGY_RANGE}
EOF

# ============================================================
# FINISHED
# ============================================================

echo "=============================================="
echo "Vina configuration files generated"
echo "=============================================="

ls -lh vina_*.txt

echo
echo "Settings:"
echo "Grid size       : ${SIZE_X} x ${SIZE_Y} x ${SIZE_Z} Å"
echo "Exhaustiveness  : ${EXHAUSTIVENESS}"
echo "Number of modes : ${NUM_MODES}"
echo "Energy range    : ${ENERGY_RANGE} kcal/mol"