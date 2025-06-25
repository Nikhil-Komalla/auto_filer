#!/bin/bash

# Get current directory name automatically
FOLDER=$(basename "$PWD")
echo "Working with folder: $FOLDER"

# Find base directory with Old_incorrect and New_H_ads
BASE="$PWD"
while [[ "$BASE" != "/" && ! (-d "$BASE/Old_incorrect" && -d "$BASE/New_H_ads") ]]; do
    BASE=$(dirname "$BASE")
done

if [[ "$BASE" == "/" ]]; then
    echo "Error: Cannot find Old_incorrect and New_H_ads directories"
    exit 1
fi

echo "Base directory: $BASE"

# Set paths
OLD="$BASE/Old_incorrect/$FOLDER"
NEW="$BASE/New_H_ads/$FOLDER"

# Create directories
mkdir -p "$NEW"/{surf+mol,surf,mol}

# Step 1: Copy main files
echo "Copying main files..."
cd "$OLD/surf+mol"
cp INCAR KPOINTS POSCAR POTCAR o_run.sh CHGCAR "$NEW/surf+mol/"

# Step 2: Distribute to surf and mol
echo "Distributing files..."
cd "$NEW/surf+mol"
cp POSCAR o_run.sh ../surf/
cp POSCAR o_run.sh ../mol/

# Step 3: Find reference (try Cu72 first, then ask)
REF="Cu72"
if [[ ! -d "$BASE/New_H_ads/$REF" && ! -d "$BASE/Old_incorrect/$REF" ]]; then
    echo "Enter reference folder name (default: Cu72):"
    read -r input
    REF=${input:-Cu72}
fi

# Step 4: Copy reference files to surf
echo "Copying reference files..."
cd "$NEW/surf"
if [[ -d "$BASE/New_H_ads/$REF/surf" ]]; then
    cp "$BASE/New_H_ads/$REF/surf"/{INCAR,KPOINTS,POTCAR,o_run.sh} .
else
    cp "$BASE/Old_incorrect/$REF/surf"/{INCAR,KPOINTS,POTCAR,o_run.sh} .
fi

# Step 5: Copy reference files to mol
cd "$NEW/mol"
if [[ -d "$BASE/New_H_ads/$REF/mol" ]]; then
    cp "$BASE/New_H_ads/$REF/mol"/{INCAR,KPOINTS,POTCAR,o_run.sh} .
else
    cp "$BASE/Old_incorrect/$REF/mol"/{INCAR,KPOINTS,POTCAR,o_run.sh} .
fi

# Step 6: Edit surf POSCAR
echo "Modifying surf POSCAR..."
cd "$NEW/surf"
cp POSCAR POSCAR.backup
sed -i '1s/H$//' POSCAR        # Remove H from end of line 1
sed -i '6s/H$//' POSCAR        # Remove H from end of line 6
sed -i '7s/1$//' POSCAR        # Remove 1 from end of line 7
sed -i '$d' POSCAR             # Remove last line

# Step 7: Edit mol POSCAR
echo "Modifying mol POSCAR..."
cd "$NEW/mol"
cp POSCAR POSCAR.backup
total_lines=$(wc -l < POSCAR)
sed -i '1s/.*/H/' POSCAR        # Replace line 1 with H
sed -i '6s/.*/H/' POSCAR        # Replace line 6 with H
sed -i '7s/.*/1/' POSCAR        # Replace line 7 with 1
# Remove lines 10 to second-to-last
if [[ $total_lines -gt 10 ]]; then
    sed -i "10,$((total_lines-1))d" POSCAR
fi

echo "Done! Files copied and modified for: $FOLDER"
echo "Location: $NEW"
