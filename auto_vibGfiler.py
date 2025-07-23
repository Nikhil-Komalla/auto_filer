#!/usr/bin/env python3

import os
import shutil
from pathlib import Path
import re

# --- Step 1: Setup paths ---
cwd = Path.cwd()
folder_name = cwd.name
ref_dir = cwd.parents[1] / folder_name

print(f"\n📂 Target directory         : {cwd}")
print(f"📁 Reference source path    : {ref_dir}")

# Check if reference directory exists
if not ref_dir.is_dir():
    print("❌ ERROR: Reference folder not found at:", ref_dir)
    exit(1)

# --- Step 2: Copy files from reference location ---
files_to_copy = ['INCAR', 'CONTCAR', 'KPOINTS', 'POTCAR', 'run.sh', 'o_run.sh']
print("\n📦 Copying files...")
for fname in files_to_copy:
    source = ref_dir / fname
    if source.exists():
        shutil.copy(source, cwd / fname)
        print(f"✓ Copied: {fname}")
    else:
        print(f"⚠️  Warning: {fname} not found in reference folder.")

# --- Step 3: Edit INCAR ---
incar_path = cwd / 'INCAR'
if incar_path.exists():
    print("\n🛠️  Editing INCAR...")
    shutil.copy(incar_path, incar_path.with_suffix(".backup"))

    with open(incar_path, 'r') as f:
        content = f.read()
    
    # FIXED: Use \g<1> syntax to avoid group reference ambiguity
    content = re.sub(r'^(\s*PREC\s*=\s*).*$', r'\g<1>Accurate', content, flags=re.MULTILINE)
    content = re.sub(r'^(\s*IBRION\s*=\s*).*$', r'\g<1>5', content, flags=re.MULTILINE)
    content = re.sub(r'^(\s*ISMEAR\s*=\s*).*$', r'\g<1>0', content, flags=re.MULTILINE)
    
    # Add missing parameters if not found
    if not re.search(r'^\s*PREC\s*=', content, re.MULTILINE):
        content += "\nPREC = Accurate\n"
    if not re.search(r'^\s*IBRION\s*=', content, re.MULTILINE):
        content += "IBRION = 5\n"
    if not re.search(r'^\s*ISMEAR\s*=', content, re.MULTILINE):
        content += "ISMEAR = 0\n"

    with open(incar_path, 'w') as f:
        f.write(content)

    print("✓ INCAR updated: PREC=Accurate, IBRION=5, ISMEAR=0")
else:
    print("⚠️  INCAR not found, skipped editing.")

# --- Step 4: Edit CONTCAR Constraints (Preserve Spacing) ---
contcar_path = cwd / 'CONTCAR'
if not contcar_path.exists():
    print("⚠️  CONTCAR not found. Skipping POSCAR creation.")
    exit(0)

print("\n🛠️  Editing CONTCAR...")
shutil.copy(contcar_path, contcar_path.with_suffix(".backup"))

with open(contcar_path, 'r') as f:
    content = f.read()

# --- Identify coordinate start ---
lines = content.split('\n')
coord_start = None
for idx, line in enumerate(lines):
    if line.strip().lower().startswith(('direct', 'cartesian')):
        coord_start = idx + 1
        break

if coord_start is None:
    print("❌ ERROR: Could not find 'Direct' or 'Cartesian' keyword in CONTCAR.")
    exit(1)

# --- Get number of atoms from line 7 ---
try:
    atom_counts = list(map(int, lines[6].split()))
    natoms = sum(atom_counts)
except Exception as e:
    print(f"❌ ERROR reading number of atoms on line 7: {e}")
    exit(1)

coord_end = coord_start + natoms

# --- Modify constraint flags in place (preserving spacing) ---
for i in range(coord_start, coord_end):
    if i < len(lines):
        line = lines[i]
        # Check if line has existing constraints (T/F pattern at end)
        if re.search(r'[TF]\s+[TF]\s+[TF]\s*$', line):
            if i == coord_end - 1:  # Last atom
                lines[i] = re.sub(r'[TF](\s+)[TF](\s+)[TF](\s*)$', r'T\g<1>T\g<2>T\g<3>', line)
            else:  # All other atoms
                lines[i] = re.sub(r'[TF](\s+)[TF](\s+)[TF](\s*)$', r'F\g<1>F\g<2>F\g<3>', line)
        else:
            # No existing constraints, add them
            if i == coord_end - 1:
                lines[i] = line.rstrip() + " T T T"
            else:
                lines[i] = line.rstrip() + " F F F"

# --- Step 5: Write POSCAR ---
poscar_path = cwd / 'POSCAR'
with open(poscar_path, 'w') as f:
    f.write('\n'.join(lines))

print("✓ Constraints modified: all atoms F F F, last T T T (spacing preserved)")
print(f"✓ Output written to: POSCAR")

print("\n✅ Script completed.")
