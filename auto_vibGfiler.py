#!/usr/bin/env python3

import os
import shutil
from pathlib import Path

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
        lines = f.readlines()

    new_lines = []
    flags = {'PREC': False, 'IBRION': False, 'ISMEAR': False}

    for line in lines:
        stripped = line.strip()
        if stripped.startswith('PREC'):
            new_lines.append("PREC = Accurate\n")
            flags['PREC'] = True
        elif stripped.startswith('IBRION'):
            new_lines.append("IBRION = 5\n")
            flags['IBRION'] = True
        elif stripped.startswith('ISMEAR'):
            new_lines.append("ISMEAR = 0\n")
            flags['ISMEAR'] = True
        else:
            new_lines.append(line)

    if not flags['PREC']:
        new_lines.append("PREC = Accurate\n")
    if not flags['IBRION']:
        new_lines.append("IBRION = 5\n")
    if not flags['ISMEAR']:
        new_lines.append("ISMEAR = 0\n")

    with open(incar_path, 'w') as f:
        f.writelines(new_lines)

    print("✓ INCAR updated: PREC=Accurate, IBRION=5, ISMEAR=0")
else:
    print("⚠️  INCAR not found, skipped editing.")

# --- Step 4: Edit CONTCAR Constraints ---
contcar_path = cwd / 'CONTCAR'
if not contcar_path.exists():
    print("⚠️  CONTCAR not found. Skipping POSCAR creation.")
    exit(0)

print("\n🛠️  Editing CONTCAR...")
shutil.copy(contcar_path, contcar_path.with_suffix(".backup"))

with open(contcar_path, 'r') as f:
    lines = f.readlines()

# --- Identify coordinate start ---
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

# --- Edit constraint lines ---
for i in range(coord_start, coord_end):
    cols = lines[i].split()
    if len(cols) >= 3:
        xyz = " ".join(cols[:3])
        if i == coord_end - 1:
            lines[i] = f"{xyz} T T T\n"
        else:
            lines[i] = f"{xyz} F F F\n"

# --- Step 5: Write POSCAR ---
poscar_path = cwd / 'POSCAR'
with open(poscar_path, 'w') as f:
    f.writelines(lines)

print("✓ Constraints applied: all atoms F F F, last T T T")
print(f"✓ Renamed output: {contcar_path.name} ➜ POSCAR")

print("\n✅ Script completed.")
