#!/usr/bin/env python3
"""
csv2xyz.py  –  convert Vicsek csv snapshots → OVITO-ready .xyz
All particles are written as species “A”.
Usage:
    python csv2xyz.py
Output:
    config_data/trial_*/xyz/config_12345.xyz
"""

import csv, pathlib, sys

SRC_DIR = pathlib.Path('data/config_data')          # root that contains trial_*/
XYZ_DIR = 'xyz'                                # sub-folder that will hold .xyz files
SPECIES = 'A'                                  # same species for every particle

def csv_to_xyz(csv_path, xyz_path):
    """Read x,y,vx,vy csv and write minimal xyz file."""
    with csv_path.open(newline='') as fc, xyz_path.open('w') as fx:
        rows = list(csv.DictReader(fc))          # first line is header
        fx.write(f'{len(rows)}\n')
        fx.write('Properties=id:I:1:pos:R:2:velo:R:2:species:S:1\n')
        for i, r in enumerate(rows):
            x, y   = float(r['x']), float(r['y'])
            vx, vy = float(r['vx']), float(r['vy'])
            fx.write(f'{i} {x} {y} {vx} {vy} {SPECIES}\n')

def main():
    if not SRC_DIR.is_dir():
        sys.exit(f'Folder {SRC_DIR} not found – run this script from the project root.')

    trials = [d for d in SRC_DIR.iterdir() if d.is_dir() and d.name.startswith('trial_')]
    if not trials:
        sys.exit('No trial_*/ directories inside config_data/')

    for trial in trials:
        csv_folder = trial
        xyz_folder = trial / XYZ_DIR
        xyz_folder.mkdir(exist_ok=True)

        csv_files = sorted(csv_folder.glob('config_*.csv'))
        if not csv_files:
            print(f'No csv files in {trial}, skipping')
            continue

        for csv_file in csv_files:
            num = csv_file.stem.split('_')[1]          # config_123.csv → 123
            xyz_file = xyz_folder / f'config_{num}.xyz'
            csv_to_xyz(csv_file, xyz_file)

        print(f'converted {len(csv_files)} snapshots → {xyz_folder}')

if __name__ == '__main__':
    main()