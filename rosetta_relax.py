import csv 
import os
import psutil
import glob
import multiprocessing
import pyrosetta
from rosetta.protocols.relax import FastRelax

pyrosetta.init("-ignore_zero_occupancy false")

def perform_relax(pdb):
    pose = pyrosetta.pose_from_pdb(pdb)

    fr = FastRelax()
    sfxn = pyrosetta.rosetta.core.scoring.get_score_function(True)

    fr.set_scorefxn(sfxn)
    fr.max_iter(100)
    fr.constrain_coords()
    fr.coord_constrain_sidechains()
    fr.apply(pose)
    final_score = sfxn(pose)
    print("this is the score:", final_score)
    return final_score

def process_directory_and_create_csv(pdb_list, output_csv):

    rows = []

    for pdb_file in pdb_list:
        score = perform_relax(pdb_file)
        rows.append([pdb_file, score])

    with open(output_csv, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Name", "Total Score (REU)"])
        csv_writer.writerows(rows)

    print(f"Results written to {output_csv}")

input_folder = "/ifs/scratch/home/eca2133/final_part/parrots_analysis/chain_A_WT_designs"
output_folder = "/ifs/scratch/home/eca2133/final_part/parrots_analysis/chain_A_WT_designs/WT_HSA_only_scores.csv"

pattern = input_folder + "/*.pdb"
list_of_binders = glob.glob(pattern)

process_directory_and_create_csv(list_of_binders, output_folder)