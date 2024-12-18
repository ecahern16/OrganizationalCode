import os
import subprocess
import pyrosetta
from rosetta.protocols.relax import FastRelax
import csv
import multiprocessing
import pandas as pd

#Computational alanine scan for the C terminal end of HSA protein

pyrosetta.init()

def perform_relax(pdb):
    pose = pyrosetta.pose_from_pdb(pdb)
    fr = FastRelax()
    sfxn = pyrosetta.get_score_function(True)

    fr.set_scorefxn(sfxn)
    fr.max_iter(100)
    fr.constrain_coords()
    fr.coord_constrain_sidechains()
    fr.apply(pose)
    final_score = sfxn(pose)
    print("this is the score:", final_score)
    return pose

def mutate_and_score(pdb: str, residue: int):
    mutater = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
    pose = pyrosetta.pose_from_pdb(pdb)
    sfxn = pyrosetta.get_score_function(True)
    print("o seq:", pose.sequence())

    mutater.set_target(residue)
    mutater.set_res_name('ALA')
    mutater.apply(pose)
    print("m seq:", pose.sequence())

    perform_relax(pose)
    final_score = sfxn(pose)
    print("this is the score:", final_score)

    return [pose.pdb_info().pose2pdb(residue), final_score]
    

def alanine_scanning(pdb: str):
    residues = range(380, 578)

    with multiprocessing.Pool() as pool:
        residue_list = pool.starmap(mutate_and_score, [(pdb, i) for i in residues])
    
    header = ["Mutated Residue", "Score (REU)"]
    with open("/ifs/scratch/home/eca2133/final_part/parrots_analysis/alanine_scanning_HSA.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(header)  # Write the header
        writer.writerows(residue_list) 
    return residue_list

