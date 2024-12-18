from math import comb
import os
from tkinter.ttk import Separator
import pyrosetta
import pandas as pd
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.protocols.simple_filters import ContactMolecularSurfaceFilter

#Get the contact molecular surface value of protein interface

pyrosetta.init()

def interface_analysis(pdb_path: os.path):

    working_pose = pyrosetta.pose_from_pdb(pdb_path)

    interface_analyzer = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover()
    interface_analyzer.apply(working_pose)

    per_residue_data = interface_analyzer.get_all_per_residue_data()
    interface_res = per_residue_data.interface_residues

    print(interface_analyzer.get_separated_interface_energy())
    print(interface_analyzer.get_complex_energy())
    print(interface_analyzer.get_crossterm_interface_energy())
    
    res_list = []
    pose_index_list = []
    for i in range(1, len(interface_res)+ 1):
        if interface_res[i] == 1:
            pose_index_list.append(i)
            pdb_index = working_pose.pdb_info().pose2pdb(i)
            res_list.append(pdb_index)

    return res_list


def calculate_contact_surface(pose: pyrosetta.Pose):
    
    results = []
    chainA = ChainSelector('A')  
    chainB = ChainSelector('B')  

    chain_A_selector = pyrosetta.rosetta.core.select.get_residue_selector_from_subset(chainA.apply(pose))
    chain_B_selector = pyrosetta.rosetta.core.select.get_residue_selector_from_subset(chainB.apply(pose))
    
    contact_surface_filter = ContactMolecularSurfaceFilter()
    
    contact_surface_filter.selector1(chain_A_selector)
    contact_surface_filter.selector2(chain_B_selector)
    
    contact_surface_filter.apolar_target(True)  # True indicates chain B is apolar target
    
    contact_surface_value = contact_surface_filter.compute(pose)
    
    return contact_surface_value

# outputs and inputs 
input_list = ["/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/design_input2/HSA_helix~00007-105mpnn_3mpnn_0-0.9.pdb", 
              "/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/design_input2/HSA_helix~00018-200_157mpnn_0mpnn-0.9.pdb"]
output_path = "/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/contact_surface_values.csv"


#Helpful syntax 
#new_method = interface_analyzer.get_all_data()
#print(new_method.complex_total_energy)
#print(new_method.complexed_interface_score)
#print(new_method.separated_total_energy)
#print(new_method.separated_interface_score)
