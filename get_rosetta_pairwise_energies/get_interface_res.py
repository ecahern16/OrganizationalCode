from math import comb
import os
import sys
from tkinter.ttk import Separator
import pyrosetta
import pandas as pd
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.protocols.simple_filters import ContactMolecularSurfaceFilter

pyrosetta.init()
def interface_analysis(pdb_path, dir_name):

    working_pose = pyrosetta.pose_from_pdb(pdb_path)

    interface_analyzer = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover()
    interface_analyzer.apply(working_pose)

    per_residue_data = interface_analyzer.get_all_per_residue_data()
    #print(interface_analyzer.get_interface_delta_sasa())
    interface_res = per_residue_data.interface_residues
    
    res_list = []
    pose_index_list = []
    for i in range(1, len(interface_res)+ 1):
        if interface_res[i] == 1:
            pose_index_list.append(i)
            pdb_chain = working_pose.pdb_info().chain(i)
            pdb_index = working_pose.pdb_info().number(i)
            res_list.append(str(pdb_chain)+ " " + str(pdb_index) + " _")
    
    print(res_list)

    list_A = [item for item in res_list if item.startswith('A')]
    print(list_A)
    list_B = [item for item in res_list if item.startswith('B')]
    print(list_B)

    with open(f'./{dir_name}/prot_faceA', 'w') as file_A:
        for item in list_A:
            file_A.write(f"{item}\n")

    with open(f'./{dir_name}/prot_faceB', 'w') as file_B:
        for item in list_B:
            file_B.write(f"{item}\n")

    return res_list

def main():
    scriptname = sys.argv[0]
    pdb = sys.argv[1]
    dir = sys.argv[2]

    print("Script name:", scriptname)
    print("pdb:", pdb)
    return pdb, dir
    
if __name__ == "__main__":
    main()

pdb = main()
my_result = interface_analysis(pdb[0], pdb[1])