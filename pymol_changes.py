import pyrosetta
import pymol
from pymol import cmd
import os
import glob

#Align designed HSA structure in complex with PD-L1 to HSA in complex with FcRn 
#Remove PD-L1 after aligned
#Remove WT HSA that was in complex with FcRn
#Save new file of designed HSA in complex with FcRn 

pyrosetta.init()

def get_new_complex(pdb_paths, fcrn_path):   
        for pdb_file in pdb_paths:
                basename = os.path.splitext(os.path.basename(pdb_file))[0]  
                print("Loading HSA_structure...")
                cmd.load(pdb_file, "HSA_structure")
                print("Loading antigen structure...")
                cmd.load(fcrn_path, "FCRN")

                cmd.select("fcrn_A", "FCRN and chain A")
                print("Selected HSA in FcRn complex.")
                cmd.align("HSA_structure", "fcrn_A")
                print("Alignment completed.")
                cmd.remove("FCRN and chain A")
                cmd.remove("HSA_structure and chain B")
                cmd.alter("HSA_structure and chain B", "chain='C'")
                cmd.create("HSA_aligned_structure", "FCRN or HSA_structure") 
                print("Merged object created.")
                cmd.save(f"FcRn_complex_{basename}.pdb", "HSA_aligned_structure")
                print("Structure saved.")
                cmd.remove("all")
                cmd.delete("all")

#Extract HSA designed protein from complex with target

def get_chainA(pdb_paths, output_path):
        for pdb_file in pdb_paths:
                basename = os.path.splitext(os.path.basename(pdb_file))[0]  
                print("Loading HSA_structure...")
                cmd.load(pdb_file, "HSA_structure")
                cmd.select("antigen", "HSA_structure and chain B")
                print("Selected antigen in HSA complex.")
                cmd.remove("antigen")
                cmd.create("HSA_only_structure", "HSA_structure and chain A")
                cmd.save(f"{output_path}/HSA_only_{basename}.pdb", "HSA_only_structure")
                cmd.remove("all")
                cmd.delete("all")

receptor_path = "/ifs/scratch/home/eca2133/4k71_2chains.pdb"

input_folder = "/ifs/scratch/home/eca2133/final_part/parrots_analysis/to_be_relaxed"
output_folder = "/ifs/scratch/home/eca2133/final_part/parrots_analysis/chain_A_WT_designs"

pattern = input_folder + "/*.pdb"
list_of_binders = glob.glob(os.path.join(input_folder, "hsa*.pdb"))

get_chainA(list_of_binders, output_folder)