import pyrosetta
import os 
import glob

def save_fasta_from_pdb(pdb_paths: list, output_folder: str):

    pyrosetta.init()

    for pdb_file in pdb_paths:
        pose = pyrosetta.pose_from_file(pdb_file)
        basename = os.path.splitext(os.path.basename(pdb_file))[0]
        fasta_filename = f"{basename}_chainA.fasta"

        num_chains = pose.num_chains()
    
        with open(os.path.join(output_folder, fasta_filename), 'w') as fasta_file:
            chain_id = 'A'
            fasta_file.write(f">Chain A {basename}\n")

        # Get the sequence of the specified chain
            chain_sequence = ""
            for i in range(1, pose.total_residue() + 1):
                if pose.pdb_info().chain(i) == chain_id:
                    chain_sequence += pose.residue(i).name1()
                    
            fasta_file.write(f"{chain_sequence}\n")
        
        print(f"FASTA file saved as {fasta_filename}")

input_folder = "/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/PDL1_PROTEIN_MPNN_ROSETTA_DESIGN_residue_test/step2_230712/pdbs/hsa_wt_pdl1_0_br_1_model/"

pattern = input_folder + "/*.pdb"
list_of_binders = glob.glob(pattern)
output_folder_name = "/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/fasta_wt_epitope1_newlist"

save_fasta_from_pdb(list_of_binders, input_folder)
