import pyrosetta
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
import os 
import shutil
import glob
from utilities import *
import multiprocessing
import psutil

pyrosetta.init()

def clash_test(pdb_list: list, firstChain: str, secondChain: str, clash_distance: float, n_thread: int):
    #added this for threads
    cpu_id:int = n_thread
    proc = psutil.Process()
    aff = proc.cpu_affinity()
    proc.cpu_affinity([cpu_id,cpu_id+1])
    aff = proc.cpu_affinity()

    pdbs_without_clashes = []

    for pdb_path in pdb_list:
        clashes = []

        file = pdb_path.split("/")[-1][0:-4]
        print("PDB PATH:", pdb_path)
        print("FILE:", file)

        working_pose: pyrosetta.Pose = pyrosetta.pose_from_file(pdb_path)
        if len(working_pose.sequence()) < 500:
            continue

        start_residue_index = working_pose.pdb_info().pdb2pose('A', 300)
        last_residue_index = working_pose.chain_end(1)

        if start_residue_index == 0 or start_residue_index > working_pose.total_residue():
            print(f"Invalid start residue index: {start_residue_index}")
            continue
        
        if last_residue_index == 0 or last_residue_index > working_pose.total_residue():
            print(f"Invalid last residue index: {last_residue_index}")
            continue

        residue_index_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        residue_index_selector.set_index_range(start_residue_index, last_residue_index)
        print("Res selector:", residue_index_selector.selection_positions(working_pose))

        selected_positions = residue_index_selector.apply(working_pose)
        selected_residues_list = [i for i in range(1, 615) if selected_positions[i]== 1]
        print(selected_residues_list)

        nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
        nbr_selector.set_focus_selector(residue_index_selector)
        nbr_selector.set_distance(clash_distance)
        nbr_selected = nbr_selector.apply(working_pose)
        selected_nbr_list = [i for i in range(1, working_pose.size()+ 1) if nbr_selected[i] == 1]
        print(selected_nbr_list)

        # Create a new list that contains elements of nbr not in residue selector
        result_list = [item for item in selected_nbr_list if item not in selected_residues_list]
        print("Result list:", result_list)
        # Create a new list that is only chain B aka greater index than end of HSA last_residue_index
        result_list2 = [item for item in result_list if item > last_residue_index and item <= working_pose.total_residue()]
        print("Result list 2:", result_list2)

        for i in range(0, len(selected_residues_list)):
            res1 = working_pose.residue(selected_residues_list[i])
                       
            # Iterate through residues in the second chain
            for j in range(0, len(result_list2)):
                res2 = working_pose.residue(result_list2[j])
                
                # Iterate through atoms in the first residue
                for atom1 in range(1, res1.natoms() + 1):
                    atom1_xyz = res1.atom(atom1).xyz()
                    
                    # Iterate through atoms in the second residue
                    for atom2 in range(1, res2.natoms() + 1):
                        atom2_xyz = res2.atom(atom2).xyz()
                        distance = atom1_xyz.distance(atom2_xyz)
                        
                        # Check for steric clash
                        if distance < 0.6 or distance == 0:
                            selected_pdb_id = working_pose.pdb_info().pose2pdb(selected_residues_list[i])
                            result_pdb_id = working_pose.pdb_info().pose2pdb(result_list2[j])             
                            clashes.append([selected_pdb_id, result_pdb_id, distance])
        print("Clashes:", clashes)
        if len(clashes) == 0:
            pdbs_without_clashes.append(pdb_path)
    print(pdbs_without_clashes)
    print("Length:", len(pdbs_without_clashes))

    for pdb_ex in pdbs_without_clashes:
        filename = os.path.basename(pdb_ex)
        destination_dir = "/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/post_clash_5000"
        destination_file = os.path.join(destination_dir, filename)
        print(f"Attempting to copy {filename} to {destination_dir}")
    
        try:
            shutil.copy(pdb_ex, destination_file)
            print(f"Copied {filename} successfully.")
        except Exception as e:
            print(f"Failed to copy {filename}: {e}")
    
    return pdbs_without_clashes

input_folder = "/ifs/scratch/home/eca2133/make_loops/connected2/HSA_pdbs"
pattern = input_folder + "/*.pdb"
list_of_binders = glob.glob(pattern)

pdb_list = ["/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/clash_problem/HSA_helix~00012-47mpnn_1mpnn_1-0.9.pdb", 
            "/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/clash_problem/HSA_helix~00064-196mpnn_3mpnn_0-0.9.pdb", 
            "/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/clash_problem/HSA_helix~00058-50mpnn_3mpnn_2-0.9.pdb", 
            "/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/clash_problem/HSA_helix~00074-156mpnn_1mpnn_0-0.9.pdb",
            "/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/design_input2/HSA_helix~00613-32mpnn_1mpnn_2-0.9.pdb", 
            "/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/design_input2/HSA_helix~00241-90mpnn_2mpnn_2-0.9.pdb"]


num_threads = 22
chunks = make_chunks(list_of_binders, num_threads)
threads = []

for thread_num in range(num_threads):
    current_thread = multiprocessing.Process(target=clash_test, args=(
        chunks[thread_num],                                                 # PDB Files
        "A",         
        "B",          
        10.0,
        thread_num                                                                                                                                                                                                                                         # ProteinMPNN?
        ))                                                                
    threads.append(current_thread)
    current_thread.start()

for t in threads:
    t.join()