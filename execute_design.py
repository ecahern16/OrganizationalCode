from datetime import datetime
import multiprocessing
import glob
import psutil
import csv
import shutil
import os
import pandas as pd
from utilities import *
import energy_methods_original

input_folder = "/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/WT_input"
directory_name = "PDL1_PROTEIN_MPNN_ROSETTA_DESIGN_residue_new_br100"
csv_suffix = "_WT_HSA_residue_new_br100.csv"

def execute_pipeline(pdb_paths: list, protein_mpnn_1_path: os.path, protein_mpnn_2_path: os.path, backrub_path: os.path, n_thread: int, # type: ignore
                     protein_mpnn_flag: bool, design_flag: bool, relaxed_flag: bool, n_trials: int, n_pass: int):

    #added this for threads
    cpu_id:int = n_thread
    proc = psutil.Process()
    aff = proc.cpu_affinity()
    proc.cpu_affinity([cpu_id,cpu_id+1])
    aff = proc.cpu_affinity()

    backrub_designs = energy_methods_original.perform_chainA_backrub(pdb_paths, backrub_path, iterations=100)

    # -------------------- First Design Step --------------------
    step_1_passed = energy_methods_original.design_round_for_WT(backrub_designs, protein_mpnn_1_path,
                                                protein_mpnn_flag, design_flag, relaxed_flag, n_thread, n_trials,
                                                n_pass)
    print('This is step_1_passed:', step_1_passed)
    step_1_list = []
    count = 0
    for file in step_1_passed:   
        print("this is the file in step_1_passed:", file)
        current_dict = {"name": file}
        current_dict= energy_methods_original.get_dgDSASA_dict(file, current_dict)
        step_1_list.append(current_dict)
        print("/nCurrent dictionary in step1:", current_dict)
        count = count + 1
        print("DEBUG:", count)
    df = pd.DataFrame(step_1_list)
    with open('scoring/step_1_design' + csv_suffix, 'a') as f:
        df.to_csv(f, index = False, header=False)

    # -------------------- Second Design Step --------------------

    step_2_passed = energy_methods_original.design_round_for_WT(step_1_passed, protein_mpnn_2_path,
                                                protein_mpnn_flag, design_flag, True, n_thread, n_trials, n_pass)
    print('This is step_2_passed:', step_2_passed)
    step_2_list = [] 
    count = 0
    for file in step_2_passed:
        print("this is the file in step_2_passed:", file)
        current_dict = {"name": file}
        current_dict= energy_methods_original.get_dgDSASA_dict(file, current_dict)
        step_2_list.append(current_dict)
        print("/nCurrent dictionary in step2:", current_dict)
        count = count + 1
        print("DEBUG:", count)
    df2 = pd.DataFrame(step_2_list)
    with open('scoring/step_2_design' + csv_suffix, 'a') as f:
        df2.to_csv(f, index = False, header=False)

    return step_1_passed
    
pyrosetta.init(
    extra_options="-mute basic -mute core -mute protocols -relax:default_repeats 1  -ignore_zero_occupancy false")

#Make csv header
make_dir("scoring")
header_flag1 = True
with open('scoring/step_1_design' + csv_suffix, 'a') as f:
        if header_flag1:
            writer = csv.writer(f)
            writer.writerow(energy_methods_original.get_dgDSASA_keys()) # write the header
            header_flag1 = False

header_flag2 = True
with open('scoring/step_2_design' + csv_suffix, 'a') as f:
    if header_flag2:
        writer = csv.writer(f)
        writer.writerow(energy_methods_original.get_dgDSASA_keys()) # write the header
        header_flag2 = False


pattern = input_folder + "/*.pdb"
list_of_binders = glob.glob(pattern)
print(list_of_binders)
num_threads = 1
chunks = make_chunks(list_of_binders,num_threads)
threads = []

make_dir(directory_name)
make_dir("temp_files")
for thread_num in range(num_threads):
    current_thread = multiprocessing.Process(target=execute_pipeline, args=(
        chunks[thread_num],                                                 # PDB Files
        directory_name + "/step1_230712",                                   # Outputs for STEP 1 design
        directory_name + "/step2_230712",                                   # Outputs for STEP 2 design
        directory_name + "/backrub_outputs",                                # Outputs for Backrub
        thread_num,                                                         # Threadnum
        True,                                                               # ProteinMPNN?
        True,                                                               # Rosetta Design?
        False,                                                              # Are inputs relaxed?
        8,                                                                  # Num of models to generate
        4))                                                                 # Num of models to pass
    threads.append(current_thread)
    current_thread.start()

for t in threads:
    t.join()
