import pyrosetta
import os
import csv
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.core.sequence import SequenceAlignment
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector


def find_matching_directories(pdb_dir, target_dir):
    pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]
    
    target_dirs = [d for d in os.listdir(target_dir) if os.path.isdir(os.path.join(target_dir, d))]

    matches = {}
    
    for pdb_file in pdb_files:
        suffix = os.path.splitext(pdb_file)[0]
        for directory in target_dirs:
            if directory.startswith(suffix):
                matches[pdb_file] = directory
                break
    return matches

pdb_directory = './design_input2/'
target_directory = './PDL1_PROTEIN_MPNN_ROSETTA_DESIGN_backrub1st_26files/step1_230712/pdbs/'

matches = find_matching_directories(pdb_directory, target_directory)


def seq_test(matches, pdb_dir, target_dir):
    pyrosetta.init()

    percentages = []

    for pdb_file, directory in matches.items():
        pdb_file_path = os.path.join(pdb_dir, pdb_file)
        print("pdb_file_path: ", pdb_file_path)
        directory2 = os.path.join(target_dir, directory)
        step1_files = [f for f in os.listdir(directory2) if f.endswith('.pdb')]
        print("step1_files: ", step1_files)
        
        for file in step1_files:
            directory_path = os.path.join(target_dir, directory + "/" + file)
            print("directory_path: ", directory_path)

            identity_count = 0
            pose1 = pyrosetta.pose_from_pdb(pdb_file_path)
            pose2 = pyrosetta.pose_from_pdb(directory_path)

            seq1 = pose1.sequence()
            print("Sequence 1:", seq1)
            seq2 = pose2.sequence()
            print("Sequence 2:", seq2)

            HSA_length= 531
            PDL1_length = 422

            mutable_region1 = seq1[HSA_length:-PDL1_length]
            print(mutable_region1)
            mutable_region2 = seq2[HSA_length:-PDL1_length]
            print(mutable_region2)

            identity_count = sum(aa1 == aa2 for aa1, aa2 in zip(mutable_region1, mutable_region2))
            print("Identity count:", identity_count)
            sequence_identity = identity_count / min(len(mutable_region1), len(mutable_region2)) * 100
            print("Seq Similarity of mutable region: " + str(sequence_identity) + "%")

            percentages.append([pdb_file_path, directory_path, sequence_identity])

    header = ["Input", "Step1", "Sequence Similarity"]
    with open('./sequence_alignment.csv', 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for row in percentages:
            writer.writerow(row)

    return sequence_identity

seq_test(matches, pdb_directory, target_directory)
