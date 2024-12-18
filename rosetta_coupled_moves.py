import os
import os.path
import subprocess
from utilities import *
import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta import init
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from rosetta.core.select.movemap import *
from pyrosetta.rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation

def perform_coupledmoves(pdb_list: list, coupledmoves_output: os.path, iterations=50):

    init()
    make_dir(coupledmoves_output)
    final: list = []
    index = 0

    for pdb_path in pdb_list:
        file = pdb_path.split("/")[-1][0:-4]
        print("PDB PATH:", pdb_path)
        print("FILE:", file)

        working_pose: pyrosetta.Pose = pyrosetta.pose_from_file(pdb_path)
        start_residue_index = working_pose.pdb_info().pdb2pose('A', 500)
        print("Index of 500:",  working_pose.pdb_info().pdb2pose('A', 500))
        last_residue_index = working_pose.chain_end(1)
        print("Last res index:", last_residue_index)

        residue_index_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        residue_index_selector.set_index_range(start_residue_index, last_residue_index)
        print("Res selector:", residue_index_selector.selection_positions(working_pose))

        nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
        nbr_selector.set_focus_selector(residue_index_selector)
        nbr_selector.set_distance(10.0)
        
        tf = TaskFactory()
        tf.push_back(operation.InitializeFromCommandline())
        prevent_repacking_rlt = operation.PreventRepackingRLT()
        prevent_subset_repacking = operation.OperateOnResidueSubset(prevent_repacking_rlt, residue_index_selector, True )
        tf.push_back(prevent_subset_repacking)

        cm_mover = pyrosetta.rosetta.protocols.coupled_moves.CoupledMovesProtocol()
        for i in range(iterations):
            cm_mover.set_main_task_factory(tf)
            cm_mover.apply(working_pose)

        output_path = os.path.join(coupledmoves_output, str(file) + "_" + str(index) + "restr_CM.pdb")
        final.append(output_path)
        print("OUTPUT PATH:", output_path)
        working_pose.dump_pdb(output_path)
        index += 1
    
    print("FINAL:", final)
    return final

input_list = ["test_designs/HSA_helix~00069-6mpnn_2mpnn_0-0.9.pdb" , "test_designs/HSA_helix~00073-99mpnn_3mpnn_1-0.9.pdb"]
backrub_designs = perform_coupledmoves(input_list, "./coupled_moved_outputs", iterations=50)
