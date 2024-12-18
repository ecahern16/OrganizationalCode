#!/bin/bash
#SBATCH --partition=glab
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem=20G
##SBATCH --gres=gpu:1
##export CUDA_VISIBLE_iDEVICES=0

conda activate newpyrosetta
pdb_file="/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/PDL1_PROTEIN_MPNN_ROSETTA_DESIGN_residue_test/step2_230712/pdbs/hsa_wt_pdl1_0_br_1_model/hsa_wt_pdl1_0_br_1_model_1_model.pdb"
new_dir="hsa_wt_pdl1_0_br_1_model_1_model"

mkdir $new_dir

python ./get_interface_res.py $pdb_file $new_dir

interface_energy=/ifs/data/glab/rosetta/main/source/bin/interface_energy.default.linuxgccrelease
database=/ifs/scratch/home/bs3281/apps/rosetta/database
pep_face="/ifs/scratch/home/eca2133/final_part/res_energy/${new_dir}/prot_faceA"
prot_face="/ifs/scratch/home/eca2133/final_part/res_energy/${new_dir}/prot_faceB"

pep_face_basename=$(basename "$pep_face")
prot_face_basename=$(basename "$prot_face")
echo "========================= $file_name, $pep_face_basename, $prot_face_basename ========================="
$interface_energy -s "$pdb_file" -face1 "$pep_face" -face2 "$prot_face" -score:hbond_bb_per_residue_energy

mv "slurm-${SLURM_JOB_ID}.out" "./${new_dir}"

input_path="/ifs/scratch/home/eca2133/final_part/res_energy/${new_dir}/slurm-${SLURM_JOB_ID}.out"

python parse_outfile.py $input_path

mv "/ifs/scratch/home/eca2133/final_part/res_energy/${new_dir}_interface_energies.csv" "./${new_dir}"