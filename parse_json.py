import json 
import os
import pandas as pd

#Parse the ColabFold output to calculate and average pLDDT value for just the C terminus of designed HSA protein

def parse_json(json_file):
    with open(json_file, 'r') as file:
        data = json.load(file)  

    plddt = data['plddt']

    sum = 0
    start = len(plddt) - 125
    for i in range(start, len(plddt)):
        sum = sum + plddt[i]
    
    print(sum)
    average = sum / 125
    print(average)

    basename = os.path.basename(json_file)
    file_name = basename[:-53]
    csv_list = []
    csv_list.append([file_name, average])

    df = pd.DataFrame(csv_list)
    with open('plddt_cterm.csv', 'a') as f:
        df.to_csv(f, index = False, header=False)

    return data

file = "/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/colabfold_output/HSA_helix~00419-128mpnn_1mpnn_3-0.9_9_br_3_model_514_559_disulfide/Chain_A_HSA_helix_00419-128mpnn_1mpnn_3-0.9_9_br_3_model_514_559_disulfide_scores_rank_001_alphafold2_ptm_model_3_seed_000.json"
parse_json(file)