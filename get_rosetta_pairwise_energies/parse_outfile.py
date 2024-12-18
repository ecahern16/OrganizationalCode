import numpy as np
import pandas as pd
import re
import os
import sys

#interface_output = '/ifs/scratch/home/eca2133/final_part/res_energy/WT_loop_check/slurm-6288.out'
#pdb_file = "hsa_wt_pdl1_0_br_1_model_4_model"
# Define a regex pattern to extract the needed parts
def parse(interface_output):
    pattern = re.compile(r'([A-Za-z]+)(\d+)_([A-Za-z])\((\d+)\) --- ([A-Za-z]+)(\d+)_([A-Za-z])\((\d+)\): ([\d\.\-e]+)')
    with open(interface_output, 'r') as file:
        lines = file.readlines()
    entry_lines = []
    current_entry = []
    recording = False
    for line in lines:
        if line.startswith('========================='):
            if current_entry:
                entry_lines.append(current_entry)
            current_entry = [line]
            recording = False
        elif line.startswith('##### PAIRWISE SHORT-RANGE ENERGIES #####'):
            recording = True
        if recording:
            current_entry.append(line)
    if current_entry:
        entry_lines.append(current_entry)
    for entry in entry_lines:
        #pdb_file = entry[0].split(' ')[1].replace(',', '')
        pdb_file = interface_output.split('/')[-2] 
        interface_energies = []
        for line in entry:
            if line.startswith('##### PAIRWISE SHORT-RANGE ENERGIES #####'):
                recording = True
                continue
            if recording:
                if line.startswith('0') or line.strip().endswith('0'):
                    continue
                interface_energies.append(line.strip())
        data = []
        for energy_line in interface_energies:
            match = pattern.search(energy_line)
            if match:
                #pep_resn, _, pep_chain, pep_resi, prot_resn, _, prot_chain, prot_resi, energy = match.groups()
                #data.append([pep_resn, pep_resi, pep_chain, prot_resn, prot_resi, prot_chain, float(energy)])
                pep_resn, pep_resi_prefix, pep_chain, pep_resi, prot_resn, prot_resi_prefix, prot_chain, prot_resi, energy = match.groups()

                pep_full_res = f"{pep_resn}{pep_resi_prefix}"  # e.g., GLU501
                prot_full_res = f"{prot_resn}{prot_resi_prefix}"  # e.g., LYS115

                data.append([pep_resn, pep_resi_prefix, pep_chain, prot_resn, prot_resi_prefix, prot_chain, float(energy)])

        # Create DataFrame
        df = pd.DataFrame(data, columns=['peptide_resn', 'peptide_resi', 'peptide_chain', 'protein_resn', 'protein_resi', 'protein_chain', 'pairwise_energy'])
        df = df.sort_values(by = 'pairwise_energy')
        # Calculate the total interface energy
        total_energy = df['pairwise_energy'].sum()
        total_row = pd.DataFrame([['TOTAL', '', '', '', '', '', total_energy]], columns=df.columns)
        # Append the total row to the DataFrame
        df = pd.concat([df, total_row], ignore_index=True)
        # Add a blank row (NaN) to the DataFrame
        blank_row = pd.DataFrame([[np.nan] * len(df.columns)], columns=df.columns)
        df = pd.concat([df, blank_row], ignore_index=True)
        # Add the interface_output path row to the DataFrame
        path_row = pd.DataFrame([[interface_output] + [''] * (len(df.columns) - 1)], columns=df.columns)
        df = pd.concat([df, path_row], ignore_index=True)
        # Save to a CSV file
        output_file = os.path.join(os.getcwd(), f'{pdb_file}_interface_energies.csv')
        df.to_csv(output_file, index=False)
        print(f'Saved: {output_file}')
#parse(interface_output, pdb_file)

def main():
    scriptname = sys.argv[0]
    path = sys.argv[1]

    print("Script name:", scriptname)
    print("Path:", path)
    return path
    
if __name__ == "__main__":
    main()

path = main()
my_result = parse(path)