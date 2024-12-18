import os
import shutil
import csv

# Define the paths
csv_file_input = './scoring/top_parrots.csv'  # Path to your CSV file
destination_folder = './post_PARROTS_3022/'   # Destination folder where files will be copied

def move_files(csv_file):
    os.makedirs(destination_folder, exist_ok=True)

    with open(csv_file, mode='r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row

        for row in reader:
            # Extract the path from the first column
            file_path = row[0]
            file_name = os.path.basename(file_path)
            destination_path = os.path.join(destination_folder, file_name)
            shutil.copy(file_path, destination_path)
            print(f"Copied {file_name} to {destination_folder}")

#FILTER STEP
def check_threshold(file_path, column1_name, column2_name, destination_dir):
    os.makedirs(destination_dir, exist_ok=True)
    with open(file_path, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file) #need to use DictReader for it to work
        next(csv_reader)  # Skip header if exists
        for row in csv_reader:
            dg_sep = float(row[column1_name])
            h_bond = float(row[column2_name])
            if dg_sep < -2.8 and h_bond > 0.07:
                source_file = row['name'] 
                filename = os.path.basename(source_file)
                destination_file = os.path.join(destination_dir, filename)
                shutil.move(source_file, destination_file)
                print("Moved {} to {}".format(source_file, destination_file))


column_name1 = 'separated_interface/dSASAx100'
column_name2 = 'hbond_energy/separated_interface'
destination = '/ifs/scratch/home/eca2133/final_part/PARROTS-Pipeline-main/post_PARROTS_3022'
csv_suffix = "_scaleup_3900files.csv"
check_threshold('scoring/step_1_design' + csv_suffix, column_name1, column_name2, destination)