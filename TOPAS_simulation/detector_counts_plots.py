import os
import re
import argparse
import matplotlib.pyplot as plt
from pydicom import dcmread
import h5py

plt.close("all")

plt.rcParams.update({
    'font.size': 14,            
    'axes.titlesize': 16,      
    'axes.labelsize': 14,       
    'xtick.labelsize': 13,      
    'ytick.labelsize': 13,     
    'legend.fontsize': 12
})

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--TopasFolderPath", nargs='?', help="Directory containing Beams folders. If you followed the workflow, it should be name as your/Directory/Protons_x.", default=".")
parser.add_argument("-r", "--RTDoseFolderPath", nargs='?', help="Directory of the DICOM Patient.", default=".")
parser.add_argument("-o", "--OutputPath", nargs='?', help="Output directory where the plots will be saved.", default=".")
parser.add_argument("-b", "--BeamNumber", type=int, nargs='+', help="Beam list", default=[1])
parser.add_argument("-f", "--SaveFigures", help="Set to True in order to save the generated Figures", default=False)
parser.add_argument("-p", "--PatientIdx", type=int, help="Patient index for the TROTS Protons Dataset (from 1 to 20).", default=1) 
parser.add_argument("-m", "--MatFilePath", help="Directory of the original .mat file named beamlist.", default=".") 

args = parser.parse_args()

SaveFigures = args.SaveFigures if type(args.SaveFigures) == bool else args.SaveFigures == "True"

#First lets link control points to energy directly from the .mat file
mat_path = os.path.join(args.MatFilePath, "beamlist.mat")
if not os.path.exists(mat_path):
    raise FileNotFoundError(f"Beamlist.mat file not found at {args.MatFilePath}")

with h5py.File(mat_path, 'r') as f:
    patient_ref = f['BeamList'][args.PatientIdx - 1, 0]
    patient_matrix = f[patient_ref][:].T

#This are the scorers names for the detector option in the code generate_scripts.py
file_types = ['nocs', 'nzps', 'nzms', 'ocs', 'zps', 'zms']
analysis_data = {ftype: {} for ftype in file_types}

for beamnumber in args.BeamNumber:
    print(f"\nPROCESSING BEAM: {beamnumber}")
    
    csv_path = os.path.join(args.TopasFolderPath, f'Beam_{beamnumber}')
    plan_path = os.path.join(args.RTDoseFolderPath, "rtplan.dcm")

    if not os.path.exists(csv_path):
        print(f"Warning: Folder {csv_path} does not exist. Skipping.")
        continue

    rtplan = dcmread(plan_path)
    beam_idx = int(beamnumber - 1)
    beam = rtplan.IonBeamSequence[beam_idx]
    beam_data_mat = patient_matrix[patient_matrix[:, 0] == beamnumber]

    cp_to_energy = {}
    current_energy = float(beam_data_mat[0, 1]) 
    mat_row_idx = 0

    for i, cp_seq in enumerate(beam.IonControlPointSequence):
        if i % 2 != 0:
            continue       
        if hasattr(cp_seq, "NumberOfScanSpotPositions") and int(cp_seq.NumberOfScanSpotPositions) > 0:
            num_spots = int(cp_seq.NumberOfScanSpotPositions)
            
            if mat_row_idx < len(beam_data_mat):
                current_energy = float(beam_data_mat[mat_row_idx, 1])
                mat_row_idx += num_spots
        
        cp_to_energy[i] = current_energy
    
    for file in sorted(os.listdir(csv_path)):
        match = re.search(r'D1(n?(?:ocs|zps|zms))_CP(\d+)\.csv', file)
        if match:
            file_type = match.group(1) 
            cp = int(match.group(2))    
            
            energy = cp_to_energy.get(cp, None)
            if energy is None:
                continue
            
            filepath = os.path.join(csv_path, file)
            
            try:
                with open(filepath, 'r') as f:
                    lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
                
                if not lines:
                    counts = 0
                elif len(lines) == 1:
                    try:
                        counts = float(lines[0].split()[0])
                    except ValueError:
                        counts = 0
                else:
                    
                    first_line_tokens = lines[0].replace(',', ' ').split()
                    if len(first_line_tokens) == 1:
                        counts = sum(float(l.replace(',', ' ').split()[0]) for l in lines)
                    else:
                        counts = len(lines)
                        
            except Exception as e:
                print(f"Error reading file {file}: {e}")
                counts = 0
            
            if energy in analysis_data[file_type]:
                analysis_data[file_type][energy] += counts
            else:
                analysis_data[file_type][energy] = counts
            
            #lets now include a global count for all the surfaces, but differentiating between proton and non-proton ancestors
            if file_type.startswith('n'):
                global_key = 'global_n'
            else:
                global_key = 'global'
            if global_key not in analysis_data:
                analysis_data[global_key] = {}
            analysis_data[global_key][energy] = analysis_data[global_key].get(energy, 0) + counts


surface_labels = {
    'ocs': 'Outer Curved Surface (OCS)',
    'zps': 'Z Plus Surface (ZPS)',
    'zms': 'Z Minus Surface (ZMS)'
}

#scatter plots for each surface with proton/non-proton differentiation
for surface_key in ['ocs', 'zps', 'zms']:
    surface_text = surface_labels.get(surface_key, surface_key.upper())
    
    plt.figure(figsize=(10, 6))
    
    #Plot proton ancestors for this surface
    proton_file_type = f'n{surface_key}'
    if proton_file_type in analysis_data and analysis_data[proton_file_type]:
        proton_energies = sorted(analysis_data[proton_file_type].keys())
        proton_counts = [analysis_data[proton_file_type][e] for e in proton_energies]
        plt.scatter(proton_energies, proton_counts, color='cyan',edgecolor='darkcyan', alpha=0.7, label='Proton Ancestor', linewidth=1.5,s=85)
    
    #Plot non-proton ancestors
    if surface_key in analysis_data and analysis_data[surface_key]:
        non_proton_energies = sorted(analysis_data[surface_key].keys())
        non_proton_counts = [analysis_data[surface_key][e] for e in non_proton_energies]
        plt.scatter(non_proton_energies, non_proton_counts,color='orchid', edgecolor='darkorchid', alpha=0.7, label='Non-Proton Ancestor', linewidth=1.5,marker="X",s=70)
    
    plt.title(f"Surface: {surface_text}", fontweight='bold')
    plt.xlabel("Energy [MeV]")
    plt.ylabel("$\gamma$ Counts")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc='best')
    plt.tight_layout()
    
    if SaveFigures:
        output_filename = f"counts_vs_energy_{surface_key}.pdf"
        output_file_path = os.path.join(args.OutputPath, output_filename)
        plt.savefig(output_file_path, dpi=300)

#Create scatter plot for total sum with proton/non-proton differentiation
plt.figure(figsize=(10, 6))

#TOTAL proton ancestors 
if 'global_n' in analysis_data and analysis_data['global_n']:
    global_n_energies = sorted(analysis_data['global_n'].keys())
    global_n_counts = [analysis_data['global_n'][e] for e in global_n_energies]
    plt.scatter(global_n_energies, global_n_counts, color='cyan', edgecolor='darkcyan', alpha=0.7, label='Proton Ancestor', linewidth=1.5,s=85)

#TOTAL non-proton ancestors
if 'global' in analysis_data and analysis_data['global']:
    global_energies = sorted(analysis_data['global'].keys())
    global_counts = [analysis_data['global'][e] for e in global_energies]
    plt.scatter(global_energies, global_counts, color='orchid', edgecolor='darkorchid', alpha=0.7, label='Non-Proton Ancestor', linewidth=1.5,marker="X",s=70)

plt.title("Total Counts Across All Surfaces", fontweight='bold')
plt.xlabel("Energy [MeV]")
plt.ylabel("$\gamma$ Counts")
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='best')
plt.tight_layout()

if SaveFigures:
    output_filename = "counts_vs_energy_total.pdf"
    output_file_path = os.path.join(args.OutputPath, output_filename)
    plt.savefig(output_file_path, dpi=300)

plt.show()
    
