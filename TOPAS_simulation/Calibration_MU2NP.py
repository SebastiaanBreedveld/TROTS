#---------------------------------------------------------------------------------
#Before running this code: read README.md
#---------------------------------------------------------------------------------

import os
import re
import argparse
import numpy as np
import nrrd
import matplotlib.pyplot as plt
from pydicom import dcmread
from lmfit import minimize, Parameters, report_fit
from numpy.polynomial.polynomial import polyval
import pandas as pd
import h5py
#plt.close("all")

plt.rcParams.update({
    'font.size': 14,            
    'axes.titlesize': 16,      
    'axes.labelsize': 14,       
    'xtick.labelsize': 13,      
    'ytick.labelsize': 13,     
    'legend.fontsize': 12})

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--TopasFolderPath", nargs='?', help="Directory containing TOPAS rtdoses.", default=".")
parser.add_argument("-r", "--RTDoseFolderPath", nargs='?', help="Directory  Directorio del Patient DICOM.", default=".")
parser.add_argument("-o", "--OutputPath", nargs='?', help="Directorio de salida.", default=".")
parser.add_argument("-s", "--StructuresPath", nargs='?', help="Directory containing structures nrrd binary labelmaps.", default=".")
parser.add_argument("-b", "--BeamNumber", type=int, nargs='+', help="Lista de números de haz.", default=[1])
parser.add_argument("-f", "--SaveFigures",help="Set to True in order to save the generated Figures",default=False)
parser.add_argument("-p", "--PatientIdx", type=int, help="Patient index for the TROTS Protons Dataset (from 1 to 20).", default=1) 
parser.add_argument("-m", "--MatFilePath", help="Directory of the original .mat file named beamlist.", default=".") 

args = parser.parse_args()

SaveFigures = args.SaveFigures if type(args.SaveFigures) == bool else args.SaveFigures == "True"

all_masks = []

for file in sorted(os.listdir(args.StructuresPath)):
    if file.endswith('.nrrd'):        
        organ_path = os.path.join(args.StructuresPath, file)        
        nrrd_data = nrrd.read(organ_path, index_order='C')
        masks = nrrd_data[0].astype(bool)
        clean_name = file.replace('1 RTSTRUCT Contours-', '')\
                         .replace('-label.nrrd', '')\
                         .replace('.nrrd', '')
                         
        all_masks.append({'organ_name': clean_name, 'mask': masks})
        print(f"Loaded organ: {clean_name}")
        
masks_org = [item['mask'] for item in all_masks]
global_mask = np.logical_or.reduce(masks_org)


def get_dvh(dosis_voxeles, dose_bins):
    hist, _ = np.histogram(dosis_voxeles, bins=dose_bins)
    cum_sum = np.cumsum(hist[::-1])[::-1]
    total_voxeles = len(dosis_voxeles)
    vols_proporcion = cum_sum / total_voxeles if total_voxeles > 0 else np.zeros_like(cum_sum)
    return np.insert(vols_proporcion * 100, 0, 100.0)


def get_topas_dvh(params, energies, topas_doses_list, dose_bins):
    polparams = [params['c0'].value, params['c1'].value, params['c2'].value]
    totalDose = np.zeros(len(topas_doses_list[0]))
    
    for energy, dose_cp_mask in zip(energies, topas_doses_list):
        scaleFactor = polyval(energy, polparams)
        totalDose += dose_cp_mask * scaleFactor
        
    return get_dvh(totalDose, dose_bins)


def residual1D(params, energies, topas_doses_list, reference_dvh, dose_bins):
    topas_dvh = get_topas_dvh(params, energies, topas_doses_list, dose_bins)
    return topas_dvh - reference_dvh

def fmt_latex(val):
    mantissa, exp = f"{val:.1e}".split('e')
    exp_clean = exp.replace('+', '') 
    return f"{mantissa} \\times 10^{{{exp_clean}}}"


mat_path = os.path.join(args.MatFilePath, "beamlist.mat")
if not os.path.exists(mat_path):
    raise FileNotFoundError(f"Beamlist.mat file not found at {args.MatFilePath}")

with h5py.File(mat_path, 'r') as f:
    patient_ref = f['BeamList'][args.PatientIdx - 1, 0]
    patient_matrix = f[patient_ref][:].T

global_sol = []

for beamnumber in args.BeamNumber:
    print(f"PROCESSING BEAM: {beamnumber}")
    
    dcmcp_path = os.path.join(args.RTDoseFolderPath, f"RTDose_Beam{beamnumber}_CPs")
    rtdose_topas_path = os.path.join(args.RTDoseFolderPath, f"RTDose_TOPAS_Beam{beamnumber}")
    plan_path = os.path.join(args.RTDoseFolderPath, "rtplan.dcm")
    beam_path = os.path.join(args.RTDoseFolderPath, f"rtdose_beam{beamnumber}.dcm")
    full_dose_path = os.path.join(args.RTDoseFolderPath, "rtdose.dcm")

    
    if not (os.path.exists(dcmcp_path) and os.path.exists(rtdose_topas_path) and os.path.exists(plan_path) and os.path.exists(beam_path)):
        print(f"Skipping Beam {beamnumber}, there are files missing.")
        continue 

    rtplan = dcmread(plan_path)
    beam_idx = int(beamnumber - 1)
    beam = rtplan.IonBeamSequence[beam_idx]
    full_dcm=dcmread(full_dose_path)
    
    full_dose = full_dcm.pixel_array * full_dcm.DoseGridScaling
    max_dose = np.max(full_dose)
    dose_bins = np.linspace(0, max_dose * 1.4, 300) 
    
    beam_data_mat = patient_matrix[patient_matrix[:, 0] == beamnumber]
    
    #We need to extract the energy before applying the RangeShifter 
    cp_to_energy = {}
    current_energy = float(beam_data_mat[0, 1]) 
    mat_row_idx = 0

    for i, cp_seq in enumerate(beam.IonControlPointSequence):
        if hasattr(cp_seq, "NumberOfScanSpotPositions") and int(cp_seq.NumberOfScanSpotPositions) > 0:
            num_spots = int(cp_seq.NumberOfScanSpotPositions)
            if mat_row_idx < len(beam_data_mat):
                current_energy = round(float(beam_data_mat[mat_row_idx, 1]), 2)
                mat_row_idx += num_spots
        cp_to_energy[i] = current_energy

    beam_dcm = dcmread(beam_path)
    beam_dose = beam_dcm.pixel_array * beam_dcm.DoseGridScaling
    max_dose_beam = np.max(beam_dose)
    dose_bins = np.linspace(0, max_dose_beam * 1.35, 300) 

    dicom_grids = {}
    cp_procesados = []
    for file in sorted(os.listdir(dcmcp_path)):
        if file.startswith(f'rtdose_beam{beamnumber}_CP'):
            match = re.search(r'CP(\d+)', file)
            if match:
                cp = int(match.group(1))
                path = os.path.join(dcmcp_path, file)
                rtdose_cp = dcmread(path)
                dicom_grids[cp] = rtdose_cp.pixel_array * rtdose_cp.DoseGridScaling
                cp_procesados.append(cp)

    topas_grids = {}
    for cp in cp_procesados:
        topas_file = f"rtdose_CP{cp}.dcm" 
        path_topas = os.path.join(rtdose_topas_path, topas_file)
        if os.path.exists(path_topas):
            rtdose_topas = dcmread(path_topas)
            topas_grids[cp] = rtdose_topas.pixel_array * rtdose_topas.DoseGridScaling
            
    #dicom dose with combined mask    
    dicom_dose = np.zeros(np.sum(global_mask))
    for cp in cp_procesados:
        dicom_dose += dicom_grids[cp][global_mask]
            
    reference_dvh = get_dvh(dicom_dose, dose_bins)

    topas_doses_list = []
    energies_list = []
    for cp in cp_procesados:
        if cp in topas_grids:
            energy = cp_to_energy.get(cp)
            if energy is not None:
                topas_doses_list.append(topas_grids[cp][global_mask])
                energies_list.append(energy)
    
    global_sol.append({
        'energies':energies_list,
        'topas doses':topas_doses_list
        })
    
df_beam={}
for i in range(len(global_sol)):
    df_beam[i]=pd.DataFrame.from_dict(global_sol[i])
sum_doses=pd.concat(df_beam).groupby('energies').sum().reset_index()

#Now that we have the complete doses together, we minimize over this global dose
energies = sum_doses['energies'].to_numpy()
topas_doses_list = sum_doses['topas doses'].to_numpy()

max_dose_plan = np.max(full_dose)
dose_bins = np.linspace(0, max_dose_plan * 1.4, 300) 

full_dose_global = full_dose[global_mask]
reference_dvh_global = get_dvh(full_dose_global, dose_bins)

params = Parameters()
params.add('c0', value=0.0, min=0.0, max=500.0, vary=True)
params.add('c1', value=8.63e3, min=0.0, vary=True)
params.add('c2', value=3.87e2, min=0.0, vary=True)

result = minimize(residual1D, params, args=(energies, topas_doses_list, reference_dvh_global, dose_bins), nan_policy='omit', method='Nelder-Mead')
report_fit(result)

opt_params = [result.params['c0'].value, result.params['c1'].value, result.params['c2'].value]

total_topas_dose_global = np.zeros(len(topas_doses_list[0]))
for energy, dose_cp_mask in zip(energies, topas_doses_list):
    scaleFactor = polyval(energy, opt_params)
    total_topas_dose_global += dose_cp_mask * scaleFactor


plt.figure(figsize=(12, 7))

table = []
for i, item in enumerate(all_masks):
    name = item['organ_name']
    mask = item['mask']
    
    sub_mask_indices = mask[global_mask]
    
    organ_dose_topas = total_topas_dose_global[sub_mask_indices]
    
    organ_dose_dicom = full_dose[mask]
    
    volume_voxel = int(np.sum(mask))
    
    if volume_voxel > 0:
        mean_dicom = np.mean(organ_dose_dicom)
        max_dicom = np.max(organ_dose_dicom)
        min_dicom = np.min(organ_dose_dicom)
        
        mean_topas = np.mean(organ_dose_topas)
        max_topas = np.max(organ_dose_topas)
        min_topas = np.min(organ_dose_topas)
    else:
        mean_dicom = max_dicom = min_dicom = 0
        mean_topas = max_topas = min_topas = 0
        
    table.append({
        'Organ': name,
        'Volume (Voxels)': volume_voxel,
        'Volume (cm^3)': volumen_voxeles*rtdose_cp.PixelSpacing[0]*rtdose_cp.PixelSpacing[1]*rtdose_cp.SliceThickness/1000,
        'Mean DICOM (Gy)': round(mean_dicom, 2),
        'Mean TOPAS (Gy)': round(mean_topas, 2),
        'Max DICOM (Gy)': round(max_dicom, 2),
        'Max TOPAS (Gy)': round(max_topas, 2),
        'Min DICOM (Gy)': round(min_dicom, 2),
        'Min TOPAS (Gy)': round(min_topas, 2)
    })
    dvh_topas = get_dvh(organ_dose_topas, dose_bins)
    dvh_dicom = get_dvh(organ_dose_dicom, dose_bins)
    
    color = plt.cm.tab10(i)
    
    plt.plot(dose_bins, dvh_dicom, linestyle='--', color=color)
    plt.plot(dose_bins, dvh_topas, label=f'{name}', linestyle='-', color=color)

df_stats = pd.DataFrame(table)
plt.title('DVH comparison',fontweight='bold')
plt.xlabel('Dose [Gy]')
plt.ylabel('Volume [%]')
plt.xlim(0, max_dose_plan * 1.35)
plt.ylim(0, 105)
plt.grid(True, linestyle=':', alpha=0.6)
leg = plt.legend(loc='best',title='--- DICOM | ── TOPAS',title_fontsize=13)
leg.get_title().set_fontweight('bold')

if SaveFigures==True:
    plt.savefig(args.OutputPath, "global_dvh.pdf")
    
stats_output_path = os.path.join(args.OutputPath, "dvh_table.csv")
df_stats.to_csv(stats_output_path, index=False, sep=';')

plt.figure(figsize=(9, 5))

energy_range = np.linspace(np.min(energies), np.max(energies), 100)
scaling_factors = polyval(energy_range, opt_params)
c0, c1, c2 = opt_params

text_polynom = f"$({fmt_latex(c2)}) \\cdot E^2 + ({fmt_latex(c1)}) \\cdot E + ({fmt_latex(c0)})$"
plt.plot(energy_range, scaling_factors, label=f'k(E)={text_polynom}', color='darkorange', linewidth=2)

plt.title('Global Calibration Curve',fontweight='bold')
plt.ticklabel_format(useMathText=True)
plt.xlabel('Energy [MeV]')
plt.ylabel('Calibration factor (k)')
plt.tight_layout()
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(loc='best')
if SaveFigures==True:
    plt.savefig(args.OutputPath, "curve.pdf")

plt.figure(figsize=(9, 5))

dvh_topas_global = get_topas_dvh(result.params, energies, topas_doses_list, dose_bins)

plt.plot(dose_bins, reference_dvh_global, label='DICOM', linestyle='--', color='teal', linewidth=2.5)
plt.plot(dose_bins, dvh_topas_global, label='TOPAS', linestyle='-', color='teal', linewidth=2)
textstr = '\n'.join((
    'Fitting method: Nelder-Mead',
    rf'$\chi^2/\nu = {result.redchi:.2f}$',
    rf'$c_0 = {fmt_latex(c0)}$',
    rf'$c_1 = {fmt_latex(c1)}$',
    rf'$c_2 = {fmt_latex(c2)}$'
))
props = dict(boxstyle='round', facecolor='lavender', alpha=0.5)
plt.text(56, 85, textstr, verticalalignment='top', bbox=props)    
   
plt.title('DVH for Virtual Organ',fontweight='bold')
plt.xlabel('Dose [Gy]')
plt.ylabel('Volume [%]')
plt.xlim(0, max_dose_plan * 1.35)
plt.ylim(0, 105)
plt.tight_layout()
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(loc='best')
plt.show()

if SaveFigures==True:
    plt.savefig(args.OutputPath, "virtual_dvh.pdf")

#Export values to calibration.txt 
calib_factor=polyval(energies, opt_params)
df_calibration=pd.DataFrame({'energies':energies,'calib factor': calib_factor})

calibration_output_path = os.path.join(args.OutputPath, "calibration_globalplan.txt")
df_calibration.to_csv(calibration_output_path,index=False, header=False, sep='\t')
