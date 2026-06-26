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
parser.add_argument("-t", "--TopasFolderPath", nargs='?', help="Directory containing TOPAS rtdoses.", default=".")
parser.add_argument("-r", "--RTDoseFolderPath", nargs='?', help="Directory containing Patient DICOM.", default=".")
parser.add_argument("-o", "--OutputPath", nargs='?', help="Output directory.", default=".")
parser.add_argument("-s", "--StructuresPath", nargs='?', help="Directory containing structures nrrd binary labelmaps.", default=".")
parser.add_argument("-b", "--BeamNumber", type=int, nargs='+', help="List of beam numbers to process.", default=[1, 2, 3])
parser.add_argument("-f", "--SaveFigures", help="Set to True in order to save the generated Figures", default=False)
parser.add_argument("-p", "--PatientIdx", type=int, help="Patient index for the TROTS Protons Dataset (from 1 to 20).", default=1) 
parser.add_argument("-m", "--MatFilePath", help="Directory of the original .mat file named beamlist.", default=".") 
parser.add_argument("-d", "--DVHCalibrated", help="Set to True to show DVH with the calibrated files", default=False)
parser.add_argument("-v", "--virtualOrgan", help="Set to True to perform the calibration over an organ computed as the sum of all the structures. By default it will be done over CTV-High", default=False)


args = parser.parse_args()

SaveFigures = args.SaveFigures if isinstance(args.SaveFigures, bool) else args.SaveFigures == "True"
DVHCalibrated = args.DVHCalibrated if isinstance(args.DVHCalibrated, bool) else args.DVHCalibrated == "True"
virtualOrgan = args.virtualOrgan if isinstance(args.virtualOrgan, bool) else args.virtualOrgan == "True"

def get_dvh(voxel_doses, dose_bins):
    hist, _ = np.histogram(voxel_doses, bins=dose_bins)
    cum_sum = np.cumsum(hist[::-1])[::-1]
    total_voxels = len(voxel_doses)
    vol_proportion = cum_sum / total_voxels if total_voxels > 0 else np.zeros_like(cum_sum)
    return np.insert(vol_proportion * 100, 0, 100.0)

def get_topas_dvh(params, energies, topas_doses_list, dose_bins):
    polparams = [params['c0'].value, params['c1'].value, params['c2'].value]
    total_dose = np.zeros(len(topas_doses_list[0]))
    
    for energy, dose_cp_mask in zip(energies, topas_doses_list):
        scale_factor = polyval(energy, polparams)
        total_dose += dose_cp_mask * scale_factor
        
    return get_dvh(total_dose, dose_bins)

def residual1D(params, energies, topas_doses_list, reference_dvh, dose_bins):
    topas_dvh = get_topas_dvh(params, energies, topas_doses_list, dose_bins)
    useful_mask = (reference_dvh > 5) & (reference_dvh < 95)
    if not np.any(useful_mask):
        return topas_dvh - reference_dvh
    return topas_dvh[useful_mask] - reference_dvh[useful_mask]

def fmt_latex(val):
    mantissa, exp = f"{val:.2e}".split('e')
    exp_clean = exp.replace('+', '') 
    return f"{mantissa} \\times 10^{{{exp_clean}}}"

def process_and_plot_individual_beam(beam_number, target_mask, opt_params, dose_bins, args, patient_matrix, max_dose_plan):
    rtdose_topas_path = os.path.join(args.RTDoseFolderPath, f"RTDose_TOPAS_Beam{beam_number}")
    dcmcp_path = os.path.join(args.RTDoseFolderPath, f"RTDose_Beam{beam_number}_CPs")
    
    if not (os.path.exists(rtdose_topas_path) and os.path.exists(dcmcp_path)):
        print(f"Skipping Beam {beam_number}: Incomplete file paths.")
        return

    rtplan = dcmread(os.path.join(args.RTDoseFolderPath, "rtplan.dcm"))
    beam_seq = rtplan.IonBeamSequence[beam_number - 1]
    beam_data = patient_matrix[patient_matrix[:, 0] == beam_number]
    
    cp_to_energy = {}
    mat_row_idx = 0
    for i, cp_seq in enumerate(beam_seq.IonControlPointSequence):
        if i % 2 != 0: 
            continue
        if hasattr(cp_seq, "NumberOfScanSpotPositions"):
            num_spots = int(cp_seq.NumberOfScanSpotPositions)
            if num_spots > 0 and mat_row_idx < len(beam_data):
                cp_to_energy[i] = float(beam_data[mat_row_idx, 1])
                mat_row_idx += num_spots

    topas_dose_target = np.zeros(np.sum(target_mask))
    dicom_dose_target = np.zeros(np.sum(target_mask))
    
    for file in sorted(os.listdir(dcmcp_path)):
        if file.endswith('.dcm'):
            match = re.search(r'CP(\d+)', file)
            if match:
                cp = int(match.group(1))
                topas_cp_file = f"rtdose_CP{cp}.dcm"
                path_topas = os.path.join(rtdose_topas_path, topas_cp_file)
                path_dicom = os.path.join(dcmcp_path, file)
                
                if os.path.exists(path_topas) and os.path.exists(path_dicom):
                    energy = cp_to_energy.get(cp, 0)
                    scale_factor = polyval(energy, opt_params)
                    
                    topas_cp_dcm = dcmread(path_topas)
                    dicom_cp_dcm = dcmread(path_dicom)
                    
                    topas_dose_vox = topas_cp_dcm.pixel_array * topas_cp_dcm.DoseGridScaling
                    dicom_dose_vox = dicom_cp_dcm.pixel_array * dicom_cp_dcm.DoseGridScaling
                    
                    topas_dose_target += topas_dose_vox[target_mask] * scale_factor
                    dicom_dose_target += dicom_dose_vox[target_mask]

    dvh_dicom = get_dvh(dicom_dose_target, dose_bins)
    dvh_topas = get_dvh(topas_dose_target, dose_bins)
    
    plt.figure(figsize=(9, 5))
    plt.plot(dose_bins, dvh_dicom, label=f'DICOM Beam {beam_number}', linestyle='-', color='blue')
    plt.plot(dose_bins, dvh_topas, label=f'TOPAS Calibrated Beam {beam_number}', linestyle='--', color='blue')
    plt.title(f'Beam {beam_number} Individual DVH Comparison (CTV High)', fontweight='bold')
    plt.xlabel('Dose [Gy]')
    plt.ylabel('Volume [%]')
    plt.xlim(0, max_dose_plan)
    plt.ylim(0, 105)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(loc='best')
    plt.tight_layout()
    
    if SaveFigures:
        plt.savefig(f'{args.OutputPath}/beam{beam_number}_individual_dvh.pdf')
    plt.show()

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

ctv_mask = None
for item in all_masks:
    if item['organ_name'].lower().strip() == 'ctv high':
        ctv_mask = item['mask']
        break

target_mask_for_individual = ctv_mask if ctv_mask is not None else global_mask

if virtualOrgan:
    print("\n>> INFO: Optimizing over the Virtual Organ (Sum of all structures).")
    opt_mask = global_mask
    opt_target_name = "Virtual Organ"
else:
    if ctv_mask is not None:
        print("\n>> INFO: Optimizing over CTV High.")
        opt_mask = ctv_mask
        opt_target_name = "CTV High"
    else:
        print("\n>> WARNING: CTV High not found! Using Virtual Organ instead.")
        opt_mask = global_mask
        opt_target_name = "Virtual Organ"


if DVHCalibrated:
    print("\n MODE: CALIBRATED DVH VALIDATION (DRAWING ONLY) ")
    
    full_dose_path = os.path.join(args.RTDoseFolderPath, "rtdose.dcm")
    if not os.path.exists(full_dose_path):
        raise FileNotFoundError(f"Missing complete reference plan 'rtdose.dcm' at {args.RTDoseFolderPath}")
        
    full_dcm = dcmread(full_dose_path)
    full_dose = full_dcm.pixel_array * full_dcm.DoseGridScaling
    max_dose_plan = np.max(full_dose)
    dose_bins = np.linspace(0, max_dose_plan * 1.35, 300)
    
    total_calibrated_topas_dose = np.zeros_like(full_dose, dtype=float)
    
    for beamnumber in args.BeamNumber:
        topas_beam_file = f"RTDOSE_TOPAS_Beam{beamnumber}.dcm"
        topas_beam_path = os.path.join(args.RTDoseFolderPath, topas_beam_file)
        
        if not os.path.exists(topas_beam_path):
            topas_beam_path = os.path.join(args.RTDoseFolderPath, f"RTDOSE_TOPAS_Beam{beamnumber}")
            
        if os.path.exists(topas_beam_path):
            print(f"Loading calibrated TOPAS beam file: {topas_beam_file}")
            beam_dcm = dcmread(topas_beam_path)
            total_calibrated_topas_dose += beam_dcm.pixel_array * beam_dcm.DoseGridScaling
        else:
            print(f"Warning: Calibrated file for Beam {beamnumber} not found at {topas_beam_path}")

    plt.figure(figsize=(12, 7))
    table = []
    
    for i, item in enumerate(all_masks):
        name = item['organ_name']
        mask = item['mask']
        
        organ_dose_topas = total_calibrated_topas_dose[mask]
        organ_dose_dicom = full_dose[mask]
        voxel_volume = int(np.sum(mask))
        # We will also compute the chi-square value for each organ to assess the goodness of fit
        # The degrees of freedom should be the number of valid voxels minus the number of fitted parameters (3)
        red_chi_sq = 0
        chi_sq=0
        valid_mask = (organ_dose_topas > 0) & (organ_dose_dicom > 0)
        valid_count = np.count_nonzero(valid_mask)
        if valid_count > 3:
            nu = valid_count 
            chi_sq = np.sum(((organ_dose_topas[valid_mask] - organ_dose_dicom[valid_mask]) ** 2) / organ_dose_dicom[valid_mask])
            red_chi_sq = np.sum(((organ_dose_topas[valid_mask] - organ_dose_dicom[valid_mask]) ** 2) / organ_dose_dicom[valid_mask]) / nu
        if voxel_volume > 0:
            mean_dicom = np.mean(organ_dose_dicom)
            max_dicom = np.max(organ_dose_dicom)
            min_dicom = np.min(organ_dose_dicom)
            
            mean_topas = np.mean(organ_dose_topas)
            max_topas = np.max(organ_dose_topas)
            min_topas = np.min(organ_dose_topas)            
            mean_error_rel = (abs(mean_topas - mean_dicom) / mean_dicom * 100) if mean_dicom > 0 else 0
            max_error_rel = (abs(max_topas - max_dicom) / max_dicom * 100) if max_dicom > 0 else 0
            min_error_rel = (abs(min_topas - min_dicom) / min_dicom * 100) if min_dicom > 0 else 0
        else:
            mean_dicom = max_dicom = min_dicom = 0
            mean_topas = max_topas = min_topas = 0
            mean_error_rel = 0
            max_error_rel = 0
            min_error_rel = 0
            
        table.append({
            'Organ': name,
            'Volume (Voxels)': voxel_volume,
            'Volume (cm^3)': voxel_volume * full_dcm.PixelSpacing[0] * full_dcm.PixelSpacing[1] * full_dcm.SliceThickness / 1000,
            'Mean DICOM (Gy)': round(mean_dicom, 2),
            'Mean TOPAS (Gy)': round(mean_topas, 2),
            'Mean Error (%)': round(mean_error_rel, 2),
            'Max DICOM (Gy)': round(max_dicom, 2),
            'Max TOPAS (Gy)': round(max_topas, 2),
            'Max Error (%)': round(max_error_rel, 2),
            'Min DICOM (Gy)': round(min_dicom, 2),
            'Min TOPAS (Gy)': round(min_topas, 2),
            'Min Error (%)': round(min_error_rel, 2),
            'Chi-square':round(chi_sq, 2),
            'nu (degrees of freedom)':nu,
            'Chi-square reduced':round(red_chi_sq, 2)
        })
        
        dvh_topas = get_dvh(organ_dose_topas, dose_bins)
        dvh_dicom = get_dvh(organ_dose_dicom, dose_bins)
        
        color = plt.cm.tab10(i)
        plt.plot(dose_bins, dvh_dicom, linestyle='--', color=color)
        plt.plot(dose_bins, dvh_topas, label=f'{name}', linestyle='-', color=color)
        
    df_stats = pd.DataFrame(table)
    plt.title('Calibrated DVH Comparison', fontweight='bold')
    plt.xlabel('Dose [Gy]')
    plt.ylabel('Volume [%]')
    plt.xlim(0, max_dose_plan * 1.1)
    plt.ylim(0, 105)
    plt.tight_layout()
    plt.grid(True, linestyle=':', alpha=0.6)
    leg = plt.legend(loc='best', title='--- DICOM | ── TOPAS', title_fontsize=13)
    leg.get_title().set_fontweight('bold')
    
    if SaveFigures:
        plt.savefig(f'{args.OutputPath}/calibrated_global_dvh.pdf')
    
    stats_output_path = os.path.join(args.OutputPath, "calibrated_dvh_table.csv")
    df_stats.to_csv(stats_output_path, index=False, sep=';')
    plt.show()
    
else:
    print("\nMODE: OPTIMIZATION & CALIBRATION ROUTINE")
    mat_path = os.path.join(args.MatFilePath, "beamlist.mat")
    if not os.path.exists(mat_path):
        raise FileNotFoundError(f"Beamlist.mat file not found at {args.MatFilePath}")

    with h5py.File(mat_path, 'r') as f:
        patient_ref = f['BeamList'][args.PatientIdx - 1, 0]
        patient_matrix = f[patient_ref][:].T

    global_sol = []
    global_cp_to_energy = {} 
    full_dose_path = os.path.join(args.RTDoseFolderPath, "rtdose.dcm")
    full_dcm = dcmread(full_dose_path)
    full_dose = full_dcm.pixel_array * full_dcm.DoseGridScaling
    
    for beamnumber in args.BeamNumber:
        print(f"Processing Control Points for BEAM: {beamnumber}")
        
        dcmcp_path = os.path.join(args.RTDoseFolderPath, f"RTDose_Beam{beamnumber}_CPs")
        rtdose_topas_path = os.path.join(args.RTDoseFolderPath, f"RTDose_TOPAS_Beam{beamnumber}")
        plan_path = os.path.join(args.RTDoseFolderPath, "rtplan.dcm")
        beam_path = os.path.join(args.RTDoseFolderPath, f"rtdose_beam{beamnumber}.dcm")

        if not (os.path.exists(dcmcp_path) and os.path.exists(rtdose_topas_path) and os.path.exists(plan_path) and os.path.exists(beam_path)):
            print(f"Skipping Beam {beamnumber}, there are files missing.")
            continue 
        global_cp_to_energy[beamnumber] = {}
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
            if hasattr(cp_seq, "NumberOfScanSpotPositions"):
                if int(cp_seq.NumberOfScanSpotPositions)<= 0:
                    print(f'Control Point {i} has an invalid number of Spots.\n')
                num_spots = int(cp_seq.NumberOfScanSpotPositions)
                
                if mat_row_idx < len(beam_data_mat):
                    current_energy = float(beam_data_mat[mat_row_idx, 1])
                    mat_row_idx += num_spots
                else:
                    print('Beam data insuficient to obtain energies.\n')
            
            cp_to_energy[i] = current_energy
            global_cp_to_energy[beamnumber][i] = current_energy

        dicom_grids = {}
        processed_cps = []
        for file in sorted(os.listdir(dcmcp_path)):
            if file.startswith(f'rtdose_beam{beamnumber}_CP'):
                match = re.search(r'CP(\d+)', file)
                if match:
                    cp = int(match.group(1))
                    path = os.path.join(dcmcp_path, file)
                    rtdose_cp = dcmread(path)
                    dicom_grids[cp] = rtdose_cp.pixel_array * rtdose_cp.DoseGridScaling
                    processed_cps.append(cp)

        topas_grids = {}
        for cp in processed_cps:
            topas_file = f"rtdose_CP{cp}.dcm" 
            path_topas = os.path.join(rtdose_topas_path, topas_file)
            if os.path.exists(path_topas):
                rtdose_topas = dcmread(path_topas)
                topas_grids[cp] = rtdose_topas.pixel_array * rtdose_topas.DoseGridScaling
                
        topas_doses_list = []
        energies_list = []
        for cp in processed_cps:
            if cp in topas_grids:
                energy = cp_to_energy.get(cp)
                if energy is not None:
                    # Apply dynamic mask depending on the parameter flag
                    topas_doses_list.append(topas_grids[cp][opt_mask])
                    energies_list.append(energy)
        
        global_sol.append({
            'energies': energies_list,
            'topas doses': topas_doses_list
        })
        
    df_beam = {}
    for i in range(len(global_sol)):
        df_beam[i] = pd.DataFrame.from_dict(global_sol[i])
    sum_doses = pd.concat(df_beam).groupby('energies').sum().reset_index()

    energies = sum_doses['energies'].to_numpy()
    topas_doses_list = sum_doses['topas doses'].to_numpy()

    max_dose_plan = np.max(full_dose)
    dose_bins = np.linspace(0, max_dose_plan * 1.35, 300) 

    full_dose_opt = full_dose[opt_mask]
    reference_dvh_opt = get_dvh(full_dose_opt, dose_bins)

    params = Parameters()
    params.add('c0', value=-4.06e6, min=-2.0e7, max=0.0)
    params.add('c1', value=93750.0, min=10000.0, max=5.0e5)
    params.add('c2', value=0.0, min=-1000.0, max=1000.0)
    
    print("\nRunning Nelder-Mead optimization...")
    result = minimize(residual1D, params, args=(energies, topas_doses_list, reference_dvh_opt, dose_bins), nan_policy='omit', method='Nelder-Mead')
    report_fit(result)

    opt_params = [result.params['c0'].value, result.params['c1'].value, result.params['c2'].value]

    total_topas_dose_global = np.zeros(np.sum(global_mask))
    
    print("Reconstructing full dose for global comparison plots...")
    total_topas_dose_reconstruct = np.zeros(np.shape(full_dose))
    for beamnumber in args.BeamNumber:
        rtdose_topas_path = os.path.join(args.RTDoseFolderPath, f"RTDose_TOPAS_Beam{beamnumber}")
        if not os.path.exists(rtdose_topas_path):
            continue
            
        for file in sorted(os.listdir(rtdose_topas_path)):
            if file.startswith('rtdose_CP'):
                match = re.search(r'CP(\d+)', file)
                if match:
                    cp = int(match.group(1))
                    energy = global_cp_to_energy[beamnumber].get(cp, 0)
                    if energy > 0: 
                        dcm_cp = dcmread(os.path.join(rtdose_topas_path, file))
                        k = polyval(energy, opt_params)
                        total_topas_dose_reconstruct += (dcm_cp.pixel_array * dcm_cp.DoseGridScaling) * k
                    else:
                        print(f"Warning: Energía no encontrada para Beam {beamnumber}, CP {cp} durante la reconstrucción.")

    plt.figure(figsize=(12, 7))
    table = []
    
    for i, item in enumerate(all_masks):
        name = item['organ_name']
        mask = item['mask']
        
        organ_dose_topas = total_topas_dose_reconstruct[mask]
        organ_dose_dicom = full_dose[mask]
        voxel_volume = int(np.sum(mask))
        
        if voxel_volume > 0:
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
            'Volume (Voxels)': voxel_volume,
            'Volume (cm^3)': voxel_volume * full_dcm.PixelSpacing[0] * full_dcm.PixelSpacing[1] * full_dcm.SliceThickness / 1000,
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
    plt.title(f'DVH Comparison (Post-Optimization on {opt_target_name})', fontweight='bold')
    plt.xlabel('Dose [Gy]')
    plt.ylabel('Volume [%]')
    plt.xlim(0, max_dose_plan * 1.35)
    plt.ylim(0, 105)
    plt.tight_layout()
    plt.grid(True, linestyle=':', alpha=0.6)
    leg = plt.legend(loc='best', title='--- DICOM | ── TOPAS', title_fontsize=13)
    leg.get_title().set_fontweight('bold')
    
    if SaveFigures:
        plt.savefig(f'{args.OutputPath}/global_dvh.pdf')
    
    stats_output_path = os.path.join(args.OutputPath, "dvh_table.csv")
    df_stats.to_csv(stats_output_path, index=False, sep= ';')

    # Calibration Curve Plot
    plt.figure(figsize=(9, 5))
    energy_range = np.linspace(np.min(energies), np.max(energies), 100)
    scaling_factors = polyval(energy_range, opt_params)
    c0, c1, c2 = opt_params

    polynomial_text = f"$({fmt_latex(c2)}) \\cdot E^2 + ({fmt_latex(c1)}) \\cdot E + ({fmt_latex(c0)})$"
    plt.plot(energy_range, scaling_factors, label=f'k(E)={polynomial_text}', color='darkorange', linewidth=2)

    plt.title(f'Calibration Curve (Optimized on {opt_target_name})', fontweight='bold')
    plt.ticklabel_format(useMathText=True)
    plt.xlabel('Energy [MeV]')
    plt.ylabel('Calibration factor (k)')
    plt.tight_layout()
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(loc='best')
    
    if SaveFigures:
        plt.savefig(f'{args.OutputPath}/curve.pdf')

    # DVH specific target fit plot
    plt.figure(figsize=(9, 5))
    dvh_topas_opt_fit = get_topas_dvh(result.params, energies, topas_doses_list, dose_bins)

    plt.plot(dose_bins, reference_dvh_opt, label='DICOM', linestyle='--', color='teal', linewidth=2.5)
    plt.plot(dose_bins, dvh_topas_opt_fit, label='TOPAS', linestyle='-', color='teal', linewidth=2)
    
    textstr = '\n'.join((
        'Fitting method: Nelder-Mead',
        rf'$\chi^2/\nu = {result.redchi:.2f}$',
        rf'$c_0 = {fmt_latex(c0)}$',
        rf'$c_1 = {fmt_latex(c1)}$',
        rf'$c_2 = {fmt_latex(c2)}$'
    ))
    props = dict(boxstyle='round', facecolor='lavender', alpha=0.5)
    plt.text(56, 85, textstr, verticalalignment='top', bbox=props)
        
    plt.title(f'DVH Performance Fit for {opt_target_name}', fontweight='bold')
    plt.xlabel('Dose [Gy]')
    plt.ylabel('Volume [%]')
    plt.xlim(0, max_dose_plan * 1.35)
    plt.ylim(0, 105)
    plt.tight_layout()
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(loc='best')
    plt.show()

    if SaveFigures:
        plt.savefig(f'{args.OutputPath}/opt_target_dvh.pdf')

    calib_factor = polyval(energies, opt_params)
    df_calibration = pd.DataFrame({'energies': energies, 'calib factor': calib_factor})
    calibration_output_path = os.path.join(args.OutputPath, "calibration_globalplan.txt")
    df_calibration.to_csv(calibration_output_path, index=False, header=False, sep='\t')

    # INDIVIDUAL BEAM DRAWING 
    print("\n>> Drawing Individual Beams ajusted to CTV High with obtained calibration...")
    for beamnumber in args.BeamNumber:
        process_and_plot_individual_beam(
            beam_number=beamnumber,
            target_mask=target_mask_for_individual, 
            opt_params=opt_params,
            dose_bins=dose_bins,
            args=args,
            patient_matrix=patient_matrix,
            max_dose_plan=max_dose_plan
        )
