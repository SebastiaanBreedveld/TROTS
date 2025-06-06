import os
import mat73
import numpy as np
import pydicom
from pydicom.dataset import Dataset 
from pydicom.uid import ExplicitVRLittleEndian
from rt_utils import RTStructBuilder
from PIL import Image, ImageDraw
import csv
import matplotlib.pyplot as plt

flipAxes = True

caseFolders = ['Prostate_CK', 'Head-and-Neck', 'Protons', 'Liver', 'Prostate_BT', 'Prostate_VMAT', 'Head-and-Neck-Alt']
patientIndexTot = 1
for folder in caseFolders:
    matfiles = [f.name for f in os.scandir('./'+folder)]
    for matFile in matfiles:
        print('Folder and file:', folder, matFile)
        patientIndexInt = matfiles.index(matFile) + 1
#        patientFolder = folder + str(patientIndexInt) 
        patientFolder = f"{folder}_{patientIndexInt:02d}"
        if patientFolder not in [f.name for f in os.scandir('./DICOMs/')]:
            os.mkdir('./DICOMs/'+patientFolder)

        mat = mat73.loadmat('./'+folder+'/'+matFile)
        if flipAxes:
            mat['patient']['CT'] = np.swapaxes(mat['patient']['CT'], 0, 1)
   
        # write planning objectives into csv
        with open("./DICOMs/"+patientFolder+'/planning.csv', 'w', newline = '') as csvfile:
            filewriter = csv.writer(csvfile)
            if folder != 'Prostate_BT':
                filewriter.writerow(['Beam angle', 'Couch angle', 'Collimator'])
                for rowIndex in range(0, len(mat['patient']['Beams']['BeamConfig'])):
                    filewriter.writerow([mat['patient']['Beams']['BeamConfig'][rowIndex]['Gantry']] + [mat['patient']['Beams']['BeamConfig'][rowIndex]['Couch']] + [mat['patient']['Beams']['BeamConfig'][rowIndex]['Collimator']])

            filewriter.writerow(['Name', 'Minimise', 'Objective', 'Sufficient', 'Weight', 'Priority', 'IsConstraint'])
            for rowIndex in range(0, len(mat['problem']['Name'])):
                filewriter.writerow([mat['problem']['Name'][rowIndex]] + [mat['problem']['Minimise'][rowIndex]] + [mat['problem']['Objective'][rowIndex]] + [mat['problem']['Sufficient'][rowIndex]] + [mat['problem']['Weight'][rowIndex]] + [mat['problem']['Priority'][rowIndex]] + [mat['problem']['IsConstraint'][rowIndex]])

        for sliceIndex in range(0, np.array(mat['patient']['CT']).shape[2]):
            print('Working on CT slice',sliceIndex,'...')
            if sliceIndex == 0:
                ReferencedSOPInstanceUID = pydicom.uid.generate_uid()
                StudyInstanceUID = pydicom.uid.generate_uid()
                SeriesInstanceUID = pydicom.uid.generate_uid()
                FORUID = pydicom.uid.generate_uid()
            meta = pydicom.Dataset()
            SOPInstanceUID = pydicom.uid.generate_uid()
            meta.MediaStorageSOPClassUID = pydicom.uid.CTImageStorage
            meta.MediaStorageSOPInstanceUID = SOPInstanceUID
            meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian  
            ds = Dataset()
            ds.file_meta = meta

            ds.PatientName = str(patientIndexInt) + '^' + folder
            ds.PatientID = "0" * (6-len(str(patientIndexTot))) + str(patientIndexTot)
            ds.PatientBirthDate = ""
            ds.PatientSex = ""
            ds.SliceThickness = mat['patient']['Resolution'][2]
            ds.SpacingBetweenSlices = ds.SliceThickness

            ds.StudyDate = "20230308"
            ds.StudyTime = "104455"
            ds.AccessionNumber = ""
            ds.Manufacturer = ""
            ds.InstitutionName = ""
            ds.ReferringPhysiciansName = ""
            ds.ImageType = ['ORIGINAL', 'PRIMARY', 'AXIAL']

            ds.ReferencedSOPClassUID = pydicom.uid.CTImageStorage
            ds.ReferencedSOPInstanceUID = ReferencedSOPInstanceUID

            ds.StudyInstanceUID = StudyInstanceUID
            ds.StudyID = "0" * (6-len(str(patientIndexTot))) + str(patientIndexTot)
            ds.SeriesInstanceUID = SeriesInstanceUID
            ds.SeriesNumber = 1

            ds.FrameOfReferenceUID = FORUID 
            ds.FrameOfReferenceIndicator = ""

            ds.SOPClassUID = pydicom.uid.CTImageStorage
            ds.SOPInstanceUID = SOPInstanceUID

            ds.Modality = "CT"

            ds.InstanceNumber = sliceIndex + 1
            ds.PatientPosition = 'HFS'
            ds.ImagePositionPatient = [mat['patient']['Offset'][0], mat['patient']['Offset'][1], mat['patient']['Offset'][2] + mat['patient']['Resolution'][2]*(np.array(mat['patient']['CT']).shape[2]-1-sliceIndex)] 
            ds.ImageOrientationPatient = [1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000] 

            ds.SamplesPerPixel = 1
            ds.PhotometricInterpretation = "MONOCHROME2"
            ds.Rows = np.array(mat['patient']['CT']).shape[0]
            ds.Columns = np.array(mat['patient']['CT']).shape[1]
            ds.PixelSpacing = [mat['patient']['Resolution'][0],mat['patient']['Resolution'][1]]

            ds.BitsAllocated = 16
            ds.BitsStored = 16
            ds.HighBit = 15
            ds.WindowCenter = "100.0"
            ds.WindowWidth = "400.0"
            ds.RescaleIntercept = "-1024.0"
            ds.RescaleSlope = "1.0"
            ds.RescaleType = "HU"
            ds.PixelRepresentation = 1

            ds.PixelData = np.array(mat['patient']['CT'][:,:,sliceIndex]+1024).tobytes() 

            ds.save_as("./DICOMs/"+patientFolder+'/CTSlice'+'0'*(3-len(str(sliceIndex)))+str(sliceIndex)+".dcm", write_like_original = False)

        # write RTStruct
        rtstruct = RTStructBuilder.create_new(dicom_series_path="./DICOMs/"+patientFolder+'/')

        for structIndex in range(0, len(mat['patient']['StructureNames'])):
            print('Working on structure',structIndex,'(',mat['patient']['StructureNames'][structIndex],')...')
            structBinMap = np.zeros(mat['patient']['CT'].shape)
            for sliceIndex in range(0, len(mat['patient']['Contours'][0])):
                if mat['patient']['Contours'][0][sliceIndex][structIndex] != None:
                    for subStructIndex in range(0, len(mat['patient']['Contours'][0][sliceIndex][structIndex])): 
                        X = [(mat['patient']['Contours'][0][sliceIndex][structIndex][subStructIndex][valueIndex][0]-mat['patient']['Offset'][0])/mat['patient']['Resolution'][0] for valueIndex in range(0, len(mat['patient']['Contours'][0][sliceIndex][structIndex][subStructIndex]))] 
                        Y = [(mat['patient']['Contours'][0][sliceIndex][structIndex][subStructIndex][valueIndex][1]-mat['patient']['Offset'][1])/mat['patient']['Resolution'][1] for valueIndex in range(0, len(mat['patient']['Contours'][0][sliceIndex][structIndex][subStructIndex]))]
                        polygon = [(Y[valueIndex], X[valueIndex]) for valueIndex in range(0, len(X))] 
                        width = structBinMap.shape[1]
                        height = structBinMap.shape[0] 
                        img = Image.new('L', (width, height), 0)
                        ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
                        mask = np.array(img)
                        structBinMap[:,:,sliceIndex] += mask
            structBinMap[np.where(structBinMap > 0)] = 1
            structBinMap = np.array(structBinMap, dtype = bool)
            rtstruct.add_roi(
              mask = structBinMap, 
              name = mat['patient']['StructureNames'][structIndex])
    
        rtstruct.save("./DICOMs/"+patientFolder+'/structs.dcm')

        RTStruct = "./DICOMs/"+patientFolder+'/structs.dcm'
        RTStructContent = pydicom.dcmread(RTStruct)

        for structIndex in range(0, len(RTStructContent.StructureSetROISequence)):
            if 'ring' in RTStructContent.StructureSetROISequence[structIndex].ROIName.lower() or 'shell' in RTStructContent.StructureSetROISequence[structIndex].ROIName.lower():
                RTStructContent.RTROIObservationsSequence[structIndex].RTROIInterpretedType = 'CONTROL'
            elif 'ptv' in RTStructContent.StructureSetROISequence[structIndex].ROIName.lower():
                RTStructContent.RTROIObservationsSequence[structIndex].RTROIInterpretedType = 'PTV'
            elif 'ctv' in RTStructContent.StructureSetROISequence[structIndex].ROIName.lower():
                RTStructContent.RTROIObservationsSequence[structIndex].RTROIInterpretedType = 'CTV'
            elif RTStructContent.StructureSetROISequence[structIndex].ROIName.lower() == 'patient' or RTStructContent.StructureSetROISequence[structIndex].ROIName.lower() == 'external':
                RTStructContent.RTROIObservationsSequence[structIndex].RTROIInterpretedType = 'EXTERNAL'
            else:
                RTStructContent.RTROIObservationsSequence[structIndex].RTROIInterpretedType = 'OAR'

        RTStructContent.save_as("./DICOMs/"+patientFolder+'/structs.dcm') 

        patientIndexTot += 1
