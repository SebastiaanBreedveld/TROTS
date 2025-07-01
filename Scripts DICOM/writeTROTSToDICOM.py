import os
import mat73
import numpy as np
import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
import csv
import seaborn as sns

caseFolders = ['Prostate_CK', 'Head-and-Neck', 'Protons', 'Liver', 'Prostate_BT', 'Prostate_VMAT', 'Head-and-Neck-Alt']
for folder in caseFolders:
    try:
        matFiles = [f.name for f in os.scandir('./'+folder)]
    except:
        continue
    for matFile in matFiles:
        print('Folder and file:', folder, matFile)
        patientIndexInt = matFiles.index(matFile) + 1
        patientFolder = f"{folder}_{patientIndexInt:02d}"
        if patientFolder not in [f.name for f in os.scandir('./DICOMs/')]:
            os.mkdir('./DICOMs/'+patientFolder)

        mat = mat73.loadmat(matFile)
        resolutionX = mat['patient']['Resolution'][0];
        resolutionY = mat['patient']['Resolution'][1];
        resolutionZ = mat['patient']['Resolution'][2];
        ctshape = mat['patient']['CT'].shape
        nRowsCT =  ctshape[1] # DICOM Rows, goes with y, is index 1 because of how matrix is stored
        nColumnsCT = ctshape[0] # DICOM Columns, goes with x, is index 0 because of how matrix is stored
        nSlicesCT = ctshape[2] # DICOM slices, goes with z
        ctsopids = {}
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

        for sliceIndex in range(0, nSlicesCT):
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

            ds.PatientName = 'TROTS' + '^' + matFile[:-4]
            ds.PatientID = 'TROTS_' + matFile[:-4]
            ds.PatientBirthDate = ""
            ds.PatientSex = ""
            ds.PatientAge = ""
            ds.SliceThickness = resolutionZ
            ds.SpacingBetweenSlices = ds.SliceThickness

            ds.StudyDate = "20230308"
            ds.StudyTime = "104455"
            ds.StudyDescription = "RT optimizers"
            ds.AccessionNumber = ""
            ds.Manufacturer = ""
            ds.InstitutionName = ""
            ds.ReferringPhysiciansName = ""
            ds.ImageType = ['ORIGINAL', 'PRIMARY', 'AXIAL']

            ds.ReferencedSOPClassUID = pydicom.uid.CTImageStorage
            ds.ReferencedSOPInstanceUID = ReferencedSOPInstanceUID

            ds.StudyInstanceUID = StudyInstanceUID
            ds.StudyID = 'TROTS'
            ds.SeriesInstanceUID = SeriesInstanceUID
            ds.SeriesNumber = 1
            ds.SeriesDescription = "TROTS anonymized " + folder

            ds.FrameOfReferenceUID = FORUID
            ds.PositionReferenceIndicator = ""

            ds.SOPClassUID = pydicom.uid.CTImageStorage
            ds.SOPInstanceUID = SOPInstanceUID
            ctsopids[sliceIndex] = SOPInstanceUID

            ds.Modality = "CT"

            ds.InstanceNumber = sliceIndex + 1
            ds.PatientPosition = 'HFS'
            ds.ImagePositionPatient = [mat['patient']['Offset'][0], mat['patient']['Offset'][1], mat['patient']['Offset'][2] + resolutionZ*sliceIndex]
            ds.ImageOrientationPatient = [1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000]

            ds.SamplesPerPixel = 1
            ds.PhotometricInterpretation = "MONOCHROME2"
            ds.Rows = nRowsCT
            ds.Columns = nColumnsCT
            ds.PixelSpacing = [resolutionX, resolutionY]

            ds.BitsAllocated = 16
            ds.BitsStored = 16
            ds.HighBit = 15
            ds.WindowCenter = "100.0"
            ds.WindowWidth = "400.0"
            ds.RescaleIntercept = "-1024.0"
            ds.RescaleSlope = "1.0"
            ds.RescaleType = "HU"
            ds.PixelRepresentation = 1

            ds.PixelData = np.array(np.swapaxes(mat['patient']['CT'][:,:,sliceIndex],0,1)+1024).tobytes()

            ds.save_as("./DICOMs/"+patientFolder+'/CTSlice'+str(sliceIndex).zfill(3)+".dcm", write_like_original = False)

        # write RTStruct
        meta = pydicom.Dataset()
        SOPInstanceUID = pydicom.uid.generate_uid()
        meta.MediaStorageSOPClassUID = pydicom.uid.RTStructureSetStorage
        meta.MediaStorageSOPInstanceUID = SOPInstanceUID
        meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
        rds = Dataset()
        rds.file_meta = meta

        rds.Modality = "RTSTRUCT"
        rds.PatientName = ds.PatientName
        rds.PatientID = ds.PatientID
        rds.PatientBirthDate = ds.PatientBirthDate
        rds.PatientSex = ds.PatientSex
        rds.PatientAge = ds.PatientAge
        rds.StudyDate = ds.StudyDate
        rds.StudyTime = ds.StudyTime
        rds.AccessionNumber = ds.AccessionNumber
        rds.Manufacturer = ds.Manufacturer
        rds.ReferringPhysicianName = ""
        rds.OperatorsName = ""
        
        rds.StudyInstanceUID = ds.StudyInstanceUID
        rds.StudyID = ds.StudyID
        rds.SeriesInstanceUID = pydicom.uid.generate_uid()
        rds.SeriesNumber = 1        
        rds.SOPClassUID = pydicom.uid.RTStructureSetStorage
        rds.SOPInstanceUID = SOPInstanceUID
        rds.FrameOfReferenceUID = ds.FrameOfReferenceUID
        rds.PositionReferenceIndicator = ds.PositionReferenceIndicator

        rds.StructureSetLabel = "Contours"
        rds.StructureSetDate = ds.StudyDate
        rds.StructureSetTime = ds.StudyTime
        
        rds.ApprovalStatus = "UNAPPROVED"

        rds.ReferencedFrameOfReferenceSequence = Sequence()
        rfor = Dataset()
        rfor.FrameOfReferenceUID = ds.FrameOfReferenceUID
        # rfor.RTReferencedStudySequence = Sequence() # Optional
        # rrs = Dataset()
        # rrs.ReferencedSOPClassUID = 
        # rfor.RTReferencedStudySequence.append(rrs)
        rds.ReferencedFrameOfReferenceSequence.append(rfor)

        rds.StructureSetROISequence = Sequence()
        rds.ROIContourSequence = Sequence()
        rds.RTROIObservationsSequence = Sequence()
        nROIs = len(mat['patient']['StructureNames'])
        COLOR_PALETTE = sns.color_palette("hls", nROIs)
        for structIndex, structName in enumerate(mat['patient']['StructureNames']):
            print('Working on structure',structIndex+1,'(',structName,')...')
            
            roi = Dataset()
            roi.ROINumber = structIndex + 1
            roi.ReferencedFrameOfReferenceUID = rds.FrameOfReferenceUID
            roi.ROIName = structName
            roi.ROIGenerationAlgorithm = ''
            rds.StructureSetROISequence.append(roi)
            
            robs = Dataset()
            robs.ObservationNumber = structIndex + 1
            robs.ReferencedROINumber = structIndex + 1
            if 'ring' in structName.lower() or 'shell' in structName.lower() or 'intermediate' in structName.lower():
                robs.RTROIInterpretedType = 'CONTROL'
            elif 'ptv' in structName.lower():
                robs.RTROIInterpretedType = 'PTV'
            elif 'ctv' in structName.lower():
                robs.RTROIInterpretedType = 'CTV'
            elif 'gtv' in structName.lower():
                robs.RTROIInterpretedType = 'GTV'
            elif structName.lower() == 'patient' or structName.lower() == 'external':
                robs.RTROIInterpretedType = 'EXTERNAL'
            else:
                robs.RTROIInterpretedType = 'OAR'
            robs.ROIInterpreter = ''
            rds.RTROIObservationsSequence.append(robs)
            
            rc = Dataset()
            rc.ROIDisplayColor = [int(c*255) for c in COLOR_PALETTE[structIndex]]
            rc.ReferencedROINumber = structIndex + 1
            rc.ContourSequence = Sequence()
            for sliceIndex in range(0, len(mat['patient']['Contours'][0])):
                if mat['patient']['Contours'][0][sliceIndex][structIndex] != None:
                    for subStructIndex in range(0, len(mat['patient']['Contours'][0][sliceIndex][structIndex])):
                        cdata = mat['patient']['Contours'][0][sliceIndex][structIndex][subStructIndex]
                        nPoints = cdata.shape[0]
                        Z = [mat['patient']['Offset'][2] + resolutionZ*sliceIndex] * nPoints
                        cdata = np.c_[cdata, Z]
                        cont = Dataset()
                        cont.ContourGeometricType = 'CLOSED_PLANAR'
                        cont.NumberOfContourPoints = nPoints
                        cont.ContourData = cdata.flatten().tolist()
                        cont.ContourImageSequence = Sequence()
                        ci = Dataset()
                        ci.ReferencedSOPClassUID = pydicom.uid.CTImageStorage
                        ci.ReferencedSOPInstanceUID = ctsopids[sliceIndex]
                        cont.ContourImageSequence.append(ci)
                        rc.ContourSequence.append(cont)
            rds.ROIContourSequence.append(rc)
            

        rds.save_as("./DICOMs/"+patientFolder+'/rtstruct.dcm', write_like_original = False)
