import argparse
import os
import mat73
import numpy as np
import copy
import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.valuerep import format_number_as_ds
from pydicom.tag import Tag
import csv
import seaborn as sns
from scipy.interpolate import griddata

# -b /opt --TreatmentMachineName CGTR_2021 --rtdose False --tuneID Spot1 --rsID RS7.4cm --keyHole True --halfGantry True --minEnergy 100 --maxEnergy 226.09
parser = argparse.ArgumentParser()
parser.add_argument("--Manufacturer", nargs='?', help="The name of the manufacturer to be saved in the DICOMs", default="")
parser.add_argument("--ManufacturerModelName", nargs='?', help="The name of the manufacturer model to be saved in the DICOMs", default="")
parser.add_argument("--InstitutionName", nargs='?', help="The name of the Institution to be saved in the DICOMs", default="")
parser.add_argument("--ReferringPhysicianName", nargs='?', help="The name of the referring physician to be saved in the DICOMs", default="")
parser.add_argument("--OperatorsName", nargs='?', help="The name of the operator to be saved in the DICOMs", default="")
parser.add_argument("--TreatmentMachineName", nargs='?', help="The name of the treatment machine name to be saved in the DICOMs", default="TROTS")
parser.add_argument("--keyHole", nargs='?', help="Whether to use keyhole technique instead of CLOSEDPLANAR_XOR", default=False)
parser.add_argument("--rtdose", nargs='?', help="Also output rtdose files", default=True)
parser.add_argument("--tuneID", nargs='?', help="ScanSpotTuneID", default='TROTS_1.0')
parser.add_argument("--rsID", nargs='?', help="Overwrite RangeShifterID", default=None)
parser.add_argument("--halfGantry", nargs='?', help="Whether to add couch angle instead of gantry angles over 180", default=False)
parser.add_argument("--minEnergy", nargs='?', help="Low threshold for supported minimum beam energy", default=0)
parser.add_argument("--maxEnergy", nargs='?', help="Low threshold for supported maximum beam energy", default=1000)
parser.add_argument("-b","--folderBasePath", nargs='?', help="The base directory in which the code is run containing all neccessary folders", default=".")
parser.add_argument("-o", "--outputPath", nargs='?', help="The output directory of the DICOM file", default="/tmp")
parser.add_argument("-n", "--DoseBeamNumber", type=list, nargs='?', help="A list of beam numbers to be calculated where a separate rtdose_<BeamNumber>.dcm is calculated, format: [BeamNumber_i, ..]. By default, all beams will be included.", default=None)
parser.add_argument("-c", "--DoseControlPoints", type=list, nargs='?', help="A list of control point numbers where a separate rtdose_<BeamNumber>_<ControlPointNumber>.dcm is calculated, format: [(BeamNumber_i,ControlPoint_i), ...]", default=[])
parser.add_argument("-s", "--DoseBeamSpots", type=list, nargs='?', help="A list of beam spot numbers where the rtdose_<BeamNumber>_<ControlPointNumber>_<BeamSpotNumber>.dcm is calculated, format: [(BeamNumber_i,ControlPoint_i,BeamSpotNumber_i), ...]", default=[])
args = parser.parse_args()

pydicom.config.settings.writing_validation_mode = pydicom.config.RAISE

caseFolders = ['Prostate_CK', 'Head-and-Neck', 'Protons', 'Liver', 'Prostate_BT', 'Prostate_VMAT', 'Head-and-Neck-Alt']
for folder in caseFolders:
    try:
        matFiles = [f.name for f in os.scandir(args.folderBasePath+"/"+folder)]
        matFiles.sort()
    except:
        print('No mat files found in folder ' + args.folderBasePath+"/" + folder + '. Skipping.')
        continue
    if(folder == 'Protons'):
        try:
            beamlistfolder = mat73.loadmat(args.folderBasePath+"/"+folder+'/BeamList.mat')
            machinedata = mat73.loadmat(args.folderBasePath+"/"+folder+'/MachineData.mat')
        except:
            print("Note: Protons folder is missing BeamList.mat and MachineData.mat, thus rtplan.dcm will not be created")

    for matFile in matFiles:
        if((matFile!='BeamList.mat') and (matFile!='MachineData.mat')):
            mat = mat73.loadmat(args.folderBasePath+"/"+folder + '/' + matFile)
        else:
            continue
        if not '01' in matFile:
            break
        print('Folder and file:', folder, matFile)
        patientFolder = matFile.split('.')[0]
        patientIndexInt = int(patientFolder.split('_')[1])
        outFolder = args.outputPath + "/DICOMs/" + folder + "/" + patientFolder + "/"
        os.makedirs(outFolder, exist_ok=True)
        resolutionX = mat['patient']['Resolution'][0];
        resolutionY = mat['patient']['Resolution'][1];
        resolutionZ = mat['patient']['Resolution'][2];
        ctshape = mat['patient']['CT'].shape
        nRowsCT =  ctshape[1] # DICOM Rows, goes with y, is index 1 because of how matrix is stored
        nColumnsCT = ctshape[0] # DICOM Columns, goes with x, is index 0 because of how matrix is stored
        nSlicesCT = ctshape[2] # DICOM slices, goes with z
        ctsopids = {}
        # write planning objectives into csv
        with open(outFolder+'planning.csv', 'w', newline = '') as csvfile:
            filewriter = csv.writer(csvfile)
            if folder != 'Prostate_BT':
                filewriter.writerow(['Beam angle', 'Couch angle', 'Collimator'])
                for rowIndex in range(0, len(mat['patient']['Beams']['BeamConfig'])):
                    filewriter.writerow([mat['patient']['Beams']['BeamConfig'][rowIndex]['Gantry']] + [mat['patient']['Beams']['BeamConfig'][rowIndex]['Couch']] + [mat['patient']['Beams']['BeamConfig'][rowIndex]['Collimator']])

            filewriter.writerow(['Name', 'DataName', 'Robust', 'Minimise', 'Objective', 'Sufficient', 'Weight', 'Priority', 'IsConstraint'])
            for rowIndex in range(0, len(mat['problem']['Name'])):
                roiName = mat['problem']['Name'][rowIndex]
                matrixIdx = int(mat['problem']['dataID'][rowIndex])-1
                dataName = mat['data']['matrix'][matrixIdx]["Name"]
                if roiName.upper() != 'MU':
                    structIdx = mat["patient"]["StructureNames"].index(roiName)
                    isRobust = int(mat["data"]["matrix"][matrixIdx]["A"].shape[0]%9 == 0 and int(mat["data"]["matrix"][matrixIdx]["A"].shape[0]/9==mat["patient"]["SampledVoxels"][structIdx].shape[1]))
                else:
                    isRobust = 0
                filewriter.writerow([roiName] + [dataName] + [isRobust] + [int(mat['problem']['Minimise'][rowIndex])] + [mat['problem']['Objective'][rowIndex]] + [mat['problem']['Sufficient'][rowIndex]] + [mat['problem']['Weight'][rowIndex]] + [int(mat['problem']['Priority'][rowIndex])] + [int(mat['problem']['IsConstraint'][rowIndex])])

        for sliceIndex in range(0, nSlicesCT):
            # print('Working on CT slice',sliceIndex,'...')
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
            ds.SliceThickness = format_number_as_ds(resolutionZ)
            ds.SpacingBetweenSlices = ds.SliceThickness

            ds.StudyDate = "20230308"
            ds.StudyTime = "104455"
            ds.StudyDescription = "RT optimizers https://sebastiaanbreedveld.nl/trots/"
            ds.AccessionNumber = ""
            ds.Manufacturer = args.Manufacturer
            ds.InstitutionName = args.InstitutionName
            ds.ReferringPhysicianName = args.ReferringPhysicianName
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
            ds.ImagePositionPatient = [format_number_as_ds(mat['patient']['Offset'][0]), format_number_as_ds(mat['patient']['Offset'][1]), format_number_as_ds(mat['patient']['Offset'][2] + resolutionZ*sliceIndex)]
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

            ds.save_as(outFolder+'CT_'+str(sliceIndex).zfill(3)+".dcm", write_like_original = False)

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
        rds.Manufacturer = args.Manufacturer
        rds.ReferringPhysicianName = args.ReferringPhysicianName
        rds.OperatorsName = args.OperatorsName

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
        rfor.RTReferencedStudySequence = Sequence() # Optional
        rrs = Dataset()
        rrs.ReferencedSOPClassUID = "1.2.840.10008.3.1.2.3.1" # Detached Study Management SOP Class (Retired)
        rrs.ReferencedSOPInstanceUID = ds.StudyInstanceUID
        rrs.RTReferencedSeriesSequence = Sequence()
        ss = Dataset()
        ss.SeriesInstanceUID = SeriesInstanceUID
        ss.ContourImageSequence = Sequence()
        for ctsopi in ctsopids.values():
            ci = Dataset()
            ci.ReferencedSOPClassUID = pydicom.uid.CTImageStorage
            ci.ReferencedSOPInstanceUID = ctsopi
            ss.ContourImageSequence.append(ci)
        rrs.RTReferencedSeriesSequence.append(ss)
        rfor.RTReferencedStudySequence.append(rrs)
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
                    keyHole = args.keyHole if type(args.keyHole)==bool else args.keyHole=='True'# True is needed by RayStation since CLOSEDPLANAR_XOR is not supported there. Both options lead to the same artefacts in Slicer3D: https://discourse.slicer.org/t/closedplanar-xor-visualization-artefact-with-holes-2d/43589/3
                    if keyHole:
                        singleContours = mat['patient']['Contours'][0][sliceIndex][structIndex]
                        lSc = len(singleContours) - 1
                        for sc in range(lSc):
                            singleContours[-1] = np.vstack([singleContours[-1], singleContours[lSc - 1 - sc][-1]])
                        cdata = np.vstack(singleContours)
                        if (cdata[0] == cdata[-1]).all():
                            cdata = cdata[:-1]
                        nPoints = cdata.shape[0]
                        Z = [mat['patient']['Offset'][2] + resolutionZ*sliceIndex] * nPoints
                        cdata = np.c_[cdata, Z]
                        cont = Dataset()
                        cont.ContourGeometricType = 'CLOSED_PLANAR' # TODO inner/outer? https://dicom.innolitics.com/ciods/rt-structure-set/roi-contour/30060039/30060040/30060050
                        cont.NumberOfContourPoints = nPoints
                        cont.ContourData = [format_number_as_ds(cp) for cp in cdata.flatten()]
                        cont.ContourImageSequence = Sequence()
                        ci = Dataset()
                        ci.ReferencedSOPClassUID = pydicom.uid.CTImageStorage
                        ci.ReferencedSOPInstanceUID = ctsopids[sliceIndex]
                        cont.ContourImageSequence.append(ci)
                        rc.ContourSequence.append(cont)
                    else:
                        for subStructIndex in range(0, len(mat['patient']['Contours'][0][sliceIndex][structIndex])):
                            cdata = mat['patient']['Contours'][0][sliceIndex][structIndex][subStructIndex]
                            # last Point is always repeated, remove it, since we define CLOSED_PLANAR
                            if (cdata[0] == cdata[-1]).all():
                                cdata = cdata[:-1]
                            nPoints = cdata.shape[0]
                            Z = [mat['patient']['Offset'][2] + resolutionZ*sliceIndex] * nPoints
                            cdata = np.c_[cdata, Z]
                            cont = Dataset()
                            cont.ContourGeometricType = 'CLOSEDPLANAR_XOR' # TODO inner/outer? https://dicom.innolitics.com/ciods/rt-structure-set/roi-contour/30060039/30060040/30060050
                            cont.NumberOfContourPoints = nPoints
                            cont.ContourData = [format_number_as_ds(cp) for cp in cdata.flatten()]
                            cont.ContourImageSequence = Sequence()
                            ci = Dataset()
                            ci.ReferencedSOPClassUID = pydicom.uid.CTImageStorage
                            ci.ReferencedSOPInstanceUID = ctsopids[sliceIndex]
                            cont.ContourImageSequence.append(ci)
                            rc.ContourSequence.append(cont)
            rds.ROIContourSequence.append(rc)

        rds.save_as(outFolder+'rtstruct.dcm', write_like_original = False)

        if((folder == 'Protons') and ('beamlistfolder' in locals()) and ('machinedata' in locals())):
            # write rtplan
            print('Working on rtplan...')
            beamSigmaEnergy = machinedata["BeamInfo"]["EnergyRangeList"]
            beamSigmas = machinedata["BeamInfo"]["Sigma"]
            SAD = [machinedata["BeamInfo"]["SAD"][0],machinedata["BeamInfo"]["SAD"][1]]
            meta = pydicom.Dataset()
            SOPInstanceUID = pydicom.uid.generate_uid()
            meta.MediaStorageSOPClassUID = pydicom.uid.RTIonPlanStorage

            meta.MediaStorageSOPInstanceUID = SOPInstanceUID
            meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian

            rtds = Dataset()
            rtds.file_meta = meta
            rtds.SOPClassUID = pydicom.uid.RTIonPlanStorage
            rtds.SOPInstanceUID = SOPInstanceUID
            rtds.StudyDate = rds.StudyDate
            rtds.SeriesDate = rds.StudyDate
            rtds.SeriesDescription = "solutionX"
            rtds.StudyTime = rds.StudyTime
            rtds.AccessionNumber = rds.AccessionNumber
            rtds.Modality = "RTPLAN"
            rtds.Manufacturer = args.Manufacturer
            # rtds.InstitutionName = args.InstitutionName # Optional
            rtds.SpecificCharacterSet='ISO_IR 100' # Optional
            rtds.InstanceCreationDate='20250708' # Optional
            rtds.InstanceCreationTime='112329.000000' # Optional
            rtds.SeriesTime='112329.000000' # Optional
            rtds.SoftwareVersions='pydicom' # Optional
            rtds.ReferringPhysicianName = args.ReferringPhysicianName
            rtds.OperatorsName = args.OperatorsName
            rtds.ManufacturerModelName    = args.ManufacturerModelName
            rtds.PatientName = rds.PatientName
            rtds.PatientID = rds.PatientID
            rtds.PatientBirthDate = rds.PatientBirthDate
            rtds.PatientSex = rds.PatientSex
            rtds.StudyInstanceUID = rds.StudyInstanceUID
            rtds.SeriesInstanceUID = rds.StudyInstanceUID
            rtds.StudyID = rds.StudyID
            rtds.SeriesNumber = 1
            # rtds.InstanceNumber = 1 # Optional
            rtds.FrameOfReferenceUID = rds.FrameOfReferenceUID
            rtds.PositionReferenceIndicator = ''
            rtds.RTPlanName = "TROTS" # Optional
            rtds.RTPlanLabel = "66/54 Gy"
            rtds.RTPlanDate = ""
            rtds.RTPlanTime = ""
            rtds.RTPlanGeometry = "PATIENT"
            rtds.TreatmentProtocols = "PencilBeamScanning" # Optional
            rtds.PlanIntent = "CURATIVE" # Optional
            rtds.ApprovalStatus = "UNAPPROVED"
            rtds.PrescriptionDescription = "IMPT_conv2" # Optional

            rtds.PatientSetupSequence = Sequence()
            patientSetup = Dataset()
            patientSetup.PatientPosition = mat["patient"]["PatientPosition"]
            patientSetup.PatientSetupNumber = 1
            # patientSetup.PatientSetupLabel = "Standard" # Optional
            rtds.PatientSetupSequence.append(patientSetup)

            patientIndex = patientIndexInt-1
            rtds.IonBeamSequence = Sequence()
            currentbeamlist = []
            beaminfo = {}
            beaminfo["BeamNumber"] = 1
            beaminfo["FileBeamNumber"] = int(beamlistfolder["BeamList"][patientIndex][0][0])
            beaminfo["FileIndexStart"] = 0
            beaminfo["FileIndexEnd"] = 0
            if(beamlistfolder["BeamList"][patientIndex][0][4]==0):
                beaminfo["RangeShifters"] = []
            else:
                beaminfo["RangeShifters"] = [beamlistfolder["BeamList"][patientIndex][0][4]]
            beaminfo["ConstantRangeShifter"] = True
            beaminfo["FinalCumulativeMetersetWeight"] = 0
            beaminfo["ControlPoints"] = []
            controlpointinfo = {}
            controlpointinfo["ControlPointNumber"] = 0
            controlpointinfo["BeamEnergy"] = beamlistfolder["BeamList"][patientIndex][0][1]
            controlpointinfo["MetersetWeights"] = []
            controlpointinfo["ScanSpotPositions"] = []
            controlpointinfo["CumulativeMetersetWeight"] = 0
            controlpointinfo["RangeShifter"] = beamlistfolder["BeamList"][patientIndex][0][4]
            controlpointinfo["FileIndexStart"] = 0
            controlpointinfo["FileIndexEnd"] = 0
            for rowindex in range(len(beamlistfolder["BeamList"][patientIndex])):
                row = beamlistfolder["BeamList"][patientIndex][rowindex]
                if((row[4]!=0) and (row[4] not in beaminfo["RangeShifters"])):
                    beaminfo["RangeShifters"].append(row[4])
                if(beaminfo["FileBeamNumber"] ==row[0]):
                    if((controlpointinfo["BeamEnergy"]==row[1]) and (controlpointinfo["RangeShifter"]==row[4])):
                        controlpointinfo["ScanSpotPositions"].extend(beamlistfolder["BeamList"][patientIndex][rowindex][2:4])
                        controlpointinfo["MetersetWeights"].append(mat["solutionX"][rowindex])
                    else:
                        controlpointinfo["FileIndexEnd"] = rowindex
                        beaminfo["ControlPoints"].append(copy.deepcopy(controlpointinfo))
                        controlpointinfo["FileIndexStart"] = controlpointinfo["FileIndexEnd"]
                        controlpointinfo["CumulativeMetersetWeight"] += sum(controlpointinfo["MetersetWeights"])
                        controlpointinfo["ControlPointNumber"] += 1
                        controlpointinfo["BeamEnergy"] = row[1]
                        controlpointinfo["MetersetWeights"] = [mat["solutionX"][rowindex]]
                        controlpointinfo["ScanSpotPositions"] = beamlistfolder["BeamList"][patientIndex][rowindex][2:4].tolist()
                        if(controlpointinfo["RangeShifter"]!=row[4]):
                            controlpointinfo["RangeShifter"] = row[4]
                            beaminfo["ConstantRangeShifter"] = False
                else:
                    controlpointinfo["FileIndexEnd"] = rowindex
                    beaminfo["FinalCumulativeMetersetWeight"] = controlpointinfo["CumulativeMetersetWeight"] + sum(controlpointinfo["MetersetWeights"])
                    beaminfo["ControlPoints"].append(copy.deepcopy(controlpointinfo))
                    beaminfo["FileIndexEnd"] = rowindex
                    currentbeamlist.append(copy.deepcopy(beaminfo))
                    beaminfo["FileIndexStart"] = beaminfo["FileIndexEnd"]
                    beaminfo["BeamNumber"] = beaminfo["BeamNumber"] + 1
                    beaminfo["FileBeamNumber"] = int(row[0])
                    beaminfo["ConstantRangeShifter"] = True
                    if(beamlistfolder["BeamList"][patientIndex][0][4]==0):
                        beaminfo["RangeShifters"] = []
                    else:
                        beaminfo["RangeShifters"] = [beamlistfolder["BeamList"][patientIndex][0][4]]
                    beaminfo["FinalCumulativeMetersetWeight"] = 0
                    beaminfo["ControlPoints"] = []
                    controlpointinfo["CumulativeMetersetWeight"] = 0

                    controlpointinfo["ControlPointNumber"] = 0
                    controlpointinfo["BeamEnergy"] = row[1]
                    controlpointinfo["MetersetWeights"] = [mat["solutionX"][rowindex]]
                    controlpointinfo["ScanSpotPositions"] = beamlistfolder["BeamList"][patientIndex][rowindex][2:4].tolist()
                    controlpointinfo["RangeShifter"] = row[4]
                    controlpointinfo["FileIndexStart"] = controlpointinfo["FileIndexEnd"]
            beaminfo["FinalCumulativeMetersetWeight"] = controlpointinfo["CumulativeMetersetWeight"] + sum(controlpointinfo["MetersetWeights"])
            controlpointinfo["FileIndexEnd"] = len(beamlistfolder["BeamList"][patientIndex])
            beaminfo["ControlPoints"].append(copy.deepcopy(controlpointinfo))
            beaminfo["FileIndexEnd"] = len(beamlistfolder["BeamList"][patientIndex])
            currentbeamlist.append(beaminfo)

            structToROINumber = {}
            for struct in rds.StructureSetROISequence:
                structToROINumber[struct.ROIName] = struct.ROINumber

            probleminfo = {}
            for index in range(len(mat["problem"]["Name"])): #all elements of mat["problem"] are the same length
                # All problem info is grouped according to the unique key = (dataID, weight, is constraint)
                if(mat["problem"]["Name"][index].lower()=="mu"):
                    continue

                # If two rows have the same key they are merged one being the upper and other being the lower bound
                key = (int(mat["problem"]["dataID"][index]),float(mat["problem"]["Weight"][index]),bool(mat["problem"]["IsConstraint"][index]))

                roiName = mat['problem']['Name'][index]
                matrixIdx = int(mat['problem']['dataID'][index])-1
                dataName = mat['data']['matrix'][matrixIdx]["Name"]
                structIdx = mat["patient"]["StructureNames"].index(roiName)
                isRobust = int(mat["data"]["matrix"][matrixIdx]["A"].shape[0]%9 == 0 and int(mat["data"]["matrix"][matrixIdx]["A"].shape[0]/9==mat["patient"]["SampledVoxels"][structIdx].shape[1]))
                # print(roiName, dataName, isRobust, mat["problem"]["IsConstraint"][index])
                if mat["problem"]["IsConstraint"][index] == 0:
                    if mat["problem"]["Objective"][index] != mat["problem"]["Sufficient"][index]:
                        print('Warning: sufficient objective',mat["problem"]["Sufficient"][index],'was specified but does not match with objective for',dataName)

                if(key not in probleminfo):
                    constraint = {}
                    constraint["Name"] = dataName
                    constraint["roiName"] = roiName
                    constraint["IsRobust"] = isRobust
                    if not mat["problem"]["IsConstraint"][index]:
                        constraint['Sufficient'] = mat["problem"]['Sufficient'][index]

                    if (("gtv" in constraint["Name"].lower()) or ("ptv" in constraint["Name"].lower()) or ("ctv" in constraint["Name"].lower())):
                        constraint["type"] = "TARGET"
                        if(mat["problem"]["Minimise"][index]==1):
                            constraint["Max"] = mat["problem"]["Objective"][index]
                            constraint["Min"] = ""
                        else:
                            constraint["Max"] = ""
                            constraint["Min"] = mat["problem"]["Objective"][index]
                    else:
                        # If it isn't a TARGET it must be an ORGAN_AT_RISK
                        constraint["type"] = "ORGAN_AT_RISK"
                        # ORGAN_AT_RISK can only have an upper bound
                        assert(mat["problem"]["Minimise"][index]==1)
                        constraint["Max"] = mat["problem"]["Objective"][index]
                    probleminfo[key] = constraint
                else:
                    # The constraints must correspond to the same structure
                    assert(probleminfo[key]["Name"] == dataName)
                    assert(probleminfo[key]["roiName"] == roiName)
                    assert(probleminfo[key]["IsRobust"] == isRobust)
                    # Only TARGET can have both upper and lower bound
                    assert(probleminfo[key]["type"]=="TARGET") # having upper and lower only make sense for target
                    # Only one upper and one lower can be defined for one key in a TARGET
                    if(mat["problem"]["Minimise"][index]==1):
                        if (probleminfo[key]["Max"]!=""):
                            print('Warning, found a duplicate max constraint,',dataName,'taking the minimum') # happens for Protons16.mat
                            probleminfo[key]["Max"] = min(mat["problem"]["Objective"][index], probleminfo[key]["Max"])
                        else:
                            probleminfo[key]["Max"] = mat["problem"]["Objective"][index]
                    else:
                        assert(probleminfo[key]["Min"]=="")
                        probleminfo[key]["Min"] = mat["problem"]["Objective"][index]

            rtds.DoseReferenceSequence = Sequence()
            for dosenumber,key in enumerate(probleminfo):
                doseReference = Dataset()
                doseReference.ReferencedROINumber = structToROINumber[probleminfo[key]["roiName"]]
                doseReference.DoseReferenceNumber = dosenumber + 1
                doseReference.DoseReferenceUID = pydicom.uid.generate_uid()
                doseReference.DoseReferenceStructureType = "VOLUME"
                doseReference.DoseReferenceDescription = ("Constraint: " if(key[2]) else "Objective: ") + probleminfo[key]["Name"]
                if(probleminfo[key]["IsRobust"]):
                    doseReference.DoseReferenceDescription += " (robust)"
                doseReference.DoseReferenceType = probleminfo[key]["type"]
                doseReference.ConstraintWeight = format_number_as_ds(float(key[1]))
                if(probleminfo[key]["type"] == "TARGET"):
                    if(probleminfo[key]["Min"]!=""):
                        doseReference.TargetMinimumDose = format_number_as_ds(float(probleminfo[key]["Min"]))
                    if(probleminfo[key]["Max"]!=""):
                        doseReference.TargetMaximumDose = format_number_as_ds(float(probleminfo[key]["Max"]))
                    sufficients = mat['problem']['Sufficient'][mat['problem']['Name'].index(probleminfo[key]["roiName"])]
                    if sufficients is not None and type(sufficients) == np.ndarray:
                        doseReference.TargetPrescriptionDose = format_number_as_ds(float(sufficients))
                else:
                    doseReference.OrganAtRiskMaximumDose = format_number_as_ds(float(probleminfo[key]["Max"]))
                    if args.TreatmentMachineName == 'CGTR_2021': #
                        continue
                rtds.DoseReferenceSequence.append(doseReference)

            rtds.FractionGroupSequence = Sequence()
            fractionds = Dataset()
            fractionds.FractionGroupNumber = 1
            fractionds.NumberOfFractionsPlanned = 1
            fractionds.NumberOfBeams = len(currentbeamlist)
            fractionds.NumberOfBrachyApplicationSetups = 0
            fractionds.ReferencedBeamSequence = Sequence()
            for beaminfo in currentbeamlist:
                be = Dataset()
                be.BeamMeterset = format_number_as_ds(beaminfo["FinalCumulativeMetersetWeight"])
                be.ReferencedBeamNumber = beaminfo["BeamNumber"]
                fractionds.ReferencedBeamSequence.append(be)
            # fractionds.ReferencedDoseReferenceSequence = Sequence() # Optional
            # rdr = Dataset()
            # rdr.ReferencedDoseReferenceNumber=1
            # fractionds.ReferencedDoseReferenceSequence.append(rdr)
            rtds.FractionGroupSequence.append(fractionds)

            rtds.ReferencedStructureSetSequence = Sequence()
            ref_rds = Dataset()
            ref_rds.ReferencedSOPClassUID = rds.SOPClassUID
            ref_rds.ReferencedSOPInstanceUID = rds.SOPInstanceUID
            rtds.ReferencedStructureSetSequence.append(ref_rds)

            totalMetersetWeightOfBeams = 0
            for beaminfo in currentbeamlist:
                be = Dataset()
                be.Manufacturer = args.Manufacturer # Optional
                # be.InstitutionName = '' # Optional # Optional
                # be.InstitutionAddress = args.InstitutionName # Optional
                # be.ManufacturerModelName = args.ManufacturerModelName # Optional
                # be.ToleranceTableNumber = '0' # Optional
                be.TreatmentMachineName = args.TreatmentMachineName
                be.PrimaryDosimeterUnit = "MU"
                be.BeamNumber             = beaminfo["BeamNumber"]
                be.BeamName               = str(beaminfo["BeamNumber"])
                be.BeamDescription        = '' # Optional
                be.BeamType               = 'STATIC'
                be.RadiationType          = "PROTON"
                be.TreatmentDeliveryType  = 'TREATMENT'
                be.NumberOfWedges         = 0
                be.NumberOfCompensators   = 0
                be.NumberOfBoli           = 0
                be.NumberOfBlocks         = 0
                be.FinalCumulativeMetersetWeight = format_number_as_ds(beaminfo["FinalCumulativeMetersetWeight"])
                be.ScanMode                   = 'MODULATED'
                be.VirtualSourceAxisDistances = SAD
                if(len(beaminfo["RangeShifters"]) ==0):
                    be.NumberOfRangeShifters = 0
                else:
                    numberOfRangeShifters = len(beaminfo["RangeShifters"])
                    be.NumberOfRangeShifters = numberOfRangeShifters
                    be.RangeShifterSequence = Sequence()
                    for RSindex in range(numberOfRangeShifters):
                        rsDataset = Dataset()
                        rsDataset.AccessoryCode = "Undefined Accessory Code"
                        rsDataset.RangeShifterNumber = RSindex
                        rsDataset.RangeShifterID = "Rs " + str(beaminfo["RangeShifters"][RSindex]) + " mm"
                        if not args.rsID == None:
                            rsDataset.RangeShifterID = args.rsID
                        rsDataset.RangeShifterType = "BINARY"
                        be.RangeShifterSequence.append(rsDataset)
                # Optional snout sequence:
                sn = Dataset()
                sn.SnoutID = 'Snout'
                be.SnoutSequence = Sequence()
                be.SnoutSequence.append(sn)
                be.NumberOfControlPoints  = 2*len(beaminfo["ControlPoints"])
                be.ReferencedPatientSetupNumber = 1 # Required for RS
                be.IonControlPointSequence = Sequence()
                be.NumberOfLateralSpreadingDevices = 0
                be.NumberOfRangeModulators = 0
                be.PatientSupportType = 'TABLE'
                be.PatientSupportID = 'TABLE' # Optional

                MetersetWeightTolerance = 1e-8
                totalMetersetWeightOfControlPoints = 0
                halfGantry = args.halfGantry if type(args.halfGantry)==bool else args.halfGantry=='True'
                for controlpointinfo in beaminfo["ControlPoints"]:
                    icpoi = Dataset()
                    icpoi.NominalBeamEnergyUnit = 'MEV'
                    icpoi.ControlPointIndex = 2*controlpointinfo["ControlPointNumber"]
                    icpoi.NominalBeamEnergy = format_number_as_ds(min(max(float(args.minEnergy),controlpointinfo["BeamEnergy"]), float(args.maxEnergy)))
                    icpoi.CumulativeMetersetWeight = format_number_as_ds(float(controlpointinfo["CumulativeMetersetWeight"]))
                    icpoi.GantryAngle = mat['patient']['Beams']['BeamConfig'][beaminfo["FileBeamNumber"]-1]['Gantry']
                    icpoi.GantryRotationDirection = 'NONE'
                    icpoi.BeamLimitingDeviceAngle = mat['patient']['Beams']['BeamConfig'][beaminfo["FileBeamNumber"]-1]['Collimator']
                    icpoi.BeamLimitingDeviceRotationDirection = 'NONE'
                    icpoi.PatientSupportAngle = mat['patient']['Beams']['BeamConfig'][beaminfo["FileBeamNumber"]-1]['Couch']
                    icpoi.PatientSupportRotationDirection = 'NONE'
                    icpoi.TableTopVerticalPosition     = None
                    icpoi.TableTopLongitudinalPosition = None
                    icpoi.TableTopLateralPosition      = None
                    icpoi.IsocenterPosition = [format_number_as_ds(isp) for isp in mat["patient"]["Isocentre"].flatten()]
                    icpoi.TableTopPitchAngle = 0
                    icpoi.TableTopPitchRotationDirection = 'NONE'
                    icpoi.TableTopRollAngle = 0
                    icpoi.TableTopRollRotationDirection = 'NONE'
                    icpoi.GantryPitchAngle = 0
                    icpoi.GantryPitchRotationDirection = 'NONE'
                    if halfGantry and icpoi.GantryAngle > 180:
                        icpoi.PatientSupportAngle += 180
                        if icpoi.PatientSupportAngle > 180: # move it to +/-180
                            icpoi.PatientSupportAngle -= 360
                        icpoi.GantryAngle = 360 - icpoi.GantryAngle
                    icpoi.SnoutPosition = 0 # or None if you remove the optional SnoutSequence
                    icpoi.ScanSpotTuneID = args.tuneID
                    icpoi.NumberOfScanSpotPositions = len(controlpointinfo["MetersetWeights"])
                    icpoi.ScanSpotPositionMap = controlpointinfo["ScanSpotPositions"]
                    icpoi.ScanSpotMetersetWeights = controlpointinfo["MetersetWeights"]
                    sigma1 = np.interp(controlpointinfo["BeamEnergy"],beamSigmaEnergy[:][0], beamSigmas[:][0])
                    sigma2 = np.interp(controlpointinfo["BeamEnergy"],beamSigmaEnergy[:][1], beamSigmas[:][1])
                    icpoi.ScanningSpotSize = [sigma1,sigma2]
                    icpoi.NumberOfPaintings = 1
                    if((beaminfo["ConstantRangeShifter"]==False) or ((controlpointinfo["ControlPointNumber"]==0) and (controlpointinfo["RangeShifter"]!=0))):
                        icpoi.RangeShifterSettingsSequence = Sequence()
                        for RSindex in range(numberOfRangeShifters):
                           rsSettings = Dataset()
                           rsSettings.RangeShifterSetting = "OUT" if(controlpointinfo["RangeShifter"]==0) else "IN"
                           rsSettings.RangeShifterWaterEquivalentThickness = beaminfo["RangeShifters"][RSindex]
                           rsSettings.ReferencedRangeShifterNumber = RSindex
                        icpoi.RangeShifterSettingsSequence.append(rsSettings)
                    assert(abs(icpoi.CumulativeMetersetWeight - totalMetersetWeightOfControlPoints) < MetersetWeightTolerance)
                    if(icpoi.NumberOfScanSpotPositions != 1):
                        totalMetersetWeightOfControlPoints += sum(icpoi.ScanSpotMetersetWeights)
                    else:
                        totalMetersetWeightOfControlPoints += icpoi.ScanSpotMetersetWeights
                    # icpoi.ReferencedDoseReferenceSequence = Sequence() # Optional
                    # for rdn in range(len(...))
                        # rdrv = Dataset()
                        # rdrv.CumulativeDoseReferenceCoefficient = None
                        # rdrv.ReferencedDoseReferenceNumber = rdn
                        # icpoi.ReferencedDoseReferenceSequence.append(rdrv)
                    be.IonControlPointSequence.append(icpoi)

                    icpoi = Dataset()
                    icpoi.NominalBeamEnergyUnit = 'MEV'
                    icpoi.ControlPointIndex = 2*controlpointinfo["ControlPointNumber"] + 1
                    icpoi.NominalBeamEnergy = format_number_as_ds(min(max(float(args.minEnergy),controlpointinfo["BeamEnergy"]), float(args.maxEnergy)))
                    icpoi.CumulativeMetersetWeight = format_number_as_ds(controlpointinfo["CumulativeMetersetWeight"]+ sum(controlpointinfo["MetersetWeights"]))
                    icpoi.ScanSpotTuneID = args.tuneID
                    icpoi.NumberOfScanSpotPositions = len(controlpointinfo["MetersetWeights"])
                    icpoi.ScanSpotPositionMap = controlpointinfo["ScanSpotPositions"]
                    icpoi.ScanSpotMetersetWeights = [0.0 for i in range(len(controlpointinfo["MetersetWeights"]))]
                    icpoi.ScanningSpotSize = [sigma1,sigma2]
                    icpoi.NumberOfPaintings = 1
                    if(beaminfo["ConstantRangeShifter"]==False):
                        icpoi.RangeShifterSettingsSequence = Sequence()
                        rsSettings = Dataset()
                        if(controlpointinfo["RangeShifter"]==0):
                            rsSettings.RangeShifterSetting = "OUT"
                            rsSettings.RangeShifterWaterEquivalentThickness = beaminfo["RangeShifters"][0]
                            rsSettings.ReferencedRangeShifterNumber = 0
                        else:
                            rsSettings.RangeShifterSetting = "IN"
                            rsSettings.RangeShifterWaterEquivalentThickness = controlpointinfo["RangeShifter"]
                            rsSettings.ReferencedRangeShifterNumber = beaminfo["RangeShifters"].index(controlpointinfo["RangeShifter"])
                        icpoi.RangeShifterSettingsSequence.append(rsSettings)
                    assert(abs(icpoi.CumulativeMetersetWeight - totalMetersetWeightOfControlPoints) < MetersetWeightTolerance)
                    if(icpoi.NumberOfScanSpotPositions != 1):
                        totalMetersetWeightOfControlPoints += sum(icpoi.ScanSpotMetersetWeights)
                    else:
                        totalMetersetWeightOfControlPoints += icpoi.ScanSpotMetersetWeights
                    be.IonControlPointSequence.append(icpoi)
                assert(abs(be.FinalCumulativeMetersetWeight - icpoi.CumulativeMetersetWeight) < MetersetWeightTolerance)
                assert(abs(be.FinalCumulativeMetersetWeight - totalMetersetWeightOfControlPoints) < MetersetWeightTolerance)
                rtds.IonBeamSequence.append(be)
                totalMetersetWeightOfBeams += be.FinalCumulativeMetersetWeight
            assert(abs(totalMetersetWeightOfBeams - sum(mat["solutionX"])) < MetersetWeightTolerance)
            rtds.save_as(outFolder+'rtplan.dcm', write_like_original = False)

            # write rtdose
            if not args.rtdose or (type(args.rtdose) == str and args.rtdose=='False'):
                continue
            print('Working on rtdose...')
            meta = pydicom.Dataset()
            SOPInstanceUID = pydicom.uid.generate_uid()
            meta.MediaStorageSOPClassUID = pydicom.uid.RTDoseStorage
            meta.MediaStorageSOPInstanceUID = SOPInstanceUID
            meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian

            doseds = Dataset()
            doseds.file_meta = meta
            doseds.SOPClassUID = pydicom.uid.RTDoseStorage
            doseds.SOPInstanceUID = SOPInstanceUID
            doseds.StudyDate = rds.StudyDate
            doseds.SeriesDate = rds.StudyDate
            doseds.SeriesDescription = "solutionX Plan"
            doseds.StudyTime = rds.StudyTime
            doseds.AccessionNumber = rds.AccessionNumber
            doseds.Modality = "RTDOSE"
            doseds.Manufacturer = args.Manufacturer
            doseds.InstitutionName = args.InstitutionName
            doseds.ReferringPhysicianName = args.ReferringPhysicianName
            doseds.OperatorsName = args.OperatorsName
            doseds.ManufacturerModelName    = args.ManufacturerModelName
            doseds.PatientName = rds.PatientName
            doseds.PatientID = rds.PatientID
            doseds.PatientBirthDate = rds.PatientBirthDate
            doseds.PatientSex = rds.PatientSex
            doseds.SliceThickness = ds.SliceThickness
            doseds.StudyInstanceUID = rds.StudyInstanceUID
            doseds.SeriesInstanceUID = pydicom.uid.generate_uid()
            doseds.StudyID = rds.StudyID
            doseds.SeriesNumber = 1
            doseds.InstanceNumber = 1
            db=mat["patient"]["DoseBox"]
            resolution = mat["patient"]["Resolution"]
            poffset = mat["patient"]["Offset"]
            doseds.ImagePositionPatient = [format_number_as_ds(poffset[i] + (db[i][0] - 1)*resolution[i]) for i in range(3)]
            doseds.ImageOrientationPatient = ds.ImageOrientationPatient
            doseds.FrameOfReferenceUID = ds.FrameOfReferenceUID
            doseds.SamplesPerPixel = ds.SamplesPerPixel
            doseds.PhotometricInterpretation = ds.PhotometricInterpretation
            doseds.FrameIncrementPointer = Tag(0x3004,0x000c)
            doseds.Rows = int(db[1][1] - db[1][0] + 1)
            doseds.Columns = int(db[0][1] - db[0][0] + 1)
            doseds.NumberOfFrames = db[2][1] - db[2][0] + 1
            doseds.GridFrameOffsetVector = [format_number_as_ds(doseds.ImagePositionPatient[2] + dz*resolution[2]) for dz in range(doseds.NumberOfFrames)]
            doseds.PixelSpacing = ds.PixelSpacing
            doseds.BitsAllocated = ds.BitsAllocated
            doseds.BitsStored = ds.BitsStored
            doseds.HighBit = ds.HighBit
            doseds.PixelRepresentation = 0
            doseds.DoseUnits = "GY"
            doseds.DoseType = "PHYSICAL"
            doseds.DoseSummationType = "PLAN"
            doseds.DoseGridScaling = 1
            doseds.PositionReferenceIndicator = ds.PositionReferenceIndicator

            doseds.ReferencedRTPlanSequence = Sequence()
            seqrt = Dataset()
            seqrt.ReferencedSOPClassUID = rtds.SOPClassUID
            seqrt.ReferencedSOPInstanceUID = rtds.SOPInstanceUID
            doseds.ReferencedRTPlanSequence.append(seqrt)

            SampledIndices = [np.array([[db[0,0], db[1,0], db[2,0]]]),
                              np.array([[db[0,1], db[1,0], db[2,0]]]),
                              np.array([[db[0,1], db[1,1], db[2,0]]]),
                              np.array([[db[0,0], db[1,1], db[2,0]]]),
                              np.array([[db[0,0], db[1,0], db[2,1]]]),
                              np.array([[db[0,1], db[1,1], db[2,1]]]),
                              np.array([[db[0,1], db[1,0], db[2,1]]]),
                              np.array([[db[0,0], db[1,1], db[2,1]]])]
            SampledDoses = [0.0 for i in range(len(SampledIndices))]

            for matrix in mat["data"]["matrix"]:
                if((matrix["A"].shape[0] > 1) and (np.any(matrix["A"].nonzero()[0]<0)==False)):
                    if(matrix["Name"] in mat["patient"]["StructureNames"]):
                        structIdx = mat["patient"]["StructureNames"].index(matrix["Name"])
                    else:
                        continue

                    if(matrix["A"].shape[0]==mat["patient"]["SampledVoxels"][structIdx].shape[1]):
                        SampledDosePart = matrix["A"]*mat["solutionX"] + matrix["b"]
                    elif (((matrix["A"].shape[0] % 9)==0) and (matrix["A"].shape[0]/9==mat["patient"]["SampledVoxels"][structIdx].shape[1])):
                        SampledDosePart = matrix["A"][:int(matrix["A"].shape[0]/9)]*mat["solutionX"] + matrix["b"]
                    else:
                        raise ValueError("")

                    SampledIndices.extend(mat["patient"]["SampledVoxels"][structIdx].T)
                    SampledDoses.extend(SampledDosePart)
            SampledIndices = np.vstack(SampledIndices)
            SampledDoses = np.array(SampledDoses)

            IndicesToCalculate = np.mgrid[db[0][0]:db[0][1] + 1, db[1][0]:db[1][1] + 1, db[2][0]:db[2][1] + 1]
            IndicesToCalculate = np.rollaxis(IndicesToCalculate, 0, 4)
            subct = mat["patient"]["CT"][db[0][0]-1:db[0][1], db[1][0]-1:db[1][1], db[2][0]-1:db[2][1]]
            mask = subct > -1024
            dose = np.zeros_like(subct)
            dose[mask] = griddata(SampledIndices, SampledDoses, IndicesToCalculate[mask])
            doseds.PixelData = np.swapaxes(dose, 2, 0).flatten().tobytes()

            doseds.save_as(outFolder + 'rtdose.dcm', write_like_original = False)

            beamnrs = [beaminfo["BeamNumber"] for beaminfo in currentbeamlist]
            if args.DoseBeamNumber is None:
                args.DoseBeamNumber = beamnrs
            for beamNumber in args.DoseBeamNumber:
                if not beamNumber in beamnrs:
                    print('Wrong beam number', beamNumber)
                    continue
                print('Working on rtdose, Beam:'+str(beamNumber)+'...')
                bdoseds = copy.deepcopy(doseds)
                bdoseds.SOPInstanceUID = pydicom.uid.generate_uid()
                bdoseds.SeriesInstanceUID = pydicom.uid.generate_uid()
                bdoseds.SeriesDescription = "solutionX Beam" + str(beamNumber)
                bdoseds.DoseSummationType = "BEAM"
                assert(len(bdoseds.ReferencedRTPlanSequence) == 1)
                bdoseds.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence = Sequence()
                rfgs = Dataset()
                rfgs.ReferencedBeamSequence = Sequence()
                rfgs.ReferencedBeamSequence.append(Dataset())
                rfgs.ReferencedBeamSequence[0].ReferencedBeamNumber = beamNumber
                rfgs.ReferencedFractionGroupNumber = 1
                bdoseds.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence.append(rfgs)

                idxStart = currentbeamlist[beamnrs.index(beamNumber)]["FileIndexStart"]
                idxEnd = currentbeamlist[beamnrs.index(beamNumber)]["FileIndexEnd"]
                SampledIndices = [np.array([[db[0,0], db[1,0], db[2,0]]]),
                                  np.array([[db[0,1], db[1,0], db[2,0]]]),
                                  np.array([[db[0,1], db[1,1], db[2,0]]]),
                                  np.array([[db[0,0], db[1,1], db[2,0]]]),
                                  np.array([[db[0,0], db[1,0], db[2,1]]]),
                                  np.array([[db[0,1], db[1,1], db[2,1]]]),
                                  np.array([[db[0,1], db[1,0], db[2,1]]]),
                                  np.array([[db[0,0], db[1,1], db[2,1]]])]
                SampledDoses = [0.0 for i in range(len(SampledIndices))]
                for matrix in mat["data"]["matrix"]:
                    if((matrix["A"].shape[0] > 1) and (np.any(matrix["A"].nonzero()[0]<0)==False)):
                        if(matrix["Name"] in mat["patient"]["StructureNames"]):
                            structIdx = mat["patient"]["StructureNames"].index(matrix["Name"])
                        else:
                            continue

                        if(matrix["A"].shape[0]==mat["patient"]["SampledVoxels"][structIdx].shape[1]):
                            SampledDosePart = matrix["A"][:,idxStart:idxEnd]*mat["solutionX"][idxStart:idxEnd] + matrix["b"]
                        elif (((matrix["A"].shape[0] % 9)==0) and (matrix["A"].shape[0]/9==mat["patient"]["SampledVoxels"][structIdx].shape[1])):
                            SampledDosePart = matrix["A"][:int(matrix["A"].shape[0]/9)][:,idxStart:idxEnd]*mat["solutionX"][idxStart:idxEnd] + matrix["b"]
                        else:
                            raise ValueError("")

                        SampledIndices.extend(mat["patient"]["SampledVoxels"][structIdx].T)
                        SampledDoses.extend(SampledDosePart)
                SampledIndices = np.vstack(SampledIndices)
                SampledDoses = np.array(SampledDoses)

                dose = np.zeros_like(subct)
                dose[mask] = griddata(SampledIndices, SampledDoses, IndicesToCalculate[mask])
                bdoseds.PixelData = np.swapaxes(dose, 2, 0).flatten().tobytes()

                bdoseds.save_as(outFolder+'rtdose_beam'+str(beamNumber)+'.dcm', write_like_original = False)

            for controlPointNumber in args.DoseControlPoints:
                print('Working on rtdose, Beam:'+str(controlPointNumber[0])+' ControlPoint:'+str(controlPointNumber[1])+'...')
                cpdoseds = copy.deepcopy(doseds)
                cpdoseds.SOPInstanceUID = pydicom.uid.generate_uid()
                cpdoseds.SeriesInstanceUID = pydicom.uid.generate_uid()
                cpdoseds.SeriesDescription = "solutionX Beam" + str(controlPointNumber[0]) + " CP" + str(controlPointNumber[1])
                cpdoseds.DoseSummationType = "CONTROL_POINT"
                if(int(controlPointNumber[1]) % 2 == 0):
                    print("NOTE: this control point contains no dose, so skipping")
                    continue
                idxStart = currentbeamlist[controlPointNumber[0]-1]["ControlPoints"][int((controlPointNumber[1]-1)/2)]["FileIndexStart"]
                idxEnd = currentbeamlist[controlPointNumber[0]-1]["ControlPoints"][int((controlPointNumber[1]-1)/2)]["FileIndexEnd"]

                assert(len(cpdoseds.ReferencedRTPlanSequence) == 1)
                cpdoseds.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence = Sequence()
                rfgs = Dataset()
                rfgs.ReferencedBeamSequence = Sequence()
                rfgs.ReferencedBeamSequence.append(Dataset())
                rfgs.ReferencedBeamSequence[0].ReferencedBeamNumber = controlPointNumber[0]
                rcps = Sequence()
                rcpp = Dataset()
                rcpp.ReferencedStartControlPointIndex = controlPointNumber[1]
                rcpp.ReferencedStopControlPointIndex = controlPointNumber[1]
                rcps.append(rcpp)
                rfgs.ReferencedBeamSequence[0].ReferencedControlPointSequence = rcps
                rfgs.ReferencedFractionGroupNumber = 1
                cpdoseds.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence.append(rfgs)

                SampledIndices = [np.array([[db[0,0], db[1,0], db[2,0]]]),
                                  np.array([[db[0,1], db[1,0], db[2,0]]]),
                                  np.array([[db[0,1], db[1,1], db[2,0]]]),
                                  np.array([[db[0,0], db[1,1], db[2,0]]]),
                                  np.array([[db[0,0], db[1,0], db[2,1]]]),
                                  np.array([[db[0,1], db[1,1], db[2,1]]]),
                                  np.array([[db[0,1], db[1,0], db[2,1]]]),
                                  np.array([[db[0,0], db[1,1], db[2,1]]])]
                SampledDoses = [0.0 for i in range(len(SampledIndices))]
                for matrix in mat["data"]["matrix"]:
                    if((matrix["A"].shape[0] > 1) and (np.any(matrix["A"].nonzero()[0]<0)==False)):
                        if(matrix["Name"] in mat["patient"]["StructureNames"]):
                            structIdx = mat["patient"]["StructureNames"].index(matrix["Name"])
                        else:
                            continue

                        if(matrix["A"].shape[0]==mat["patient"]["SampledVoxels"][structIdx].shape[1]):
                            SampledDosePart = matrix["A"][:,idxStart:idxEnd]*mat["solutionX"][idxStart:idxEnd] + matrix["b"]
                        elif (((matrix["A"].shape[0] % 9)==0) and (matrix["A"].shape[0]/9==mat["patient"]["SampledVoxels"][structIdx].shape[1])):
                            SampledDosePart = matrix["A"][:int(matrix["A"].shape[0]/9)][:,idxStart:idxEnd]*mat["solutionX"][idxStart:idxEnd] + matrix["b"]
                        else:
                            raise ValueError("")

                        SampledIndices.extend(mat["patient"]["SampledVoxels"][structIdx].T)
                        SampledDoses.extend(SampledDosePart)
                SampledIndices = np.vstack(SampledIndices)
                SampledDoses = np.array(SampledDoses)

                dose = np.zeros_like(subct)
                dose[mask] = griddata(SampledIndices, SampledDoses, IndicesToCalculate[mask])
                cpdoseds.PixelData = np.swapaxes(dose, 2, 0).flatten().tobytes()

                cpdoseds.save_as(outFolder+'rtdose_beam'+str(controlPointNumber[0])+'_CP'+str(controlPointNumber[1])+'.dcm', write_like_original = False)

        for beamSpotNumber in args.DoseBeamSpots:
                print('Working on rtdose, Beam:'+str(beamSpotNumber[0])+' ControlPoint:'+str(beamSpotNumber[1])+' BeamSpot:'+str(beamSpotNumber[2])+'...')

                bsdoseds = copy.deepcopy(doseds)
                bsdoseds.SOPInstanceUID = pydicom.uid.generate_uid()
                bsdoseds.SeriesInstanceUID = pydicom.uid.generate_uid()
                bsdoseds.SeriesDescription = "solutionX Beam" + str(beamSpotNumber[0]) + " CP" + str(beamSpotNumber[1]) + " SP" + str(beamSpotNumber[2])
                bsdoseds.DoseSummationType = "CONTROL_POINT"
                assert(len(bsdoseds.ReferencedRTPlanSequence) == 1)
                bsdoseds.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence = Sequence()
                rfgs = Dataset()
                rfgs.ReferencedBeamSequence = Sequence()
                rfgs.ReferencedBeamSequence.append(Dataset())
                rfgs.ReferencedBeamSequence[0].ReferencedBeamNumber = beamSpotNumber[0]
                rcps = Sequence()
                rcpp = Dataset()
                rcpp.ReferencedStartControlPointIndex = beamSpotNumber[1]
                rcpp.ReferencedStopControlPointIndex = beamSpotNumber[1]
                rcps.append(rcpp)
                rfgs.ReferencedBeamSequence[0].ReferencedControlPointSequence = rcps
                rfgs.ReferencedFractionGroupNumber = 1
                bsdoseds.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence.append(rfgs)

                if(int(beamSpotNumber[1]) % 2 == 0):
                    print("NOTE: this control point contains no dose, so skipping.")
                    continue

                idx = currentbeamlist[beamSpotNumber[0]-1]["ControlPoints"][int((beamSpotNumber[1]-1)/2)]["FileIndexStart"] + beamSpotNumber[2] - 1

                SampledIndices = [np.array([[db[0,0], db[1,0], db[2,0]]]),
                                  np.array([[db[0,1], db[1,0], db[2,0]]]),
                                  np.array([[db[0,1], db[1,1], db[2,0]]]),
                                  np.array([[db[0,0], db[1,1], db[2,0]]]),
                                  np.array([[db[0,0], db[1,0], db[2,1]]]),
                                  np.array([[db[0,1], db[1,1], db[2,1]]]),
                                  np.array([[db[0,1], db[1,0], db[2,1]]]),
                                  np.array([[db[0,0], db[1,1], db[2,1]]])]
                SampledDoses = [0.0 for i in range(len(SampledIndices))]
                for matrix in mat["data"]["matrix"]:
                    if((matrix["A"].shape[0] > 1) and (np.any(matrix["A"].nonzero()[0]<0)==False)):
                        if(matrix["Name"] in mat["patient"]["StructureNames"]):
                            structIdx = mat["patient"]["StructureNames"].index(matrix["Name"])
                        else:
                            continue

                        if(matrix["A"].shape[0]==mat["patient"]["SampledVoxels"][structIdx].shape[1]):
                            SampledDosePart = matrix["A"][:,idx]*mat["solutionX"][idx] + matrix["b"]
                        elif (((matrix["A"].shape[0] % 9)==0) and (matrix["A"].shape[0]/9==mat["patient"]["SampledVoxels"][structIdx].shape[1])):
                            SampledDosePart = matrix["A"][:int(matrix["A"].shape[0]/9)][:,idx]*mat["solutionX"][idx] + matrix["b"]
                        else:
                            raise ValueError("")

                        SampledIndices.extend(mat["patient"]["SampledVoxels"][structIdx].T)
                        SampledDoses.extend(SampledDosePart)
                SampledIndices = np.vstack(SampledIndices)
                SampledDoses = np.array(SampledDoses)

                dose = np.zeros_like(subct)
                dose[mask] = griddata(SampledIndices, SampledDoses, IndicesToCalculate[mask])
                bsdoseds.PixelData = np.swapaxes(dose, 2, 0).flatten().tobytes()

                bsdoseds.save_as(outFolder+'rtdose_beam'+str(beamSpotNumber[0])+'_CP'+str(beamSpotNumber[1])+'_SP'+str(beamSpotNumber[2])+'.dcm', write_like_original = False)

        print('DICOM files writen to ' + outFolder)
