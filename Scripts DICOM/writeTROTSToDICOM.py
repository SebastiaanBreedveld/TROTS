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

def ConvertCTIndexToRealSpace(CTPoints, resolution, offset):
    return resolution*(CTPoints - 1) + offset # -1 because matlab indexing starts at 1

pydicom.config.settings.writing_validation_mode = pydicom.config.RAISE

parser = argparse.ArgumentParser()
parser.add_argument("-b","--folderBasePath", nargs='?', help="The base directory in which the code is run containing all neccessary folders", default=".")
parser.add_argument("-o", "--outputPath", nargs='?', help="The output directory of the DICOM file", default="/tmp")
parser.add_argument("-n", "--DoseBeamNumber", type=list, nargs='?', help="A list of numbers to be calculated where a seperate rtdose_<BeamNumber>.dcm is calculated,, format: [BeamNumber_i, ..]", default=[])
parser.add_argument("-c", "--DoseControlPoints", type=list, nargs='?', help="A list of control point numbers where a seperate rtdose_<BeamNumber>_<ControlPointNumber>.dcm is calculated, format: [(BeamNumber_i,ControlPoint_i), ...]", default=[])
parser.add_argument("-s", "--DoseBeamSpots", type=list, nargs='?', help="A list of beam spot numbers where the rtdose_<BeamNumber>_<ConrolPointNumber>_<BeamSpotNumber>.dcm is calculated, format: [(BeamNumber_i,ControlPoint_i,BeamSpotNumber_i), ...]", default=[])
args = parser.parse_args()

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
        print('Folder and file:', folder, matFile)
        patientFolder = matFile.split('.')[0]
        patientIndexInt = int(patientFolder.split('_')[1])
        os.makedirs(args.outputPath + "/DICOMs/" + patientFolder, exist_ok=True)
        resolutionX = mat['patient']['Resolution'][0];
        resolutionY = mat['patient']['Resolution'][1];
        resolutionZ = mat['patient']['Resolution'][2];
        ctshape = mat['patient']['CT'].shape
        nRowsCT =  ctshape[1] # DICOM Rows, goes with y, is index 1 because of how matrix is stored
        nColumnsCT = ctshape[0] # DICOM Columns, goes with x, is index 0 because of how matrix is stored
        nSlicesCT = ctshape[2] # DICOM slices, goes with z
        ctsopids = {}
        # write planning objectives into csv
        with open(args.outputPath + "/DICOMs/"+patientFolder+'/planning.csv', 'w', newline = '') as csvfile:
            filewriter = csv.writer(csvfile)
            if folder != 'Prostate_BT':
                filewriter.writerow(['Beam angle', 'Couch angle', 'Collimator'])
                for rowIndex in range(0, len(mat['patient']['Beams']['BeamConfig'])):
                    filewriter.writerow([mat['patient']['Beams']['BeamConfig'][rowIndex]['Gantry']] + [mat['patient']['Beams']['BeamConfig'][rowIndex]['Couch']] + [mat['patient']['Beams']['BeamConfig'][rowIndex]['Collimator']])

            filewriter.writerow(['Name', 'Minimise', 'Objective', 'Sufficient', 'Weight', 'Priority', 'IsConstraint'])
            for rowIndex in range(0, len(mat['problem']['Name'])):
                filewriter.writerow([mat['problem']['Name'][rowIndex]] + [mat['problem']['Minimise'][rowIndex]] + [mat['problem']['Objective'][rowIndex]] + [mat['problem']['Sufficient'][rowIndex]] + [mat['problem']['Weight'][rowIndex]] + [mat['problem']['Priority'][rowIndex]] + [mat['problem']['IsConstraint'][rowIndex]])

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
            ds.StudyDescription = "RT optimizers"
            ds.AccessionNumber = ""
            ds.Manufacturer = ""
            ds.InstitutionName = ""
            ds.ReferringPhysicianName = ""
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

            ds.save_as(args.outputPath + "/DICOMs/"+patientFolder+'/CTSlice'+str(sliceIndex).zfill(3)+".dcm", write_like_original = False)

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
                    keyHole = False # Both options lead to the same artefacts in Slicer3D: https://discourse.slicer.org/t/closedplanar-xor-visualization-artefact-with-holes-2d/43589/3
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

        rds.save_as(args.outputPath + "/DICOMs/"+patientFolder+'/rtstruct.dcm', write_like_original = False)

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
            rtds.Manufacturer = ds.Manufacturer
            rtds.InstitutionName = ""
            rtds.ReferringPhysicianName = ""
            rtds.OperatorsName = ""
            rtds.ManufacturerModelName    = ""
            rtds.PatientName = rds.PatientName
            rtds.PatientID = rds.PatientID
            rtds.PatientBirthDate = rds.PatientBirthDate
            rtds.PatientSex = rds.PatientSex
            rtds.StudyInstanceUID = rds.StudyInstanceUID
            rtds.SeriesInstanceUID = rds.StudyInstanceUID
            rtds.StudyID = rds.StudyID
            rtds.SeriesNumber = 1
            rtds.InstanceNumber = 1
            rtds.FrameOfReferenceUID = rds.FrameOfReferenceUID
            rtds.PositionReferenceIndicator = ''
            rtds.RTPlanLabel = "TROTS solutionX"
            rtds.RTPlanDate = ""
            rtds.RTPlanTime = ""
            rtds.RTPlanGeometry = "PATIENT"

            rtds.PatientSetupSequence = Sequence()
            patientSetup = Dataset()
            patientSetup.PatientPosition = mat["patient"]["PatientPosition"]
            patientSetup.PatientSetupNumber = 1
            patientSetup.PatientSetupLabel = "Standard"
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
                # All problem info is grouped accordinding to the unique key = (Structure name, weight, is constraint)
                # If two rows have the same key they are merged one being the upper and other being the lower bound
                key = (mat["problem"]["Name"][index],float(mat["problem"]["Weight"][index]),bool(mat["problem"]["IsConstraint"][index]))
                if(key[0].lower()=="mu"):
                    continue
                if(key not in probleminfo):
                    constraint = {}
                    if (("gtv" in key[0].lower()) or ("ptv" in key[0].lower()) or ("ctv" in key[0].lower())):
                        constraint["type"] = "TARGET"
                        if(mat["problem"]["Minimise"][index]==1):
                            constraint["Min"] = mat["problem"]["Objective"][index]
                            constraint["Max"] = ""
                        else:
                            constraint["Min"] = ""
                            constraint["Max"] = mat["problem"]["Objective"][index]
                    else:
                        # If it isn't a TARGET it must be a ORGAN_AT_RISK
                        constraint["type"] = "ORGAN_AT_RISK"
                        # ORGAN_AT_RISK can only have an upper bound
                        assert(mat["problem"]["Minimise"][index]==1)
                        constraint["Max"] = mat["problem"]["Objective"][index]
                    probleminfo[key] = constraint
                else:
                    # Only TARGET can have both upper and lower bound
                    assert(probleminfo[key]["type"]=="TARGET") # having upper and lower only make sense for target
                    # Only one upper and one lower can be defined for one key in a TARGET
                    if(mat["problem"]["Minimise"]==1):
                        assert(probleminfo[key]["Max"]=="")
                        probleminfo[key]["Max"] = mat["problem"]["Objective"][index]
                    else:
                        assert(probleminfo[key]["Min"]=="")
                        probleminfo[key]["Min"] = mat["problem"]["Objective"][index]

            rtds.DoseReferenceSequence = Sequence()
            for dosenumber,key in enumerate(probleminfo):
                doseReference = Dataset()
                doseReference.ReferencedROINumber = structToROINumber[key[0]]
                doseReference.DoseReferenceNumber = dosenumber + 1
                doseReference.DoseReferenceUID = pydicom.uid.generate_uid()
                doseReference.DoseReferenceStructureType = "VOLUME"
                doseReference.DoseReferenceDescription = key[0] + (" Constraint" if(key[2]) else " Objective")
                doseReference.DoseReferenceType = probleminfo[key]["type"]
                doseReference.ConstraintWeight = format_number_as_ds(float(key[1]))
                if(probleminfo[key]["type"] == "TARGET"):
                    if(probleminfo[key]["Min"]!=""):
                        doseReference.TargetMinimumDose = format_number_as_ds(float(probleminfo[key]["Min"]))
                    if(probleminfo[key]["Max"]!=""):
                        doseReference.TargetMaximumDose = format_number_as_ds(float(probleminfo[key]["Max"]))
                else:
                    doseReference.OrganAtRiskMaximumDose = format_number_as_ds(float(probleminfo[key]["Max"]))
                rtds.DoseReferenceSequence.append(doseReference)

            rtds.FractionGroupSequence = Sequence()
            fractionds = Dataset()
            fractionds.FractionGroupNumber = 1
            fractionds.NumberOfFractionsPlanned = 1
            fractionds.NumberOfBeams = len(currentbeamlist)
            fractionds.NumberOfBrachyApplicationSetups = 0
            fractionds.BeamSequence = Sequence()
            for beaminfo in currentbeamlist:
                be = Dataset()
                be.BeamMeterset = format_number_as_ds(beaminfo["FinalCumulativeMetersetWeight"])

                be.BeamNumber = beaminfo["BeamNumber"]
                fractionds.BeamSequence.append(be)
            rtds.FractionGroupSequence.append(fractionds)

            rtds.ReferencedStructureSetSequence = Sequence()
            ref_rds = Dataset()
            ref_rds.ReferencedSOPClassUID = rds.SOPClassUID
            ref_rds.ReferencedSOPInstanceUID = rds.SOPInstanceUID
            rtds.ReferencedStructureSetSequence.append(ref_rds)

            totalMetersetWeightOfBeams = 0
            for beaminfo in currentbeamlist:
                be = Dataset()
                be.Manufacturer = ""
                be.InstitutionName  = ""
                be.ManufacturerModelName  = ""
                be.TreatmentMachineName   = ""
                be.InstitutionAddress = ""
                be.PrimaryDosimeterUnit = "MU"
                be.BeamNumber             = beaminfo["BeamNumber"]
                be.BeamName               = str(beaminfo["BeamNumber"])
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
                        rsDataset.RangeShifterType = "BINARY"
                        be.RangeShifterSequence.append(rsDataset)
                be.NumberOfControlPoints  = 2*len(beaminfo["ControlPoints"])
                be.IonControlPointSequence = Sequence()
                be.NumberOfLateralSpreadingDevices = 0
                be.NumberOfRangeModulators = 0
                be.PatientSupportType = 'TABLE'

                MetersetWeightTolerance = 1e-8
                totalMetersetWeightOfControlPoints = 0
                for controlpointinfo in beaminfo["ControlPoints"]:
                    icpoi = Dataset()
                    icpoi.NominalBeamEnergyUnit = 'MEV'
                    icpoi.ControlPointIndex = 2*controlpointinfo["ControlPointNumber"]
                    icpoi.NominalBeamEnergy = format_number_as_ds(controlpointinfo["BeamEnergy"])
                    icpoi.CumulativeMetersetWeight = format_number_as_ds(float(controlpointinfo["CumulativeMetersetWeight"]))
                    icpoi.GantryAngle = mat['patient']['Beams']['BeamConfig'][beaminfo["FileBeamNumber"]-1]['Gantry']
                    icpoi.GantryRotationDirection = 'NONE'
                    icpoi.BeamLimitingDeviceAngle = mat['patient']['Beams']['BeamConfig'][beaminfo["FileBeamNumber"]-1]['Collimator']
                    icpoi.BeamLimitingDeviceRotationDirection = 'NONE'
                    icpoi.PatientSupportAngle = mat['patient']['Beams']['BeamConfig'][beaminfo["FileBeamNumber"]-1]['Couch']
                    icpoi.PatientSupportRotationDirection = 'NONE'
                    icpoi.TableTopVerticalPosition     = 0
                    icpoi.TableTopLongitudinalPosition = 0
                    icpoi.TableTopLateralPosition      = 0
                    icpoi.IsocenterPosition = [format_number_as_ds(isp) for isp in mat["patient"]["Isocentre"].flatten()]
                    icpoi.TableTopPitchAngle = 0
                    icpoi.TableTopPitchRotationDirection = 'NONE'
                    icpoi.TableTopRollAngle = 0
                    icpoi.TableTopRollRotationDirection = 'NONE'
                    icpoi.GantryPitchAngle = 0
                    icpoi.GantryPitchRotationDirection = 'NONE'
                    icpoi.SnoutPosition = None
                    icpoi.ScanSpotTuneID = 'TROTS_1.0'
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
                    be.IonControlPointSequence.append(icpoi)

                    icpoi = Dataset()
                    icpoi.NominalBeamEnergyUnit = 'MEV'
                    icpoi.ControlPointIndex = 2*controlpointinfo["ControlPointNumber"] + 1
                    icpoi.NominalBeamEnergy = format_number_as_ds(controlpointinfo["BeamEnergy"])
                    icpoi.CumulativeMetersetWeight = format_number_as_ds(controlpointinfo["CumulativeMetersetWeight"]+ sum(controlpointinfo["MetersetWeights"]))
                    icpoi.ScanSpotTuneID = 'TROTS_1.0'
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
                be.PatientSetupNumber = 1
                be.ToleranceTableNumber = 0
                assert(abs(be.FinalCumulativeMetersetWeight - icpoi.CumulativeMetersetWeight) < MetersetWeightTolerance)
                assert(abs(be.FinalCumulativeMetersetWeight - totalMetersetWeightOfControlPoints) < MetersetWeightTolerance)
                rtds.IonBeamSequence.append(be)
                totalMetersetWeightOfBeams += be.FinalCumulativeMetersetWeight
            assert(abs(totalMetersetWeightOfBeams - sum(mat["solutionX"])) < MetersetWeightTolerance)
            rtds.save_as(args.outputPath + "/DICOMs/"+patientFolder+'/rtplan.dcm', write_like_original = False)

            # write rtplan
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
            doseds.SeriesDescription = "solutionX"
            doseds.StudyTime = rds.StudyTime
            doseds.AccessionNumber = rds.AccessionNumber
            doseds.Modality = "RTDOSE"
            doseds.Manufacturer = ds.Manufacturer
            doseds.InstitutionName = ""
            doseds.ReferringPhysicianName = ""
            doseds.OperatorsName = ""
            doseds.ManufacturerModelName    = ""
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
            doseds.ImagePositionPatient = [format_number_as_ds(mat['patient']['Offset'][0]), format_number_as_ds(mat['patient']['Offset'][1]), format_number_as_ds(mat['patient']['Offset'][2])]
            doseds.ImageOrientationPatient = ds.ImageOrientationPatient
            doseds.FrameOfReferenceUID = ds.FrameOfReferenceUID
            doseds.SamplesPerPixel = ds.SamplesPerPixel
            doseds.PhotometricInterpretation = ds.PhotometricInterpretation
            doseds.FrameIncrementPointer = Tag(0x3004,0x000c)
            doseds.Rows = ds.Rows
            doseds.Columns = ds.Columns
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
            
            resolution = mat["patient"]["Resolution"]
            offset = mat["patient"]["Offset"]
            
            CTIndexDoseEdges = np.stack([mat["patient"]["DoseBox"][:, 0] - 1, mat["patient"]["DoseBox"][:, 1] + 1])
            DoseEdges = ConvertCTIndexToRealSpace(CTIndexDoseEdges, resolution, offset)
            SampledPoints = [np.array([[DoseEdges[0,0], DoseEdges[0,1],DoseEdges[0,2]]]),
                             np.array([[DoseEdges[1,0], DoseEdges[0,1],DoseEdges[0,2]]]),
                             np.array([[DoseEdges[1,0], DoseEdges[1,1],DoseEdges[0,2]]]),
                             np.array([[DoseEdges[1,0], DoseEdges[1,1],DoseEdges[1,2]]]),
                             np.array([[DoseEdges[0,0], DoseEdges[0,1],DoseEdges[1,2]]]),
                             np.array([[DoseEdges[0,0], DoseEdges[1,1],DoseEdges[1,2]]]),
                             np.array([[DoseEdges[0,0], DoseEdges[1,1],DoseEdges[0,2]]]),
                             np.array([[DoseEdges[1,0], DoseEdges[0,1],DoseEdges[1,2]]])]
            SampledDoses = [0.0 for i in range(len(SampledPoints))]
            
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
            
                    SampledPoints.append(ConvertCTIndexToRealSpace(mat["patient"]["SampledVoxels"][structIdx].T, resolution, offset))
                    SampledDoses.extend(SampledDosePart)
            SampledPoints = np.vstack(SampledPoints)
            SampledDoses = np.array(SampledDoses)
            
            mask = np.zeros_like(mat["patient"]["CT"])
            db=mat["patient"]["DoseBox"]
            mask[db[0, 0]-1:db[0, 1],db[1, 0]-1:db[1, 1],db[2, 0]-1:db[2, 1]] = True
            mask[np.where(mat["patient"]["CT"]<=-1024)] = False
            
            CTIndicesToCalculate = np.vstack(np.where(mask)).T + 1 # matlab indices start at 1
            PointsToCalculate = ConvertCTIndexToRealSpace(CTIndicesToCalculate, resolution, offset)

            dose = np.zeros_like(mat["patient"]["CT"])
            dose[mask==1] = griddata(SampledPoints, SampledDoses, PointsToCalculate, method="linear")
            doseds.PixelData = np.moveaxis(np.moveaxis(dose, 2, 0), 2, 1).flatten().tobytes()

            doseds.NumberOfFrames = dose.shape[2]
            doseds.GridFrameOffsetVector = [format_number_as_ds(ConvertCTIndexToRealSpace(np.array([0,0,i+1]), resolution, offset)[2]) for i in range(dose.shape[2])]

            doseds.save_as(args.outputPath + "/DICOMs/"+patientFolder+'/rtdose.dcm', write_like_original = False)
            
            for beamNumber in args.DoseBeamNumber:
                print('Working on rtplan, Beam:'+str(beamNumber)+'...')
                bdoseds = copy.deepcopy(doseds)
                bdoseds.DoseSummationType = "BEAM"
                idxStart = currentbeamlist[beamNumber - 1]["FileIndexStart"]
                idxEnd = currentbeamlist[beamNumber - 1]["FileIndexEnd"]
                
                SampledPoints = [np.array([[DoseEdges[0,0], DoseEdges[0,1],DoseEdges[0,2]]]),
                             np.array([[DoseEdges[1,0], DoseEdges[0,1],DoseEdges[0,2]]]),
                             np.array([[DoseEdges[1,0], DoseEdges[1,1],DoseEdges[0,2]]]),
                             np.array([[DoseEdges[1,0], DoseEdges[1,1],DoseEdges[1,2]]]),
                             np.array([[DoseEdges[0,0], DoseEdges[0,1],DoseEdges[1,2]]]),
                             np.array([[DoseEdges[0,0], DoseEdges[1,1],DoseEdges[1,2]]]),
                             np.array([[DoseEdges[0,0], DoseEdges[1,1],DoseEdges[0,2]]]),
                             np.array([[DoseEdges[1,0], DoseEdges[0,1],DoseEdges[1,2]]])]
                SampledDoses = [0.0 for i in range(len(SampledPoints))]
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

                        SampledPoints.append(ConvertCTIndexToRealSpace(mat["patient"]["SampledVoxels"][structIdx].T, resolution, offset))
                        SampledDoses.extend(SampledDosePart)
                SampledPoints = np.vstack(SampledPoints)
                SampledDoses = np.array(SampledDoses)
                dose = np.zeros_like(mat["patient"]["CT"])
                dose[mask==1] = griddata(SampledPoints, SampledDoses, PointsToCalculate, method="linear")
                bdoseds.PixelData = np.moveaxis(np.moveaxis(dose, 2, 0), 2, 1).flatten().tobytes()

                bdoseds.save_as(args.outputPath + "/DICOMs/"+patientFolder+'/rtdose_'+str(beamNumber)+'.dcm', write_like_original = False)
            
            for controlPointNumber in args.DoseControlPoints:
                print('Working on rtplan, Beam:'+str(controlPointNumber[0])+' ControlPoint:'+str(controlPointNumber[1])+'...')
                cpdoseds = copy.deepcopy(doseds)
                cpdoseds.DoseSummationType = "CONTROL_POINT"
                if(int(controlPointNumber[1]) % 2 == 0):
                    print("NOTE: this control point contains no dose")
                    dose = np.zeros_like(mat["patient"]["CT"])
                else:
                    idxStart = currentbeamlist[controlPointNumber[0]-1]["ControlPoints"][int((controlPointNumber[1]-1)/2)]["FileIndexStart"]
                    idxEnd = currentbeamlist[controlPointNumber[0]-1]["ControlPoints"][int((controlPointNumber[1]-1)/2)]["FileIndexEnd"]

                    SampledPoints = [np.array([[DoseEdges[0,0], DoseEdges[0,1],DoseEdges[0,2]]]),
                                 np.array([[DoseEdges[1,0], DoseEdges[0,1],DoseEdges[0,2]]]),
                                 np.array([[DoseEdges[1,0], DoseEdges[1,1],DoseEdges[0,2]]]),
                                 np.array([[DoseEdges[1,0], DoseEdges[1,1],DoseEdges[1,2]]]),
                                 np.array([[DoseEdges[0,0], DoseEdges[0,1],DoseEdges[1,2]]]),
                                 np.array([[DoseEdges[0,0], DoseEdges[1,1],DoseEdges[1,2]]]),
                                 np.array([[DoseEdges[0,0], DoseEdges[1,1],DoseEdges[0,2]]]),
                                 np.array([[DoseEdges[1,0], DoseEdges[0,1],DoseEdges[1,2]]])]
                    SampledDoses = [0.0 for i in range(len(SampledPoints))]
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

                            SampledPoints.append(ConvertCTIndexToRealSpace(mat["patient"]["SampledVoxels"][structIdx].T, resolution, offset))
                            SampledDoses.extend(SampledDosePart)
                    SampledPoints = np.vstack(SampledPoints)
                    SampledDoses = np.array(SampledDoses)
                    dose = np.zeros_like(mat["patient"]["CT"])
                    dose[mask==1] = griddata(SampledPoints, SampledDoses, PointsToCalculate, method="linear")
                    cpdoseds.PixelData = np.moveaxis(np.moveaxis(dose, 2, 0), 2, 1).flatten().tobytes()

                cpdoseds.save_as(args.outputPath + "/DICOMs/"+patientFolder+'/rtdose_'+str(controlPointNumber[0])+'_'+str(controlPointNumber[1])+'.dcm', write_like_original = False)

        for beamSpotNumber in args.DoseBeamSpots:
                print('Working on rtplan, Beam:'+str(beamSpotNumber[0])+' ControlPoint:'+str(beamSpotNumber[1])+' BeamSpot:'+str(beamSpotNumber[2])+'...')
                bsdoseds = copy.deepcopy(doseds)
                bsdoseds.DoseSummationType = "CONTROL_POINT"
                if(int(beamSpotNumber[1]) % 2 == 0):
                    print("NOTE: this control point contains no dose")
                    dose = np.zeros_like(mat["patient"]["CT"])
                else:
                    idx = currentbeamlist[beamSpotNumber[0]-1]["ControlPoints"][int((beamSpotNumber[1]-1)/2)]["FileIndexStart"] + beamSpotNumber[2] - 1

                    SampledPoints = [np.array([[DoseEdges[0,0], DoseEdges[0,1],DoseEdges[0,2]]]),
                                 np.array([[DoseEdges[1,0], DoseEdges[0,1],DoseEdges[0,2]]]),
                                 np.array([[DoseEdges[1,0], DoseEdges[1,1],DoseEdges[0,2]]]),
                                 np.array([[DoseEdges[1,0], DoseEdges[1,1],DoseEdges[1,2]]]),
                                 np.array([[DoseEdges[0,0], DoseEdges[0,1],DoseEdges[1,2]]]),
                                 np.array([[DoseEdges[0,0], DoseEdges[1,1],DoseEdges[1,2]]]),
                                 np.array([[DoseEdges[0,0], DoseEdges[1,1],DoseEdges[0,2]]]),
                                 np.array([[DoseEdges[1,0], DoseEdges[0,1],DoseEdges[1,2]]])]
                    SampledDoses = [0.0 for i in range(len(SampledPoints))]
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

                            SampledPoints.append(ConvertCTIndexToRealSpace(mat["patient"]["SampledVoxels"][structIdx].T, resolution, offset))
                            SampledDoses.extend(SampledDosePart)
                    SampledPoints = np.vstack(SampledPoints)
                    SampledDoses = np.array(SampledDoses)
                    dose = np.zeros_like(mat["patient"]["CT"])
                    dose[mask==1] = griddata(SampledPoints, SampledDoses, PointsToCalculate, method="linear")
                    bsdoseds.PixelData = np.moveaxis(np.moveaxis(dose, 2, 0), 2, 1).flatten().tobytes()

                bsdoseds.save_as(args.outputPath + "/DICOMs/"+patientFolder+'/rtdose_'+str(beamSpotNumber[0])+'_'+str(beamSpotNumber[1])+'_'+str(beamSpotNumber[2])+'.dcm', write_like_original = False)
        
        print('DICOM files writen to ' + args.outputPath + "/DICOMs/"+patientFolder)
