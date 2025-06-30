import os
import mat73
import numpy as np
import copy
import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.valuerep import format_number_as_ds
import csv
import seaborn as sns

pydicom.config.settings.writing_validation_mode = pydicom.config.RAISE

caseFolders = ['Prostate_CK', 'Head-and-Neck', 'Protons', 'Liver', 'Prostate_BT', 'Prostate_VMAT', 'Head-and-Neck-Alt']
for folder in caseFolders:
    try:
        matFiles = [f.name for f in os.scandir('./'+folder)]
    except:
        continue
    for matFile in matFiles:
        print('Folder and file:', folder, matFile)
        patientFolder = matFile.split('.')[0]
        patientIndexInt = int(patientFolder.split('_')[1])
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

        rds.save_as("./DICOMs/"+patientFolder+'/rtstruct.dcm', write_like_original = False)
        
        if(folder == 'Protons'):
            beamlistfolder = mat73.loadmat('./'+folder+'/BeamList.mat')
            
            machinedata = mat73.loadmat('./'+folder+'/MachineData.mat')
            beamSigmas = {}
            for i in range(machinedata["BeamInfo"]["Sigma"].shape[0]):
                if(machinedata["BeamInfo"]["EnergyRangeList"][i][0]<machinedata["BeamInfo"]["EnergyRangeList"][i][1]):
                    beamSigmas[tuple([machinedata["BeamInfo"]["EnergyRangeList"][i][0],machinedata["BeamInfo"]["EnergyRangeList"][i][1]])] = machinedata["BeamInfo"]["Sigma"][i]
                else:
                    beamSigmas[tuple([machinedata["BeamInfo"]["EnergyRangeList"][i][1],machinedata["BeamInfo"]["EnergyRangeList"][i][0]])] = machinedata["BeamInfo"]["Sigma"][i]
            SAD = [machinedata["BeamInfo"]["SAD"][0],machinedata["BeamInfo"]["SAD"][1]]
            
            # write rtplan
            print('Working on rtplan...')
            meta = pydicom.Dataset()
            SOPInstanceUID = pydicom.uid.generate_uid()
            meta.MediaStorageSOPClassUID = pydicom.uid.RTPlanStorage
            meta.MediaStorageSOPInstanceUID = SOPInstanceUID
            meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian  
            
            ds = Dataset()
            ds.file_meta = meta
            ds.SOPClassUID = pydicom.uid.RTPlanStorage
            ds.SOPInstanceUID = SOPInstanceUID
            ds.StudyDate = "20230308"
            ds.SeriesDate = "20230308"
            ds.StudyTime = "104455"
            ds.SeriesTime = "104455"
            ds.AccessionNumber = ""
            ds.Modality = "RTPLAN"
            ds.Manufacturer = ""
            ds.InstitutionName = ""
            ds.ReferringPhysicianName = ""
            ds.OperatorsName = ""
            ds.ManufacturerModelName    = ""
            ds.PatientName = str(patientIndexInt) + '^' + folder
            ds.PatientID = "0" * (6-len(str(patientIndexTot))) + str(patientIndexTot)
            ds.PatientBirthDate = ""
            ds.PatientSex = ""
            ds.PatientIdentityRemoved = "YES"
            ds.PatientAge = ""
            ds.SoftwareVersions         = "pydicom3.0.1"
            ds.StudyInstanceUID = StudyInstanceUID
            ds.SeriesInstanceUID = SeriesInstanceUID
            ds.StudyID = "0" * (6-len(str(patientIndexTot))) + str(patientIndexTot)
            ds.SeriesNumber = 1
            ds.InstanceNumber = 1 
            ds.FrameOfReferenceUID = FORUID 
            ds.PositionReferenceIndicator = ''
            ds.RTPlanDate = ""
            ds.RTPlanTime = ""
            ds.RTPlanGeometry = "PATIENT"
            
            ds.PatientSetupSequence = Sequence()
            patientSetup = Dataset()
            patientSetup.PatientPosition = mat["patient"]["PatientPosition"]
            patientSetup.PatientSetupNumber = 1
            patientSetup.PatientSetupLabel = "Standard"
            ds.PatientSetupSequence.append(patientSetup)
            
            ds.IonBeamSequence = Sequence()
            currentbeamlist = []
            beaminfo = {}
            beaminfo["BeamNumber"] = 1
            beaminfo["FileBeamNumber"] = int(beamlistfolder["BeamList"][int(patientFolder[-2:])-1][0][0])
            beaminfo["RangeShifter"] = beamlistfolder["BeamList"][int(patientFolder[-2:])-1][0][4]
            beaminfo["FinalCumulativeMetersetWeight"] = 0
            beaminfo["ControlPoints"] = []
            controlpointinfo = {}
            controlpointinfo["ControlPointNumber"] = 0
            controlpointinfo["BeamEnergy"] = beamlistfolder["BeamList"][int(patientFolder[-2:])-1][0][1]
            controlpointinfo["MetersetWeights"] = []
            controlpointinfo["ScanSpotPositions"] = []
            controlpointinfo["CumulativeMetersetWeight"] = 0
            for rowindex in range(len(beamlistfolder["BeamList"][int(patientFolder[-2:])-1])):
                row = beamlistfolder["BeamList"][int(patientFolder[-2:])-1][rowindex]
                if(beaminfo["FileBeamNumber"] ==row[0] and beaminfo["RangeShifter"] == row[4]):
                    if(controlpointinfo["BeamEnergy"]==row[1]):
                        controlpointinfo["ScanSpotPositions"].extend(beamlistfolder["BeamList"][int(patientFolder[-2:])-1][rowindex][2:4])
                        controlpointinfo["MetersetWeights"].append(mat["solutionX"][rowindex])
                    else:
                        beaminfo["ControlPoints"].append(copy.deepcopy(controlpointinfo))
                        controlpointinfo["CumulativeMetersetWeight"] += sum(controlpointinfo["MetersetWeights"])
                        controlpointinfo["ControlPointNumber"] += 1
                        controlpointinfo["BeamEnergy"] = row[1]
                        controlpointinfo["MetersetWeights"] = [mat["solutionX"][rowindex]]
                        controlpointinfo["ScanSpotPositions"] = beamlistfolder["BeamList"][int(patientFolder[-2:])-1][rowindex][2:4].tolist()
                else:
                    beaminfo["FinalCumulativeMetersetWeight"] = controlpointinfo["CumulativeMetersetWeight"] + sum(controlpointinfo["MetersetWeights"])
                    beaminfo["ControlPoints"].append(copy.deepcopy(controlpointinfo))
                    currentbeamlist.append(copy.deepcopy(beaminfo))
                    beaminfo["BeamNumber"] = beaminfo["BeamNumber"] + 1
                    beaminfo["FileBeamNumber"] = int(row[0])
                    beaminfo["RangeShifter"] = row[4]
                    beaminfo["FinalCumulativeMetersetWeight"] = 0
                    beaminfo["ControlPoints"] = []
                    controlpointinfo["CumulativeMetersetWeight"] = 0
                    controlpointinfo["ControlPointNumber"] += 1
                    controlpointinfo["BeamEnergy"] = row[1]
                    controlpointinfo["MetersetWeights"] = [mat["solutionX"][rowindex]]
                    controlpointinfo["ScanSpotPositions"] = beamlistfolder["BeamList"][int(patientFolder[-2:])-1][rowindex][2:4].tolist()
            currentbeamlist.append(beaminfo) 
            
            ds.FractionGroupSequence = Sequence()
            fractionds = Dataset()
            fractionds.FractionGroupNumber = 1
            fractionds.NumberOfFractionsPlanned = 1
            fractionds.NumberOfBeams = len(currentbeamlist)
            fractionds.NumberOfBrachyApplicationSetups = 0
            fractionds.BeamSequence = Sequence()
            for beaminfo in currentbeamlist:
                be = Dataset()
                be.BeamMeterset = beaminfo["FinalCumulativeMetersetWeight"]
                be.BeamNumber = beaminfo["BeamNumber"]
                fractionds.BeamSequence.append(be)
            ds.FractionGroupSequence.append(fractionds)
            
            for beaminfo in currentbeamlist:
                be = Dataset()
                be.Manufacturer = ""
                be.InstitutionName  = ""
                be.ManufacturerModelName  = ""
                be.TreatmentMachineName   = ""
                be.InstitutionAddress = ""
                be.PrimaryDosimeterUnit = "MU"
                be.BeamNumber             = beaminfo["BeamNumber"]
                be.BeamName               = ""
                be.BeamType               = 'STATIC'
                be.RadiationType          = "PROTON"
                be.TreatmentDeliveryType  = 'TREATMENT'
                be.NumberOfWedges         = 0
                be.NumberOfCompensators   = 0
                be.NumberOfBoli           = 0
                be.NumberOfBlocks         = 0
                be.FinalCumulativeMetersetWeight = beaminfo["FinalCumulativeMetersetWeight"]
                be.ScanMode                   = 'MODULATED'
                be.SourceAxisDistance = SAD
                if(beaminfo["RangeShifter"] ==0):
                    be.NumberOfRangeShifters = 0
                else:
                    be.NumberOfRangeShifters = 1
                    be.RangeShifterSequence = Sequence()
                    rsDataset = Dataset()
                    rsDataset.AccessoryCode = "Undefined Accessory Code"
                    rsDataset.RangeShifterNumber = 1
                    rsDataset.RangeShifterID = "Rs " + str(beaminfo["RangeShifter"]) + " mm"
                    rsDataset.RangeShifterType = "BINARY"
                    be.RangeShifterSequence.append(rsDataset)
                be.NumberOfControlPoints  = 2*len(beaminfo["ControlPoints"])
                be.IonControlPointSequence = Sequence()
                be.NumberOfLateralSpreadingDevices = 0
                be.NumberOfRangeModulators = 0
                be.PatientSupportType = 'TABLE'
                
                for controlpointinfo in beaminfo["ControlPoints"]:
                    icpoi = Dataset()
                    icpoi.NominalBeamEnergyUnit = 'MEV'
                    icpoi.ControlPointIndex = 2*controlpointinfo["ControlPointNumber"]
                    icpoi.NominalBeamEnergy = controlpointinfo["BeamEnergy"]
                    icpoi.CumulativeMetersetWeight = controlpointinfo["CumulativeMetersetWeight"]
                    icpoi.GantryAngle = mat['patient']['Beams']['BeamConfig'][beaminfo["FileBeamNumber"]-1]['Gantry']
                    icpoi.GantryRotationDirection = 'NONE'
                    icpoi.BeamLimitingDeviceAngle = mat['patient']['Beams']['BeamConfig'][beaminfo["FileBeamNumber"]-1]['Collimator']
                    icpoi.BeamLimitingDeviceRotationDirection = 'NONE'
                    icpoi.PatientSupportAngle = mat['patient']['Beams']['BeamConfig'][beaminfo["FileBeamNumber"]-1]['Couch'] # I think...
                    icpoi.PatientSupportRotationDirection = 'NONE'
                    
                    icpoi.TableTopVerticalPosition     = 0
                    icpoi.TableTopLongitudinalPosition = 0
                    icpoi.TableTopLateralPosition      = 0
                    icpoi.IsocenterPosition = mat["patient"]["Isocentre"].flatten().tolist()
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
                    for (low,high),sigma in beamSigmas.items():
                        if((controlpointinfo["BeamEnergy"]>=low) and (controlpointinfo["BeamEnergy"]<=high)):
                            icpoi.ScanningSpotSize = [sigma[0],sigma[1]]
                            break
                    icpoi.NumberOfPaintings = 1
                    be.IonControlPointSequence.append(icpoi)
                    
                    icpoi = Dataset()
                    icpoi.NominalBeamEnergyUnit = 'MEV'
                    icpoi.ControlPointIndex = 2*controlpointinfo["ControlPointNumber"] + 1
                    icpoi.NominalBeamEnergy = controlpointinfo["BeamEnergy"]
                    icpoi.CumulativeMetersetWeight = controlpointinfo["CumulativeMetersetWeight"]+ sum(controlpointinfo["MetersetWeights"])
                    icpoi.ScanSpotTuneID = 'TROTS_1.0'
                    icpoi.NumberOfScanSpotPositions = len(controlpointinfo["MetersetWeights"])
                    icpoi.ScanSpotPositionMap = controlpointinfo["ScanSpotPositions"]
                    icpoi.ScanSpotMetersetWeights = [0.0 for i in range(len(controlpointinfo["MetersetWeights"]))]
                    for (low,high),sigma in beamSigmas.items():
                        if((controlpointinfo["BeamEnergy"]>=low) and (controlpointinfo["BeamEnergy"]<=high)):
                            icpoi.ScanningSpotSize = [sigma[0],sigma[1]]
                            break
                    icpoi.NumberOfPaintings = 1
                    be.IonControlPointSequence.append(icpoi)
                be.PatientSetupNumber = 1
                be.ToleranceTableNumber = 0 
                ds.IonBeamSequence.append(be)
            
