import os
import numpy as np
import pydicom
from pydicom.valuerep import format_number_as_ds
from pydicom.tag import Tag
import pydicom.uid
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--BinFile",nargs='?',help="Full path of the .bin file outputted by TOPAS dose scorer.",default="")
parser.add_argument("--DICOMDirectory", nargs='?',help="Directory containing TROTS Protons_nn CT Files.",default="")
parser.add_argument("--OutputDirectory", nargs='?',help="Directory where converted RTDOSE will be saved. By default, it will be the same as BinFile directory. Set to DICOMDir if you want to be in the same place as CT",default=None)
parser.add_argument("-d", "--ParticlesperHistory", nargs='?', help="Decimation factor of particles in RTplan used when TOPAS rtion simulation was run, i.e. by how much the number of histories is scaled down to have a reasonable computation time", default=1)
parser.add_argument("--DoseSummationType", nargs='?', help="Type of dose summation: PLAN, BEAM, CONTROL_POINT... ", default="PLAN")

args = parser.parse_args()

if args.OutputDirectory is None:
    args.OutputDirectory=os.path.dirname(os.path.abspath(args.BinFile))

pph=float(args.ParticlesperHistory)

ct_files = sorted([
    os.path.join(args.DICOMDirectory, f)
    for f in os.listdir(args.DICOMDirectory)
    if f.endswith(".dcm") and f.startswith("CT_")
])
ct_slices = [pydicom.dcmread(f) for f in ct_files]
ct_slices.sort(key=lambda x: float(x.ImagePositionPatient[2]))
ct0 = ct_slices[0]
dx,dy=ct0.PixelSpacing
dz=ct0.SliceThickness
nx=ct0.Columns
ny=ct0.Rows
nz=len(ct_slices)
dose = np.fromfile(args.BinFile, dtype="<d")
dose *= pph
expected = nx * ny * nz

if len(dose) != expected:
    raise ValueError(f"Wrong .bin size: expected {expected}, found {len(dose)}. Maybe you did not run writeTROTStoDICOM with doseBoxLikeCT set ?")
dose = dose.reshape((nz, ny, nx))


meta = pydicom.Dataset()
SOPInstanceUID = pydicom.uid.generate_uid()
meta.MediaStorageSOPClassUID = pydicom.uid.RTDoseStorage
meta.MediaStorageSOPInstanceUID = SOPInstanceUID
meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian

doseds = Dataset()
doseds.file_meta = meta
doseds.SOPClassUID = pydicom.uid.RTDoseStorage
doseds.SOPInstanceUID = SOPInstanceUID
doseds.StudyDate = ct0.StudyDate
doseds.SeriesDate = ct0.StudyDate
doseds.SeriesDescription = f"TOPAS {args.DoseSummationType}"
doseds.StudyTime = ct0.StudyTime
doseds.Modality = "RTDOSE"

doseds.PatientName = ct0.PatientName
doseds.PatientID = ct0.PatientID

doseds.SliceThickness = str(dz)
doseds.StudyInstanceUID = ct0.StudyInstanceUID
doseds.SeriesInstanceUID = pydicom.uid.generate_uid()
doseds.SeriesNumber = 1
doseds.InstanceNumber = 1

doseds.ImagePositionPatient = ct0.ImagePositionPatient
doseds.ImageOrientationPatient = ct0.ImageOrientationPatient
doseds.FrameOfReferenceUID = ct0.FrameOfReferenceUID
doseds.SamplesPerPixel = 1
doseds.PhotometricInterpretation = "MONOCHROME2"
doseds.FrameIncrementPointer = Tag(0x3004, 0x000C)

doseds.Rows = ny
doseds.Columns = nx
doseds.NumberOfFrames = nz

# Case a of DICOM C.8.8.3.2 (Offsets relativos en Z)
doseds.GridFrameOffsetVector = [format_number_as_ds(i * dz) for i in range(nz)]
doseds.PixelSpacing = [str(dx), str(dy)]

doseds.BitsAllocated = 16
doseds.BitsStored = 16
doseds.HighBit = 15
doseds.PixelRepresentation = 0
doseds.DoseUnits = "GY"
doseds.DoseType = "PHYSICAL"
doseds.DoseSummationType = args.DoseSummationType


max_dose = np.max(dose)
scaling = max_dose / 65535 if max_dose > 0 else 1.0
doseds.DoseGridScaling = format_number_as_ds(scaling)

dose_uint16 = (dose / scaling).astype(np.uint16)
doseds.PixelData = dose_uint16.flatten().tobytes()

doseds.ReferencedRTPlanSequence = Sequence()
seqrt = Dataset()
seqrt.ReferencedSOPClassUID = pydicom.uid.RTDoseStorage 
seqrt.ReferencedSOPInstanceUID = SOPInstanceUID
doseds.ReferencedRTPlanSequence.append(seqrt)

inputFileNameWithoutExtension = os.path.splitext(os.path.basename(os.path.abspath(args.BinFile)))[0]
outputFileName = args.OutputDirectory+ '/' + inputFileNameWithoutExtension  + '.dcm'
print(f"RTDOSE saved in: {outputFileName}")
if int(pydicom.__version_info__[0]) >= 3:
    doseds.save_as(outputFileName, enforce_file_format=True)
else:
    doseds.save_as(outputFileName, write_like_original = False)
