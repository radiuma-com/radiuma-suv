# =============================================================================
# radiuma-suv: PET TO SUVbw CONVERSION
# =============================================================================
#
# This program converts PET images into body-weight normalized Standardized
# Uptake Value (SUVbw) images using PET pixel values and DICOM metadata.
# It supports conventional single-frame and Enhanced multi-frame PET images.
#
# SUPPORTED INPUTS AND OUTPUTS:
#
# 1. Input:  NumPy array + corresponding PET DICOM metadata dictionary
#    Output: Floating-point SUVbw NumPy array + metadata dictionary (in-memory)
#
# 2. Input:  Single-frame PET DICOM file
#    Output: New single-frame SUVbw DICOM file
#
# 3. Input:  Multi-frame Enhanced PET DICOM file
#    Output: New multi-frame Enhanced PET SUVbw DICOM file
#
# 4. Input:  Folder containing single-frame PET DICOM images
#    Output: Folder containing corresponding SUVbw DICOM images
#
# 5. Input:  Folder containing multi-frame Enhanced PET DICOM images
#    Output: Folder containing corresponding Enhanced PET SUVbw DICOM images
#
# Multiple independent PET series are processed separately.
# Original DICOM files are not overwritten.
#
# PUBLIC API:
#
# radiuma_suv(input_array, input_metadata, parameter_estimation="OFF")
#     Converts an in-memory PET NumPy array to SUVbw.
#
# radiuma_suv_from_input(image_input, metadata_json=None,
#                        parameter_estimation="OFF", output_path=None)
#     Converts supported PET inputs and produces corresponding SUVbw outputs.
#
# RadiumaSUVLoader(parameter_estimation="OFF", pet_metadata=None,
#                  output_path=None)
#     Provides an object-oriented interface using .convert(image_input).
#
# PARAMETER ESTIMATION:
#     OFF: Use existing DICOM metadata without filling missing parameters.
#     ON:  Fill eligible missing metadata using configured default values.
#
# DICOM outputs use calibrated 16-bit pixels to represent SUVbw values.
# NumPy outputs preserve floating-point SUVbw values in memory.

# ============================================================
# radiuma-suv - INSTALL AND IMPORT REQUIRED LIBRARIES
# ============================================================

# Install required third-party libraries
# Dependencies are installed through pyproject.toml (not at module import time).

# Import required libraries
import copy
import json
import math
import warnings
import glob
import logging
import os

from datetime import datetime, timedelta
from types import SimpleNamespace
from collections import defaultdict

import numpy as np
import pydicom
import SimpleITK as sitk

from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.datadict import tag_for_keyword, keyword_for_tag
from pydicom.tag import Tag



# Import-time diagnostic removed for library compatibility.
# Configure logger
logger = logging.getLogger("Dev_logger")

# ============================================================================
# SECTION 2 - User Contorolled Parameters
# ============================================================================
# Automatically detects single DICOM files and DICOM directories.
image_input = None  # Set a valid path only when executing this module directly.
parameter_estimation = "OFF"  # "ON": fill eligible missing DICOM values; "OFF": preserve original metadata.
output_path = None  # Optional output location for DICOM input. None value save output in the code directory
# Import-time diagnostic removed for library compatibility.


# ============================================================================
# Section 3 - CONFIGURED DEFAULT PARAMETERS
# ============================================================================
PatientWeight = "75.5"  # Patient weight (kg); patient-specific  # Changing is not Recommended
PatientAge = "045Y"  # Age (DICOM AS format); patient-specific
PatientSize = "1.75"  # Height (m); patient-specific  # Changing is not Recommended
PatientSex = "M"  # DICOM M/F/O; patient-specific  # Changing is not Recommended
RadionuclideTotalDose = "370000000"  # Administered dose; conventional PET Bq, Enhanced PET MBq; must be verified against source IOD  # Changing is not Recommended
RadionuclideHalfLife = "6586.2"  # Radionuclide half-life (s); must match actual tracer  # Changing is not Recommended
RadiopharmaceuticalStartTime = "093015"  # Administration start time (TM, HHMMSS); deprecated, prefer full DateTime  # Changing is not Recommended
RadiopharmaceuticalStartDateTime = "20261010093015"  # Administration start datetime (DT)  # Changing is not Recommended
RadionuclidePositronFraction = "0.967"  # Radionuclide positron fraction [0,1]
RadiopharmaceuticalSpecificActivity = "1000000000"  # Specific activity in Bq/micromole
RadiopharmaceuticalVolume = "5.0"  # Administered radiopharmaceutical volume (cm^3)
RadiopharmaceuticalStopTime = "093115"  # Administration stop time (TM); deprecated, prefer full DateTime
RadiopharmaceuticalStopDateTime = None  # Standard DICOM DT (0018,1079); use actual recorded stop datetime
ActualFrameDuration = "180000"  # Acquisition-frame duration (ms)  # Changing is not Recommended
FrameReferenceTime = "90000"  # Frame reference time (ms)  # Changing is not Recommended
AcquisitionDate = "20261010"  # Acquisition date (DA YYYYMMDD)  # Changing is not Recommended
AcquisitionTime = "100000"  # Acquisition time (TM HHMMSS)  # Changing is not Recommended
SeriesDate = "20261010"  # Series date (DA YYYYMMDD)  # Changing is not Recommended
SeriesTime = "095500"  # Series time (TM HHMMSS)  # Changing is not Recommended
FrameTime = "1000"  # Frame time (ms); conditional on acquisition type
TriggerTime = "0"  # Trigger time (ms); conditional on acquisition type
NominalInterval = "800"  # Nominal interval (ms); conditional on gated acquisition
IntervalsAcquired = "1000"  # Count of intervals acquired; conditional
IntervalsRejected = "10"  # Count of intervals rejected; conditional
NumberOfTimeSlices = "1"  # Number of time slices; conditional
NumberOfTimeSlots = "1"  # Number of time slots; conditional
Units = "BQML"  # PET pixel-value units (e.g. BQML, GML, CM2ML, CNTS)  # Changing is not Recommended
SUVType = "BW"  # SUV normalization type (e.g. BW)  # Changing is not Recommended
RescaleSlope = "1.0"  # Stored-pixel rescale slope (DS)  # Changing is not Recommended
RescaleIntercept = "0.0"  # Stored-pixel rescale intercept (DS)  # Changing is not Recommended
PhilipsSUVScaleFactor = None  # Philips vendor-private factor; not a public DICOM keyword  # Changing is not Recommended
PhilipsActivityConcentrationScaleFactor = None  # Philips vendor-private activity scaling; not a public DICOM keyword  # Changing is not Recommended
DecayCorrection = "START"  # Decay correction reference (START, ADMIN, NONE)  # Changing is not Recommended
DecayFactor = "1.2"  # Decay factor actually applied
DoseCalibrationFactor = "1.0"  # Dose calibration factor actually applied
SliceSensitivityFactor = "1.0"  # Slice sensitivity correction factor actually applied
ScatterFractionFactor = "0.25"  # Scatter fraction correction factor actually applied
DeadTimeFactor = "1.05"  # Detector dead-time factor actually applied
PrimaryPromptsCountsAccumulated = "10000000"  # Primary accumulated prompt counts
SecondaryCountsAccumulated = "1500000"  # Secondary accumulated counts
CorrectedImage = ["ATTN", "SCAT", "DECY", "RAN"]  # Applied correction codes (DICOM multi-valued CS)
RandomsCorrectionMethod = "DLYD"  # Randoms correction method (e.g. DLYD)
AttenuationCorrectionMethod = "CT-based attenuation correction"  # Attenuation correction method (descriptive text)
ScatterCorrectionMethod = "Model-based scatter correction"  # Scatter correction method (descriptive text)
ReconstructionMethod = "OSEM"  # Reconstruction method (descriptive text)
Rows = 256  # Image rows (US integer); must match NumPy array
Columns = 256  # Image columns (US integer); must match NumPy array
PixelSpacing = [4.0, 4.0]  # Row/column spacing [mm, mm]
SliceThickness = "3.0"  # Slice thickness (mm)
SpacingBetweenSlices = "3.0"  # Slice spacing (mm)
ImagePositionPatient = [0.0, 0.0, 0.0]  # Image position [x,y,z] mm; frame/slice specific
ImageOrientationPatient = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]  # Row/column direction cosines [6 values]
SeriesType = ["STATIC", "IMAGE"]  # PET series type (multi-valued CS)
CountsSource = "EMISSION"  # PET counts source (EMISSION or TRANSMISSION)
ReconstructionDiameter = "600.0"  # Reconstruction diameter (mm)
PETFrameCorrectionFactorsSequence = None  # Enhanced PET frame corrections (DICOM SQ); supply verified structure
DecayCorrectionDateTime = "20261010100000"  # Enhanced PET decay correction reference (DT)  # Changing is not Recommended
DecayCorrected = "YES"  # Enhanced PET decay correction flag YES/NO  # Changing is not Recommended
AttenuationCorrected = "YES"  # Enhanced PET attenuation correction flag YES/NO
ScatterCorrected = "YES"  # Enhanced PET scatter correction flag YES/NO
RandomsCorrected = "YES"  # Enhanced PET randoms correction flag YES/NO
erode_voxels = "0"  # Legacy application-only setting; NOT a DICOM attribute

_MANUAL_FIELDS = (
    "PatientWeight", "ActualFrameDuration", "FrameReferenceTime",
    "AcquisitionDate", "AcquisitionTime", "SeriesDate", "SeriesTime",
    "Units", "SUVType", "DecayCorrection", "RescaleSlope", "RescaleIntercept",
    "PatientAge",
    "PatientSize",
    "PatientSex",
    "DecayFactor",
    "DoseCalibrationFactor",
    "FrameTime",
    "TriggerTime",
    "NominalInterval",
    "IntervalsAcquired",
    "IntervalsRejected",
    "NumberOfTimeSlices",
    "NumberOfTimeSlots",
    "SliceSensitivityFactor",
    "ScatterFractionFactor",
    "DeadTimeFactor",
    "PrimaryPromptsCountsAccumulated",
    "SecondaryCountsAccumulated",
    "CorrectedImage",
    "RandomsCorrectionMethod",
    "AttenuationCorrectionMethod",
    "ScatterCorrectionMethod",
    "ReconstructionMethod",
    "SeriesType",
    "CountsSource",
    "ReconstructionDiameter",
    "PixelSpacing",
    "SliceThickness",
    "SpacingBetweenSlices",
    "ImagePositionPatient",
    "ImageOrientationPatient",
    "Rows",
    "Columns",
    "DecayCorrectionDateTime",
    "DecayCorrected",
    "AttenuationCorrected",
    "ScatterCorrected",
    "RandomsCorrected",
)
_MANUAL_RADIOPHARM_FIELDS = (
    "RadionuclideTotalDose", "RadionuclideHalfLife",
    "RadiopharmaceuticalStartTime", "RadiopharmaceuticalStartDateTime",
    "RadionuclidePositronFraction",
    "RadiopharmaceuticalSpecificActivity",
    "RadiopharmaceuticalStopTime", "RadiopharmaceuticalStopDateTime",
    "RadiopharmaceuticalVolume",
)
ENHANCED_PET_SOP_UID = "1.2.840.10008.5.1.4.1.1.130"

# Original scientific configuration defaults; applied only when parameter_estimation value is ON.
_DEFAULT_PARAMETERS = {name: copy.deepcopy(globals()[name]) for name in (
    'PatientWeight',
    'PatientAge',
    'PatientSize',
    'PatientSex',
    'RadionuclideTotalDose',
    'RadionuclideHalfLife',
    'RadiopharmaceuticalStartTime',
    'RadiopharmaceuticalStartDateTime',
    'RadionuclidePositronFraction',
    'RadiopharmaceuticalSpecificActivity',
    'RadiopharmaceuticalVolume',
    'RadiopharmaceuticalStopTime',
    'RadiopharmaceuticalStopDateTime',
    'ActualFrameDuration',
    'FrameReferenceTime',
    'AcquisitionDate',
    'AcquisitionTime',
    'SeriesDate',
    'SeriesTime',
    'FrameTime',
    'TriggerTime',
    'NominalInterval',
    'IntervalsAcquired',
    'IntervalsRejected',
    'NumberOfTimeSlices',
    'NumberOfTimeSlots',
    'Units',
    'SUVType',
    'RescaleSlope',
    'RescaleIntercept',
    'PhilipsSUVScaleFactor',
    'PhilipsActivityConcentrationScaleFactor',
    'DecayCorrection',
    'DecayFactor',
    'DoseCalibrationFactor',
    'SliceSensitivityFactor',
    'ScatterFractionFactor',
    'DeadTimeFactor',
    'PrimaryPromptsCountsAccumulated',
    'SecondaryCountsAccumulated',
    'CorrectedImage',
    'RandomsCorrectionMethod',
    'AttenuationCorrectionMethod',
    'ScatterCorrectionMethod',
    'ReconstructionMethod',
    'Rows',
    'Columns',
    'PixelSpacing',
    'SliceThickness',
    'SpacingBetweenSlices',
    'ImagePositionPatient',
    'ImageOrientationPatient',
    'SeriesType',
    'CountsSource',
    'ReconstructionDiameter',
    'PETFrameCorrectionFactorsSequence',
    'DecayCorrectionDateTime',
    'DecayCorrected',
    'AttenuationCorrected',
    'ScatterCorrected',
    'RandomsCorrected',
)}
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# Section 4 - INPUT DETECTION
# ============================================================================
def _radiuma_image_groups(directory):
    """Find independent PET inputs in a folder, preserving each DICOM series.

    DICOM slices sharing a SeriesInstanceUID represent one image volume.
    NIfTI, NRRD, and NHRD files are not accepted as PET image inputs.
    """
    from collections import defaultdict
    dicom_series = defaultdict(list)
    standalone = []
    for root, dirs, names in os.walk(directory):
        dirs.sort()
        for name in sorted(names):
            path = os.path.join(root, name)
            lower = name.lower()
            if lower.endswith((".nii", ".nii.gz", ".nrrd", ".nhdr")):
                continue
            # DICOM series can have either a conventional or no extension.
            try:
                header = pydicom.dcmread(path, stop_before_pixels=True)
            except Exception:
                continue
            if str(getattr(header, "Modality", "")) != "PT":
                continue
            series_uid = str(getattr(header, "SeriesInstanceUID", "") or "")
            # Without a UID, independent files must not be merged by guesswork.
            key = (series_uid, root) if series_uid else ("single", path)
            dicom_series[key].append(path)
    groups = [("file", path) for path in standalone]
    for key, filenames in sorted(dicom_series.items(), key=lambda item: str(item[0])):
        if len(filenames) > 1:
            # Match the file ordering used by SimpleITK when possible.
            reader = sitk.ImageSeriesReader()
            series_uid = key[0]
            if series_uid != "single":
                try:
                    gdcm_names = list(reader.GetGDCMSeriesFileNames(key[1], series_uid))
                    if gdcm_names and set(gdcm_names) == set(filenames):
                        filenames = gdcm_names
                except Exception:
                    pass
        groups.append(("dicom", filenames))
    return groups
def _suv_detect_source_groups(image_input):
    """Return independently processable PET DICOM groups under a folder."""
    if os.path.isfile(image_input):
        return [[os.fspath(image_input)]]
    if not os.path.isdir(image_input):
        raise FileNotFoundError(image_input)
    groups = _radiuma_image_groups(os.fspath(image_input))
    if not groups:
        raise ValueError('No PET DICOM series found')
    if any(kind != 'dicom' for kind, _ in groups):
        raise ValueError('Unsupported PET file group')
    return [list(paths) for _, paths in groups]
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# Section 5 - PET DATA EXTRACTION
# ============================================================================
def _pet_header_dictionary(ds):
    """
    DICOM keywords required by radiuma-suv, recursively including sequences.
    """
    result = {}
    for elem in ds:
        key = elem.keyword
        if elem.VR == "SQ" and key:
            result[key] = [_pet_header_dictionary(item) for item in elem.value]
        elif key and elem.VR != "SQ":
            result[key] = _python_metadata(elem.value)
        elif elem.tag in (Tag(0x7053, 0x1000), Tag(0x7053, 0x1009)):
            result[f"{int(elem.tag):08X}"] = _python_metadata(elem.value)
    return result
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# Section 6 - METADATA VALIDATION AND PARAMETER ESTIMATION
# ============================================================================
'''
This section validates PET DICOM metadata and optionally fills missing fields using predefined configuration values
when parameter_estimation="ON". It checks numerical values, dates, times, image geometry, and radiopharmaceutical information, 
including manufacturer-specific Philips scaling factors. Existing metadata values are preserved, and changes are applied to a 
copy of the original DICOM dataset. When parameter_estimation="OFF", the dataset is returned unchanged.
'''
def safe_float(x, default=None):
    try:
        # If the value is missing or empty, return the default
        if x is None or x == "":
            return default
        # Try converting to float
        return float(x)
    except Exception:
        # If conversion fails, return the default
        return default
def _dicom_value_missing(value):
    return value is None or (isinstance(value, str) and not value.strip()) or (hasattr(value, "__len__") and not isinstance(value, (str, bytes)) and len(value) == 0)
def _dicom_find_value(ds, keyword):
    """Locate a field in the original header, including nested DICOM sequences."""
    from pydicom.dataset import Dataset
    if not isinstance(ds, Dataset):
        return None
    try:
        element = ds.data_element(keyword)
    except KeyError:
        element = None
    if element is not None and not _dicom_value_missing(element.value):
        return element.value
    for elem in ds:
        if elem.VR == "SQ":
            for item in elem.value:
                found = _dicom_find_value(item, keyword)
                if not _dicom_value_missing(found):
                    return found
    return None
def apply_configured_pet_metadata(ds, estimation, settings):
    """In ON mode, fill missing DICOM metadata from configured values."""
    def invalid_metadata(keyword, value):
        if _dicom_value_missing(value):
            return True
        positive = {
            "PatientWeight", "PatientSize", "RadionuclideTotalDose",
            "RadionuclideHalfLife", "ActualFrameDuration", "FrameTime",
            "NominalInterval", "ReconstructionDiameter", "SliceThickness",
            "SpacingBetweenSlices", "Rows", "Columns", "NumberOfTimeSlices",
            "NumberOfTimeSlots", "RadiopharmaceuticalSpecificActivity",
            "RadiopharmaceuticalVolume", "DoseCalibrationFactor", "DecayFactor",
            "SliceSensitivityFactor", "DeadTimeFactor",
        }
        nonnegative = {
            "FrameReferenceTime", "TriggerTime", "IntervalsAcquired",
            "IntervalsRejected", "PrimaryPromptsCountsAccumulated",
            "SecondaryCountsAccumulated", "ScatterFractionFactor",
        }
        if keyword in positive or keyword in nonnegative or keyword in {
            "RescaleSlope", "RescaleIntercept", "RadionuclidePositronFraction"
        }:
            number = safe_float(value)
            if number is None or not math.isfinite(number):
                return True
            if keyword in positive:
                return number <= 0
            if keyword in nonnegative:
                return number < 0
            if keyword == "RescaleSlope":
                return number == 0
            if keyword == "RadionuclidePositronFraction":
                return not (0 < number <= 1)
            return False
        permitted = {
            "Units": {"BQML", "GML", "CM2ML", "CNTS"},
            "SUVType": {"BW", "LBM", "LBMJAMES128", "LBMJANMA", "IBW", "BSA"},
            "DecayCorrection": {"START", "ADMIN", "NONE"},
            "PatientSex": {"M", "F", "O"},
            "DecayCorrected": {"YES", "NO"},
            "AttenuationCorrected": {"YES", "NO"},
            "ScatterCorrected": {"YES", "NO"},
            "RandomsCorrected": {"YES", "NO"},
            "CountsSource": {"EMISSION", "TRANSMISSION"},
        }
        if keyword in permitted:
            return str(value).strip().upper() not in permitted[keyword]
        if keyword == "PatientAge":
            import re
            return re.fullmatch(r"\d{3}[DWMY]", str(value).strip()) is None
        if keyword in ("AcquisitionDate", "SeriesDate"):
            return parse_dicom_da(value) is None
        if keyword in ("AcquisitionTime", "SeriesTime",
                       "RadiopharmaceuticalStartTime", "RadiopharmaceuticalStopTime"):
            parsed = parse_dicom_tm(value)
            return (parsed is None or not (0 <= parsed[0] < 24 and
                    0 <= parsed[1] < 60 and 0 <= parsed[2] < 60))
        if keyword in ("RadiopharmaceuticalStartDateTime", "RadiopharmaceuticalStopDateTime", "DecayCorrectionDateTime"):
            return parse_dicom_dt(value) is None
        if keyword in ("PixelSpacing", "ImagePositionPatient", "ImageOrientationPatient"):
            try:
                numbers = [float(x) for x in (value.split("\\") if isinstance(value, str) else value)]
                required = {"PixelSpacing": 2, "ImagePositionPatient": 3,
                            "ImageOrientationPatient": 6}[keyword]
                return (len(numbers) != required or
                        any(not math.isfinite(x) for x in numbers) or
                        (keyword == "PixelSpacing" and any(x <= 0 for x in numbers)))
            except (TypeError, ValueError):
                return True
        if keyword in ("CorrectedImage", "SeriesType"):
            items = value.split("\\") if isinstance(value, str) else value
            return not items or any(not str(x).strip() for x in items)
        return False
    if estimation not in ("ON", "OFF"):
        raise ValueError('parameter_estimation must be "ON" or "OFF"')
    if estimation == "OFF":
        return ds
    # Original headers are never modified on disk.
    ds = copy.deepcopy(ds)
    for keyword in _MANUAL_FIELDS:
        value = settings[keyword]
        if not invalid_metadata(keyword, value) and _dicom_value_missing(_dicom_find_value(ds, keyword)):
            if keyword in ("PixelSpacing", "ImagePositionPatient", "ImageOrientationPatient"):
                value = [float(part.strip()) for part in (value.split("\\") if isinstance(value, str) else value)]
            elif keyword in ("CorrectedImage", "SeriesType"):
                value = [part.strip() for part in (value.split("\\") if isinstance(value, str) else value)]
            elif keyword in ("Rows", "Columns", "IntervalsAcquired", "IntervalsRejected",
                             "NumberOfTimeSlices", "NumberOfTimeSlots"):
                value = int(value)
            setattr(ds, keyword, value)
    supplied_radio = {keyword: settings[keyword] for keyword in _MANUAL_RADIOPHARM_FIELDS
                      if not invalid_metadata(keyword, settings[keyword]) and
                      _dicom_value_missing(_dicom_find_value(ds, keyword))}
    if supplied_radio:
        from pydicom.dataset import Dataset
        from pydicom.sequence import Sequence
        seq = getattr(ds, "RadiopharmaceuticalInformationSequence", None)
        if not seq:
            ds.RadiopharmaceuticalInformationSequence = Sequence([Dataset()])
        for keyword, value in supplied_radio.items():
            setattr(ds.RadiopharmaceuticalInformationSequence[0], keyword, value)
    # Supply valid configured Philips private scale factors only where absent.
    if "PHILIPS" in str(getattr(ds, "Manufacturer", "") or "").upper() and str(getattr(ds, "Units", "") or "").upper() == "CNTS":
        for tag, configured in (((0x7053, 0x1000), settings["PhilipsSUVScaleFactor"]),
                                ((0x7053, 0x1009), settings["PhilipsActivityConcentrationScaleFactor"])):
            original = ds.get(tag)
            original_number = safe_float(original.value if original is not None else None)
            configured_number = safe_float(configured)
            if (original_number is None
                    and configured_number is not None and math.isfinite(configured_number) and configured_number > 0):
                ds.add_new(tag, "DS", str(configured_number))
    if (settings["PETFrameCorrectionFactorsSequence"] is not None and
            _dicom_value_missing(_dicom_find_value(ds, "PETFrameCorrectionFactorsSequence"))):
        from pydicom.sequence import Sequence
        if not isinstance(settings["PETFrameCorrectionFactorsSequence"], Sequence):
            raise TypeError("PETFrameCorrectionFactorsSequence must be a pydicom Sequence")
        ds.PETFrameCorrectionFactorsSequence = copy.deepcopy(settings["PETFrameCorrectionFactorsSequence"])
    return ds
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SECTION 7 - DICOM DATE AND TIME PARSING
# ============================================================================
'''
This section parses and converts standard and private DICOM date, time, and datetime fields into Python-compatible formats, 
supporting fractional seconds and multiple datetime formats. It provides standardized temporal information for PET acquisition 
timing and radiopharmaceutical decay calculations.
'''
def parse_dicom_tm(tm_str):
    # Return None if empty
    if not tm_str:
        return None
    # Convert to string and remove surrounding spaces
    s = str(tm_str).strip()
    # Default fractional seconds = 0
    frac = 0
    # If there is a decimal point, extract fractional seconds
    if "." in s:
        _, frac_str = s.split(".", 1)
        # Pad or trim to 6 digits for microseconds
        frac = int((frac_str + "000000")[:6])
    # Take the part before decimal and pad to HHMMSS
    base = s.split(".")[0].ljust(6, "0")
    try:
        # Parse hours
        h = int(base[0:2])
        # Parse minutes
        m = int(base[2:4])
        # Parse seconds
        sec = int(base[4:6])
        # Return time components plus microseconds
        return h, m, sec, frac
    except Exception:
        # Return None if parsing fails
        return None
def parse_dicom_da(da_str):
    # Return None if empty
    if not da_str:
        return None
    try:
        # Parse date in YYYYMMDD format
        return datetime.strptime(str(da_str), "%Y%m%d")
    except Exception:
        # Return None if invalid
        return None
def combine_da_tm(date_str, time_str):
    # Parse the date
    d = parse_dicom_da(date_str)
    # Parse the time
    t = parse_dicom_tm(time_str)
    # If either failed, return None
    if d is None or t is None:
        return None
    # Replace time part of the date object and return full datetime
    return d.replace(hour=t[0], minute=t[1], second=t[2], microsecond=t[3])
def parse_dicom_dt(dt_str):
    # Return None if empty
    if not dt_str:
        return None
    # Convert to clean string
    s = str(dt_str).strip()
    # Remove timezone offset if present, because strptime below does not use it
    for sign in ["+", "-"]:
        if len(s) > 14 and sign in s[14:]:
            s = s.split(sign, 1)[0]
            break
    # Try several possible DICOM datetime formats
    for fmt in ["%Y%m%d%H%M%S.%f", "%Y%m%d%H%M%S", "%Y%m%d%H%M", "%Y%m%d"]:
        try:
            return datetime.strptime(s, fmt)
        except Exception:
            pass
    # If none worked, return None
    return None
def parse_private_dt(ds, tag):
    elem = ds.get(tag)
    if elem is None:
        return None
    value = getattr(elem, "value", None)
    if isinstance(value, bytes):
        value = value.decode("ascii", errors="ignore").strip("\x00 ")
    return parse_dicom_dt(value)
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SECTION 9 - PATIENT AND RADIOPHARMACEUTICAL DATA
# ============================================================================
'''
This section extracts and validates essential patient and radiopharmaceutical information from PET DICOM metadata, including patient weight,
height, sex, injected dose, and radionuclide half-life. It also determines radiopharmaceutical administration and scan reference times, 
accounting for manufacturer-specific metadata, timing inconsistencies, and midnight transitions to support accurate SUV calculation and 
radioactive decay correction.
'''
def get_radiopharm_item(ds):
    # Get RadiopharmaceuticalInformationSequence
    seq = getattr(ds, "RadiopharmaceuticalInformationSequence", None)
    # Fail if missing
    if not seq or len(seq) == 0:
        raise ValueError("Missing RadiopharmaceuticalInformationSequence")
    # Use the first item
    return seq[0]
def get_weight_kg(ds):
    # Read patient weight
    w = safe_float(getattr(ds, "PatientWeight", None))
    # Fail if missing or invalid
    if w is None or w <= 0:
        raise ValueError("Missing or zero PatientWeight")
    # Convert grams to kg if it looks too large
    if w >= 1000:
        w /= 1000.0
    # Return kg
    return w
def get_weight_g(ds):
    # Convert kg to grams
    return get_weight_kg(ds) * 1000.0
def get_height_m(ds):
    # Read patient height from PatientSize
    h = safe_float(getattr(ds, "PatientSize", None))
    # Fail if invalid
    if h is None or h <= 0:
        raise ValueError("Missing or invalid PatientSize")
    return h
def get_height_cm(ds):
    # Convert meters to cm
    return get_height_m(ds) * 100.0
def get_patient_sex(ds):
    # Reject absent or invalid values rather than disguising them as "O".
    sex = str(getattr(ds, "PatientSex", "") or "").strip().upper()
    if sex not in ("M", "F", "O"):
        raise ValueError("Missing or invalid PatientSex")
    return sex
def get_total_dose_bq(ds):
    # Get radiopharmaceutical metadata block
    rp = get_radiopharm_item(ds)
    # Read total dose
    dose = safe_float(getattr(rp, "RadionuclideTotalDose", None))
    # Fail if missing or invalid
    if dose is None or dose <= 0:
        raise ValueError("Missing or invalid RadionuclideTotalDose")
    # If it looks too small, assume MBq and convert to Bq
    if dose < 1e4:
        dose *= 1e6
    return dose
def get_half_life_s(ds):
    # Get radiopharmaceutical metadata block
    rp = get_radiopharm_item(ds)
    # Read half-life
    hl = safe_float(getattr(rp, "RadionuclideHalfLife", None))
    # Fail if missing or invalid
    if hl is None or hl <= 0:
        raise ValueError("Missing or invalid RadionuclideHalfLife")
    return hl
def get_radiopharm_start_datetime(ds, reference_dt=None):
    rp = get_radiopharm_item(ds)
    half_life_s = get_half_life_s(ds)
    full_dt = parse_dicom_dt(getattr(rp, "RadiopharmaceuticalStartDateTime", None))
    start_time = getattr(rp, "RadiopharmaceuticalStartTime", None)
    if reference_dt is None:
        reference_dt = parse_dicom_dt(getattr(ds, "AcquisitionDateTime", None))
        if reference_dt is None:
            reference_dt = combine_da_tm(getattr(ds, "AcquisitionDate", None),
                                          getattr(ds, "AcquisitionTime", None))
    if reference_dt is None:
        raise ValueError("Missing decay-correction reference datetime")
    if full_dt is not None:
        delta_s = (reference_dt - full_dt).total_seconds()
        if -3600.0 <= delta_s < 2.0 * half_life_s:
            return full_dt
        if half_life_s >= 41400.0:
            raise ValueError("Unreliable administration date for long-lived radionuclide")
        warnings.warn("Administration date is inconsistent; using reference date and administration time")
        start_time = full_dt.strftime("%H%M%S.%f")
    else:
        if not start_time:
            raise ValueError("Missing radiopharmaceutical start date/time")
        if half_life_s >= 41400.0:
            raise ValueError("Full administration datetime required for long-lived radionuclide")
        warnings.warn("Only administration time available; reference date is assumed")
    t = parse_dicom_tm(start_time)
    if t is None:
        raise ValueError("Invalid radiopharmaceutical administration time")
    start_dt = reference_dt.replace(hour=t[0], minute=t[1], second=t[2], microsecond=t[3])
    if (reference_dt - start_dt).total_seconds() < -3600.0:
        start_dt -= timedelta(days=1)
        warnings.warn("Administration occurred before midnight; using previous day")
    return start_dt
def get_reference_time_start(ds):
    manufacturer = str(getattr(ds, "Manufacturer", "") or "").upper()
    is_ge = "GE" in manufacturer or "GEMS" in manufacturer
    is_siemens = "SIEMENS" in manufacturer
    # The private scan-start tags are valid for their own manufacturers only.
    if is_ge:
        ge_dt = parse_private_dt(ds, (0x0009, 0x100D))
        if ge_dt is not None:
            return ge_dt
    if is_siemens:
        siemens_dt = parse_private_dt(ds, (0x0071, 0x1022))
        if siemens_dt is not None:
            return siemens_dt
    acq_dt = parse_dicom_dt(getattr(ds, "AcquisitionDateTime", None))
    if acq_dt is None:
        acq_dt = combine_da_tm(getattr(ds, "AcquisitionDate", None),
                                getattr(ds, "AcquisitionTime", None))
    if acq_dt is None:
        raise ValueError("Missing acquisition datetime")
    series_dt = combine_da_tm(getattr(ds, "SeriesDate", None),
                               getattr(ds, "SeriesTime", None))
    if series_dt is not None and series_dt == acq_dt:
        return acq_dt
    frame_ref_ms = safe_float(getattr(ds, "FrameReferenceTime", None))
    if frame_ref_ms is None or frame_ref_ms < 0:
        raise ValueError("Missing or invalid FrameReferenceTime")
    frame_ref_s = frame_ref_ms / 1000.0
    # GE's acquisition timestamp needs no average-count-rate offset.
    if is_ge:
        return acq_dt - timedelta(seconds=frame_ref_s)
    frame_duration_ms = safe_float(getattr(ds, "ActualFrameDuration", None))
    if frame_duration_ms is None or frame_duration_ms <= 0:
        raise ValueError("Missing or invalid ActualFrameDuration")
    half_life_s = get_half_life_s(ds)
    tave_s = average_count_rate_time_s(half_life_s, frame_duration_ms / 1000.0)
    return acq_dt + timedelta(seconds=tave_s - frame_ref_s)
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SUBSECTION 10 - SPATIAL ORDERING AND SERIES VALIDATION
# ============================================================================
'''
This subsection sorts PET DICOM slices according to their spatial position and orientation to ensure correct anatomical 
ordering. It also verifies consistency of slice orientation and essential metadata, including image units, decay correction, 
patient weight, height, and sex, to prevent processing of inconsistent PET series.
'''
def sort_pet_slices(datasets):
    def key(ds):
        orientation = getattr(ds, "ImageOrientationPatient", None)
        position = getattr(ds, "ImagePositionPatient", None)
        if orientation is None or len(orientation) != 6 or position is None or len(position) != 3:
            raise ValueError("Missing DICOM geometry; PET slice order cannot be verified")
        row_dir = np.asarray(orientation[:3], dtype=float)
        col_dir = np.asarray(orientation[3:], dtype=float)
        normal = np.cross(row_dir, col_dir)
        return (float(np.dot(np.asarray(position, dtype=float), normal)),
                safe_float(getattr(ds, "InstanceNumber", None), 0))
    if datasets:
        first_orientation = np.asarray(datasets[0].ImageOrientationPatient, dtype=float)
        if any(not np.allclose(np.asarray(ds.ImageOrientationPatient, dtype=float),
                                   first_orientation, atol=1e-5)
               for ds in datasets):
            raise ValueError("Inconsistent PET slice orientation")
    return sorted(datasets, key=key)
def assert_series_consistency(datasets):
    # Attributes that should not vary slice to slice
    attrs_to_check = [
        "Units",
        "DecayCorrection",
        "PatientWeight",
        "PatientSize",
        "PatientSex",
    ]
    # Check each attribute
    for attr in attrs_to_check:
        # Collect values as strings
        vals = [str(getattr(ds, attr, None)) for ds in datasets]
        # If more than one unique value exists, fail
        if len(set(vals)) > 1:
            raise ValueError(f"Inconsistent series attribute: {attr} -> {set(vals)}")
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SUBSECTION 11 - DICOM METADATA DICTIONARY CONVERSION
# ============================================================================
'''
This subsection converts DICOM metadata between Python dictionaries and pydicom datasets, supporting 
standard and private tags, nested sequences, and data validation. It also preserves original metadata 
while inserting only missing fields, ensuring existing information remains unchanged.
'''
def _key_tag(key):
    if isinstance(key, (tuple, list)) and len(key) == 2:
        return Tag(int(key[0]), int(key[1]))
    if isinstance(key, int):
        return Tag(key)
    if not isinstance(key, str):
        raise TypeError("DICOM dictionary keys must be keywords or tag identifiers")
    if tag_for_keyword(key) is not None:
        return tag_for_keyword(key)
    compact = key.replace("(", "").replace(")", "").replace(",", "").replace(" ", "").replace("0x", "")
    if len(compact) == 8 and all(c in "0123456789abcdefABCDEF" for c in compact):
        return Tag(int(compact, 16))
    raise ValueError(f"Unknown DICOM keyword or tag: {key!r}")
def _dicom_dataset(metadata):
    if not isinstance(metadata, dict):
        raise TypeError("DICOM metadata must be a Python dictionary")
    ds = Dataset()
    for key, value in metadata.items():
        tag = _key_tag(key)
        if isinstance(value, dict):
            raise TypeError(f"Metadata {key!r} must not be a nested dictionary outside a sequence")
        if isinstance(value, (list, tuple)) and any(isinstance(x, dict) for x in value):
            if not all(isinstance(x, dict) for x in value):
                raise TypeError(f"Invalid DICOM sequence for {key!r}")
            seq = Sequence([_dicom_dataset(x) for x in value])
            if tag_for_keyword(str(key)) == tag:
                setattr(ds, str(key), seq)
            else:
                ds.add_new(tag, "SQ", seq)
        elif tag_for_keyword(str(key)) == tag:
            setattr(ds, str(key), copy.deepcopy(value))
        else:
            if isinstance(value, pydicom.dataelem.DataElement):
                ds.add(copy.deepcopy(value))
            else:
                # Philips CNTS private scaling values are DS.
                if tag in (Tag(0x7053, 0x1000), Tag(0x7053, 0x1009)):
                    ds.add_new(tag, "DS", str(value))
                else:
                    raise ValueError(f"Unknown value representation for private tag {key!r}")
    return ds
def _python_metadata(value):
    if isinstance(value, Dataset):
        result = {}
        for item in value:
            key = item.keyword if item.keyword else f"{int(item.tag):08X}"
            result[key] = _python_metadata(item.value)
        return result
    if isinstance(value, (Sequence, list, tuple, pydicom.multival.MultiValue)):
        return [_python_metadata(x) for x in value]
    if isinstance(value, pydicom.dataelem.DataElement):
        return _python_metadata(value.value)
    if isinstance(value, (bytes, str, int, float, bool)) or value is None:
        return copy.deepcopy(value)
    return copy.deepcopy(value)
def _metadata_with_insertions(original, ds):
    """Copy source metadata and insert only newly supplied DICOM fields."""
    def insert_missing(destination, source):
        for key, value in source.items():
            if key not in destination:
                destination[key] = copy.deepcopy(value)
            elif isinstance(destination[key], list) and isinstance(value, list):
                if (len(destination[key]) == len(value)
                        and all(isinstance(x, dict) for x in destination[key])
                        and all(isinstance(x, dict) for x in value)):
                    for dst, src in zip(destination[key], value):
                        insert_missing(dst, src)
    result = copy.deepcopy(original)
    insert_missing(result, _python_metadata(ds))
    return result
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SECTION 12 - CONVENTIONAL PET METADATA PREPARATION
# ============================================================================
'''
This section prepares and validates slice-level metadata for conventional PET images, ensuring consistency
between the NumPy array dimensions and DICOM headers. It applies optional metadata completion, verifies series 
consistency, extracts rescale slopes and intercepts for SUV conversion, and preserves the original metadata
while incorporating missing fields when parameter estimation is enabled.
'''
def _prepare_conventional(array, input_metadata, estimation, settings):
    if "slices" in input_metadata:
        slice_dicts = input_metadata["slices"]
        if not isinstance(slice_dicts, list) or len(slice_dicts) != array.shape[0]:
            raise ValueError("metadata['slices'] must be a list matching PET slice count")
        datasets = [_dicom_dataset(d) for d in slice_dicts]
    else:
        if array.shape[0] != 1:
            raise ValueError("Conventional 3D PET requires metadata['slices'] with one header per slice")
        slice_dicts = [input_metadata]
        datasets = [_dicom_dataset(input_metadata)]
    applied = [apply_configured_pet_metadata(ds, estimation, settings) for ds in datasets]
    for index, ds in enumerate(applied):
        if getattr(ds, "Rows", None) is not None and int(ds.Rows) != array.shape[1]:
            raise ValueError(f"Slice {index}: metadata Rows disagrees with input array")
        if getattr(ds, "Columns", None) is not None and int(ds.Columns) != array.shape[2]:
            raise ValueError(f"Slice {index}: metadata Columns disagrees with input array")
    assert_series_consistency(applied)
    # The input NumPy array must already be in the same order as metadata['slices'].
    slopes = []
    intercepts = []
    for ds in applied:
        slope = safe_float(getattr(ds, "RescaleSlope", None))
        intercept = safe_float(getattr(ds, "RescaleIntercept", None))
        if slope is None or intercept is None:
            raise ValueError("Missing RescaleSlope or RescaleIntercept")
        slopes.append(slope)
        intercepts.append(intercept)
    result_metadata = copy.deepcopy(input_metadata)
    if estimation == "ON":
        if "slices" in result_metadata:
            result_metadata["slices"] = [
                _metadata_with_insertions(original, ds)
                for original, ds in zip(slice_dicts, applied)
            ]
        else:
            result_metadata = _metadata_with_insertions(input_metadata, applied[0])
    return applied, np.asarray(slopes, dtype=np.float64), np.asarray(intercepts, dtype=np.float64), result_metadata
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# section 13 - SUVbw CALCULATION ENGINE
# ============================================================================

# ============================================================================
# SUBSECTION 13-1: FRAME-TIMING CALCULATION
# ============================================================================
'''
This subsection calculates the effective average count-rate time within a PET acquisition frame using the 
radionuclide half-life and frame duration, accounting for radioactive decay to support accurate PET scan 
timing and decay correction.
''' 
def average_count_rate_time_s(half_life_s, frame_duration_s):
    # Decay constant
    lam = math.log(2.0) / half_life_s
    # Return average count rate time
    return (1.0 / lam) * math.log((lam * frame_duration_s) / (1.0 - math.exp(-lam * frame_duration_s)))
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SUBSECTION 13-2: BODY COMPOSITION NORMALIZATION
# ============================================================================
'''
This subsection calculates body composition parameters, including lean body mass (LBM) using the James128,
Morgan, and Janma formulas, ideal body weight (IBW), and body surface area (BSA). These calculations use patient weight,
height, and sex to provide normalization factors for PET standardized uptake value (SUV) calculations.

'''
def lbm_james128_kg(weight_kg, height_cm, sex):
    male = 1.10 * weight_kg - 128.0 * (weight_kg / height_cm) ** 2
    female = 1.07 * weight_kg - 148.0 * (weight_kg / height_cm) ** 2
    if sex == "M":
        return male
    if sex == "F":
        return female
    return 0.5 * (male + female)
def lbm_morgan_kg(weight_kg, height_cm, sex):
    male = 1.10 * weight_kg - 120.0 * (weight_kg / height_cm) ** 2
    female = 1.07 * weight_kg - 148.0 * (weight_kg / height_cm) ** 2
    if sex == "M":
        return male
    if sex == "F":
        return female
    return 0.5 * (male + female)
def lbm_janma_kg(weight_kg, height_cm, sex):
    bmi = weight_kg / ((height_cm / 100.0) ** 2)
    male = 9270.0 * weight_kg / (6680.0 + 216.0 * bmi)
    female = 9270.0 * weight_kg / (8780.0 + 244.0 * bmi)
    if sex == "M":
        return male
    if sex == "F":
        return female
    return 0.5 * (male + female)
def ibw_kg(height_cm, sex):
    male = 48.0 + 1.06 * (height_cm - 152.0)
    female = 45.5 + 0.91 * (height_cm - 152.0)
    if sex == "M":
        return male
    if sex == "F":
        return female
    return 0.5 * (male + female)
def bsa_cm2(weight_kg, height_cm):
    bsa_m2 = 0.007184 * (height_cm ** 0.725) * (weight_kg ** 0.425)
    return bsa_m2 * 10000.0
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SUBSECTION 13-3: CONVENTIONAL PET INTENSITY RESCALING
# ============================================================================
'''
This function applies slice-specific DICOM rescale slopes and intercepts to a PET image 
volume, converting stored pixel values into rescaled values using a linear transformation. 
It returns the rescaled volume as a 64-bit floating-point NumPy array for subsequent PET 
processing and SUV calculations.
'''
def apply_rescale(volume, slopes, intercepts):
    # Create empty output array
    out = np.empty_like(volume, dtype=np.float64)
    # Process each slice
    for z in range(volume.shape[0]):
        # Apply slope and intercept
        out[z] = volume[z] * slopes[z] + intercepts[z]
    # Return rescaled values
    return out
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SUBSECTION 13-4: CONVENTIONAL PET SUVbw UNIT CONVERSION
# ============================================================================
'''
This subsection converts conventional PET image values into body-weight-normalized standardized 
uptake values (SUVbw) based on DICOM units (BQML, GML, CM2ML, and CNTS). It applies rescale factors, 
radioactive decay correction, patient-specific normalization, and Philips-specific scaling factors while 
validating metadata and rejecting unsupported or invalid conversion conditions.
'''
def convert_bqml_to_suvbw(datasets, values):
    out = np.empty_like(values, dtype=np.float64)
    # Use each frame's dose, half-life, and reference/administration times.
    for z, ds in enumerate(datasets):
        weight_g = get_weight_g(ds)
        dose_bq_admin = get_total_dose_bq(ds)
        half_life_s = get_half_life_s(ds)
        lam = math.log(2.0) / half_life_s
        decay_correction = str(getattr(ds, "DecayCorrection", "") or "").upper().strip()
        if decay_correction == "ADMIN":
            dose_at_ref = dose_bq_admin
        elif decay_correction == "START":
            tref = get_reference_time_start(ds)
            tadm = get_radiopharm_start_datetime(ds, reference_dt=tref)
            dose_at_ref = dose_bq_admin * math.exp(-lam * (tref - tadm).total_seconds())
        elif decay_correction == "NONE":
            acq_dt = parse_dicom_dt(getattr(ds, "AcquisitionDateTime", None))
            if acq_dt is None:
                acq_dt = combine_da_tm(getattr(ds, "AcquisitionDate", None),
                                       getattr(ds, "AcquisitionTime", None))
            if acq_dt is None:
                raise ValueError("Missing acquisition time for DecayCorrection=NONE")
            frame_duration_ms = safe_float(getattr(ds, "ActualFrameDuration", None))
            if frame_duration_ms is None or frame_duration_ms <= 0:
                raise ValueError("Missing or invalid ActualFrameDuration for DecayCorrection=NONE")
            tave_s = average_count_rate_time_s(half_life_s, frame_duration_ms / 1000.0)
            tref = acq_dt + timedelta(seconds=tave_s)
            tadm = get_radiopharm_start_datetime(ds, reference_dt=tref)
            dose_at_ref = dose_bq_admin * math.exp(-lam * (tref - tadm).total_seconds())
        else:
            raise NotImplementedError(f"Unsupported DecayCorrection: {decay_correction}")
        if dose_at_ref <= 0 or not math.isfinite(dose_at_ref):
            raise ValueError("Invalid decay-corrected radionuclide dose")
        out[z] = values[z] * (weight_g / dose_at_ref)
    return out
def convert_gml_to_suvbw(ds0, values):
    suv_type = str(getattr(ds0, "SUVType", "") or "").upper().strip()
    if suv_type == "":
        raise ValueError("Missing SUVType in DICOM header")
    if suv_type == "BW":
        return values
    weight_kg = get_weight_kg(ds0)
    height_cm = get_height_cm(ds0)
    sex = get_patient_sex(ds0)
    if suv_type == "LBM":
        factor = lbm_morgan_kg(weight_kg, height_cm, sex)
    elif suv_type == "LBMJAMES128":
        factor = lbm_james128_kg(weight_kg, height_cm, sex)
    elif suv_type == "LBMJANMA":
        factor = lbm_janma_kg(weight_kg, height_cm, sex)
    elif suv_type == "IBW":
        factor = ibw_kg(height_cm, sex)
    else:
        raise NotImplementedError(f"Unsupported SUVType in GML: {suv_type}")
    if factor <= 0:
        raise ValueError(f"Computed non-positive normalization factor for SUVType={suv_type}")
    return values * (weight_kg / factor)
def convert_cm2ml_to_suvbw(ds0, values):
    suv_type = str(getattr(ds0, "SUVType", "") or "").upper().strip()
    # Allow missing SUVType, but if present it should be BSA
    if suv_type not in ["", "BSA"]:
        raise ValueError(f"Unexpected SUVType for CM2ML: {suv_type}")
    weight_kg = get_weight_kg(ds0)
    weight_g = weight_kg * 1000.0
    height_cm = get_height_cm(ds0)
    bsa = bsa_cm2(weight_kg, height_cm)
    if bsa <= 0:
        raise ValueError("Computed non-positive BSA")
    return values * (weight_g / bsa)
def convert_cnts_to_suvbw(datasets, values):
    out = np.empty_like(values, dtype=np.float64)
    for z, ds in enumerate(datasets):
        if "PHILIPS" not in str(getattr(ds, "Manufacturer", "") or "").upper():
            raise ValueError("CNTS scaling requires Philips manufacturer")
        suv_sf_tag = ds.get((0x7053, 0x1000))
        ac_sf_tag = ds.get((0x7053, 0x1009))
        suv_sf = safe_float(suv_sf_tag.value if suv_sf_tag is not None else None)
        ac_sf = safe_float(ac_sf_tag.value if ac_sf_tag is not None else None)
        # Prefer activity concentration when available (IBSI v3).
        if ac_sf is not None and ac_sf != 0:
            out[z] = convert_bqml_to_suvbw([ds], values[z:z+1] * ac_sf)[0]
        elif suv_sf is not None and suv_sf != 0:
            suv_type = str(getattr(ds, "SUVType", "") or "").upper().strip()
            if suv_type not in ("", "BW"):
                raise ValueError("Philips direct SUV factor does not specify SUVbw")
            out[z] = values[z] * suv_sf
        else:
            raise ValueError("CNTS without valid Philips scale factors")
    return out
def convert_to_suvbw(datasets, volume, slopes, intercepts):
    ds0 = datasets[0]
    if np.any(np.abs(intercepts) > 1e-12):
        raise ValueError("RescaleIntercept is not zero. Excluding dataset.")
    values = apply_rescale(volume, slopes, intercepts)
    units = str(getattr(ds0, "Units", "") or "").upper().strip()
    if units == "BQML":
        return convert_bqml_to_suvbw(datasets, values)
    if units == "GML":
        return convert_gml_to_suvbw(ds0, values)
    if units == "CM2ML":
        return convert_cm2ml_to_suvbw(ds0, values)
    if units == "CNTS":
        return convert_cnts_to_suvbw(datasets, values)
    raise NotImplementedError(f"Unsupported Units: {units}")
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SUBSECTION 13-5: ENHANCED PET FRAME METADATA AND CALIBRATION
# ============================================================================
def enhanced_group_attribute(dataset, per_frame, shared, attribute):
  
    '''
    This subsection extracts frame-specific and shared metadata from Enhanced PET DICOM images, identifies measurement 
    units, and applies real-world value mappings using linear rescale factors or lookup tables (LUTs). It prioritizes 
    appropriate SUV mappings, supports fallback pixel transformations, and validates calibration parameters to prepare 
    accurate quantitative PET values for subsequent SUVbw conversion.
    
    '''
    for group in (per_frame, shared):
        if group is None:
            continue
        direct = getattr(group, attribute, None)
        if direct is not None and direct != "":
            return direct
        for elem in group:
            if getattr(elem, "VR", None) != "SQ":
                continue
            for item in elem.value:
                found = getattr(item, attribute, None)
                if found is not None and found != "":
                    return found
    return getattr(dataset, attribute, None)
def enhanced_world_mappings(dataset, per_frame, shared):
    """Per-frame RWVM takes precedence over shared and root RWVM."""
    for group in (per_frame, shared, dataset):
        if group is None:
            continue
        mapping = getattr(group, "RealWorldValueMappingSequence", None)
        if mapping:
            return mapping
    return None
def enhanced_unit(mapping):
    units = getattr(mapping, "MeasurementUnitsCodeSequence", None)
    if not units:
        return None
    code_item = units[0]
    value = (getattr(code_item, "CodeValue", None)
             or getattr(code_item, "LongCodeValue", None)
             or getattr(code_item, "URNCodeValue", None))
    if not value:
        return None
    code = str(value).replace(" ", "").lower()
    if code in ("bq/ml", "bqml"):
        return "BQML", ""
    if code in ("g/ml{suvbw}", "g/ml{suvbodyweight}"):
        return "GML", "BW"
    if code in ("g/ml{suvlbm}", "g/ml{suvlbmmorgan}"):
        return "GML", "LBM"
    if code == "g/ml{suvlbmjames128}":
        return "GML", "LBMJAMES128"
    if code in ("g/ml{suvlbmjanma}", "g/ml{suvlbmjanmahasatian}"):
        return "GML", "LBMJANMA"
    if code == "g/ml{suvibw}":
        return "GML", "IBW"
    if code in ("cm2/ml{suvbsa}", "cm^2/ml{suvbsa}"):
        return "CM2ML", "BSA"
    return None
def enhanced_apply_world_mapping(pixel, mapping):
    """Apply the RWVM slope/intercept or a complete real-world-value LUT."""
    slope = safe_float(getattr(mapping, "RealWorldValueSlope", None))
    intercept = safe_float(getattr(mapping, "RealWorldValueIntercept", None))
    if slope is not None and intercept is not None:
        if not math.isfinite(slope) or not math.isfinite(intercept):
            raise ValueError("Invalid Enhanced PET real-world mapping parameters")
        return pixel * slope + intercept
    lut = getattr(mapping, "RealWorldValueLUTData", None)
    if lut is None:
        raise ValueError("Enhanced PET mapping lacks slope/intercept or LUT")
    first = getattr(mapping, "RealWorldValueFirstValueMapped", None)
    last = getattr(mapping, "RealWorldValueLastValueMapped", None)
    if first is None or last is None:
        first = getattr(mapping, "DoubleFloatRealWorldValueFirstValueMapped", None)
        last = getattr(mapping, "DoubleFloatRealWorldValueLastValueMapped", None)
    if first is None or last is None:
        raise ValueError("Enhanced PET LUT missing mapped-value bounds")
    first = float(first)
    last = float(last)
    lut_values = np.asarray(lut, dtype=np.float64).ravel()
    if (not np.isfinite(first) or not np.isfinite(last)
            or last < first or not float(first).is_integer()
            or lut_values.size != int(last - first + 1)
            or not np.all(np.isfinite(lut_values))):
        raise ValueError("Invalid Enhanced PET real-world-value LUT")
    indices = pixel - first
    if (np.any(indices < 0) or np.any(indices >= lut_values.size)
            or not np.all(indices == np.rint(indices))):
        raise ValueError("Enhanced PET voxel outside real-world-value LUT range")
    return lut_values[indices.astype(np.int64)]
def enhanced_mapping_for_frame(ds, per_frame, shared):
    mapping_seq = enhanced_world_mappings(ds, per_frame, shared)
    if mapping_seq:
        # A higher-priority mapping must have sufficient scaling information.
        valid = []
        for m in mapping_seq:
            unit = enhanced_unit(m)
            slope = safe_float(getattr(m, "RealWorldValueSlope", None))
            intercept = safe_float(getattr(m, "RealWorldValueIntercept", None))
            has_linear = (slope is not None and intercept is not None
                          and math.isfinite(slope) and math.isfinite(intercept))
            has_lut = (getattr(m, "RealWorldValueLUTData", None) is not None
                       and (getattr(m, "RealWorldValueFirstValueMapped", None) is not None
                            or getattr(m, "DoubleFloatRealWorldValueFirstValueMapped", None) is not None))
            if unit is not None and (has_linear or has_lut):
                valid.append((unit, m))
        # Prefer direct SUVbw, then another SUV type, then activity concentration.
        for wanted in ("GML_BW", "GML", "CM2ML", "BQML"):
            for unit, m in valid:
                if unit is None:
                    continue
                kind, suv_type = unit
                if (wanted == "GML_BW" and kind == "GML" and suv_type == "BW"
                        or wanted == "GML" and kind == "GML" and suv_type != "BW"
                        or wanted == kind and wanted in ("CM2ML", "BQML")):
                    return unit, m
    # IBSI fallback: Enhanced PET Pixel Value Transformation Sequence.
    for group in (per_frame, shared):
        seq = getattr(group, "PixelValueTransformationSequence", None) if group is not None else None
        if not seq:
            continue
        item = seq[0]
        kind = str(getattr(item, "RescaleType", "") or "").upper().strip()
        if kind not in ("BQML", "GML", "CM2ML"):
            continue
        slope = safe_float(getattr(item, "RescaleSlope", None))
        intercept = safe_float(getattr(item, "RescaleIntercept", None))
        if slope is None or intercept is None:
            continue
        suv_type = (str(getattr(ds, "SUVType", "") or "").upper().strip()
                    if kind == "GML" else "BSA" if kind == "CM2ML" else "")
        if kind == "GML" and suv_type == "":
            raise ValueError("Missing SUVType in Enhanced PET DICOM header")
        return (kind, suv_type), SimpleNamespace(RealWorldValueSlope=slope,
                                                RealWorldValueIntercept=intercept)
    raise ValueError("Enhanced PET lacks a usable real-world mapping or pixel transformation")
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SUBSECTION 13-6: ENHANCED PET SUVbw CONVERSION
# ============================================================================
def load_enhanced_pet_suvbw(pixels, source):
    """Return spatially sorted frame metadata, SUVbw volume, and original object.

    DICOM Per-frame Functional Groups Sequence items index the raw frames.
    The original object is retained for writing an Enhanced PET DICOM output.
    """
    if str(getattr(source, "SOPClassUID", "")) != ENHANCED_PET_SOP_UID:
        raise ValueError("Expected Enhanced PET SOP Class")
    pixels = np.asarray(pixels, dtype=np.float64)
    if pixels.ndim == 2:
        pixels = pixels[np.newaxis, :, :]
    if pixels.ndim != 3:
        raise ValueError("Enhanced PET pixel array must be [frame, row, column]")
    nframes = int(getattr(source, "NumberOfFrames", pixels.shape[0]))
    if nframes != pixels.shape[0]:
        raise ValueError("Enhanced PET NumberOfFrames does not match PixelData")
    per_frames = getattr(source, "PerFrameFunctionalGroupsSequence", None)
    if not per_frames or len(per_frames) != nframes:
        raise ValueError("Missing or inconsistent Per-Frame Functional Groups Sequence")
    shared_frames = getattr(source, "SharedFunctionalGroupsSequence", None)
    shared = shared_frames[0] if shared_frames else None
    metadata = []
    suv_frames = []
    for index in range(nframes):
        pf = per_frames[index]
        orient = enhanced_group_attribute(source, pf, shared, "ImageOrientationPatient")
        position = enhanced_group_attribute(source, pf, shared, "ImagePositionPatient")
        spacing = enhanced_group_attribute(source, pf, shared, "PixelSpacing")
        if orient is None or position is None or spacing is None:
            raise ValueError(f"Enhanced frame {index}: missing image geometry")
        unit_and_type, mapping = enhanced_mapping_for_frame(source, pf, shared)
        kind, suv_type = unit_and_type
        real_values = enhanced_apply_world_mapping(pixels[index], mapping)
        frame = SimpleNamespace(
            _frame_index=index,
            ImageOrientationPatient=orient,
            ImagePositionPatient=position,
            PixelSpacing=spacing,
            InstanceNumber=index + 1,
            PatientWeight=getattr(source, "PatientWeight", None),
            PatientSize=getattr(source, "PatientSize", None),
            PatientSex=getattr(source, "PatientSex", None),
            Units=kind, SUVType=suv_type,
            RadiopharmaceuticalInformationSequence=enhanced_group_attribute(
                source, pf, shared, "RadiopharmaceuticalInformationSequence"),
        )
        if kind == "BQML":
            # Enhanced PET uses full administration datetime, never time-only fallback.
            rp = get_radiopharm_item(frame)
            if parse_dicom_dt(getattr(rp, "RadiopharmaceuticalStartDateTime", None)) is None:
                raise ValueError("Enhanced PET missing full RadiopharmaceuticalStartDateTime")
            corrected = str(enhanced_group_attribute(source, pf, shared,
                                                    "DecayCorrected") or "").upper().strip()
            if corrected not in ("YES", "NO"):
                raise ValueError("Enhanced PET missing or invalid DecayCorrected")
            if corrected == "YES":
                reference = parse_dicom_dt(enhanced_group_attribute(
                    source, pf, shared, "DecayCorrectionDateTime"))
                if reference is None:
                    raise ValueError("Enhanced PET missing DecayCorrectionDateTime")
            else:
                reference = parse_dicom_dt(enhanced_group_attribute(
                    source, pf, shared, "FrameReferenceDateTime"))
                if reference is None:
                    start = parse_dicom_dt(enhanced_group_attribute(
                        source, pf, shared, "FrameAcquisitionDateTime"))
                    duration = safe_float(enhanced_group_attribute(
                        source, pf, shared, "FrameAcquisitionDuration"))
                    if start is None or duration is None or duration <= 0:
                        raise ValueError("Enhanced PET lacks measurement datetime and frame duration")
                    reference = start + timedelta(seconds=average_count_rate_time_s(
                        get_half_life_s(frame), duration / 1000.0))
            frame.DecayCorrection = "START" if corrected == "YES" else "NONE"
            admin = get_radiopharm_start_datetime(frame, reference_dt=reference)
            dose = get_total_dose_bq(frame)
            half_life = get_half_life_s(frame)
            dose_ref = dose * math.exp(-math.log(2.0) / half_life *
                                      (reference - admin).total_seconds())
            if not math.isfinite(dose_ref) or dose_ref <= 0:
                raise ValueError("Enhanced PET invalid decay-corrected dose")
            result = real_values * (get_weight_g(frame) / dose_ref)
        elif kind == "GML":
            frame.DecayCorrection = str(getattr(source, "DecayCorrection", "") or "")
            result = convert_gml_to_suvbw(frame, real_values)
        elif kind == "CM2ML":
            frame.DecayCorrection = str(getattr(source, "DecayCorrection", "") or "")
            result = convert_cm2ml_to_suvbw(frame, real_values)
        else:
            raise ValueError(f"Unsupported Enhanced PET units: {kind}")
        if not np.all(np.isfinite(result)):
            raise ValueError("Non-finite SUVbw in Enhanced PET frame")
        metadata.append(frame)
        suv_frames.append(result)
    # Geometry sorting must reorder pixel frames and metadata together.
    sorted_metadata = sort_pet_slices(metadata)
    volume = np.stack([suv_frames[frame._frame_index] for frame in sorted_metadata], axis=0)
    return sorted_metadata, volume, source
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# Section 14 - PUBLIC API
# ============================================================================
def radiuma_suv(input_array, input_metadata, parameter_estimation="OFF"):
    """Convert in-memory PET stored pixels to SUVbw, returning (array, metadata).

    Metadata OFF: returned dictionary is a deep copy of the input dictionary.
    Metadata ON: missing DICOM values are filled from existing configured defaults
    and the fields actually filled are returned.
    Output array contains floating-point SUVbw, not uint16-encoded pixels.
    """
    if parameter_estimation not in ("ON", "OFF"):
        raise ValueError('parameter_estimation must be "ON" or "OFF"')
    if not isinstance(input_array, np.ndarray):
        raise TypeError("input_array must be a numpy.ndarray")
    if input_array.ndim not in (2, 3) or not np.issubdtype(input_array.dtype, np.number):
        raise ValueError("PET input must be a numeric 2D or 3D NumPy array")
    if not np.all(np.isfinite(input_array)):
        raise ValueError("PET input array contains non-finite values")
    if not isinstance(input_metadata, dict):
        raise TypeError("input_metadata must be a DICOM metadata dictionary")
    settings = copy.deepcopy(_DEFAULT_PARAMETERS)
    array = input_array.astype(np.float64)
    if array.ndim == 2:
        array = array[np.newaxis, :, :]
    if "slices" not in input_metadata:
        ds = _dicom_dataset(input_metadata)
        if str(getattr(ds, "SOPClassUID", "")) == ENHANCED_PET_SOP_UID:
            prepared = apply_configured_pet_metadata(ds, parameter_estimation, settings)
            _, suv, _ = load_enhanced_pet_suvbw(array, prepared)
            output_metadata = (copy.deepcopy(input_metadata) if parameter_estimation == "OFF"
                               else _metadata_with_insertions(input_metadata, prepared))
            return (suv[0] if input_array.ndim == 2 else suv), output_metadata
    datasets, slopes, intercepts, output_metadata = _prepare_conventional(
        array, input_metadata, parameter_estimation, settings
    )
    suv = convert_to_suvbw(datasets, array, slopes, intercepts)
    if not np.all(np.isfinite(suv)):
        raise ValueError("SUVbw conversion produced non-finite pixel values")
    return (suv[0] if input_array.ndim == 2 else suv), output_metadata
def main(input_array, input_metadata, parameter_estimation="OFF"):
    """Optional in-memory launcher; equivalent to calling radiuma_suv()."""
    return radiuma_suv(input_array, input_metadata, parameter_estimation)
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# Section 15 - OUTPUT GENERATION
# ============================================================================
def _new_suv_output_directory(input_directory):
    """Choose a new sibling _processed path without overwriting prior results."""
    source = os.path.abspath(os.fspath(input_directory))
    if os.path.isfile(source):
        base = os.path.join(os.path.dirname(source),
                            os.path.splitext(os.path.basename(source))[0] + "_processed.dcm")
        stem, ext = os.path.splitext(base)
    else:
        base = source.rstrip(os.sep) + "_processed"
        stem, ext = base, ""
    target = base
    index = 2
    while os.path.exists(target):
        target = f"{stem}_{index}{ext}"
        index += 1
    return target
# Import-time diagnostic removed for library compatibility.

# ==========================================================================
# Section 16 - TYPE-PRESERVING SUVbw PUBLIC INPUT / OUTPUT API
# ==========================================================================
# The original quantitative radiuma_suv() engine and image reading functions
# above are preserved. Only the outer input/output interface is replaced.
# DICOM cannot generally store arbitrary floating-point SUVbw pixels within
# the original PET Image SOP Class. Output uses signed 16-bit quantized pixels
# and an explicit SUVbw calibration (rescale for conventional PET; RWVM for
# Enhanced PET). The on-disk values are reproducible within half a quantization
# step; the unchanged in-memory engine computes float64 SUVbw.

# ============================================================================
# SUBSECTION 16-1: OUTPUT DICOM METADATA PREPARATION
# ============================================================================
def _suv_prepare_output_header(ds, estimation):
    """Fill missing source tags only in ON mode, never modify source on disk."""
    if estimation not in ("ON", "OFF"):
        raise ValueError('parameter_estimation must be "ON" or "OFF"')
    return apply_configured_pet_metadata(ds, estimation, copy.deepcopy(_DEFAULT_PARAMETERS))

# ============================================================================
# SUBSECTION 16-2: QUANTITATIVE SUVbw PIXEL ENCODING
# ============================================================================
def _suv_int_pixels(suv):
    """Reversibly calibratable DICOM integer representation of float SUVbw."""
    arr = np.asarray(suv, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise ValueError("SUVbw output contains non-finite values")
    if arr.ndim not in (2, 3):
        raise ValueError("Expected a single 2D frame or an Enhanced 3D frame stack")
    # Use 32-bit signed, preserving dynamic range and minimizing quantization.
    safe_limit = 32760
    largest = float(np.max(np.abs(arr))) if arr.size else 0.0
    # A positive slope even for an all-zero image.
    slope = largest / safe_limit if largest > 0 else 1.0 / safe_limit
    encoded = np.rint(arr / slope).astype('<i2')
    if not np.allclose(encoded.astype(np.float64) * slope, arr,
                       rtol=0, atol=slope * 0.50001):
        raise ValueError("DICOM SUVbw quantization validation failed")
    return encoded, slope
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SUBSECTION 16-2: DERIVED DICOM IDENTIFIERS AND PIXEL STORAGE
# ============================================================================
def _suv_new_identity(ds, output_series_uid):
    """Create derived-series/object identity and document SUV derivation."""
    from pydicom.uid import generate_uid
    ds.SeriesInstanceUID = output_series_uid
    ds.SOPInstanceUID = generate_uid()
    if hasattr(ds, 'file_meta') and ds.file_meta is not None:
        ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
        ds.file_meta.MediaStorageSOPClassUID = ds.SOPClassUID
    original_types = [str(v) for v in getattr(ds, 'ImageType', [])]
    ds.ImageType = ['DERIVED', 'PRIMARY'] + original_types[2:]
    ds.DerivationDescription = (
        'radiuma_suv: body-weight normalized standardized uptake value '
        '(SUVbw) calculated from original stored PET pixels and DICOM metadata. '
        'Pixel values integer-encoded with DICOM calibration.'
    )
    # Original series identifiers and patient identifiers remain on the
    # input files; output files are new derived SOP instances.
def _suv_set_pixel_storage(ds, integer_pixels):
    from pydicom.uid import ExplicitVRLittleEndian
    # Standard uncompressed integer pixel representation.
    ds.Rows = int(integer_pixels.shape[-2])
    ds.Columns = int(integer_pixels.shape[-1])
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 1
    ds.PixelData = np.ascontiguousarray(integer_pixels, dtype='<i2').tobytes()
    ds['PixelData'].is_undefined_length = False
    if getattr(ds, 'file_meta', None) is None:
        from pydicom.dataset import FileMetaDataset
        ds.file_meta = FileMetaDataset()
    ds.file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    ds.is_little_endian = True
    ds.is_implicit_VR = False
    # Encapsulated compressed-frame helpers no longer correspond to PixelData.
    for keyword in ('ExtendedOffsetTable', 'ExtendedOffsetTableLengths',
                    'FloatPixelData', 'DoubleFloatPixelData'):
        if keyword in ds:
            del ds[keyword]
    if integer_pixels.ndim == 3:
        ds.NumberOfFrames = str(integer_pixels.shape[0])
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SUBSECTION 16-4: CONVENTIONAL PET DICOM OUTPUT
# ============================================================================
def _suv_write_conventional(src_ds, suv_frame, dst, estimation, series_uid):
    ds = _suv_prepare_output_header(src_ds, estimation)
    integer, slope = _suv_int_pixels(suv_frame)
    _suv_set_pixel_storage(ds, integer)
    ds.RescaleSlope = str(slope)
    ds.RescaleIntercept = '0'
    ds.Units = 'GML'
    ds.SUVType = 'BW'
    _suv_new_identity(ds, series_uid)
    ds.save_as(dst, enforce_file_format=True)
    _suv_verify_saved(dst, suv_frame, slope)
    return dst
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SUBSECTION 16-5: ENHANCED PET DICOM MAPPING HELPERS
# ============================================================================
def _suv_remove_world_mappings(group):
    if group is not None and 'RealWorldValueMappingSequence' in group:
        del group.RealWorldValueMappingSequence
def _suv_set_identity_pvts(group):
    """Prevent an original transformation from rescaling SUV pixels again."""
    if group is None:
        return
    seq = getattr(group, 'PixelValueTransformationSequence', None)
    if seq:
        for item in seq:
            item.RescaleSlope = '1'
            item.RescaleIntercept = '0'
            item.RescaleType = 'GML'
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SUBSECTION 16-6: ENHANCED PET DICOM OUTPUT
# ============================================================================
def _suv_write_enhanced(src_ds, suv_sorted, dst, estimation, series_uid):
    """Preserve the original frame count/order and publish per-frame SUV RWVM."""
    from pydicom.dataset import Dataset
    from pydicom.sequence import Sequence
    ds = _suv_prepare_output_header(src_ds, estimation)
    per_frames = getattr(ds, 'PerFrameFunctionalGroupsSequence', None)
    if not per_frames:
        raise ValueError('Enhanced PET lacks PerFrameFunctionalGroupsSequence')
    shared_seq = getattr(ds, 'SharedFunctionalGroupsSequence', None)
    shared = shared_seq[0] if shared_seq else None
    positions = []
    orientations = []
    for frame in per_frames:
        pos = enhanced_group_attribute(ds, frame, shared, 'ImagePositionPatient')
        ori = enhanced_group_attribute(ds, frame, shared, 'ImageOrientationPatient')
        if pos is None or ori is None:
            raise ValueError('Missing Enhanced PET spatial geometry')
        positions.append(np.asarray(pos, dtype=float))
        orientations.append(np.asarray(ori, dtype=float))
    first_orientation = orientations[0]
    if any(not np.allclose(o, first_orientation, atol=1e-5)
           for o in orientations):
        raise ValueError('Enhanced PET orientations disagree')
    normal = np.cross(first_orientation[:3], first_orientation[3:])
    if np.linalg.norm(normal) == 0:
        raise ValueError('Invalid Enhanced PET orientation')
    ordered_indices = sorted(range(len(per_frames)),
                             key=lambda i: float(np.dot(positions[i], normal)))
    if len(ordered_indices) != suv_sorted.shape[0]:
        raise ValueError('Frame count does not match SUV conversion')
    # radiuma_suv's Enhanced PET engine sorts by spatial coordinate. Undo
    # that sort to preserve the input DICOM's PerFrameFunctionalGroups order.
    original_order = [None] * len(ordered_indices)
    for sorted_idx, input_idx in enumerate(ordered_indices):
        original_order[input_idx] = suv_sorted[sorted_idx]
    originals = np.stack(original_order, axis=0)
    integer, slope = _suv_int_pixels(originals)
    _suv_set_pixel_storage(ds, integer)
    # Remove old real-world mappings at every inheritance level and replace
    # with an unambiguous per-frame SUVbw mapping.
    _suv_remove_world_mappings(ds)
    _suv_remove_world_mappings(shared)
    _suv_set_identity_pvts(shared)
    for frame in per_frames:
        _suv_remove_world_mappings(frame)
        _suv_set_identity_pvts(frame)
        item = Dataset()
        item.RealWorldValueSlope = float(slope)
        item.RealWorldValueIntercept = 0.0
        item.RealWorldValueFirstValueMapped = int(integer.min())
        item.RealWorldValueLastValueMapped = int(integer.max())
        unit = Dataset()
        unit.CodeValue = 'g/ml{SUVbw}'
        unit.CodingSchemeDesignator = 'UCUM'
        unit.CodeMeaning = 'Standardized Uptake Value body weight'
        item.MeasurementUnitsCodeSequence = Sequence([unit])
        frame.RealWorldValueMappingSequence = Sequence([item])
    ds.Units = 'GML'
    ds.SUVType = 'BW'
    ds.RescaleSlope = '1'
    ds.RescaleIntercept = '0'
    _suv_new_identity(ds, series_uid)
    ds.save_as(dst, enforce_file_format=True)
    _suv_verify_saved(dst, originals, slope)
    return dst
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SUBSECTION 16-7: OUTPUT DICOM VALIDATION
# ============================================================================
def _suv_verify_saved(path, expected, slope):
    """Read the written DICOM back to catch corrupted pixel encoding."""
    ds = pydicom.dcmread(path)
    decoded = np.asarray(ds.pixel_array)
    target = np.asarray(expected)
    if decoded.shape != target.shape:
        raise ValueError(f'Written DICOM shape mismatch: {path}')
    if not np.all(np.abs(decoded.astype(np.float64) * slope - target) <= slope * .50001):
        raise ValueError(f'Written DICOM SUVbw values do not match expected: {path}')
# Import-time diagnostic removed for library compatibility.

# ============================================================================
# SUBSECTION 16-8: DICOM SERIES OUTPUT WRITER
# ============================================================================
def _suv_processed_filename(source_path):
    """Keep the source filename and append _processed before the DICOM extension."""
    filename = os.path.basename(os.fspath(source_path))
    stem, extension = os.path.splitext(filename)
    if extension.lower() == '.dcm':
        return stem + '_processed' + extension
    return filename + '_processed.dcm'

def _suv_load_and_write_group(paths, destination, estimation):
    """One source DICOM PET volume/series to matching derived DICOM objects."""
    from pydicom.uid import generate_uid
    os.makedirs(destination, exist_ok=False)
    input_headers = [pydicom.dcmread(p) for p in paths]
    if not input_headers or any(str(getattr(ds, 'Modality', '')) != 'PT'
                                for ds in input_headers):
        raise ValueError('Every input must be a PET DICOM object (Modality=PT)')
    enhanced = [str(getattr(ds, 'SOPClassUID', '')) == ENHANCED_PET_SOP_UID
                for ds in input_headers]
    if any(enhanced) and not all(enhanced):
        raise ValueError('Cannot combine Enhanced and conventional PET in one series')
    series_uid = generate_uid()
    output_paths = []
    if all(enhanced):
        for index, (path, ds) in enumerate(zip(paths, input_headers), start=1):
            header = _pet_header_dictionary(ds)
            suv, _ = radiuma_suv(np.asarray(ds.pixel_array), header,
                                 parameter_estimation=estimation)
            if suv.ndim == 2:
                suv = suv[np.newaxis, :, :]
            # Preserve one output DICOM per source multi-frame DICOM object.
            dst = os.path.join(destination, _suv_processed_filename(path))
            output_paths.append(_suv_write_enhanced(ds, suv, dst,
                                                     estimation, series_uid))
    else:
        if any(int(getattr(ds, 'NumberOfFrames', 1) or 1) != 1
               for ds in input_headers):
            raise ValueError('Multi-frame non-Enhanced PET is unsupported')
        for ds in input_headers:
            if not hasattr(ds, 'ImagePositionPatient') or not hasattr(ds, 'ImageOrientationPatient'):
                raise ValueError('Conventional series geometry is required')
        source_paths_by_header = {id(ds): path for ds, path in zip(input_headers, paths)}
        ordered = sort_pet_slices(input_headers)
        raw = np.stack([np.asarray(ds.pixel_array) for ds in ordered])
        headers = {'slices': [_pet_header_dictionary(ds) for ds in ordered]}
        suv, _ = radiuma_suv(raw, headers, parameter_estimation=estimation)
        for index, (ds, suv_frame) in enumerate(zip(ordered, suv), start=1):
            dst = os.path.join(destination, _suv_processed_filename(source_paths_by_header[id(ds)]))
            output_paths.append(_suv_write_conventional(
                ds, suv_frame, dst, estimation, series_uid))
    print(f"[radiuma-suv] PET DICOM series conversion completed successfully. Files saved: {len(output_paths)}, Output directory: {destination}")
    return output_paths

# ============================================================================
# SUBSECTION 16-9: FORMAT-PRESERVING PUBLIC INPUT OUTPUT API
# ============================================================================
def radiuma_suv_from_input(image_input, metadata_json=None, parameter_estimation='OFF',
                           output_path=None):
    """Convert PET to SUVbw, preserving the input *representation*.

    Input type                         Returned output
    -------------------------------    -------------------------------
    NumPy ndarray + DICOM metadata     (float64 SUVbw ndarray, metadata dict)
    Single-frame PET DICOM file        (new DICOM file path, output header dict)
    Multi-frame Enhanced PET file      (new DICOM file path, output header dict)
    Folder of conventional PET         (new DICOM folder path, list of headers)
    Folder of Enhanced PET files       (new DICOM folder path, list of headers)
    Folder with independent series     (new root folder path, nested headers)

    Source DICOMs are never overwritten. The original patient/study/geometry
    content is copied. Output header changes required for valid derived SUVbw
    include pixel encoding, calibration, derived image/series UID and provenance.
    Missing parameters are added from defaults only when estimation is ON.

    DICOM pixel values are quantized to signed 16-bit integer data; use the
    output DICOM RescaleSlope (conventional) or RealWorldValueMappingSequence
    (Enhanced PET) to retrieve physical SUVbw values.
    """
    if parameter_estimation not in ('OFF', 'ON'):
        raise ValueError('parameter_estimation must be "ON" or "OFF"')
    if isinstance(image_input, np.ndarray):
        if metadata_json is None:
            raise ValueError('NumPy PET input requires associated DICOM metadata_json')
        if output_path is not None:
            raise ValueError('NumPy input returns array and metadata in memory; output_path is for DICOM')
        return radiuma_suv(image_input, metadata_json,
                           parameter_estimation=parameter_estimation)
    if not isinstance(image_input, (str, os.PathLike)):
        raise TypeError('image_input must be a NumPy array or DICOM file/folder path')
    src = os.fspath(image_input)
    if not (os.path.isfile(src) or os.path.isdir(src)):
        raise FileNotFoundError(src)
    if os.path.isfile(src) and src.lower().endswith(('.nii', '.nii.gz', '.nrrd', '.nhdr', '.npy')):
        raise ValueError('PET path inputs must be DICOM; for NumPy pass an array and metadata')
    source_groups = _suv_detect_source_groups(src)
    if os.path.isfile(src) and len(source_groups) != 1:
        raise ValueError('A single file cannot contain multiple source groups')
    if output_path is None:
        output_path = _new_suv_output_directory(src)
    out_root = os.path.abspath(os.fspath(output_path))
    if output_path is not None:
        if os.path.isfile(src):
            parent, filename = os.path.split(out_root)
            stem, extension = os.path.splitext(filename)
            if extension.lower() == '.dcm':
                if not stem.endswith('_processed'):
                    out_root = os.path.join(parent, stem + '_processed' + extension)
            elif not filename.endswith('_processed'):
                out_root += '_processed'
        elif not os.path.basename(out_root.rstrip(os.sep)).endswith('_processed'):
            out_root = out_root.rstrip(os.sep) + '_processed' 
    if os.path.exists(out_root):
        raise FileExistsError(f'Refusing to overwrite existing output: {out_root}')
    if os.path.isfile(src):
        # Preserve the single-file output representation.
        if not out_root.lower().endswith('.dcm'):
            out_root += '.dcm'
        if os.path.exists(out_root):
            raise FileExistsError(out_root)
        from pydicom.uid import generate_uid
        ds = pydicom.dcmread(src)
        if str(getattr(ds, 'Modality', '')) != 'PT':
            raise ValueError('Input must be PET DICOM (Modality=PT)')
        if str(getattr(ds, 'SOPClassUID', '')) == ENHANCED_PET_SOP_UID:
            suv, _ = radiuma_suv(np.asarray(ds.pixel_array),
                                 _pet_header_dictionary(ds),
                                 parameter_estimation=parameter_estimation)
            if suv.ndim == 2:
                suv = suv[np.newaxis, ...]
            _suv_write_enhanced(ds, suv, out_root, parameter_estimation, generate_uid())
        else:
            if int(getattr(ds, 'NumberOfFrames', 1) or 1) != 1:
                raise ValueError('Non-Enhanced multi-frame PET is unsupported')
            pixels = np.asarray(ds.pixel_array)
            suv, _ = radiuma_suv(pixels, _pet_header_dictionary(ds),
                                 parameter_estimation=parameter_estimation)
            _suv_write_conventional(ds, suv, out_root,
                                    parameter_estimation, generate_uid())
        return out_root, _pet_header_dictionary(pydicom.dcmread(out_root,
                                                                stop_before_pixels=True))
    # Directory input returns a directory of the same number/types of images.
    # Check input/output separation before writing anything.
    if os.path.commonpath([os.path.abspath(src), out_root]) == os.path.abspath(src):
        raise ValueError('Output directory must not be nested within input folder')
    # Preserve relative directory paths for input folders at any depth.
    from collections import Counter, defaultdict
    relative_folders = []
    for files in source_groups:
        parents = {os.path.dirname(os.path.abspath(path)) for path in files}
        if len(parents) != 1:
            raise ValueError('A source DICOM series spans multiple folders')
        relative = os.path.relpath(parents.pop(), os.path.abspath(src))
        if relative.startswith(os.pardir):
            raise ValueError('A source DICOM file is outside the input directory')
        relative_folders.append(relative)

    try:
        folder_counts = Counter(relative_folders)
        folder_indices = defaultdict(int)
        results = []
        for files, relative_folder in zip(source_groups, relative_folders):
            if relative_folder == os.curdir:
                target = out_root
            else:
                target = os.path.join(out_root, *(
                    part + '_processed' for part in relative_folder.split(os.sep)))
            if folder_counts[relative_folder] > 1:
                folder_indices[relative_folder] += 1
                target = os.path.join(
                    target, f'series_{folder_indices[relative_folder]:04d}_processed')
            # _suv_load_and_write_group creates its destination, so a
            # single-series root requires the existing output root as parent.
            if target == out_root and not os.path.exists(target):
                pass
            output_files = _suv_load_and_write_group(
                files, target, parameter_estimation)
            results.append([_pet_header_dictionary(pydicom.dcmread(
                p, stop_before_pixels=True)) for p in output_files])
        if len(results) == 1:
            return out_root, results[0]
        return out_root, results
    except Exception:
        # Do not leave partial derived series that might be mistaken for valid.
        import shutil
        shutil.rmtree(out_root, ignore_errors=True)
        print("[radiuma-suv] NumPy PET-to-SUVbw conversion completed successfully (in-memory output).")
        raise

# ============================================================================
# SUBSECTION 16-10: OBJECT-ORIENTED PUBLIC WRAPPER
# ============================================================================
class RadiumaSUVLoader:
    """convenience wrapper for representation-preserving output."""
    def __init__(self, parameter_estimation='OFF', pet_metadata=None,
                 output_path=None):
        self.parameter_estimation = parameter_estimation
        self.pet_metadata = pet_metadata
        self.output_path = output_path
    def convert(self, image_input):
        return radiuma_suv_from_input(
            image_input,
            parameter_estimation=self.parameter_estimation,
            metadata_json=self.pet_metadata,
            output_path=self.output_path,
        )

# ============================================================================
# SUBSECTION 16-11: EXECUTE PET DICOM SUVbw CONVERSION
# ============================================================================
if __name__ == "__main__":
    from .cli import main as cli_main
    raise SystemExit(cli_main())


