# radiuma-suv: PET to SUVbw Conversion – Python Script and Library

[PyPi Link](https://pypi.org/project/radiuma-suv/)

[GitHub Link](https://github.com/radiuma-com/radiuma-suv) 

**radiuma-suv** is the proposed library name for the PET-to-SUVbw conversion engine in `src/radiuma-suv/core.py`. The code converts PET pixels to **body-weight-normalized standardized uptake values (SUVbw)** using DICOM acquisition and radiopharmaceutical metadata. It supports conventional single-frame PET, Enhanced multi-frame PET, NumPy arrays with matching DICOM metadata, and folders containing one or several PET series. In-memory results are floating-point arrays; DICOM results are newly written, calibrated derived images.

**Code Overview:** `radiuma-suv` is a Python-based PET-to-SUVbw conversion tool that converts PET image pixel values into body-weight-normalized Standardized Uptake Values (SUVbw) using DICOM metadata. It supports conventional single-frame PET, Enhanced multi-frame PET, NumPy arrays, and DICOM directories containing multiple PET series.

The code provides three main interfaces:

- **`radiuma-suv(input_array, input_metadata, parameter_estimation="OFF")`**: Converts a 2D/3D PET NumPy array into a floating-point SUVbw array and returns the corresponding metadata dictionary.
- **`radiuma-suv_from_input(image_input, parameter_estimation="OFF", output_path=None)`**: Automatically detects supported PET inputs and produces SUVbw outputs in the corresponding format (NumPy array or derived DICOM file/directory).
- **`RadiumaSUVLoader(parameter_estimation="OFF", pet_metadata=None, output_path=None)`**: Provides an object-oriented interface for PET-to-SUVbw conversion through the `.convert(image_input)` method.

**Parameter Estimation:** `"OFF"` uses existing DICOM metadata without filling missing fields, while `"ON"` fills eligible missing metadata using predefined configuration values. Original DICOM files are not overwritten.

**Scientific caution:** `parameter_estimation="ON"` can fill missing metadata with notebook-configured values, including patient weight, radionuclide activity, timing, and other acquisition details. Those values must not be assumed correct for another patient. Prefer `"OFF"` with complete verified source DICOM metadata. No clinical validation is established by the notebook.

## 🔍 Table of Contents

- [🎥 YouTube Tutorials](#-youtube-tutorials)
- [🧩 IBSI Standardization & Colibration](#-ibsi-standardization)
- [✨ Key Features](#-key-features)
- [📚 Library Usage](#-library-usage)
  - [📂 Single File Processing](#-single-file-processing)
  - [🧠 In-Memory Array Processing](#-in-memory-array-processing)
  - [🔧 Advanced Configuration](#-advanced-configuration)
- [📥 Installation](#-installation)
- [📂 Data Structure Requirements](#-data-structure-requirements)
- [📋 radiuma-suv Parameters Reference](#-radiuma-suv-parameters-reference)
- [📊 Parameter Compatibility](#-parameter-compatibility)
- [📚 API Reference](#-api-reference)
- [📊 Output Structure](#-output-structure)
- [📁 Supported File Formats](#-supported-file-formats)
- [🎯 Library Examples](#-library-examples)
- [⚡ Performance Tips](#-performance-tips)
- [💻 Integration of radiuma-suv to GUI](#-integration-of-radiuma-suv-to-gui)
- [❓ Troubleshooting](#-troubleshooting)
- [🕒 Version History](#-version-history)
- [📬 Contact](#-contact)
- [👥 Authors](#-authors)
- [📚 Citation](#-citation)
- [🙏 Acknowledgment](#acknowledgment)
- [📜 License](#-license)

## 🎥 YouTube Tutorials

[Watch our tutorial videos on our YouTube channel](https://www.youtube.com/@Radiuma-software)
Video 1 (coming soon)


## 🧩 IBSI Standardization & Colibration

We have [three versions of the software/code](https://github.com/radiuma-com/radiuma-suv/tree/main/IBSI-Calibration/IBSI-Validated-Code-Versions) for [IBSI Compliance-SUV Calibrations](https://theibsi.github.io/ibsi-suv/).
-	[Version 1](https://github.com/radiuma-com/radiuma-suv/blob/main/IBSI-Calibration/IBSI-Validated-Code-Versions/Version1.ipynb): before the start of IBSI-SUV
-	[Version 2](https://github.com/radiuma-com/radiuma-suv/blob/main/IBSI-Calibration/IBSI-Validated-Code-Versions/Version2.ipynb): Calibrated with the first standardization attempt and with the first 17 DROs.  
-	[Version 3](https://github.com/radiuma-com/radiuma-suv/blob/main/IBSI-Calibration/IBSI-Validated-Code-Versions/Version3.ipynb): Calibrated with the second standardization attempt and with all 58 DROs.
IBSI calibration datasets are available for [download](https://github.com/oncoray/suv_computation)

Report of Calibrations

**Calibration 1 (Just the First 17 DROs)**

- Use Code Version 2 and the first 17 DROs.
  - [Result 1, called Version 2_Calibration 1.csv](https://github.com/radiuma-com/radiuma-suv/blob/main/IBSI-Calibration/IBSI-Validated-Results/Result%201%2C%20called%20Version%202_Callibration%201.csv)
    - **Note:** We already submitted this result, but we include it in this collection for your easy access.
          

**Calibration 2 (All 58 DROs)**

- Use Code Version 3 and all 58 DROs. 
  - [Result 2, called Version 3_Callibration 2.csv](https://github.com/radiuma-com/radiuma-suv/blob/main/IBSI-Calibration/IBSI-Validated-Results/Result%201%2C%20called%20Version%202_Callibration%201.csv)”

 - Use Code Version 2 and all 58 DROs. 
   - [Result 3, called Version 2_Callibration 2.csv](https://github.com/radiuma-com/radiuma-suv/blob/main/IBSI-Calibration/IBSI-Validated-Results/Result%202%2C%20called%20Version%203_Callibration%202.csv)”

 - Use Code Version 1 and all 58 DROs. 
   - Result 4, called “[Version 1_Callibration 2.csv](https://github.com/radiuma-com/radiuma-suv/blob/main/IBSI-Calibration/IBSI-Validated-Results/Result%203%2C%20called%20Version%202_Callibration%202.csv)”

**Calibration 3 ( All 16 Real Patients: Final Validation Phase)**
 - Use Code Version 3 and all 16 real patients (Validation phase). 
   - [Result 5, called Version 3_Callibration 3.csv](https://github.com/radiuma-com/radiuma-suv/blob/main/IBSI-Calibration/IBSI-Validated-Results/Result%204%2C%20called%20Version%201_Callibration%202.csv)”

 - Use Code Version 2 and all 16 real patients (Validation phase). 
   - [Result 6, called “Version 2_Callibration 3.csv](https://github.com/radiuma-com/radiuma-suv/blob/main/IBSI-Calibration/IBSI-Validated-Results/Result%206%2C%20called%20Version%202_Callibration%203.csv)

 - Use Code Version 1 and all 16 real patients (Validation phase). 
   - [Result 7, called “Version 1_Callibration 3.csv](https://github.com/radiuma-com/radiuma-suv/blob/main/IBSI-Calibration/IBSI-Validated-Results/Result%207%2C%20called%20Version%201_Callibration%203.csv)


## ✨ Key Features

- **PET-to-SUVbw conversion** using PET image pixels and matching DICOM metadata.
- **Conventional PET DICOM**: single-frame PET images and conventional PET slice series.
- **Enhanced PET DICOM**: multi-frame images with functional group metadata.
- **NumPy arrays**: convert in memory using a corresponding DICOM metadata dictionary.
- **Hierarchical PET directories**: recursively identify PET DICOM groups and process independent series separately.
- **Type-preserving API**: NumPy input returns a NumPy array; DICOM input creates derived DICOM outputs.
- **Optional metadata completion**: `parameter_estimation="ON"` fills selected missing values using notebook settings.
- **Derived DICOM output**: new DICOM identifiers, calibrated signed 16-bit pixels, and output validation against computed SUVbw values.
- **Source protection**: output paths are checked to avoid overwriting original DICOM files.
- **Object-oriented convenience wrapper**: `RadiumaSUVLoader`.

These are implementation capabilities, not independent quantitative or clinical validation claims.

### 📦 Library Installation via pip

After publishing the intended package under the distribution name `radiuma-suv`:

```bash
python -m pip install radiuma-suv
```

For local development in the directory containing `pyproject.toml`:

```bash
python -m pip install -e .
```

## 📚 Library Usage

### 💻 Script Execution

The `radiuma-suv` code supports two execution modes: **Standalone Mode** and **Embedded Mode**.

#### 1. Standalone Mode

Run the Python script directly to convert PET DICOM images into SUVbw.

**Step 1: Choose the input path and options at the CLI (do not edit the installed library):**

```bash
radiuma-suv PATH_TO_PET_DICOM_FILE_OR_FOLDER --parameter-estimation OFF
```

**Step 2: Alternatively execute using Python:**

```bash
python -m radiuma-suv PATH_TO_PET_DICOM_FILE_OR_FOLDER --parameter-estimation OFF
```

**Output:** Creates calibrated SUVbw DICOM files in a new `_processed` location without overwriting original images.

#### 2. Embedded Mode

Import the conversion functions into another Python application, pipeline, or GUI.

**Option A: DICOM file or directory**

```python
from radiuma-suv import radiuma-suv_from_input

output, metadata = radiuma-suv_from_input(
    image_input="PET_Input",
    parameter_estimation="OFF",
    output_path="SUV_Output"
)
```

**Option B: NumPy array with DICOM metadata**

```python
from radiuma-suv import radiuma-suv

suv_array, metadata = radiuma-suv(
    input_array=pet_array,
    input_metadata=pet_metadata,
    parameter_estimation="OFF"
)
```

**Option C: Object-oriented integration**

```python
from radiuma-suv import RadiumaSUVLoader

converter = RadiumaSUVLoader(parameter_estimation="OFF")
output, metadata = converter.convert("PET_Input")
```

**Note:** DICOM inputs automatically provide metadata through their headers. NumPy inputs require a corresponding DICOM metadata dictionary.
`parameter_estimation` controls metadata replacement:

| Option | Behavior |
|---|---|
| `"OFF"` | Does not fill missing parameters from the configured notebook defaults. |
| `"ON"` |  Fills eligible missing metadata fields using configured defaults; existing values are not generally replaced. |

The notebook has configurable values for patient weight, age, height, sex; administered dose, half-life and radiopharmaceutical timing; acquisition and series times; PET units; SUV type; decay correction; rescale slope and intercept; image geometry; and other fields. These settings are **not additional keyword arguments** of `radiuma-suv_from_input`.

The source implements unit conversion paths involving `BQML`, `GML`, `CM2ML`, and `CNTS`. Successful processing depends on the necessary metadata and calibration factors. The `BQML` converter handles `DecayCorrection` values `ADMIN`, `START`, and `NONE` with different reference-time logic. The source also includes body-composition helper functions used for certain `GML` SUV normalization variants.


## 📥 Installation

Dependencies imported by the notebook are `numpy`, `pydicom`, and `SimpleITK`; other referenced Python modules are standard-library modules. In the original notebook the installation line is `%pip install numpy pydicom SimpleITK`; for a published library, declare these in `pyproject.toml` instead.

### 🌐 GitHub Installation

https://git.talkfluent.ir/admin/radiuma-suv/src/branch/main

```bash
git clone https://github.com/YOUR_USERNAME/radiuma-suv.git
cd radiuma-suv
python -m pip install -e .
```

### macOS/Linux Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

### Windows Setup

In PowerShell inside the project directory:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e .
```

### 💻 Python Script - Command Line Interface (CLI)

The notebook contains an example `if __name__ == "__main__":` block but **does not implement an argument-parsing CLI or a `radiuma-suv` terminal command**. Create a short Python script when you want command-line execution:

```python
# convert_pet.py
from radiuma-suv import radiuma-suv_from_input
path, headers = radiuma-suv_from_input("pet_scan.dcm", parameter_estimation="OFF")
print(path)
```

```bash
python convert_pet.py
```

### 📦 Library Installation via pip

```bash
python -m pip install radiuma-suv
```

This will install the intended project **only after it is actually published** under that name on PyPI.

## 📂 Data Structure Requirements

The notebook processes PET DICOM input **without masks**. Examples of appropriate layouts:

**Single PET image:**

```text
PET/
└── scan.dcm
```

**Conventional PET image series:**

```text
PET_Series/
├── slice_001.dcm
├── slice_002.dcm
└── slice_003.dcm
```

**Hierarchical directory:**

```text
Patients/
├── patient_001/
│   └── PET/
│       ├── slice_001.dcm
│       └── slice_002.dcm
└── patient_002/
    └── PET/
        └── enhanced_pet.dcm
```

The source searches directories recursively, reads DICOM headers, selects `Modality=PT`, and groups series by their `SeriesInstanceUID` and folder. A DICOM series spanning different source folders is rejected by the writer. The corresponding output folders preserve relative structure with `_processed` naming.

**Output must not be inside the source directory.** Existing destinations are not overwritten. Non-DICOM files may be skipped during directory discovery, but `.nii`, `.nii.gz`, `.nrrd`, and `.nhdr` are not valid radiuma-suv PET inputs. For NumPy, pass an in-memory `ndarray` and DICOM metadata; a `.npy` *path* is not accepted by `radiuma-suv_from_input`.

## 📋 radiuma-suv Parameters Reference

### Public function parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `input_array` | `numpy.ndarray` | Required | PET array passed to `radiuma-suv`. |
| `input_metadata` | `dict` | Required | Corresponding PET DICOM metadata dictionary. |
| `image_input` | `numpy.ndarray`, string, path-like | Required | NumPy array, PET DICOM file, or PET DICOM directory. |
| `parameter_estimation` | `str` | `"OFF"` | `"ON"` or `"OFF"` only. |
| `output_path` | string / path-like / `None` | `None` | Destination for DICOM outputs; not valid with NumPy input. |
| `pet_metadata` | `dict` or `None` | `None` | `RadiumaSUVLoader` metadata for a NumPy input. |

### Selected configured DICOM fields

These are **example literal settings present in the notebook**, not a universal dataset or recommended patient parameters.

| Metadata field | Notebook setting | Meaning |
|---|---|---|
| `PatientWeight` | `"75.5"` | Patient weight, kg |
| `PatientAge` | `"045Y"` | Age in DICOM AS format |
| `PatientSize` | `"1.75"` | Height, m |
| `PatientSex` | `"M"` | DICOM sex code |
| `RadionuclideTotalDose` | `"370000000"` | Radionuclide activity; representation-specific unit handling must be verified |
| `RadionuclideHalfLife` | `"6586.2"` | Half-life, s |
| `RadiopharmaceuticalStartDateTime` | `"20261010093015"` | Radiopharmaceutical administration start |
| `ActualFrameDuration` | `"180000"` | Frame duration, ms |
| `FrameReferenceTime` | `"90000"` | Frame reference time, ms |
| `AcquisitionDate` | `"20261010"` | DICOM DA |
| `AcquisitionTime` | `"100000"` | DICOM TM |
| `Units` | `"BQML"` | PET units | BQML, GML, CM2ML, CNTS
| `DecayCorrection` | `"START"` | Decay reference |
| `RescaleSlope` | `"1.0"` | Pixel calibration slope |
| `RescaleIntercept` | `"0.0"` | Pixel calibration intercept |
| `Rows`, `Columns` | `256`, `256` | Example image size |
| `PixelSpacing` | `[4.0, 4.0]` | Row/column spacing, mm |
| `SliceThickness` | `"3.0"` | Slice thickness, mm |
| `ImagePositionPatient` | `[0.0, 0.0, 0.0]` | Image coordinates, mm |
| `ImageOrientationPatient` | `[1, 0, 0, 0, 1, 0]` | Orientation cosines |

The source defines many additional correction, geometry, timing, and vendor-related fields. Consult the notebook's **User Controlled Parameters** section for the complete set. These values should not be inferred or substituted without acquisition-specific evidence.

## 📊 Parameter Compatibility

| Feature / argument | NumPy + metadata | Single PET DICOM | DICOM directory |
|---|---:|---:|---:|
| `parameter_estimation` | ✓ | ✓ | ✓ |
| `output_path` | Not allowed | ✓ | ✓ |
| Floating-point array return | ✓ | No | No |
| Derived DICOM output | No | ✓ | ✓ |
| Independent series grouping | No | No | ✓ |


## 📚 API Reference

### `radiuma-suv(input_array, input_metadata, parameter_estimation="OFF")`

Converts an in-memory PET NumPy array to SUVbw. Returns `(suv_array, output_metadata)`.

### `radiuma-suv(input_array, input_metadata, parameter_estimation="ON")` 

Converts an in-memory PET NumPy array to SUVbw, it replace missed parameters in input image DICOM header with default vlues defined in code. Returns `(suv_array, output_metadata)`.

### `radiuma-suv_from_input(image_input, parameter_estimation="OFF", output_path=None)`

Detects input representation and returns a matching result: a NumPy array or a new DICOM output file/folder path, plus corresponding metadata. DICOM source objects are not deliberately overwritten.

### `radiuma-suv_from_input(image_input, parameter_estimation="ON", output_path=None)`
Detects input representation and returns a matching result: a new DICOM output file/folder path, it replace missed parameters in input image DICOM header with default vlues defined in code. DICOM source objects are not deliberately overwritten.

```python
from radiuma-suv import RadiumaSUVLoader
converter = RadiumaSUVLoader(parameter_estimation="OFF", output_path="pet_output.dcm")
output, metadata = converter.convert("pet_input.dcm")
```

## 📊 Output Structure

| Input | First return value | Second return value |
|---|---|---|
| NumPy array | Floating-point SUVbw NumPy array | Metadata dictionary |
| Single PET DICOM file | New DICOM file path | Output header dictionary |
| PET DICOM folder, one series | Output directory path | List of header dictionaries |
| Folder, multiple series | Output root directory path | Nested header lists |

**DICOM representation:** The source encodes SUVbw into calibrated **signed 16-bit pixels**. For conventional PET, the output `RescaleSlope` is needed to recover actual SUVbw. For Enhanced PET, use its real-world value mapping. The writer checks the saved file's pixel shape and reconstruction against calculated values within quantization tolerance. This is not an independent reference SUV validation.


## 📁 Supported File Formats

| Format | Supported | Notes |
|---|---|---|
| Conventional PET DICOM `.dcm` | Yes | Single frame |
| Enhanced PET DICOM `.dcm` | Yes | Multi-frame Enhanced PET SOP class |
| Extensionless DICOM in a folder | Detected if readable | Must identify as PET (`Modality=PT`) |
| PET DICOM directory | Yes | Recursive series grouping |
| `.npy` path | No | Load with `np.load()` and pass array |

## 🎯 Library Examples

### Example 1: Single DICOM

```python
from radiuma-suv import radiuma-suv_from_input
result, meta = radiuma-suv_from_input("pet.dcm", parameter_estimation="OFF")
print(result)
```
```python
from radiuma-suv import radiuma-suv_from_input
result, meta = radiuma-suv_from_input("pet.dcm", parameter_estimation="ON")
print(result)
```

### Example 2: Multiple PET folders

```python
from radiuma-suv import radiuma-suv_from_input
out, headers = radiuma-suv_from_input(
    "PET_Input", parameter_estimation="OFF", output_path="SUV_Output"
)
```
```python
from radiuma-suv import radiuma-suv_from_input
out, headers = radiuma-suv_from_input(
    "PET_Input", parameter_estimation="ON", output_path="SUV_Output"
)
```

### Example 3: NumPy data

```python
import numpy as np
from radiuma-suv import radiuma-suv
pixels = np.load("pet_pixels.npy")
# Load the genuine DICOM metadata for the same scan.
suv, meta = radiuma-suv(pixels, matching_pet_dicom_metadata, "OFF")
```
```python
import numpy as np
from radiuma-suv import radiuma-suv
pixels = np.load("pet_pixels.npy")
# Load the genuine DICOM metadata for the same scan.
suv, meta = radiuma-suv(pixels, matching_pet_dicom_metadata, "ON")
```
## ⚡ Performance Tips

1. Use the folder interface when processing multiple PET DICOM series.
2. Keep output outside the input directory and ensure the output path does not exist.
3. Avoid making unnecessary copies of large arrays when working with NumPy.
4. Use real acquisition-specific metadata for correct dose and decay calculations.
5. Validate numerical values after conversion and DICOM quantization.
6. The notebook contains no benchmarking results, so no speed claims are made here.



## ❓ Troubleshooting

| Message / problem | Recommended check |
|---|---|
| `ModuleNotFoundError: No module named 'radiuma-suv'` | Confirm the package was created, installed, and imported from the active Python environment. |
| `parameter_estimation must be "ON" or "OFF"` | Use the exact supported values. |
| `Input must be PET DICOM (Modality=PT)` | Confirm the file really contains PET DICOM. |
| `Non-Enhanced multi-frame PET is unsupported` | Input format is outside the implemented multi-frame path. |
| `No PET DICOM series found` | Verify DICOM readability, modality, and folder contents. |
| `Output directory must not be nested within input folder` | Choose an independent output directory. |
| `Refusing to overwrite existing output` | Choose another destination. |
| Missing activity, weight, half-life, or timing | Review source metadata; do not blindly apply configured defaults. |
| Unexpected SUVbw values | Recheck units, calibration, corrected dose, decay time, and actual acquisition data. |

### ❓ Get Help

- **Documentation:** Refer to this README for installation, configuration, supported formats, and usage examples.
- **Execution:** See the **Script Execution** section for standalone and embedded usage.
- **Troubleshooting:** Check the **Troubleshooting** section for common errors and solutions.
- **Issues and Bug Reports:** [GitHub Issues](https://github.com/radiuma-com/radiuma-suv/issues).
- **Technical Support:** Contact the development team using the details in the **Contact** or **Maintenance** section. 

## 🕒 Version History


```text
v1.0.1 
- Conventional PET DICOM to SUVbw
- Enhanced multi-frame PET to SUVbw
- NumPy input and metadata output
- Recursive DICOM series processing
```

## 📬 Contact

radiuma-suv is available **free of charge**.
For access, questions, or feedback:

**Dr. Mohammad R. Salmanpour (Team Lead)**  
📧[msalman@bccrc.ca](mailto:msalman@bccrc.ca) | [m.salmanpoor66@gmail.com](mailto:m.salmanpoor66@gmail.com), [m.salmanpour@ubc.ca](mailto:m.salmanpour@ubc.ca)

---
## 🛠️ Maintenance

For technical support and maintenance inquiries, please contact:

**Dr. Mohammad R. Salmanpour (Team Lead)**  
 msalman@bccrc.ca – m.salmanpoor66@gmail.com – m.salmanpour@ubc.ca

**Ramin Afrah (Maintainer)**  
 afrah.ramin@gmail.com
## 👥 Authors

- **Dr. Mohammad R. Salmanpour (Team Lead, Fund Provider, Evaluator, Medical Imaging Expert, Backend Development, Code Refactoring, Debugging, Library Management, IBSI Standardization, and Activation of the radiuma-suv library)** – [msalman@bccrc.ca](mailto:msalman@bccrc.ca), [m.salmanpoor66@gmail.com](mailto:m.salmanpoor66@gmail.com), [m.salmanpour@ubc.ca](mailto:m.salmanpour@ubc.ca)
- **Ramin Afrah (Development, Code Refactoring, Debugging, and Library Management)** - [afrah.ramin@gmail.com](mailto:afrah.ramin@gmail.com)
- **Dr. Mehrdad Oveisi (Evaluator, Software Engineer, and Advisor)** 
- **Dr. Arman Rahmim (Fund Provider, Medical Imaging Expert, Evaluator, and Advisor)** 

## 📚 Citation


```bibtex
Coming soon
```

## 📜 License

This open-source software is released under the **MIT License**, which grants permission to use, modify, and distribute it for any purpose, including research or commercial use, without requiring modified versions to be shared as open source. See the [LICENSE](LICENSE) file for details.

## Support


- **Issues**: [GitHub Issues](https://github.com/radiuma-com/radiuma-suv/issues)
- **Documentation**: This README and the included guides
- **Examples**: See `examples/basic_usage.py`


## Acknowledgment

- [🔬 **Qu**antitative **R**adiomolecular **I**maging and **T**herapy (Qurit) Lab, University of British Columbia, Vancouver, BC, Canada](https://www.qurit.ca)  
- [🏥 BC Cancer Research Institute, Department of Basic and Translational Research, Vancouver, BC, Canada](https://www.bccrc.ca/)  
- [💻 **Vir**tual **Collab**oration (VirCollab) Group, Vancouver, BC, Canada](https://vircollab.com/#/) 
- [🏭 **Tec**hnological **Vi**rtual **Co**llaboration **Corp**oration (TECVICO Corp.), Vancouver, BC, Canada](https://www.tecvico.com)  
