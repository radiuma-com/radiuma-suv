"""Basic demonstration (requires your own de-identified PET DICOM input)."""
from pathlib import Path
from radiuma_suv import radiuma_suv_from_input

PET_INPUT = Path("/path/to/your/PET_input.dcm")

if __name__ == "__main__":
    if not PET_INPUT.exists():
        raise SystemExit("Edit PET_INPUT to point to your own verified PET DICOM input.")
    output, _ = radiuma_suv_from_input(PET_INPUT, parameter_estimation="OFF")
    print("Derived SUVbw DICOM:", output)
