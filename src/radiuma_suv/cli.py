"""CLI for converting a DICOM PET file or directory to SUVbw."""
import argparse
from .core import radiuma_suv_from_input

def main(argv=None):
    parser = argparse.ArgumentParser(description="Convert PET DICOM to derived SUVbw DICOM.")
    parser.add_argument("input", help="PET DICOM file or directory")
    parser.add_argument("--output", default=None, help="Output file or directory (must not exist)")
    parser.add_argument("--parameter-estimation", choices=("ON", "OFF"), default="OFF",
                        help="OFF (recommended) avoids filling missing patient metadata")
    args = parser.parse_args(argv)
    result, _ = radiuma_suv_from_input(args.input, parameter_estimation=args.parameter_estimation,
                                       output_path=args.output)
    print(result)
    return 0
