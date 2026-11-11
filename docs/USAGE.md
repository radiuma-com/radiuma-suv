# Usage notes

## Input assumptions

- PET DICOM modality is `PT`.
- Conventional multi-slice input requires consistent series metadata and verifiable geometry.
- Enhanced PET processing requires per-frame geometry and suitable quantitative mappings.
- For a NumPy input, use `radiuma_suv(pixels, matching_header_dict)`; 3D conventional arrays require `input_metadata={"slices": [header_per_slice, ...]}`.
- Multiple PET series in a directory are processed separately.

## Command line

```bash
radiuma-suv INPUT_PATH [--output OUTPUT_PATH] [--parameter-estimation {OFF,ON}]
```

Output is never allowed to overwrite the source. Supplying `--output` may still result in an `_processed` suffix.

## Safety

Test conversion and quantitative accuracy on validated non-identifiable reference PET data before research or clinical interpretation. `ON` uses example default values, **not** measured values; verify all dose, weight, temporal and unit metadata.

## Exact public API naming

Install the **`radiuma-suv`** distribution and invoke its **`radiuma-suv`** CLI. Python code imports the **`radiuma_suv`** module, which exports **`radiuma_suv`**, **`radiuma_suv_from_input`**, and **`RadiumaSUVLoader`**. the array converter uses `input_metadata`, and the loader constructor uses `pet_metadata`; these distinct parameter names are retained for backward compatibility.
