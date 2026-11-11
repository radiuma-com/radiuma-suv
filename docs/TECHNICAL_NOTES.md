# Technical notes and limitations

This repository packages the user-provided conversion engine. Changes to the conversion calculations were intentionally avoided.

## Packaging changes

- Removed notebook-only `%pip` command. Runtime requirements are in `pyproject.toml`.
- Removed top-level debugging `print` calls that evaluated undefined local variables during import.
- Replaced machine-specific demo input with a `radiuma-suv` command-line entry point.
- Added a minimal automated test suite and GitHub Actions workflow.

## Validation status

**No independent numerical reference validation is claimed.** Unit tests establish limited software behavior only, not quantitative or clinical correctness. Investigate any unexpectedly large, negative, or non-finite SUVbw. DICOM outputs quantize computed values to signed 16-bit pixels.

## Scientific concerns requiring review

- The original source heuristically treats some smaller injected dose values as MBq and high patient-weight values as grams. Verify modality-specific units against original acquisition data.
- The optional metadata estimation settings contain a fixed example patient's weight, dates, administered activity, and other parameters. Never use those placeholders as real patient data.
- DICOM transfer syntaxes and vendor-private scaling require validation with representative inputs.
- Enhanced PET RWVM and functional-group behavior need interoperability testing with real-world images.
- DICOM identity and derivation metadata should be checked against the intended institution's imaging and research data handling requirements.
- No PET data are bundled because they may contain protected health information.

## Naming contract and historical files

- The installed distribution and executable command are `radiuma-suv`.
- Python package imports use `radiuma_suv` because hyphens are not valid Python identifiers.
- The exported API names are exactly `radiuma_suv`, `radiuma_suv_from_input`, and `RadiumaSUVLoader`.
- `radiuma_suv` takes `input_metadata`; `radiuma_suv_from_input` and `RadiumaSUVLoader` takes `pet_metadata`. These are deliberately kept as implemented to avoid breaking callers.
- `ORIGINAL_README.md` and `ORIGINAL_SOURCE.py.txt` are historical references with naming comments/examples normalized for this release; they are not verbatim copies and are not executed. Use the current top-level README and `src/radiuma_suv/core.py` for current behavior.
