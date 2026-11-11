"""Limited API behavior regression tests; not quantitative PET validation."""
import numpy as np
import pytest
from radiuma_suv import RadiumaSUVLoader, radiuma_suv, radiuma_suv_from_input
from radiuma_suv.core import _suv_int_pixels

def _gml_header():
    return {"Units": "GML", "SUVType": "BW", "PatientWeight": "70",
            "RescaleSlope": "1", "RescaleIntercept": "0", "Rows": 2, "Columns": 2}

def test_gml_bw_identity_in_memory():
    pixels = np.array([[0, 1], [2, 3]], dtype=np.int16)
    result, meta = radiuma_suv(pixels, _gml_header(), "OFF")
    np.testing.assert_array_equal(result, pixels.astype(float))
    assert meta == _gml_header()

def test_gml_rescale_applied():
    header = _gml_header()
    header["RescaleSlope"] = "2"
    result, _ = radiuma_suv(np.ones((2, 2)), header, "OFF")
    np.testing.assert_array_equal(result, np.full((2, 2), 2.0))

def test_invalid_parameter_estimation_rejected():
    with pytest.raises(ValueError, match="parameter_estimation"):
        radiuma_suv(np.ones((2,2)), _gml_header(), "INVALID")

def test_numpy_input_needs_metadata():
    with pytest.raises(ValueError, match="metadata_json"):
        radiuma_suv_from_input(np.ones((2,2)))

def test_loader_matches_function():
    loader = RadiumaSUVLoader(parameter_estimation="OFF", pet_metadata=_gml_header())
    actual, _ = loader.convert(np.ones((2,2)))
    np.testing.assert_array_equal(actual, np.ones((2,2)))

def test_int_pixel_encoding_within_half_step():
    sample = np.array([[0.0, 1.0], [2.7, 10.0]])
    stored, slope = _suv_int_pixels(sample)
    assert stored.dtype == np.dtype("<i2")
    assert np.all(np.abs(stored.astype(float) * slope - sample) <= slope * .50001)
