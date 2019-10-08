import helpers
import numpy as np

def test_divide_into_equal_passes_aligned():
    assert np.allclose(helpers.divide_into_equal_passes(0, 10, 2),
                       [0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
                       
def test_divide_into_equal_passes_misaligned():
    assert np.allclose(helpers.divide_into_equal_passes(0, 10, 2.2),
                       [0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
                       
def test_divide_into_equal_passes_non_zero_start():
    assert np.allclose(helpers.divide_into_equal_passes(2, 12, 2.2),
                       [2.0, 4.0, 6.0, 8.0, 10.0, 12.0])
                       
def test_divide_into_equal_passes_backwards():
    assert np.allclose(helpers.divide_into_equal_passes(12, 2, 2.2),
                       [12.0, 10.0, 8.0, 6.0, 4.0, 2.0])

def test_equally_spaced_points_no_edge():
    assert np.allclose(helpers.center_equally_spaced_points_in_range(0, 10, 2),
                       [0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
                       
def test_equally_spaced_points_non_zero_start():
    assert np.allclose(helpers.center_equally_spaced_points_in_range(2, 10, 2),
                       [2.0, 4.0, 6.0, 8.0, 10.0])
    
def test_equally_spaced_points_half_edge():
    assert np.allclose(helpers.center_equally_spaced_points_in_range(0, 10, 2, 0.5),
                       [1.0, 3.0, 5.0, 7.0, 9.0])
    
def test_equally_spaced_points_multi_edge():
    assert np.allclose(helpers.center_equally_spaced_points_in_range(0, 10, 2, 3.5), 
                       [4.0, 6.0])
    
def test_equally_spaced_points_backwards():
    assert np.allclose(helpers.center_equally_spaced_points_in_range(15, 5, 3, 1.0),
                       [13., 10.,  7.])
    
def test_equally_spaced_points_negative():
    assert np.allclose(helpers.center_equally_spaced_points_in_range(-15, -5, 3, 1.0),
                       [-13., -10.,  -7.])
                       
def test_divide_with_clearout_depths_no_clearout_when_max_equals():
    assert np.allclose(helpers.divide_with_clearout_depths(0.5, 0.25, 0.5),
                       [0.25, 0.5])
                       
def test_divide_with_clearout_depths_single_clearout():
    assert np.allclose(helpers.divide_with_clearout_depths(0.75, 0.125, 0.5),
                       [0.125, 0.25 , 0.375, 0.5,
                        0.25, 0.5, 
                        0.625, 0.75])
                       
def test_divide_with_clearout_depths_multiple_clearouts():
    assert (helpers.divide_with_clearout_depths(2.0, 0.125, 0.5) ==
                       [0.125, 0.25, 0.375, 0.5, 
                        0.25, 0.5, # first clearout
                        0.625, 0.75, 0.875, 1.0,
                        0.5, 1.0,  # second clearout
                        1.125, 1.25, 1.375, 1.5,
                        0.75, 1.5, # third clearout
                        1.625, 1.75, 1.875, 2.0])

    