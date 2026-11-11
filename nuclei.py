"""Nucleus detection on one 2D image (pure numpy / scikit-image, no file I/O).

Otsu threshold on a Gaussian-smoothed image, a small opening, then a
distance-transform watershed to split touching nuclei. Nuclei that touch the
image border or are smaller than ``min_area_m2`` are dropped. Same method and
defaults as ia_PoC_002, so results are comparable.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import label, regionprops
from skimage.morphology import disk, opening
from skimage.segmentation import clear_border, watershed


@dataclass(frozen=True)
class Nucleus:
    id: int          # label in the segmentation
    row: float       # centroid, image pixels
    col: float
    area_m2: float


def find_nuclei(img: np.ndarray, pixel_size_m: float, *,
                smooth_px: float = 2.0, open_px: int = 3,
                min_peak_distance_m: float = 5e-6,
                min_area_m2: float = 25e-12) -> tuple[list[Nucleus], np.ndarray]:
    """Detect nuclei in a 2D image.

    Args:
        img: (Y, X) image, any numeric dtype.
        pixel_size_m: square pixel size of ``img`` (m).
        smooth_px: Gaussian sigma before the threshold (pixels).
        open_px: radius of the opening that removes small specks (pixels).
        min_peak_distance_m: closest two watershed seeds may be (m).
        min_area_m2: smallest nucleus kept (m²).

    Returns:
        (nuclei in label order, label image). An image with no foreground
        (e.g. a flat field) gives no nuclei, not an error.
    """
    arr = np.asarray(img, dtype=np.float32)
    smoothed = gaussian(arr, sigma=smooth_px)
    if float(smoothed.max()) == float(smoothed.min()):
        return [], np.zeros(arr.shape, dtype=np.int32)
    binary = opening(smoothed > threshold_otsu(smoothed), disk(open_px))

    distance = ndi.distance_transform_edt(binary)
    min_dist_px = max(1, round(min_peak_distance_m / pixel_size_m))
    peaks = peak_local_max(distance, min_distance=min_dist_px, labels=binary)
    seeds = np.zeros(distance.shape, dtype=bool)
    seeds[tuple(peaks.T)] = True
    segmented = clear_border(watershed(-distance, label(seeds), mask=binary))

    px_area = pixel_size_m * pixel_size_m
    nuclei = [Nucleus(int(r.label), float(r.centroid[0]), float(r.centroid[1]), r.area * px_area)
              for r in regionprops(segmented) if r.area * px_area >= min_area_m2]
    return nuclei, segmented
