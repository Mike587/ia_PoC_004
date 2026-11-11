"""find_nuclei on synthetic images, and the whole script on a real PoC overview
(skipped when the ia_PoC_002 example CZIs are not on this machine)."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from nuclei import find_nuclei  # noqa: E402

PX = 0.35e-6                     # about the PoC overview pixel (20x/0.95 + 2x camera)
EXAMPLE = Path(r"C:\Users\zeiss\Mike\Image_Analysis\ia_PoC_002\experiment_pos_000077.czi")


def disc_image(centres_px, radius_px=12, shape=(300, 400), noise=2.0, seed=0):
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    img = 100 + rng.normal(0, noise, shape)
    for r, c in centres_px:
        img[(yy - r) ** 2 + (xx - c) ** 2 <= radius_px ** 2] += 1000
    return img.astype(np.uint16)


def test_finds_separate_discs_at_their_centres():
    centres = [(80, 100), (150, 250), (220, 120)]
    nuclei, seg = find_nuclei(disc_image(centres), PX)
    assert len(nuclei) == 3
    found = sorted((round(n.row), round(n.col)) for n in nuclei)
    assert all(abs(fr - r) <= 1 and abs(fc - c) <= 1 for (fr, fc), (r, c) in zip(found, sorted(centres)))
    area = np.pi * 12 ** 2 * PX ** 2
    assert all(n.area_m2 == pytest.approx(area, rel=0.15) for n in nuclei)
    assert seg.shape == (300, 400)


def test_border_touching_and_small_objects_are_dropped():
    img = disc_image([(150, 200), (5, 200)])            # second one touches the top edge
    img[60:63, 60:63] += 1000                           # a 1 um speck
    nuclei, _ = find_nuclei(img, PX)
    assert [(round(n.row), round(n.col)) for n in nuclei] == [(150, 200)]


def test_touching_pair_is_split():
    nuclei, _ = find_nuclei(disc_image([(150, 180), (150, 200)]), PX)   # discs overlap by 4 px
    assert len(nuclei) == 2


def test_flat_image_has_no_nuclei():
    nuclei, seg = find_nuclei(np.full((100, 100), 500, dtype=np.uint16), PX)
    assert nuclei == [] and not seg.any()


@pytest.mark.skipif(not EXAMPLE.exists(), reason="ia_PoC_002 example CZI not on this machine")
def test_script_on_a_real_overview(tmp_path):
    p = subprocess.run([sys.executable, str(ROOT / "analyze_nuclei.py"), str(EXAMPLE), str(tmp_path),
                        "--prefix", "T"], capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr
    res = json.loads((tmp_path / "T_targets.json").read_text())
    assert res["found"] and res["n_nuclei"] == len(res["targets"]) > 10
    cx, cy = res["image_centre_m"]
    for t in res["targets"]:                            # every target inside the field of view
        assert abs(t["x_m"] - cx) < 4096 * PX / 2 and abs(t["y_m"] - cy) < 3008 * PX / 2
    assert "provenance" in res
    assert (tmp_path / "T_nuclei.png").exists()
