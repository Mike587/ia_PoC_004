# ia_PoC_004

Image analysis for the SmartMic proof of concept
(`SmartMic/projects/smartmic_poc`): find nuclei in a widefield DAPI overview
and return them as stage positions for the detailed (confocal) pass.

Built on [ia_common](https://github.com/Mike587/ia_common) (CZI reading,
pixel ↔ stage, logging, the SmartMic script contract), pinned to a commit in
`pixi.toml`. Replaces ia_PoC_002 (bioio-based, own CZI reading).

## Run

```
pixi run python analyze_nuclei.py <overview.czi> <output_dir> --prefix TAG [--channel DAPI]
pixi run test
```

SmartMic runs it through `MS_image_analysis.run_analysis_json(..., suffix="targets")`.

## Result: `<TAG>_targets.json`

```json
{"found": true, "n_nuclei": 421, "truncated": false, "channel": "DAPI",
 "image_centre_m": [0.09225, 0.06475],
 "targets": [{"id": 3, "x_m": 0.0921, "y_m": 0.0646, "area_m2": 6.1e-11, "row": 812.4, "col": 230.9}],
 "provenance": {"script_commit": "...", "ia_common_commit": "...", "options": {}}}
```

- `x_m` / `y_m`: stage position of the nucleus centroid in metres, from
  ia_common's pixel frame (`SingleImage.frame.to_stage_m`).
- Nothing found is a valid result (`found: false`, exit 0). A broken analysis
  exits 1 and writes no result file.
- QC image: `<TAG>_nuclei.png` (kept nuclei marked with their id, scale bar).

## Method (`nuclei.py`)

Gaussian smoothing (σ 2 px), Otsu threshold, opening (radius 3 px),
distance-transform watershed (seeds ≥ 5 µm apart), border-touching nuclei and
nuclei < 25 µm² dropped. The same method and defaults as ia_PoC_002.

## Checked against ia_PoC_002 (2026-11-11)

On the two ia_PoC_002 example overviews (20x/0.95 + 2x camera, 0.347 µm pixel):

| File | ia_PoC_002 | ia_PoC_004 | Position difference |
|---|---|---|---|
| `Experiment-2124.czi` | 421 nuclei | 421, same labels and areas | +0.17 µm in X and Y, constant |
| `experiment_pos_000077.czi` | 548 | 548, same labels and areas | +0.26 / +0.12 µm, constant |

The constant 0.17 µm is half a pixel: ia_common places a pixel at its centre,
ia_PoC_002 used `col - width/2`. The extra 0.09 µm on the second file comes
from the pixel frame vs ia_PoC_002's planned scene centre. Both are below the
stage repeatability.
