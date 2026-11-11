"""
analyze_nuclei.py - find nuclei in a SmartMic PoC overview and return them as
stage positions for the detailed pass.

    pixi run python analyze_nuclei.py <overview.czi> <output_dir> --prefix TAG
        [--channel DAPI] [--max-nuclei 1000] [--min-area-um2 25]

SmartMic (MS_image_analysis.run_analysis_json, suffix "targets") runs this and
reads <output_dir>/<TAG>_targets.json if the exit code is 0. The fixed part
(arguments, log, result file, exit code, provenance) is
ia_common.cli.run_analysis_script; the detection is nuclei.find_nuclei.

Writes to <output_dir>:
    <TAG>_targets.json   found, n_nuclei, truncated, targets[] (+ provenance)
    <TAG>_targets.log    run log (also on stdout -> SmartMic's run log)
    <TAG>_nuclei.png     QC: the channel, each kept nucleus marked with its id

Each target: id, x_m / y_m (stage position of the centroid, metres, from
ia_common's pixel frame), area_m2, row / col (centroid, image pixels).
"""

import numpy as np

from ia_common.cli import run_analysis_script
from ia_common.czi import read_single
from ia_common.display import save_png

from nuclei import find_nuclei

SUFFIX = "targets"


def add_args(parser):
    parser.add_argument("--channel", default="DAPI", help="channel name (default DAPI)")
    parser.add_argument("--max-nuclei", type=int, default=1000,
                        help="keep at most this many nuclei (label order)")
    parser.add_argument("--min-area-um2", type=float, default=25.0, help="smallest nucleus kept")
    parser.add_argument("--min-peak-distance-um", type=float, default=5.0,
                        help="closest two watershed seeds may be")


def analyse(run):
    a = run.args
    im = read_single(run.image)
    if a.channel not in im.channel_names:
        raise ValueError(f"channel {a.channel!r} not in {im.channel_names}")
    img = im.data[im.channel_names.index(a.channel)].max(axis=0)          # (Y, X), max over Z
    cx, cy = im.stage_xy_m
    run.log.info(f"{im.data.shape} (C, Z, Y, X), pixel {im.pixel_size_m * 1e6:.4f} um, "
                 f"centre ({cx * 1e6:.1f}, {cy * 1e6:.1f}) um, placed by {im.frame_source}")

    nuclei, _ = find_nuclei(img, im.pixel_size_m,
                            min_peak_distance_m=a.min_peak_distance_um * 1e-6,
                            min_area_m2=a.min_area_um2 * 1e-12)
    kept = nuclei[:a.max_nuclei]
    truncated = len(nuclei) > len(kept)
    run.log.info(f"{len(nuclei)} nuclei (border-touching and < {a.min_area_um2:g} um2 dropped)"
                 + (f"; kept the first {len(kept)}" if truncated else ""))

    targets = []
    for n in kept:
        x, y = im.frame.to_stage_m(n.row, n.col)
        targets.append({"id": n.id, "x_m": float(x), "y_m": float(y), "area_m2": n.area_m2,
                        "row": n.row, "col": n.col})

    def overlay(ax):
        for n in kept:
            ax.plot(n.col, n.row, "+", color="orange", ms=8)
            ax.text(n.col + 6, n.row, str(n.id), color="orange", fontsize=6, va="center")
    run.qc(save_png, run.output_path("nuclei.png"), np.asarray(img), im.pixel_size_m,
           overlay=overlay, title=f"{len(kept)} nuclei")

    return {"found": bool(kept), "n_nuclei": len(nuclei), "truncated": truncated,
            "channel": a.channel, "image_centre_m": [cx, cy], "targets": targets}


if __name__ == "__main__":
    raise SystemExit(run_analysis_script(analyse, SUFFIX, add_args=add_args))
