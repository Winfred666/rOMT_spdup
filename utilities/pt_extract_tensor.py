#!/usr/bin/env python3
"""Extract selected volume from a .pt file into raw binary.

Expected input structure:
  obj['samples'][0]['x'] with shape (2, D, H, W)

This utility is intended to be called from MATLAB (`pt2mat.m`) via `system(...)`
so it does not depend on MATLAB's `pyenv` CPython compatibility.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import torch


def parse_pair_from_filename(path: Path) -> tuple[int, int] | None:
    m = re.match(r"^(\d+)_(\d+)_samples$", path.stem)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def parse_index_list(txt: str, name: str) -> np.ndarray:
    parts = [p.strip() for p in txt.split(",") if p.strip()]
    if not parts:
        raise ValueError(f"{name} is empty")
    vals = np.asarray([int(p) for p in parts], dtype=np.int64)
    if np.any(vals < 1):
        raise ValueError(f"{name} must use 1-based positive indices")
    return vals - 1  # convert to 0-based for numpy


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract x tensor from a .pt sample file")
    parser.add_argument("--input", required=True, help="Path to .pt file")
    parser.add_argument("--out-bin", required=True, help="Output binary file path (float64)")
    parser.add_argument("--x1", required=True, help="Comma-separated 1-based indices for dimension D")
    parser.add_argument("--x2", required=True, help="Comma-separated 1-based indices for dimension H")
    parser.add_argument("--x3", required=True, help="Comma-separated 1-based indices for dimension W")
    parser.add_argument("--frame", type=int, default=None, help="Optional absolute frame id to select sample index")
    args = parser.parse_args()

    in_path = Path(args.input)
    if not in_path.is_file():
        raise FileNotFoundError(f"Input .pt file not found: {in_path}")

    obj = torch.load(str(in_path), map_location="cpu")
    if not isinstance(obj, dict):
        raise TypeError("Top-level .pt object is not a dict")

    samples = obj.get("samples", None)
    if not isinstance(samples, (list, tuple)) or len(samples) == 0:
        raise KeyError("Missing list-like key: samples with at least one element")

    s0 = samples[0]
    if not isinstance(s0, dict):
        raise TypeError("samples[0] must be a dict")

    x = s0.get("x", None)
    if x is None:
        raise KeyError("Missing key: samples[0]['x']")

    if hasattr(x, "detach"):
        x = x.detach().cpu().numpy()
    else:
        x = np.asarray(x)

    x = np.asarray(x)
    if x.ndim != 4:
        raise ValueError(f"Expected x with 4 dims (2,D,H,W), got shape {x.shape}")
    if x.shape[0] != 2:
        raise ValueError(f"Expected first dim size 2 for frame pair, got {x.shape[0]}")

    sample_idx = 0
    frame_pair = parse_pair_from_filename(in_path)
    if args.frame is not None:
        if frame_pair is None:
            raise ValueError(
                f"Frame was provided ({args.frame}) but filename does not encode frame pair as <f1>_<f2>_samples.pt: {in_path.name}"
            )
        if args.frame == frame_pair[0]:
            sample_idx = 0
        elif args.frame == frame_pair[1]:
            sample_idx = 1
        else:
            raise ValueError(
                f"Requested frame {args.frame} does not match file frame pair {frame_pair} in {in_path}"
            )

    x_sel = x[sample_idx]
    if x_sel.ndim != 3:
        raise ValueError(f"Selected sample must be 3D (D,H,W), got shape {x_sel.shape}")

    i1 = parse_index_list(args.x1, "x1")
    i2 = parse_index_list(args.x2, "x2")
    i3 = parse_index_list(args.x3, "x3")

    d, h, w = x_sel.shape
    if i1.max(initial=-1) >= d or i2.max(initial=-1) >= h or i3.max(initial=-1) >= w:
        raise ValueError(
            f"Requested indices out of bounds: max requested [D,H,W]=[{i1.max()+1},{i2.max()+1},{i3.max()+1}], "
            f"available=[{d},{h},{w}]"
        )

    # Preserve index order exactly as provided.
    vol = x_sel[np.ix_(i1, i2, i3)]
    x64 = np.asarray(vol, dtype=np.float64)

    out_bin = Path(args.out_bin)
    out_bin.parent.mkdir(parents=True, exist_ok=True)

    x64.tofile(out_bin)


if __name__ == "__main__":
    main()
