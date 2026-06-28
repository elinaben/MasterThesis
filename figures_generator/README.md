# Figures generator

Standalone scripts that build composite figures for the thesis from raw image
data, so the figures can be regenerated reproducibly.

## Scripts

### `make_jnk_figure.py`
Builds `figs/June26_JNK_figure.png` (JNK activation: Ctrl vs proNGF).

- **Panel A** — p-JNK and Total JNK Western blots, cropped at full resolution
  from `ImagesTif/ELINA 2026-06-26 10h55m32s(IRDye 800CW).tif`. That scan holds
  two membranes: p-JNK = left (smaller), Total JNK = right (larger); left lane =
  Ctrl, right lane = proNGF. The `54/46 kDa` markers label the JNK p54/p46
  isoforms. Crop boxes are set in `BOXES` near the top of the script.
- **Panel B** — p-JNK / Total JNK ratio bars (values in `RATIOS`, from the
  quantification in `ImagesCSV/June26_*_Results.csv`).

Run from the repo root:

```bash
python3 figures_generator/make_jnk_figure.py
```

Requires `numpy`, `pillow`, and `matplotlib`.
