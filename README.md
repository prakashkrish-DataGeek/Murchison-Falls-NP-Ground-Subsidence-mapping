# 🛰️ Sentinel-1 SAR Subsidence Monitor — Murchison Falls National Park

> 4-year Sentinel-1 C-band SAR backscatter anomaly time-series mapping
> land surface change and subsidence proxy signals across Murchison Falls
> National Park, Uganda — 48 monthly frames, May 2022 → April 2026.

[![GEE](https://img.shields.io/badge/Google%20Earth%20Engine-4285F4?style=flat&logo=google&logoColor=white)](https://earthengine.google.com)
[![Python](https://img.shields.io/badge/Python-3.11-blue?style=flat&logo=python)](https://python.org)
[![Sentinel-1](https://img.shields.io/badge/Sentinel--1-SAR-blueviolet?style=flat)](https://sentinel.esa.int/web/sentinel/missions/sentinel-1)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Live Map](https://img.shields.io/badge/Live%20Map-View%20Online-orange)](https://prakashkrish-DataGeek.github.io/Murchison-Falls-NP-Ground-Subsidence-mapping/)

---

## 🗺 Live Project Page

**[→ View Interactive Time-Slider Map](https://prakashkrish-DataGeek.github.io/Murchison-Falls-NP-Ground-Subsidence-mapping/)**

---

## Overview

Murchison Falls National Park (~**3,840 km²**) sits on Pleistocene lacustrine
sediments at the edge of the Albertine Rift — one of East Africa's most
tectonically active zones. The Albert Nile delta, wetland margins, and areas
adjacent to oil exploration corridors are at elevated risk of **sediment
compaction, hydraulic subsidence and surface deformation**.

This project builds an **open-source, fully reproducible** 4-year subsidence
monitoring pipeline using free Sentinel-1 SAR data processed entirely in
Google Earth Engine and Python, with a self-contained interactive HTML output.

### Key Outputs

| Output | Description |
|--------|-------------|
| Monthly SAR Anomaly Grid | Per-pixel seasonal-baseline-removed VV backscatter anomaly (1 km, 48 months) |
| Z-score Intensity Map | Normalised \|Z\| surface-change signal → [0, 1] per pixel per month |
| HeatMapWithTime | Self-contained Folium time-slider HTML — 48 monthly frames |
| Subsidence Proxy Raster | Full-period temporal std-dev + Z-score anomaly GeoTIFF (20m) |
| Hotspot Mask | Binary raster of pixels with temporal StdDev > 2.5 dB |

---

## Anomaly Method

```
Baseline VV(pixel, cal_month) = mean VV for the same calendar month
                                  across the first 12 months of record

Anomaly(t)   = VV(t) − Baseline(pixel, cal_month(t))

Z-score(t)   = Anomaly(t) / σ_pixel

Intensity(t) = clip(|Z-score(t)|, 0, 3) / 3     → [0, 1]
```

Removing the per-calendar-month seasonal baseline filters Uganda's double-peaked
wet seasons (March–May, October–November) so the residual signal reflects genuine
multi-year surface change rather than vegetation moisture cycles.

---

## Repository Structure

```
Murchison-Falls-NP-Ground-Subsidence-mapping/
├── gee_scripts/
│   └── mfnp_s1_subsidence_export.js    ← S1 monthly composites + subsidence proxy
├── python/
│   ├── mfnp_01_process_sar.py          ← Anomaly, Z-score, intensity normalisation
│   └── mfnp_02_timeslider_map.py       ← Folium HeatMapWithTime builder
├── web/
│   ├── index.html                      ← GitHub Pages project showcase
│   └── mfnp_subsidence_timeslider.html ← Generated time-slider map (commit this)
├── data/
│   └── .gitkeep                        ← GEE CSV exports go here (gitignored)
├── requirements.txt
└── README.md
```

---

## Quick Start

### Step 1 — GEE Export (run in [code.earthengine.google.com](https://code.earthengine.google.com))

1. Open `gee_scripts/mfnp_s1_subsidence_export.js`
2. Copy-paste into a new GEE Script and click **Run**
3. In the **Tasks** tab click **Run** next to `MFNP_S1_Monthly_VV_2022_2026`
4. Download the exported CSV from Google Drive → `data/`

### Step 2 — Python Pipeline

```bash
# 1. Clone repo
git clone https://github.com/prakashkrish-DataGeek/Murchison-Falls-NP-Ground-Subsidence-mapping.git
cd Murchison-Falls-NP-Ground-Subsidence-mapping

# 2. Install dependencies
pip3 install -r requirements.txt

# 3. Process real SAR data
python3 python/mfnp_01_process_sar.py --input data/MFNP_S1_Monthly_VV_2022_2026.csv

# OR run demo mode (no GEE required — synthetic but realistic data)
python3 python/mfnp_01_process_sar.py --demo

# 4. Build the time-slider map
python3 python/mfnp_02_timeslider_map.py
# → Generates: mfnp_subsidence_timeslider.html
```

### Step 3 — Publish

```bash
git add mfnp_subsidence_timeslider.html
git commit -m "Update time-slider with latest Sentinel-1 data"
git push
```

GitHub Pages redeploys automatically within ~2 minutes.

---

## Data Sources

| Dataset | Source | Resolution | GEE Collection ID |
|---------|--------|------------|-------------------|
| Sentinel-1 GRD | ESA Copernicus | 20m | `COPERNICUS/S1_GRD` |
| NASADEM Elevation | NASA JPL | 30m | `NASA/NASADEM_HGT/001` |
| Basemap tiles | CartoDB / Esri | — | Via Folium / Leaflet |

**Sentinel-1 filter settings:** IW mode · VV polarisation · Ascending pass · Monthly median composite · 1 km sample grid

---

## True InSAR Note

GEE does not support phase-based InSAR processing. For millimetre-precision
vertical displacement over MFNP, use:

- **ASF HyP3 SBAS** (cloud-based, free for research): https://hyp3-docs.asf.alaska.edu/
- **MintPy** time-series InSAR: https://mintpy.readthedocs.io/
- **ESA SNAP + StaMPS** with Sentinel-1 SLC data from [ASF Data Search](https://search.asf.alaska.edu/)
- **Copernicus EGMS** covers Europe only — not applicable to Uganda

This project uses SAR backscatter temporal variance as a **free, globally
available proxy** suitable for hotspot screening and 4-year trend monitoring.

---

## Results Summary

- **Monitoring period:** May 2022 → April 2026 (48 monthly composites)
- **Grid coverage:** ~6,400 analysis points at 1 km spacing across MFNP
- **Mean VV backscatter:** −13 to −9 dB (tropical savanna / riverine forest)
- **Highest-change zones:** Albert Nile delta, Pakuba riverbank, Buligi peninsula
- **Year-on-year drift:** −0.018 dB (2023) → −0.048 dB (2024) → −0.076 dB (2025) → −0.090 dB (2026)
- **High-anomaly pixels (|Z| > 1.5σ):** ~16% of grid — concentrated in wetland and Nile corridor zones

---

## Author

**Prakash Krishnamachari**
- GitHub: [@prakashkrish-DataGeek](https://github.com/prakashkrish-DataGeek)
- LinkedIn: [linkedin.com/in/prakashkrishnamachari](https://linkedin.com/in/prakashkrishnamachari)
- Email: prakash.krishnamachari@gmail.com

---

## License

MIT License — see [LICENSE](LICENSE) for details.
Satellite imagery data is subject to ESA Copernicus and respective agency terms of use.

---

*Built with ❤️ for open Earth observation and conservation science*
