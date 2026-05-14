"""
mfnp_01_process_sar.py
======================
Process Sentinel-1 SAR monthly VV backscatter data for
Murchison Falls National Park (MFNP) subsidence monitoring.

Workflow:
  1. Load GEE-exported monthly VV CSV  (or generate synthetic demo data)
  2. Filter NoData / fill values
  3. Compute per-pixel seasonal baseline (first 12 months, calendar-month matched)
  4. Compute anomaly dB = VV(t) − baseline_same_calendar_month(pixel)
  5. Compute per-pixel Z-score anomaly
  6. Normalise |Z-score| → intensity [0, 1] (capped at 3σ) for HeatMapWithTime
  7. Export processed CSV → data/mfnp_sar_processed.csv

Input (GEE export):
    data/MFNP_S1_Monthly_VV_2022_2026.csv
    Columns: date (YYYY-MM), longitude, latitude, VV (dB), n_images

Output:
    data/mfnp_sar_processed.csv
    Columns: longitude, latitude, date, vv_db, anomaly_db,
             anomaly_zscore, intensity_norm

Usage:
    # Real GEE data:
    python python/mfnp_01_process_sar.py --input data/MFNP_S1_Monthly_VV_2022_2026.csv

    # Synthetic demo (generates realistic MFNP-like data):
    python python/mfnp_01_process_sar.py --demo

    # With custom output path:
    python python/mfnp_01_process_sar.py --demo --output data/mfnp_sar_processed.csv
"""

import os
import json
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
warnings.filterwarnings('ignore')

# ── Configuration ──────────────────────────────────────────────────────────────
DATA_DIR  = 'data'
OUT_CSV   = os.path.join(DATA_DIR, 'mfnp_sar_processed.csv')

START_DATE      = '2022-05'     # Analysis window start
END_DATE        = '2026-04'     # Analysis window end (inclusive)
BASELINE_MONTHS = 12            # First N months used as seasonal baseline
NODATA_VALUE    = -999.0        # GEE fill value when no S1 images in a month
VV_VALID_RANGE  = (-35.0, 5.0) # Plausible S1 VV range in dB (land surfaces)
GRID_SPACING    = 0.01          # ~1 km at equatorial Uganda (degrees)

# MFNP bounding box [lon_min, lat_min, lon_max, lat_max]
MFNP_BBOX = {
    'lon_min': 31.40, 'lon_max': 32.20,
    'lat_min':  1.75, 'lat_max':  2.55
}

# ── Date utilities ─────────────────────────────────────────────────────────────

def generate_date_list(start: str = START_DATE,
                       end:   str = END_DATE) -> list:
    """Return sorted list of 'YYYY-MM' strings from start to end inclusive."""
    dates   = []
    cur_yr, cur_mo = map(int, start.split('-'))
    end_yr, end_mo = map(int, end.split('-'))
    while (cur_yr, cur_mo) <= (end_yr, end_mo):
        dates.append(f'{cur_yr:04d}-{cur_mo:02d}')
        cur_mo += 1
        if cur_mo > 12:
            cur_mo = 1
            cur_yr += 1
    return dates


# ── Demo data generator ────────────────────────────────────────────────────────

def generate_demo_data(output_path: str) -> pd.DataFrame:
    """
    Generate synthetic but physically realistic Sentinel-1 VV data for MFNP.

    Spatial structure:
      • Background: land-type gradient (savanna ≈ -13 dB, forest ≈ -10 dB,
        water/wetland ≈ -18 to -20 dB, near Lake Albert)
      • Albert Nile corridor: higher seasonal variability (flooding)
    Temporal structure:
      • Uganda double-peaked rainy seasons: long rains (Mar–May),
        short rains (Oct–Nov); dry seasons (Dec–Feb, Jun–Sep)
      • Seasonal VV amplitude: ±1.5–2.5 dB (vegetation moisture)
      • Per-pixel white noise: σ ≈ 0.4 dB
    Subsidence hotspots (5 zones with linear negative trend):
      • Pakuba area (river bank): −0.04 dB/month
      • Southern sector:          −0.06 dB/month
      • Albert Nile bank:         −0.03 dB/month
      • Northern boundary:        −0.05 dB/month
      • Buligi peninsula:         −0.04 dB/month
    """
    np.random.seed(42)
    dates   = generate_date_list()
    n_months = len(dates)

    # ── Grid ──────────────────────────────────────────────────────────────────
    lons = np.arange(MFNP_BBOX['lon_min'] + GRID_SPACING / 2,
                     MFNP_BBOX['lon_max'],
                     GRID_SPACING)
    lats = np.arange(MFNP_BBOX['lat_min'] + GRID_SPACING / 2,
                     MFNP_BBOX['lat_max'],
                     GRID_SPACING)
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    n_pts = lon_grid.size
    lons_f = lon_grid.flatten()
    lats_f = lat_grid.flatten()

    print(f"  Grid: {len(lons)} × {len(lats)} = {n_pts:,} points")
    print(f"  Months: {n_months}  ({dates[0]} → {dates[-1]})")

    # ── Spatial background VV (dB) ────────────────────────────────────────────
    lon_norm = (lons_f - MFNP_BBOX['lon_min']) / (MFNP_BBOX['lon_max'] - MFNP_BBOX['lon_min'])
    lat_norm = (lats_f - MFNP_BBOX['lat_min']) / (MFNP_BBOX['lat_max'] - MFNP_BBOX['lat_min'])

    # Base backscatter: savanna -13 dB, slight E→W moisture gradient
    bg_vv = -13.0 + 1.5 * (1 - lon_norm) - 1.0 * lat_norm

    # Lake Albert shoreline (western edge): low backscatter (water)
    lake_mask = lons_f < 31.58
    bg_vv[lake_mask] -= np.clip(4.0 * (31.58 - lons_f[lake_mask]) / 0.18, 0, 5.0)

    # Albert Nile corridor: wetter, slightly lower mean
    nile_mask = ((lons_f > 31.55) & (lons_f < 31.78) &
                 (lats_f > 2.05)  & (lats_f < 2.50))
    bg_vv[nile_mask] -= 1.5

    # Murchison Falls gorge vicinity: rocky, higher backscatter
    falls_dist = np.sqrt((lons_f - 31.682)**2 + (lats_f - 2.278)**2)
    bg_vv += 1.5 * np.exp(-falls_dist / 0.05)

    # Seasonal amplitude: larger near wetlands and Nile corridor
    seas_amp = 1.5 + 1.0 * nile_mask.astype(float)

    # ── Subsidence hotspot catalogue ──────────────────────────────────────────
    # (lon_center, lat_center, radius_deg, trend_dB_per_month)
    hotspots = [
        (31.920, 2.180, 0.08, -0.040),   # Pakuba — riverbank compaction
        (32.050, 1.985, 0.06, -0.060),   # Southern sector
        (31.700, 2.350, 0.07, -0.030),   # Albert Nile bank
        (32.120, 2.285, 0.05, -0.050),   # Northern boundary
        (31.780, 1.885, 0.06, -0.040),   # Buligi peninsula
    ]

    # ── Build long-form dataframe ─────────────────────────────────────────────
    rows = []
    for t_idx, date_str in enumerate(dates):
        mo_num = int(date_str.split('-')[1]) - 1  # 0-indexed calendar month

        # Uganda double-peaked seasonal signal (equatorial, ~2 wet seasons/year)
        seasonal = seas_amp * (
            0.65 * np.sin(2 * np.pi * mo_num / 12.0 + 0.55) +   # long rains
            0.35 * np.sin(4 * np.pi * mo_num / 12.0 + 1.10)     # short rains
        )

        noise = np.random.normal(0.0, 0.40, n_pts)
        vv    = bg_vv + seasonal + noise

        # Apply subsidence trends (cumulative linear decline per hotspot)
        for hs_lon, hs_lat, hs_r, hs_trend in hotspots:
            dist = np.sqrt((lons_f - hs_lon)**2 + (lats_f - hs_lat)**2)
            hs_m = dist < hs_r
            # Gaussian-weighted trend so edge of hotspot declines more slowly
            weight = np.exp(-0.5 * (dist[hs_m] / (hs_r * 0.6))**2)
            vv[hs_m] += hs_trend * t_idx * weight

        for i in range(n_pts):
            rows.append({
                'longitude': round(float(lons_f[i]), 4),
                'latitude':  round(float(lats_f[i]), 4),
                'date':      date_str,
                'vv_db':     round(float(vv[i]), 3),
                'n_images':  3   # typical S1 revisits per month
            })

    df = pd.DataFrame(rows)

    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
    raw_path = os.path.join(DATA_DIR, 'mfnp_sar_raw_demo.csv')
    os.makedirs(DATA_DIR, exist_ok=True)
    df.to_csv(raw_path, index=False)
    print(f"  Demo raw data saved: {raw_path}  ({len(df):,} rows)")
    return df


# ── GEE export loader ──────────────────────────────────────────────────────────

def load_gee_export(csv_path: str) -> pd.DataFrame:
    """
    Load and normalise a GEE-exported monthly S1 CSV.

    Expected columns (from mfnp_s1_subsidence_export.js, selector list):
        date        YYYY-MM
        longitude   decimal degrees (WGS84)
        latitude    decimal degrees (WGS84)
        VV          median backscatter dB
        n_images    number of S1 scenes composited

    GEE sometimes exports a '.geo' column (GeoJSON Point) instead of
    explicit lon/lat. Both formats are handled below.
    """
    df = pd.read_csv(csv_path)
    print(f"  Loaded: {csv_path}  ({len(df):,} rows, columns: {list(df.columns)})")

    # ── Parse geometry if explicit lat/lon not present ────────────────────────
    if 'longitude' not in df.columns or 'latitude' not in df.columns:
        if '.geo' in df.columns:
            def _parse_geo(g):
                try:
                    coords = json.loads(g)['coordinates']
                    return pd.Series({'longitude': coords[0], 'latitude': coords[1]})
                except Exception:
                    return pd.Series({'longitude': np.nan, 'latitude': np.nan})
            df[['longitude', 'latitude']] = df['.geo'].apply(_parse_geo)
        else:
            raise ValueError(
                "GEE CSV has no 'longitude'/'latitude' columns and no '.geo' column. "
                "Ensure the GEE export uses the 'selectors' parameter to include "
                "longitude and latitude as explicit feature properties."
            )

    # ── Normalise column names ────────────────────────────────────────────────
    df.rename(columns={'VV': 'vv_db'}, inplace=True)

    # ── Validate date format ──────────────────────────────────────────────────
    df['date'] = df['date'].astype(str).str[:7]   # keep YYYY-MM only

    print(f"  Date range: {df['date'].min()} → {df['date'].max()}")
    print(f"  Unique months: {df['date'].nunique()}")
    print(f"  Grid points (approx): {df.groupby(['longitude','latitude']).ngroups:,}")
    return df


# ── Cleaning ───────────────────────────────────────────────────────────────────

def clean_data(df: pd.DataFrame) -> pd.DataFrame:
    """Remove fill values and out-of-range VV readings."""
    n_before = len(df)

    # Drop GEE fill (-999)
    df = df[df['vv_db'] > NODATA_VALUE + 1].copy()

    # Drop physically implausible values
    df = df[df['vv_db'].between(*VV_VALID_RANGE)].copy()

    # Drop rows missing spatial or temporal key
    df.dropna(subset=['longitude', 'latitude', 'date', 'vv_db'], inplace=True)

    n_after = len(df)
    if n_before > n_after:
        print(f"  Cleaned {n_before - n_after:,} invalid rows "
              f"({(n_before - n_after)/n_before*100:.1f}%)")
    return df


# ── Anomaly computation ────────────────────────────────────────────────────────

def compute_anomalies(df: pd.DataFrame) -> pd.DataFrame:
    """
    Per-pixel, seasonally-adjusted backscatter anomaly.

    Methodology:
      1. Baseline period: first BASELINE_MONTHS months (May 2022–Apr 2023).
         Per pixel, per calendar month → mean VV (removes seasonal bias).
      2. Anomaly (dB) = VV(t) − baseline_vv_same_calendar_month
      3. Z-score = anomaly / per-pixel temporal std (robust to varying means)
      4. Intensity = |Z-score| clipped to [0, 3σ] then scaled → [0, 1]
         0 = stable / no change
         1 = strong anomaly (potential subsidence or surface disturbance)
    """
    dates          = sorted(df['date'].unique())
    baseline_dates = dates[:BASELINE_MONTHS]

    print(f"  Baseline: {baseline_dates[0]} → {baseline_dates[-1]} "
          f"({len(baseline_dates)} months)")
    print(f"  Full analysis: {dates[0]} → {dates[-1]} "
          f"({len(dates)} months)")

    df = df.copy()
    df['cal_month'] = df['date'].str[5:7].astype(int)   # 1–12

    # ── Per-pixel seasonal baseline ────────────────────────────────────────────
    baseline = (df[df['date'].isin(baseline_dates)]
                .groupby(['longitude', 'latitude', 'cal_month'])['vv_db']
                .mean()
                .reset_index()
                .rename(columns={'vv_db': 'baseline_vv'}))

    df = df.merge(baseline, on=['longitude', 'latitude', 'cal_month'], how='left')

    # Where baseline is NaN (no data in baseline window for that calendar month),
    # use the full-pixel mean as a fallback
    pixel_means = (df.groupby(['longitude', 'latitude'])['vv_db']
                   .mean()
                   .reset_index()
                   .rename(columns={'vv_db': 'pixel_mean'}))
    df = df.merge(pixel_means, on=['longitude', 'latitude'], how='left')
    df['baseline_vv'] = df['baseline_vv'].fillna(df['pixel_mean'])

    df['anomaly_db'] = df['vv_db'] - df['baseline_vv']

    # ── Per-pixel Z-score ──────────────────────────────────────────────────────
    pixel_std = (df.groupby(['longitude', 'latitude'])['anomaly_db']
                 .std()
                 .reset_index()
                 .rename(columns={'anomaly_db': 'anomaly_std'}))
    df = df.merge(pixel_std, on=['longitude', 'latitude'], how='left')
    # Clip std from below to avoid division by zero in uniform pixels
    df['anomaly_std'] = df['anomaly_std'].clip(lower=0.05)
    df['anomaly_zscore'] = df['anomaly_db'] / df['anomaly_std']

    # ── Normalised intensity → [0, 1] ─────────────────────────────────────────
    # |Z| capped at 3σ so extreme outliers don't wash out the colormap
    df['intensity_norm'] = (df['anomaly_zscore'].abs().clip(upper=3.0) / 3.0)

    return df


# ── Summary statistics ─────────────────────────────────────────────────────────

def print_summary(df: pd.DataFrame) -> None:
    """Print processing summary for QA."""
    n_pts   = df.groupby(['longitude', 'latitude']).ngroups
    n_dates = df['date'].nunique()
    pct_high = (df['intensity_norm'] > 0.5).mean() * 100.0

    print("\n" + "=" * 60)
    print("Processing Summary")
    print("=" * 60)
    print(f"  Grid points      : {n_pts:,}")
    print(f"  Time steps       : {n_dates}  ({df['date'].min()} → {df['date'].max()})")
    print(f"  Total records    : {len(df):,}")
    print()
    print(f"  VV backscatter   : {df['vv_db'].min():.1f} → {df['vv_db'].max():.1f} dB")
    print(f"  Anomaly (dB)     : {df['anomaly_db'].min():.2f} → {df['anomaly_db'].max():.2f}")
    print(f"  Z-score range    : {df['anomaly_zscore'].min():.2f} → {df['anomaly_zscore'].max():.2f}")
    print(f"  Intensity > 0.5  : {pct_high:.1f}%  (high-change pixels)")
    print()

    # Per-year mean anomaly (drift check)
    df['year'] = df['date'].str[:4]
    yr_anom = df.groupby('year')['anomaly_db'].mean()
    print("  Mean anomaly by year:")
    for yr, val in yr_anom.items():
        bar = '▓' * int(abs(val) * 4) + ('↓' if val < 0 else '↑')
        print(f"    {yr}: {val:+.3f} dB  {bar}")
    print("=" * 60)


# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Process Sentinel-1 SAR data for MFNP subsidence monitoring'
    )
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument('--input', type=str,
                   help='Path to GEE-exported CSV (MFNP_S1_Monthly_VV_2022_2026.csv)')
    g.add_argument('--demo', action='store_true',
                   help='Generate synthetic demo data (no GEE required)')
    parser.add_argument('--output', type=str, default=OUT_CSV,
                        help=f'Output processed CSV path (default: {OUT_CSV})')
    args = parser.parse_args()

    os.makedirs(DATA_DIR, exist_ok=True)
    print("=" * 60)

    if args.demo:
        print("DEMO MODE — generating synthetic Sentinel-1 data for MFNP")
        print("=" * 60)
        df_raw = generate_demo_data(args.output)
    else:
        print(f"Loading GEE export: {args.input}")
        print("=" * 60)
        df_raw = load_gee_export(args.input)

    print("\n[Clean] Removing fill values and out-of-range readings...")
    df_clean = clean_data(df_raw)

    print("\n[Anomaly] Computing seasonal-baseline anomalies and Z-scores...")
    df_proc = compute_anomalies(df_clean)

    print_summary(df_proc)

    # ── Save ──────────────────────────────────────────────────────────────────
    cols_out = [
        'longitude', 'latitude', 'date',
        'vv_db', 'anomaly_db', 'anomaly_zscore', 'intensity_norm'
    ]
    df_proc[cols_out].to_csv(args.output, index=False)
    print(f"\n✓ Processed data → {args.output}")
    print(f"\nNext step:")
    print(f"  python python/mfnp_02_timeslider_map.py --input {args.output}")


if __name__ == '__main__':
    main()
