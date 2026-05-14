"""
mfnp_02_timeslider_map.py
=========================
Generate an interactive Folium map with a 4-year monthly time slider
showing Sentinel-1 SAR backscatter anomaly (subsidence proxy) over
Murchison Falls National Park, Uganda.

Layers:
  1. HeatMapWithTime — 48 monthly SAR anomaly frames (primary layer)
  2. MFNP park boundary polygon (gold outline)
  3. Albert Nile + Victoria Nile river centrelines
  4. Key landmarks (Murchison Falls, Paraa, park gates)
  5. Multiple basemaps (Dark, Satellite, OSM)
  6. Legend panel + title bar + minimap + fullscreen control

Output:
    web/mfnp_subsidence_timeslider.html   (~self-contained, publishable)

Usage:
    python python/mfnp_02_timeslider_map.py
    python python/mfnp_02_timeslider_map.py --input data/mfnp_sar_processed.csv
    python python/mfnp_02_timeslider_map.py --input data/mfnp_sar_processed.csv \\
                                             --output web/mfnp_subsidence_timeslider.html
"""

import os
import argparse
import numpy as np
import pandas as pd
import folium
from folium import plugins
from folium.plugins import (
    HeatMapWithTime, MeasureControl, MiniMap, Fullscreen, LocateControl
)
import warnings
warnings.filterwarnings('ignore')

# ── Configuration ──────────────────────────────────────────────────────────────
DATA_DIR  = 'data'
WEB_DIR   = 'web'
OUT_HTML  = os.path.join(WEB_DIR, 'mfnp_subsidence_timeslider.html')
IN_CSV    = os.path.join(DATA_DIR, 'mfnp_sar_processed.csv')

MFNP_CENTER  = [2.15, 31.80]   # [lat, lon] map centre
MFNP_ZOOM    = 9

# ── Park boundary (simplified polygon, ~3,840 km²) ────────────────────────────
# Approximated from GADM / Protected Planet; replace with exact UNMA boundary
# for publication-grade mapping. Format: [[lat, lon], ...]
MFNP_BOUNDARY = [
    [2.515, 31.530],  # NW — Lake Albert shore
    [2.550, 31.690],
    [2.510, 31.850],
    [2.455, 31.975],
    [2.380, 32.140],
    [2.215, 32.180],
    [2.015, 32.135],
    [1.825, 32.100],
    [1.755, 31.905],
    [1.780, 31.720],
    [1.845, 31.620],
    [1.985, 31.535],
    [2.200, 31.480],
    [2.345, 31.485],
    [2.450, 31.505],
    [2.515, 31.530],  # close ring
]

# ── River centrelines (simplified, for visual reference) ─────────────────────
ALBERT_NILE = [        # flows N from Lake Albert → Sudan
    [2.515, 31.540],
    [2.415, 31.600],
    [2.310, 31.658],
    [2.200, 31.700],
    [2.100, 31.685],
    [1.990, 31.660],
    [1.870, 31.625],
]

VICTORIA_NILE = [      # E→W from Lake Victoria through Murchison Falls
    [2.230, 32.140],
    [2.195, 32.050],
    [2.160, 31.950],
    [2.130, 31.860],
    [2.100, 31.780],
    [2.278, 31.682],   # Murchison Falls gorge
    [2.250, 31.650],
    [2.210, 31.610],
    [2.155, 31.580],   # confluence with Albert Nile
]

# ── Key landmarks: [lat, lon, emoji_label, tooltip_text] ─────────────────────
LANDMARKS = [
    [2.278, 31.682, '💧', 'Murchison Falls — Victoria Nile drops 43 m through 7 m gorge'],
    [2.295, 31.660, '⛴️', 'Paraa — Park HQ & Nile ferry crossing'],
    [2.160, 32.115, '🚩', 'Karuma Gate — eastern entrance'],
    [1.812, 31.620, '🚩', 'Tangi Gate — southern entrance'],
    [2.400, 31.545, '🌊', 'Lake Albert shoreline'],
    [1.950, 31.550, '🦁', 'Buligi Game Reserve'],
]

# ── Colour gradient for HeatMapWithTime ──────────────────────────────────────
# 0.0 → 1.0: stable (blue) → moderate (cyan/green) → elevated (yellow/orange) → high (red)
SAR_GRADIENT = {
    0.00: '#000080',   # navy — no anomaly
    0.15: '#0000ff',   # blue
    0.30: '#00aaff',   # light blue
    0.45: '#00ffaa',   # cyan-green
    0.60: '#ffff00',   # yellow
    0.75: '#ff8800',   # orange
    1.00: '#ff0000',   # red — strong anomaly
}

MONTH_NAMES = [
    'Jan','Feb','Mar','Apr','May','Jun',
    'Jul','Aug','Sep','Oct','Nov','Dec'
]


# ── Data preparation ───────────────────────────────────────────────────────────

def build_heatmap_data(df: pd.DataFrame):
    """
    Convert processed SAR DataFrame to HeatMapWithTime input format.

    Returns:
        data_list  : list of n_timesteps lists, each = [[lat, lon, intensity], ...]
        date_labels: list of 'Mon YYYY' strings for the slider index
    """
    dates = sorted(df['date'].unique())
    data_list   = []
    date_labels = []

    for d in dates:
        yr, mo = d.split('-')
        date_labels.append(f"{MONTH_NAMES[int(mo)-1]} {yr}")

        monthly  = df[df['date'] == d]
        # Filter low-intensity points to keep HTML size manageable
        monthly  = monthly[monthly['intensity_norm'] > 0.04]
        pts      = monthly[['latitude', 'longitude', 'intensity_norm']].values.tolist()
        data_list.append(pts)

    return data_list, date_labels


# ── HTML overlay helpers ───────────────────────────────────────────────────────

def legend_html() -> str:
    return """
<div id="sar-legend" style="
    position: fixed;
    bottom: 60px; left: 12px;
    z-index: 1000;
    background: rgba(10, 10, 28, 0.92);
    color: #dde8f0;
    padding: 14px 16px;
    border-radius: 10px;
    border: 1px solid rgba(100, 160, 220, 0.35);
    font-family: 'Segoe UI', Arial, sans-serif;
    font-size: 11px;
    box-shadow: 0 4px 18px rgba(0,0,0,0.6);
    max-width: 230px;
    line-height: 1.5;
">
  <div style="font-size:13px;font-weight:700;color:#7ec8f4;margin-bottom:8px;">
    📡 SAR Anomaly — Subsidence Proxy
  </div>

  <!-- Gradient bar -->
  <div style="
    width:100%;height:10px;
    background:linear-gradient(to right,#000080,#0000ff,#00aaff,#00ffaa,#ffff00,#ff8800,#ff0000);
    border-radius:4px;margin-bottom:4px;">
  </div>
  <div style="display:flex;justify-content:space-between;
              font-size:9px;color:#99aabb;margin-bottom:10px;">
    <span>Stable</span><span>Moderate</span><span>High Change</span>
  </div>

  <div style="border-top:1px solid rgba(100,160,220,0.25);padding-top:8px;">
    <div>📏 <b>Sensor:</b> Sentinel-1 IW VV (ascending)</div>
    <div>📅 <b>Period:</b> May 2022 – Apr 2026</div>
    <div>📐 <b>Grid:</b> ~1 km × 1 km</div>
    <div>📊 <b>Metric:</b> |ΔVV| / σ per pixel (Z-score)</div>
  </div>

  <div style="margin-top:8px;padding-top:6px;
              border-top:1px solid rgba(100,160,220,0.25);
              font-size:9px;color:#e09090;">
    ⚠️ SAR backscatter variance proxy — not true InSAR.<br>
    High values indicate surface change signal.<br>
    For mm-precision see: ASF HyP3 / MintPy InSAR.
  </div>

  <div style="margin-top:7px;font-size:9px;color:#556677;">
    by Prakash Krishnamachari · 2025
  </div>
</div>
"""


def title_html(n_months: int, first_label: str, last_label: str) -> str:
    return f"""
<div style="
    position: fixed;
    top: 10px; left: 50%; transform: translateX(-50%);
    z-index: 1000;
    background: rgba(8, 12, 30, 0.90);
    color: #e0eef8;
    padding: 10px 22px;
    border-radius: 9px;
    border: 1px solid rgba(90, 160, 230, 0.50);
    font-family: 'Segoe UI', Arial, sans-serif;
    text-align: center;
    pointer-events: none;
    box-shadow: 0 2px 14px rgba(0,60,180,0.45);
    white-space: nowrap;
">
  <div style="font-size:15px;font-weight:700;letter-spacing:0.5px;">
    🛰️ MURCHISON FALLS NATIONAL PARK
  </div>
  <div style="font-size:10px;color:#80b8e0;margin-top:3px;">
    Sentinel-1 SAR Subsidence Monitor &nbsp;·&nbsp;
    {n_months} monthly frames &nbsp;·&nbsp; {first_label} → {last_label}
  </div>
</div>
"""


# ── Map builder ────────────────────────────────────────────────────────────────

def build_map(df: pd.DataFrame) -> folium.Map:
    """Construct and return the full Folium time-slider map."""

    # ── Base map ─────────────────────────────────────────────────────────────
    m = folium.Map(
        location=MFNP_CENTER,
        zoom_start=MFNP_ZOOM,
        tiles=None,
        control_scale=True,
        prefer_canvas=True
    )

    # ── Basemap tile layers ───────────────────────────────────────────────────
    folium.TileLayer(
        'CartoDB dark_matter',
        name='🌑 Dark Basemap (default)',
        attr='© CartoDB, © OpenStreetMap contributors',
        show=True
    ).add_to(m)

    folium.TileLayer(
        tiles=(
            'https://server.arcgisonline.com/ArcGIS/rest/services/'
            'World_Imagery/MapServer/tile/{z}/{y}/{x}'
        ),
        name='🛰️ Esri Satellite',
        attr='© Esri, Maxar, GeoEye, USGS',
        show=False
    ).add_to(m)

    folium.TileLayer(
        'OpenStreetMap',
        name='🗺️ OpenStreetMap',
        show=False
    ).add_to(m)

    # ── Layer 1: HeatMapWithTime — SAR anomaly (primary) ─────────────────────
    print("  Building HeatMapWithTime (48 monthly frames)...")
    data_list, date_labels = build_heatmap_data(df)

    HeatMapWithTime(
        data             = data_list,
        index            = date_labels,
        name             = '📡 SAR Anomaly Heatmap',
        auto_play        = False,
        max_opacity      = 0.88,
        min_opacity      = 0.00,
        radius           = 14,
        blur             = 0.85,   # 0–1 scale in folium ≥0.15
        gradient         = SAR_GRADIENT,
        display_index    = True,
        min_speed        = 0.1,
        max_speed        = 10,
        speed_step       = 0.5,
        position         = 'bottomright'
    ).add_to(m)

    # ── Layer 2: Park boundary ────────────────────────────────────────────────
    boundary_layer = folium.FeatureGroup(
        name='🟨 MFNP Park Boundary', show=True)
    folium.Polygon(
        locations = MFNP_BOUNDARY,
        color     = '#FFD700',       # gold
        weight    = 2.5,
        opacity   = 0.90,
        fill      = False,
        tooltip   = 'Murchison Falls National Park (~3,840 km²) · Uganda'
    ).add_to(boundary_layer)
    boundary_layer.add_to(m)

    # ── Layer 3: River network ────────────────────────────────────────────────
    rivers_layer = folium.FeatureGroup(
        name='🔵 Nile River Network', show=True)

    folium.PolyLine(
        ALBERT_NILE,
        color='#4fc3f7', weight=3.0, opacity=0.75,
        tooltip='Albert Nile (Nile tributary flowing north to Sudan)'
    ).add_to(rivers_layer)

    folium.PolyLine(
        VICTORIA_NILE,
        color='#29b6f6', weight=2.5, opacity=0.75,
        tooltip='Victoria Nile / Murchison Falls section'
    ).add_to(rivers_layer)

    rivers_layer.add_to(m)

    # ── Layer 4: Landmarks ────────────────────────────────────────────────────
    landmarks_layer = folium.FeatureGroup(
        name='📍 Landmarks', show=True)

    for lat, lon, emoji, tip in LANDMARKS:
        folium.Marker(
            location = [lat, lon],
            tooltip  = tip,
            popup    = folium.Popup(
                f'<div style="font-family:sans-serif;font-size:12px;">'
                f'{emoji} <b>{tip}</b></div>',
                max_width=240
            ),
            icon = folium.DivIcon(
                html=f'<div style="font-size:20px;'
                     f'text-shadow:0 0 5px #000,0 0 10px #000;">{emoji}</div>',
                icon_size   = (28, 28),
                icon_anchor = (14, 14)
            )
        ).add_to(landmarks_layer)

    landmarks_layer.add_to(m)

    # ── Plugins ──────────────────────────────────────────────────────────────
    folium.LayerControl(collapsed=False, position='topright').add_to(m)
    MiniMap(tile_layer='CartoDB dark_matter', toggle_display=True,
            position='bottomleft').add_to(m)
    MeasureControl(position='topleft',
                   primary_length_unit='kilometers').add_to(m)
    Fullscreen(position='topright').add_to(m)
    LocateControl(position='topright').add_to(m)

    # ── Overlays ──────────────────────────────────────────────────────────────
    m.get_root().html.add_child(folium.Element(legend_html()))

    dates      = sorted(df['date'].unique())
    yr0, mo0   = dates[0].split('-')
    yr1, mo1   = dates[-1].split('-')
    first_lbl  = f"{MONTH_NAMES[int(mo0)-1]} {yr0}"
    last_lbl   = f"{MONTH_NAMES[int(mo1)-1]} {yr1}"
    m.get_root().html.add_child(
        folium.Element(title_html(len(dates), first_lbl, last_lbl))
    )

    return m


# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Generate MFNP SAR subsidence time-slider Folium map'
    )
    parser.add_argument('--input',  type=str, default=IN_CSV,
                        help=f'Processed SAR CSV (default: {IN_CSV})')
    parser.add_argument('--output', type=str, default=OUT_HTML,
                        help=f'Output HTML path (default: {OUT_HTML})')
    args = parser.parse_args()

    print("=" * 60)
    print("MFNP Subsidence Time-Slider Map Builder")
    print("=" * 60)

    # ── Load data ─────────────────────────────────────────────────────────────
    if not os.path.exists(args.input):
        print(f"ERROR: Processed CSV not found at {args.input}")
        print("Run first: python python/mfnp_01_process_sar.py --demo")
        return

    print(f"\n[Load] {args.input}")
    df = pd.read_csv(args.input)
    print(f"  {len(df):,} records | "
          f"{df['date'].nunique()} months | "
          f"{df.groupby(['longitude','latitude']).ngroups:,} grid points")
    print(f"  Date range: {df['date'].min()} → {df['date'].max()}")

    # ── Build map ─────────────────────────────────────────────────────────────
    print("\n[Build] Constructing Folium map...")
    m = build_map(df)

    # ── Save ──────────────────────────────────────────────────────────────────
    os.makedirs(WEB_DIR, exist_ok=True)
    m.save(args.output)

    size_kb = os.path.getsize(args.output) / 1024
    n_frames = df['date'].nunique()
    max_pts  = max(
        len(df[(df['date'] == d) & (df['intensity_norm'] > 0.04)])
        for d in df['date'].unique()
    )

    print(f"\n✓ Map saved → {args.output}")
    print(f"  File size : {size_kb:.0f} KB")
    print(f"  Frames    : {n_frames} monthly steps")
    print(f"  Max points/frame: {max_pts:,}")
    print(f"\nTo publish on GitHub Pages:")
    print(f"  1. Copy {args.output} to your GitHub repo web/ folder")
    print(f"  2. Update web/index.html iframe src to point to this file")
    print(f"  3. Push and enable GitHub Pages from the /web folder")


if __name__ == '__main__':
    main()
