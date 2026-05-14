// =============================================================================
// Script 01: Sentinel-1 SAR Backscatter Temporal Analysis
//            Land Subsidence Proxy — Murchison Falls National Park, Uganda
// =============================================================================
// Objective : Detect areas of anomalous surface change using Sentinel-1 GRD
//             backscatter time-series over a 4-year window (May 2022 → latest).
//
//             True InSAR (phase-based) processing is NOT natively available in
//             GEE. For mm-precision vertical displacement use:
//               • Copernicus Ground Motion Service (EGMS) — Europe-focused
//               • NASA ARIA / MintPy: https://aria.jpl.nasa.gov/
//               • ESA SNAP + StaMPS / SARscape for SLC-based processing
//               • Prepare Sentinel-1 SLCs via ASF HyP3 (free SBAS/InSAR):
//                 https://hyp3-docs.asf.alaska.edu/
//             This script computes a FREE, fast backscatter-variance proxy
//             suitable for hotspot screening and 4-year trend monitoring.
//
// Outputs:
//   (a) Monthly median VV composites (48 months) → sampled CSV for Python
//   (b) Temporal statistics rasters (mean, std-dev, Z-score)
//   (c) Subsidence proxy index raster
//
// Export target: Google Drive → folder "MFNP_Subsidence"
//
// Pipeline:
//   GEE (this script) → data/MFNP_S1_Monthly_VV_2022_2026.csv
//                      → python/mfnp_01_process_sar.py
//                      → python/mfnp_02_timeslider_map.py
//                      → web/mfnp_subsidence_timeslider.html
// =============================================================================

// ── 1. STUDY AREA — MURCHISON FALLS NATIONAL PARK ────────────────────────────
// GEE Rectangle format: [west, south, east, north] (WGS84)
// MFNP: ~3,840 km², spans the Albert Nile / Victoria Nile confluence
var mfnp = ee.Geometry.Rectangle([31.40, 1.75, 32.20, 2.55]);

Map.centerObject(mfnp, 9);
Map.addLayer(mfnp, {color: 'FFD700'}, 'MFNP Study Area', true);

// ── 2. DATE RANGE — 4-year window ending at latest available data ────────────
var startDate = ee.Date('2022-05-01');
var endDate   = ee.Date('2026-05-01');   // GEE will clip to available data

print('Analysis period:', startDate.format('YYYY-MM-dd'), '→',
      endDate.format('YYYY-MM-dd'));

// ── 3. LOAD SENTINEL-1 GRD COLLECTION ────────────────────────────────────────
// IW mode (Interferometric Wide Swath), VV polarisation, ascending pass.
// NOTE: For Uganda, ascending pass overpass is ~17:00-18:00 local time.
//       Mixing ascending + descending is not recommended for subsidence
//       analysis due to different look angles (LOS geometry).
var s1_all = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(mfnp)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .select('VV');

// Separate ascending / descending for QA (use ascending for analysis)
var s1_asc  = s1_all.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));
var s1_desc = s1_all.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));

print('Total S1 images (both passes):', s1_all.size());
print('Ascending pass images:', s1_asc.size());
print('Descending pass images:', s1_desc.size());

// Primary collection: ascending pass
var s1 = s1_asc;

// ── 4. TEMPORAL STATISTICS (full period) ─────────────────────────────────────
var meanVV = s1.mean().clip(mfnp);
var stdVV  = s1.reduce(ee.Reducer.stdDev()).clip(mfnp);
var minVV  = s1.min().clip(mfnp);
var maxVV  = s1.max().clip(mfnp);

// ── 5. Z-SCORE ANOMALY MAP ────────────────────────────────────────────────────
// Z = (mean - min) / std : high values = sudden backscatter change (disturbance)
// Negative anomaly trend = possible subsidence / vegetation loss signal
var zScore = meanVV.subtract(minVV)
                   .divide(stdVV.add(1e-6))
                   .rename('z_score');

// ── 6. SUBSIDENCE PROXY INDEX ─────────────────────────────────────────────────
// High temporal std relative to mean → likely land-use change / surface disturbance
// Areas: sediment compaction (Albert Nile delta), oil exploration zones,
//        wetland drainage, forest clearance
var subsidenceProxy = stdVV.divide(meanVV.abs().add(1e-6))
                           .rename('subsidence_proxy');

// ── 7. HOTSPOT MASK (std-dev > 2.5 dB) ──────────────────────────────────────
var hotspots = stdVV.gt(2.5).selfMask().rename('hotspot');

// ── 8. VISUALISATION ─────────────────────────────────────────────────────────
var visMean  = {min: -25, max: -5,  palette: ['000000','404040','808080','cccccc','ffffff']};
var visStd   = {min: 0,   max: 6,   palette: ['0000ff','00aaff','00ff00','ffff00','ff8000','ff0000']};
var visProxy = {min: 0,   max: 0.5, palette: ['2b83ba','abdda4','ffffbf','fdae61','d7191c']};
var visHot   = {palette: ['ff0000']};

Map.addLayer(meanVV,        visMean,  'S1 Mean Backscatter VV (dB)',   false);
Map.addLayer(stdVV,         visStd,   'S1 Temporal Std-Dev (Instability)');
Map.addLayer(subsidenceProxy, visProxy, 'Subsidence Proxy Index');
Map.addLayer(hotspots,      visHot,   'Hotspots (StdDev > 2.5 dB)',   true);

// ── 9. MONTHLY TIME SERIES CHART ─────────────────────────────────────────────
// Representative point: near Murchison Falls
var fallsPoint = ee.Geometry.Point([31.682, 2.278]);

var tsChart = ui.Chart.image.series({
  imageCollection: s1,
  region:          fallsPoint,
  reducer:         ee.Reducer.mean(),
  scale:           30
}).setOptions({
  title:     'Sentinel-1 VV Backscatter — Murchison Falls (4-year)',
  hAxis:     {title: 'Date'},
  vAxis:     {title: 'Backscatter (dB)'},
  lineWidth: 2,
  pointSize: 3,
  colors:    ['#1a73e8']
});
print(tsChart);

// Albert Nile delta point (potential sediment compaction / subsidence)
var nilePoint = ee.Geometry.Point([31.65, 2.40]);
var nileChart = ui.Chart.image.series({
  imageCollection: s1,
  region:          nilePoint,
  reducer:         ee.Reducer.mean(),
  scale:           30
}).setOptions({
  title:  'Sentinel-1 VV — Albert Nile Delta (potential subsidence zone)',
  hAxis:  {title: 'Date'},
  vAxis:  {title: 'Backscatter (dB)'},
  colors: ['#e8711a']
});
print(nileChart);

// ── 10. MONTHLY COMPOSITE EXPORT (for Python time-slider) ───────────────────
// Strategy: For each of the 48 months, compute median VV composite,
// sample to 1 km grid, add date property, then merge all into one
// FeatureCollection and export as a single CSV.
//
// Output CSV columns:
//   .geo   (GeoJSON Point: {type:"Point",coordinates:[lon,lat]})
//   date   (YYYY-MM string)
//   VV     (median dB for that month)
//   n_images (number of S1 scenes averaged, for QA)
// =============================================================================

var startYear  = 2022;
var startMonth = 5;   // May
var nMonths    = 48;  // 4 years

// Build list of [year, month] pairs
var monthList = ee.List.sequence(0, nMonths - 1).map(function(offset) {
  var d = startDate.advance(ee.Number(offset), 'month');
  return d.format('YYYY-MM');
});
print('Month list (first 6):', monthList.slice(0, 6));

// For each month: filter S1, compute median, sample 1km grid
var makeMonthlyFC = function(dateStr) {
  var date    = ee.Date(ee.String(dateStr).cat('-01'));
  var nextMon = date.advance(1, 'month');

  var monthly = s1.filterDate(date, nextMon);
  var count   = monthly.size();
  var img     = ee.Algorithms.If(
    count.gt(0),
    monthly.median(),
    ee.Image.constant(-999).rename('VV')  // fill value when no images
  );
  img = ee.Image(img).clip(mfnp);

  // Add pixel lon/lat as bands for export without geometry parsing issues
  var lonlat = ee.Image.pixelLonLat();
  var imgWithCoords = ee.Image(img).addBands(lonlat);

  // Sample the image at 1 km grid
  var sampled = imgWithCoords.sample({
    region:     mfnp,
    scale:      1000,            // 1 km grid
    projection: 'EPSG:4326',
    seed:       42,
    geometries: false            // coords already in lon/latitude bands
  });

  // Rename and add date property
  sampled = sampled.map(function(f) {
    return f.set('date', dateStr)
             .set('n_images', count)
             .set('longitude', f.get('longitude'))
             .set('latitude',  f.get('latitude'))
             .set('VV',        f.get('VV'));
  });

  return sampled;
};

// Map over month list → list of FeatureCollections → merge
var allMonthly = ee.FeatureCollection(
  monthList.map(makeMonthlyFC)
).flatten();

print('Total sampled features (approx):', allMonthly.size());
print('First feature:', allMonthly.first());

// ── 11. EXPORT: MONTHLY CSV (primary output for Python pipeline) ──────────────
Export.table.toDrive({
  collection:  allMonthly,
  description: 'MFNP_S1_Monthly_VV_2022_2026',
  folder:      'MFNP_Subsidence',
  fileFormat:  'CSV',
  selectors:   ['date', 'longitude', 'latitude', 'VV', 'n_images']
});

// ── 12. EXPORT: SUBSIDENCE PROXY RASTERS ─────────────────────────────────────
// UTM Zone 36N (EPSG:32636) for metric accuracy over Uganda
var ugandaCRS = 'EPSG:32636';

Export.image.toDrive({
  image:          subsidenceProxy,
  description:    'MFNP_Subsidence_Proxy_4yr',
  folder:         'MFNP_Subsidence',
  fileNamePrefix: 'mfnp_subsidence_proxy',
  region:         mfnp,
  scale:          20,
  crs:            ugandaCRS,
  maxPixels:      1e10
});

Export.image.toDrive({
  image:          stdVV.rename('VV_temporal_std'),
  description:    'MFNP_S1_Temporal_StdDev',
  folder:         'MFNP_Subsidence',
  fileNamePrefix: 'mfnp_s1_temporal_std',
  region:         mfnp,
  scale:          20,
  crs:            ugandaCRS,
  maxPixels:      1e10
});

Export.image.toDrive({
  image:          hotspots,
  description:    'MFNP_Subsidence_Hotspots',
  folder:         'MFNP_Subsidence',
  fileNamePrefix: 'mfnp_subsidence_hotspots',
  region:         mfnp,
  scale:          20,
  crs:            ugandaCRS,
  maxPixels:      1e10
});

// ── 13. QA STATISTICS ─────────────────────────────────────────────────────────
var stats = subsidenceProxy.reduceRegion({
  reducer:   ee.Reducer.percentile([10, 25, 50, 75, 90, 95]),
  geometry:  mfnp,
  scale:     100,
  maxPixels: 1e9
});
print('Subsidence Proxy Percentiles:', stats);

var meanStats = meanVV.reduceRegion({
  reducer:   ee.Reducer.mean().combine(ee.Reducer.stdDev(), '', true),
  geometry:  mfnp,
  scale:     100,
  maxPixels: 1e9
});
print('Mean VV Statistics:', meanStats);

// ── NOTE ON TRUE InSAR ────────────────────────────────────────────────────────
// For mm-precision subsidence over MFNP use:
// 1. Download Sentinel-1 SLC from Copernicus Open Access Hub (ASF Mirror)
// 2. Process with ESA HyP3 SBAS (cloud-based, free for research):
//    https://hyp3-docs.asf.alaska.edu/
// 3. Or use MintPy for time-series InSAR:
//    https://mintpy.readthedocs.io/
// 4. EGMS covers Europe only — not applicable to Uganda.
// This backscatter proxy approach is free, fast, globally available,
// and suitable for 1-5 dB level surface change screening.
