//// ----------------- Define Regions of Interest ----------------- ////

// Study Regions (buffer by 10km to create representative regional polygons)
var kharkiv = ee.Geometry.Point([36.25, 49.99]).buffer(10000);
var donetsk = ee.Geometry.Point([37.80, 48.01]).buffer(10000);
var luhansk = ee.Geometry.Point([39.30, 48.57]).buffer(10000);
var mariupol = ee.Geometry.Point([37.55, 47.09]).buffer(10000);
var belgorod = ee.Geometry.Point([36.58, 50.59]).buffer(10000);

// Control Regions
var lviv = ee.Geometry.Point([24.03, 49.83]).buffer(10000);
var st_petersburg = ee.Geometry.Point([30.32, 59.93]).buffer(10000);

// Feature Collections for DiD aggregation
var studyRegions = ee.FeatureCollection([
  ee.Feature(kharkiv, {name: 'Kharkiv'}),
  ee.Feature(donetsk, {name: 'Donetsk'}),
  ee.Feature(luhansk, {name: 'Luhansk'}),
  ee.Feature(mariupol, {name: 'Mariupol'}),
  ee.Feature(belgorod, {name: 'Belgorod'})
]);

var controlRegions = ee.FeatureCollection([
  ee.Feature(lviv, {name: 'Lviv'}),
  ee.Feature(st_petersburg, {name: 'St. Petersburg'})
]);

var allRegions = studyRegions.merge(controlRegions);

//// ----------------- Cloud Masking & Scaling Function ----------------- ////

// Robust cloud masking and scaling for Landsat 8 (C02, T1_L2)
function maskL8sr(image) {
  var qa = image.select('QA_PIXEL');

  // Bitmasks
  var fill = qa.bitwiseAnd(1 << 0).eq(0); 
  var dilatedCloud = qa.bitwiseAnd(1 << 1).eq(0);
  var cirrus = qa.bitwiseAnd(1 << 2).eq(0);
  var cloud = qa.bitwiseAnd(1 << 3).eq(0);
  var cloudShadow = qa.bitwiseAnd(1 << 4).eq(0);
  var snow = qa.bitwiseAnd(1 << 5).eq(0);

  // Combine masks
  var mask = fill.and(dilatedCloud).and(cirrus).and(cloud).and(cloudShadow).and(snow);

  // Apply scaling factors for surface reflectance
  var opticalBands = image.select('SR_B.')
      .multiply(0.0000275)
      .add(-0.2);

  return image
      .addBands(opticalBands, null, true)
      .updateMask(mask)
      .copyProperties(image, ["system:time_start"]);
}

//// ----------------- Compute Spectral Indices ----------------- ////

function addIndices(image) {
  // NDVI: (NIR - Red) / (NIR + Red) -> (B5 - B4) / (B5 + B4)
  var ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
  
  // NDBI: (SWIR1 - NIR) / (SWIR1 + NIR) -> (B6 - B5) / (B6 + B5)
  var ndbi = image.normalizedDifference(['SR_B6', 'SR_B5']).rename('NDBI');
  
  return image.addBands([ndvi, ndbi]);
}

//// ----------------- Processing Baseline vs Conflict ----------------- ////

var l8_collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(allRegions.geometry());

// Baseline Period: 2013
var baseline_2013 = l8_collection
  .filterDate('2013-04-01', '2013-12-31')
  .map(maskL8sr)
  .map(addIndices)
  .median()
  .clip(allRegions.geometry());

// Conflict Period: 2024
var conflict_2024 = l8_collection
  .filterDate('2024-01-01', '2024-12-31')
  .map(maskL8sr)
  .map(addIndices)
  .median()
  .clip(allRegions.geometry());

// Compute the chronological difference (2024 - 2013)
// Positive values = Increase in index; Negative values = Decrease in index
var difference_img = conflict_2024.subtract(baseline_2013);

//// ----------------- Compute Difference in Differences ----------------- ////

// Function to calculate mean change for a given collection
function getMeanChange(featureCollection, indexBand) {
  var reduced = difference_img.select(indexBand).reduceRegions({
    collection: featureCollection,
    reducer: ee.Reducer.mean(),
    scale: 30, // Landsat resolution
    tileScale: 4
  });
  // Aggregate the mean across all regions in the group
  return reduced.aggregate_mean('mean');
}

// Extract values for NDVI (Vegetation)
var study_ndvi_change = getMeanChange(studyRegions, 'NDVI');
var control_ndvi_change = getMeanChange(controlRegions, 'NDVI');
var did_ndvi = ee.Number(study_ndvi_change).subtract(ee.Number(control_ndvi_change));

// Extract values for NDBI (Urban Structure)
var study_ndbi_change = getMeanChange(studyRegions, 'NDBI');
var control_ndbi_change = getMeanChange(controlRegions, 'NDBI');
var did_ndbi = ee.Number(study_ndbi_change).subtract(ee.Number(control_ndbi_change));

// Print results to the Console
print('--- Difference-in-Differences (DiD) Results ---');
print('1. Vegetation (NDVI) Change');
print('Study Regions NDVI Change (Mean):', study_ndvi_change);
print('Control Regions NDVI Change (Mean):', control_ndvi_change);
print('DiD Estimator (NDVI isolated impact):', did_ndvi);

print('2. Urban Built-up (NDBI) Change');
print('Study Regions NDBI Change (Mean):', study_ndbi_change);
print('Control Regions NDBI Change (Mean):', control_ndbi_change);
print('DiD Estimator (NDBI isolated impact):', did_ndbi);

//// ----------------- Visualizatios ----------------- ////

// Center map on Kharkiv for primary viewing
Map.setCenter(36.25, 49.99, 10);

// Map styles
var diff_vis_veg = {
  bands: ['NDVI'],
  min: -0.2,
  max: 0.2,
  palette: ['red', 'white', 'green'] // Red = vegetation loss, Green = gain
};

var diff_vis_urban = {
  bands: ['NDBI'],
  min: -0.15,
  max: 0.15,
  palette: ['blue', 'white', 'orange'] // Blue = decreased urban signature, Orange = increased
};

// Add layers
Map.addLayer(baseline_2013, {bands: ['SR_B4', 'SR_B3', 'SR_B2'], min: 0, max: 0.15}, 'Baseline RGB (2013)', false);
Map.addLayer(conflict_2024, {bands: ['SR_B4', 'SR_B3', 'SR_B2'], min: 0, max: 0.15}, 'Conflict RGB (2024)', false);

Map.addLayer(difference_img, diff_vis_veg, 'NDVI Change (2024-2013)', true);
Map.addLayer(difference_img, diff_vis_urban, 'NDBI Change (2024-2013)', true);

// Add study/control boundaries for context
var emptyStyle = {color: '000000', fillColor: '00000000'};
Map.addLayer(studyRegions.style(emptyStyle), {}, 'Study Regions');
Map.addLayer(controlRegions.style({color: '888888', fillColor: '00000000'}), {}, 'Control Regions');