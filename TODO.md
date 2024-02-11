# Plan

## 1. Tag

-- [v] `Tag` are used to add a name and an optional *namespaced unique* `Id`,
     eg namespace + code to an element (datum, CRS, ...).
     When carrying an `Tag`, an element `PartialEq::eq` implementation
     must consider tag first.

## 2. Geodesic elements

-- [v] ellipsoid: semi major and semi minor axes length in meters,
     inverse flattening (maybe positive infinity if sphere)
-- [v] prime meridian: greenwich longitude in radians
-- [v] datum transformation: N variants (datum shift, Helmert7params,...)
-- [v] geodetic datum:
     tag + ellipsoid + prime meridian + to_reference (tag + datum transformation)

## 3. CRS

-- [v] geocentric axes: 1 variant. Skipped as only 1 variant without parameters
     (only consider meter as unit)
-- [v] geodetic datum + geocentric axes -> Geocentric CRS

-- [v] geodetic axes: N variants (EastNorthUp, NorthEastUp, EastNorth, NorthEast...)
-- [v] geodetic datum + any geodetic axes variant -> Geographic CRS
-- [v] geodetic datum + any projected axes variant + projection -> Projected CRS


## 4. Coordinates Operations

Provide traits for coordinates operations.

Goals:

- handle 1-off transformation
- handle transformation of packed coordinates sequence
- allow building chain of transformation at runtime based on CRS elements
- allow to determine at runtime of a transformation chain is invertible.

-- [v] `trait DynOperation`: 2-way coordinates operation: fwd() (mandatory),
     is_invertible(), bwd() (if is_invertible),
-- [v] `trait Operation`: 1-way direct coordinates operation: apply()
     + apply_seq() with default implementation,
-- [v] `struct Fwd<T: DynOperatiton>: Operation`: select direct transformation.
-- [v] `struct Bwd<T: DynOperation>: Operation`: select reverse transformation.
-- DynOperations:
--- [v] Normalization
--- [v] geographic to geocentric conversion,
--- [ ] Projections
--- [v] GeocentricTranslation (3-parameters),
--- [v] Helmert 7-parameter transformation,
--- [ ] Molodensky-Badekas 10-parameter transformation,
--- [ ] Abridged Molodensky transformation
--- [ ] Geographic Offsets

## 5. CRS transformations

Goals:

- [ ] handle hub transformation using WGS84 geocentric as a hub.
- [ ] handle transformation from EPSG database

-- [ ] TransformationPath: direct transformation + optional inverse transformation 
-- [ ] Provider trait: transform(from: CRS, to: CRS) returns an optional transformation path.
