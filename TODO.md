# List of things to (maybe) do

## TODO

- [v] Define a 'Float' type alias to 'f64' in 'crate::math'
- [v] Define a ~~'Real'~~ 'Length' type using a 'Float' in 'crate::cs::r1'
- [v] Use 'Float' wherever we're using 'f64': 
    - 'Length'
    - 'Angle'
- [v] Use 'Length' to define coordinates types in 'cs::cartesian': 
    - 'Length'
- [v] Define type alias to 'Length' for 'Height' coordinates types in 'cs::geodetic'

- [v] Remove DatumTransformation from 'GeodeticDatum' 

- [ ] Rework geodesic unit test
- [ ] Use Clenshaw summation
- [ ] Write solve_inverse for Karney solver
- [ ] Full integration testing of geodesic

## Maybe DO

- [ ] Move ellipsoids, prime meridian and datum constants into a provider module (maybe a feature)
- [ ] Decide whether to provide well-known ellipsoids, prime meridians, datums (and datum toWGS84) via constants ? via files ? 
      behing a feature ? 
- [ ] Add some *provider* features (core, epsg...) that provides ellipsoid, prime meridian, datum (datum transformation ?) 


