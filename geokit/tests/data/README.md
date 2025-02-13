# Geodesy test data

## xyz_to_llh_10_000.txt

Contains 10_000 test cases for geocentric to geodetic coordinates conversion.

Each line contains geodentric coordinates x, y, z followed 
by the cooresponding geodetic coordinates lon, lat, height.

Geocentric coordinates x, y and z are generated from a uniform sample on S2,
then scaled using WGS84 a and b. 
Geodetic coordinates lon, lat and height are converted from geocentric coordinates
using the following proj command:
```bash
cct -I +proj=cart +ellps=WGS84
```

## GeodTest-short.dat.gz 

Obtained from https://sourceforge.net/projects/geographiclib/files/testdata/GeodTest.dat.gz

Contains 10_000 geodesics test cases using WGS84 ellipsoid.

Each line gives 10 space separated numbers:
- [0]: latitude at point 1, lat1 (degrees, exact)
- [1]: longitude at point 1, lon1 (degrees, always 0)
- [2]: azimuth at point 1, azi1 (clockwise from north in degrees, exact)
- [3]: latitude at point 2, lat2 (degrees, accurate to 10−18 deg)
- [4]: longitude at point 2, lon2 (degrees, accurate to 10−18 deg)
- [5]: azimuth at point 2, azi2 (degrees, accurate to 10−18 deg)
- [6]: geodesic distance from point 1 to point 2, s12 (meters, exact)
- [7]: arc distance on the auxiliary sphere, a12 (degrees, accurate to 10−18 deg)
- [8]: reduced length of the geodesic, m12 (meters, accurate to 0.1 pm)
- [9]: the area under the geodesic, S12 (m2, accurate to 1 mm2)

