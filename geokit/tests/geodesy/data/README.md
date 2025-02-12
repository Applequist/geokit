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

