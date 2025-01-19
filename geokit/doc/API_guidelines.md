# API Guidelines

## Units, values and ranges

So we *could* 

1. Use newtype for:
- Angles in radians produced by multiplying a f64 value by a unit value, eg '45. * Deg',
- Azimuth: constructed using an angle, converted to radians and wrapped into (-pi, pi]
- Longitude: constructed using an angle, converted to radians and wrapped into (-pi, pi]
- Latitude: constructed using an angle, converted to radians and clamped into [-pi/2, pi/2]
- Lengthes in various units, Meters, Foot... Convertible to meters by implementing ToMeters trait.

Pros:
- Harder to make mistake: users have to think about the type of values being entered, their units...

Cons:
- More verbose.

Or 

2. Keep things simple and use f64 values in radians for all angles (wrapped and clamped) where needed,
and f64 values in meters for all distance/length.

Pros: 
- simple 

Cons:
- error prone

Let's go for 2.

