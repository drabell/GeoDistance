# GeoDistance | Java Class Lib
#### Contains 4 static Methods to calculate the great circle (orthodromic) distance between two geo-points on Earth specified by their coordinates in decimal format (Latitude, Longitude), e.g.
* John F. Kennedy International Airport, USA (JFK): {40.641766,-73.780968}
* London Heathrow Airport, UK (LHR):	{51.470020	-0.454295}
----------------------------------------------------------------------------
Sample calculated great-circle (orthodromic) distance between two geo-points:
* JFK {40.641766,-73.780968}
* LHR {51.470020,-0.454295}
 
| Method                      | Distance, km    | Distance, mi    | Note
|:----------------------------|:----------------|:----------------|:------------------|
| Haversine                   | 5,540.175419    | 3,442.505405    | high accuracy     |
| Spherical Law of Cosines    | 5,540.175419    | 3,442.505405    | high accuracy     |
| Spherical Earth Projection  | 5,784.908563    | 3,594.575531    | lower accuracy    |
| Inverse Vincenty            | 5,555.065686    | 3,451.757788    | highest accuracy  |
| Expected value              | ~ 5554.5        | ~ 3451.4        | approx.           |

#### Theory
The calculation of the great-circle (orthodromic) distance between two geographic points on the Earth's surface is one of the core problems in Geographic Information Systems (GIS). Although this may appear to be a trivial geodesic task, it requires a non-trivial algorithmic implementation. If the problem were confined to plane geometry, the Pythagorean theorem would provide a simple solution. However, actual GIS computations involve three-dimensional models—specifically spherical or ellipsoidal representations of the Earth—which demand more elaborate approaches.

The spherical model introduces a systemic error, with a typical margin of about 0.1% to 0.5%. This level of accuracy is generally acceptable for many commercial-grade, general-purpose applications. In contrast, the more accurate ellipsoidal model of the Earth (as implemented, for example, in the inverse Vincenty algorithm) theoretically reduces the error margin to fractions of a millimeter, though at the cost of significantly increased computational complexity.
##### Systemic Error Margin
The Earth is an oblate spheroid, flattened at the poles and bulging at the equator. Its equatorial radius is approximately 6,378 km, while its polar radius is about 6,357 km. This difference of roughly 21 km means that the "true" radius varies with latitude.

The spherical Earth model is an approximation, and its primary systemic error arises from ignoring the Earth's oblateness. The error margin in orthodromic distance calculations using the spherical model depends on both the computational algorithm employed and the fixed radius assumed. For the commonly used mean Earth radius of 6,370 km, the maximum relative error is typically less than 0.56%, corresponding to about 21 km for very long distances, depending on the specific points and formula applied.
####  Spherical Earth Math Model/Algorithms
* Haversine
* Spherical Law of Cosines
* Spherical Earth Projection
####  Ellipsoidal Earth Math Model/Algorithm
* Inverse Vincenty formula
*** 
