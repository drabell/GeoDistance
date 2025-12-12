# GeoDistance | Java Class Lib
#### Contains 4 static Methods to calculate the great circle (orthodromic) distance between two geo-points on Earth specified by their coordinates in decimal format (Latitude, Longitude), e.g.
* John F. Kennedy International Airport (JFK): {40.641766,-73.780968}
* Los Angeles International Airport (LAX): {33.942791,-118.410042}
----------------------------------------------------------------------------
Sample calculated great-circle (orthodromic) distance between two geo-points:
* JFK {40.641766,-73.780968}
* LHR {51.470020,-0.454295}

| Method                      | Distance, km/mi       | Note               |
|:----------------------------|:----------------------|:-------------------|
| Haversine                   | 5540.175419079547 km  | high accuracy      |
| Spherical Law of Cosines    | 5540.175419079548 km  | high accuracy      |
| Inverse Vincenty            | 5555.065686009474 km  | highest accuracy   |
| Spherical Earth Projection  | 5784.908563389233 km  | lower accuracy     |
| Expected value              | ~ 5554.5 km           |                    |
| Haversine                   | 3442.5054053574295 mi | high accuracy      |
| Spherical Law of Cosines    | 3442.5054053574304 mi | high accuracy      |
| Inverse Vincenty            | 3451.7577882724104 mi | highest accuracy   |
| Spherical Earth Projection  | 3594.5755310171303 mi | lower accuracy     |
| Expected value              | ~ 3451.4 mi           |                    |

----------------------------------------------------------------------------
