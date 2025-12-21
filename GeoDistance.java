/******************************************************************************
 Module         : GeoDistance.java |Class Lib
 Description    : Methods to calculate the great circle (orthodromic) distance
                : between two geo-points on Earth
 Version        : 21.1.001
 ------------------------------------------------------------------------------
 Copyright      :  2011-2025 Alexander Bell
-------------------------------------------------------------------------------
 DISCLAIMER     : This Module is provided on AS IS basis without any warranty.
                : The user assumes the entire risk as to the accuracy and the
                : use of this module. In no event shall the author be liable
                : for any damages arising out of the use of or inability
                : to use this module.
 TERMS OF USE   : This module is copyrighted.
                : Please keep the Copyright notice intact.
 ****************************************************************************/
import java.lang.Math;

/// <summary>
/// Class GeoDistance contains four static methods to calculate the
/// great-circle (orthodromic) distance between two geo-points on Earth
/// specified by coordinates in decimal format (Latitude, Longitude), e.g.
/// John F. Kennedy International Airport (JFK): {40.641766,-73.780968},
/// Los Angeles International Airport (LAX): {33.942791,-118.410042}
/// Sample output:
/// =======================================================================
/// Great-circle (orthodromic) distance between two geo-points:
/// JFK {40.641766,-73.780968} to LHR {51.470020,-0.454295}
/// km --------------------------------------------------------------------
/// Haversine					: 5540.175419079547 (high accuracy)
/// Spherical Law of Cosines	: 5540.175419079548 (high accuracy)
/// Inverse Vincenty			: 5555.065686009474 (highest accuracy)
/// Spherical Earth Projection	: 5784.908563389233 (lower accuracy)
/// Expected value              :~5554.5 km
/// miles -----------------------------------------------------------------
/// Haversine					: 3442.5054053574295 (high accuracy)
/// Spherical Law of Cosines	: 3442.5054053574304 (high accuracy)
/// Inverse Vincenty			: 3451.7577882724104 (highest accuracy)
/// Spherical Earth Projection	: 3594.5755310171303 (lower accuracy)
/// Expected value              :~3451.4 miles
/// =======================================================================
/// </summary>
public class GeoDistance {
    //region public enum
    public enum Units { SI, US } // SI: km, US: miles
    //endregion

    //region private const
    // Earth mean radius, km
    private static final double meanRadius = 6371.009;
    // Conversion factor mile to km
    private static final double mi2km = 1.609344;
    //endregion

    //region Haversine algorithm **********************************************************
    /// <summary>
    /// Haversine algorithm implemented in this method enables
    /// the high-accuracy geodesic calculation of the great-circle
    /// (aka orthodromic) distance (km/miles) between two geographic
    /// points on the Earth's surface.
    /// </summary>
    /// <param name="Lat1">double: 1st point Latitude</param>
    /// <param name="Lon1">double: 1st point Longitude</param>
    /// <param name="Lat2">double: 2nd point Latitude</param>
    /// <param name="Lon2">double: 2nd point Longitude</param>
    /// <returns>double: distance, km/miles</returns>
    public static double Haversine(double Lat1, double Lon1,
                                   double Lat2, double Lon2,
                                   Units Unit) {
        try  {
            // convert coordinates' latitude to radians
            double φ1 = Math.toRadians(Lat1);
            double φ2 = Math.toRadians(Lat2);

             double _a = Math.sin((φ2 - φ1)/2);
            _a *= _a; // square

            double _b = Math.sin(Math.toRadians((Lon2 - Lon1)/2));
            _b *= _b * Math.cos(φ1) * Math.cos(φ2);

            // central angle θ, aka arc segment angular distance
            double θ = 2* Math.asin(Math.sqrt(_a + _b));

            // great-circle (orthodromic) distance, km/miles
            return θ * (Unit == Units.SI? 1:1/ mi2km) * meanRadius;
        }
        catch( Exception e) { return -1; } // indicates error
    }
    //endregion

    //region Spherical Law of Cosines *****************************************************
    /// <summary>
    /// Spherical Law of Cosines (SLC) algorithm implemented in this
    /// method enables the high-accuracy geodesic calculation of the
    /// great-circle (aka orthodromic) distance (km/miles) between
    /// two geographic points on the Earth's surface.
    /// Note: results are very close to the Haversine formula, which is
    /// generally preferred for numerical stability with small distances.
    /// </summary>
    /// <param name="Lat1">double: 1st point Latitude</param>
    /// <param name="Lon1">double: 1st point Longitude</param>
    /// <param name="Lat2">double: 2nd point Latitude</param>
    /// <param name="Lon2">double: 2nd point Longitude</param>
    /// <returns>double: distance, km/miles</returns>
    public static double SLC(double Lat1, double Lon1,
                             double Lat2, double Lon2,
                             Units Unit) {
        try  {
            // convert coordinates to radians
            double φ1 = Math.toRadians(Lat1); // Lat1;
            double φ2 = Math.toRadians(Lat2); // Lat2;
            double Δλ = Math.toRadians(Lon1-Lon2); // delta Lon;

            // central angle θ, aka arc segment angular distance
            double θ = Math.acos(Math.sin(φ1) * Math.sin(φ2) +
                    Math.cos(φ1) * Math.cos(φ2) * Math.cos(Δλ));

            // great-circle (orthodromic) distance, km/miles
            return θ * (Unit == Units.SI? 1:1/ mi2km) * meanRadius;
        }
        catch( Exception e) { return -1; } // indicates error
    }
    //endregion

    //region Vincenty algorithm (ellipsoid) ***********************************************
    /// <summary>
    /// Inverse Vincenty (ellipsoid) algorithm implemented in this method enables
    /// the very high-accuracy geodesic calculation of the great-circle  (orthodromic)
    /// distance (km/miles) between two geographic points on the Earth's surface.
    /// Notes -----------------------------------------------------------------------------
    /// Inverse Vincenty (ellipsoid) algorithm provides the highest accuracy among
    /// the common spherical/ellipsoidal computational methods, but it is not a closed-form.
    /// This inverse solution (distance and bearings between two points) is an efficient
    /// iterative algorithm with nested expressions well-suited for
    /// the software implementation.
    /// Regarding the accuracy and robustness of the method, see the practical notes below.
    /// - Convergence:
    /// The inverse method can fail near antipodal points.
    /// Use a max-iteration guard and a small epsilon; if it fails, fall back
    /// to a more robust geodesic algorithm.
    /// - Precision:
    /// Double precision is sufficient; avoid premature rounding of inputs.
    /// Keep lat/lon in radians for the loop.
    /// - Model choice:
    /// WGS84 is standard. For a different datum (e.g., GRS80), set 𝑎 and 𝑓 accordingly.
    /// - Outputs:
    /// Besides distance, this method can return initial/final bearings.
    /// - AI vibe coding:
    /// This Inverse Vincenty geodesic algorithm was implemented in AI-assisted
    /// pair programming (vibe coding) interactive session with AI Copilot.
    /// -----------------------------------------------------------------------------------
    /// </summary>
    /// <returns>double: orthodromic distance, km/miles</returns>
    public static double Vincenty(double Lat1, double Lon1,
                                  double Lat2, double Lon2,
                                  Units Unit) {
        // WGS84 constants
        double a = 6378137.0; // Earth equatorial radius, m
        double f = 1.0 / 298.257223563;
        double b = a * (1.0 - f);
        try {
            // Convert to radians
            double φ1 = Math.toRadians(Lat1), φ2 = Math.toRadians(Lat2);
            double Δλ = Math.toRadians(Lon2 - Lon1);

            // Reduced latitudes
            double U1 = Math.atan((1 - f) * Math.tan(φ1));
            double U2 = Math.atan((1 - f) * Math.tan(φ2));

            double sinU1 = Math.sin(U1), cosU1 = Math.cos(U1);
            double sinU2 = Math.sin(U2), cosU2 = Math.cos(U2);

            double λ = Δλ;
            double λPrev;
            double iterLimit = 100;
            double ε = 1e-12;

            double sinσ, cosσ, σ, sinα, cos2α, cos2σm;
            double u2, A, B, Δσ;

            do {
                double sinλ = Math.sin(λ), cosλ = Math.cos(λ);
                double term1 = cosU2 * sinλ;
                double term2 = cosU1 * sinU2 - sinU1 * cosU2 * cosλ;

                sinσ = Math.sqrt(term1 * term1 + term2 * term2);
                if (sinσ == 0.0) return 0.0; // coincident points

                cosσ = sinU1 * sinU2 + cosU1 * cosU2 * cosλ;
                σ = Math.atan2(sinσ, cosσ);

                sinα = (cosU1 * cosU2 * sinλ) / sinσ;
                double sin2α = sinα * sinα;
                cos2α = 1.0 - sin2α;

                if (cos2α != 0.0) cos2σm = cosσ - (2.0 * sinU1 * sinU2) / cos2α;
                else  cos2σm = 0.0; // equatorial line

                u2 = (cos2α * (a * a - b * b)) / (b * b);

                A = 1.0 + (u2 / 16384.0) * (4096.0 + u2 * (-768.0 + u2 *
                        (320.0 - 175.0 * u2)));
                B = (u2 / 1024.0) * (256.0 + u2 * (-128.0 + u2 * (74.0 - 47.0 * u2)));

                double cos2σm2 = cos2σm * cos2σm;
                Δσ = B * sinσ * (cos2σm + (B / 4.0) * (cosσ * (-1.0 + 2.0 * cos2σm2)
                        - (B / 6.0) * cos2σm * (-3.0 + 4.0 * sinσ * sinσ) *
                        (-3.0 + 4.0 * cos2σm2)));

                double C = (f / 16.0) * cos2α * (4.0 + f * (4.0 - 3.0 * cos2α));

                λPrev = λ;
                λ = Δλ + (1.0 - C) * f * sinα *
                        (σ + C * sinσ * (cos2σm + C * cosσ * (-1.0 + 2.0 * cos2σm2)));

                if (Math.abs(λ - λPrev) < ε) break;
            } while (--iterLimit > 0);

            // If not converged, you may fall back to a robust solver
            if (iterLimit == 0) throw new ArithmeticException("No Convergence");

            double s = b * A * (σ - Δσ);

            // Optional: initial α1 and final α2 bearings calculation
            double α1 = Math.atan2(cosU2 * Math.sin(λ),
                    cosU1 * sinU2 - sinU1 * cosU2 * Math.cos(λ));
            double α2 = Math.atan2(cosU1 * Math.sin(λ),
                    -sinU1 * cosU2 + cosU1 * sinU2 * Math.cos(λ));

            // orthodromic distance, km/miles
            return s * (Unit == Units.SI? 1:1/ mi2km)/1000;
        }
        catch( Exception e) { return -1; } //indicates error
    }
    //endregion

    //region Spherical Earth Projection (SEP) *********************************************
    /// <summary>
    /// Spherical Earth Projection (SEP) to a plane formula
    /// implemented in this method enables the calculation
    /// of a great-circle (orthodromic) distance(km/miles) between two
    /// geographic points on the Earth using Pythagorean Theorem:
    /// Central Angle = Sqrt((φ2 - φ1)^2 + (Cos((φ1 + φ2)/2) * (Lon2 - Lon1))^2)
    /// Note: this is a relatively low accuracy computation approach
    /// suitable for small distances (e.g., within a city or small region);
    /// it is shown mostly for a didactic purpose. For higher accuracy over
    /// longer distances, use either Haversine, or Spherical Law of Cosines,
    /// or Inverse Vincenty methods (the latter provides the highest accuracy).
    /// </summary>
    /// <param name="Lat1">double: 1st point Latitude</param>
    /// <param name="Lon1">double: 1st point Longitude</param>
    /// <param name="Lat2">double: 2nd point Latitude</param>
    /// <param name="Lon2">double: 2nd point Longitude</param>
    /// <returns>double: distance, km/miles</returns>
    public static double SEP(double Lat1, double Lon1,
                             double Lat2, double Lon2,
                             Units UnitSys){
        try {
            // convert to radians
            double φ1 = Math.toRadians(Lat1);
            double φ2 = Math.toRadians(Lat2);
            double Δφ = Math.toRadians(Lat2 - Lat1);

            double _a = Math.toRadians((Lon2 - Lon1) * Math.cos((φ1 + φ2) / 2));

            // central angle θ (arc segment angular distance)
            double θ = Math.sqrt(_a * _a + Δφ * Δφ);

            // great-circle (orthodromic) distance, km/miles
            return θ * (UnitSys == Units.SI ? 1 : 1 / mi2km) * meanRadius;
        }
        catch( Exception e) { return -1; } // indicates error
    }
    //endregion
}
