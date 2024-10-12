#pragma once
#include <string>
#include <iostream>

class CelestialObject {
public:
    // Member variables
    double cx, cy, cz; // [x, y, z] coordinates (meters)
    double vx, vy, vz; // [x, y, z] velocity (meters per second)
    double mass; // mass (kilograms)
    double radius; // radius (meters)
    std::string name; // name of the object
    // Constructor
    constexpr CelestialObject(double x, double y, double z, double vx, double vy, double vz, double m, double r,
                              std::string n) : cx(x), cy(y), cz(z), vx(vx), vy(vy), vz(vz), mass(m), radius(r),
                                               name(n) {
    }

    // Apply acceleration
    constexpr void applyAcceleration(double v1, double v2, double v3, double time) {
        vx += v1 / time;
        vy += v2 / time;
        vz += v3 / time;
    }

    // Apply velocity to update coordinates
    constexpr void applyVelocity(double time) {
        cx += vx / time;
        cy += vy / time;
        cz += vz / time;
    }
};

std::ostream &operator<<(std::ostream &lhs, const CelestialObject &sol) {
    lhs << "CelestialObject: " << sol.name << "("
            << sol.cx << ", " << sol.cy << ", " << sol.cz << ", "
            << sol.vx << ", " << sol.vy << ", " << sol.vz << ", "
            << sol.mass << ", " << sol.radius << ", " << sol.name << ")";
    return lhs; // Return the modified stream to allow chaining
}
