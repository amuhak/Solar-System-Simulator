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
    CelestialObject(double x, double y, double z, double vx, double vy, double vz, double m, double r,
                    std::string n) : cx(x), cy(y), cz(z), vx(vx), vy(vy), vz(vz), mass(m), radius(r), name(n) {
    }

    // Apply acceleration
    void applyAcceleration(double v1, double v2, double v3, double time) {
        vx += v1 / time;
        vy += v2 / time;
        vz += v3 / time;
    }

    // Apply velocity to update coordinates
    void applyVelocity(double time) {
        cx += vx / time;
        cy += vy / time;
        cz += vz / time;
    }
};
