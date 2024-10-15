#pragma once
#include <cmath>
#include <string>
#include <iostream>

template<typename T>
class vec3 {
public:
    T x, y, z;

    T operator [](int i) const {
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }

    vec3<T> operator +(const vec3<T> &v) const {
        return {x + v.x, y + v.y, z + v.z};
    }

    vec3<T> operator +(const T &v) const {
        return {x + v, y + v, z + v};
    }

    vec3<T> operator -(const vec3<T> &v) const {
        return {x - v.x, y - v.y, z - v.z};
    }

    vec3<T> operator -(const T &v) const {
        return {x - v, y - v, z - v};
    }

    vec3<T> operator *(const vec3<T> &v) const {
        return {x * v.x, y * v.y, z * v.z};
    }

    vec3<T> operator *(const T &v) const {
        return {x * v, y * v, z * v};
    }

    vec3<T> operator /(const vec3<T> &v) const {
        return {x / v.x, y / v.y, z / v.z};
    }

    vec3<T> operator /(const T &v) const {
        return {x / v, y / v, z / v};
    }

    T norm() {
        return std::sqrt(std::abs(x * x + y * y + z * z));
    }

    vec3<T> cross(const vec3<T> &vec3) {
        return {y * vec3.z - z * vec3.y, z * vec3.x - x * vec3.z, x * vec3.y - y * vec3.x};
    }

    vec3<T> round(int precision) {
        return {
            std::round(x * std::pow(10, precision)) / std::pow(10, precision),
            std::round(y * std::pow(10, precision)) / std::pow(10, precision),
            std::round(z * std::pow(10, precision)) / std::pow(10, precision)
        };
    }
};


class CelestialObject {
public:
    // Member variables
    vec3<long double> coordinates; // [x, y, z] coordinates (meters)
    vec3<long double> velocity; // [x, y, z] velocity (meters per second)
    long double mass; // mass (kilograms)
    long double radius; // radius (meters)
    std::string name; // name of the object
    // Constructor
    CelestialObject(long double x, long double y, long double z, long double vx, long double vy, long double vz,
                    long double m, long double r,
                    std::string n) : mass(m), radius(r), name(n) {
        coordinates = {x, y, z};
        velocity = {vx, vy, vz};
    }

    // Apply acceleration
    void applyAcceleration(vec3<long double> otherV, long double time) {
        velocity = velocity + otherV / time;
    }

    // Apply velocity to update coordinates
    void applyVelocity(long double time) {
        coordinates = coordinates + velocity / time;
    }
};

std::ostream &operator<<(std::ostream &lhs, const CelestialObject &sol) {
    lhs << "CelestialObject: " << sol.name << "("
            << sol.coordinates.x << ", " << sol.coordinates.y << ", " << sol.coordinates.z << ", "
            << sol.velocity.x << ", " << sol.velocity.y << ", " << sol.velocity.z << ", "
            << sol.mass << ", " << sol.radius << ", " << sol.name << ")";
    return lhs; // Return the modified stream to allow chaining
}

template<typename T>
T distance_calc(vec3<T> pos_one, vec3<T> pos_two) {
    return std::sqrt(
        std::pow(pos_two.x - pos_one.x, 2) +
        std::pow(pos_two.y - pos_one.y, 2) +
        std::pow(pos_two.z - pos_one.z, 2)
    );
}
