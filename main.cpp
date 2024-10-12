#include <algorithm>
#include <cstdio>
#include <fstream>
#include <span>
#include <vector>

#include "celestial_objects.h"

constexpr u_int32_t starting_time_UNIX = 1704067200;
constexpr u_int32_t ending_time_UNIX = 2020291200;
constexpr double ticks_per_second = 1. / 1600;
constexpr u_int32_t ticks = (ending_time_UNIX - starting_time_UNIX) * ticks_per_second;
std::vector<std::string> eclipse_dates;
double current_time = starting_time_UNIX;
std::vector<std::pair<int, int> > celestial_object_pairs;


std::pair<vec3<double>, vec3<double> > calculate_acceleration(CelestialObject &body_one, CelestialObject &body_two) {
    vec3<double> direction_vector = body_two.coordinates - body_one.coordinates;
    double distance = direction_vector.norm();
    double force_magnitude = 6.67430E-11 * body_one.mass * body_two.mass / (distance * distance);
    vec3<double> force_vector = direction_vector / distance * force_magnitude;
    vec3<double> body_one_acceleration = force_vector / body_one.mass;
    vec3<double> body_two_acceleration = (force_vector * -1) / body_two.mass;
    return {body_one_acceleration, body_two_acceleration};
}


std::tuple<int, vec3<double>, int, vec3<double> > calculate_acceleration_for_pair(int a, int b) {
    CelestialObject body_one = celestial_objects[a];
    CelestialObject body_two = celestial_objects[b];
    auto [body_one_acceleration, body_two_acceleration] = calculate_acceleration(body_one, body_two);
    return {a, body_one_acceleration, b, body_two_acceleration};
}

void iterate_simulation(std::vector<CelestialObject> &celestial_objects, double tps) {
    vec3<double> direction_vector = celestial_objects[3].coordinates - celestial_objects[0].coordinates;
    direction_vector /= direction_vector.norm();
    int alignment_threshold = 1000;
    vec3<double> sol_to_luna_vector = celestial_objects[4].coordinates - celestial_objects[0].coordinates;
    double perpendicular_distance = (sol_to_luna_vector.cross(direction_vector)).norm();
    bool is_aligned = perpendicular_distance <= alignment_threshold + celestial_objects[4].radius + celestial_objects[3]
                      .radius;

    if (is_aligned) {
        std::string eclipse_type = "Solar";
        if ((celestial_objects[0].coordinates - celestial_objects[4].coordinates).norm() >
            (celestial_objects[0].coordinates - celestial_objects[3].coordinates).norm()) {
            eclipse_type = "Lunar";
        }
        eclipse_dates.emplace_back(eclipse_type + " Eclipse on: " + std::to_string(current_time));
    }
    std::vector<vec3<double> > body_cumulative_acceleration(celestial_objects.size(), {0, 0, 0});
    std::vector<std::tuple<int, vec3<double>, int, vec3<double> > > body_accelerations;
    for (auto &pair: celestial_object_pairs) {
        body_accelerations.push_back(calculate_acceleration_for_pair(pair.first, pair.second));
    }

    for (auto &[body_one_index, body_one_accel, body_two_index, body_two_acc]: body_accelerations) {
        body_cumulative_acceleration[body_one_index] = body_cumulative_acceleration[body_one_index] + body_one_accel;
        body_cumulative_acceleration[body_two_index] = body_cumulative_acceleration[body_two_index] + body_two_acc;
    }

    for (int i = 0; i < celestial_objects.size(); ++i) {
        celestial_objects[i].applyAcceleration(body_cumulative_acceleration[i], tps);
        celestial_objects[i].applyVelocity(tps);
    }

    double total_mass = 0;
    vec3<double> barycenter = {0, 0, 0};
    for (auto &obj: celestial_objects) {
        double mass = obj.mass;
        total_mass += mass;
        barycenter = barycenter + obj.coordinates * mass;
    }
    barycenter = barycenter / total_mass;
    for (auto &obj: celestial_objects) {
        obj.coordinates = (obj.coordinates - barycenter).round(4);
        obj.velocity = obj.velocity.round(4);
    }
}

int main() {
    for (int i = 0; i < celestial_objects.size(); ++i) {
        for (int j = i + 1; j < celestial_objects.size(); ++j) {
            celestial_object_pairs.emplace_back(i, j);
        }
    }

        for (int count = 0; count < ticks; ++count) {
        iterate_simulation(celestial_objects, ticks_per_second);
        double percent_finished = (count + 1.) / ticks * 100;
        std::cout << "\rProgress: " << percent_finished << "%";
        std::cout << std::flush;
        current_time += 1. / ticks_per_second;
    }

    std::cout << "Final UNIX timestamp: " << current_time << std::endl;
    std::sort(eclipse_dates.begin(), eclipse_dates.end());
    std::ofstream myfile;
    myfile.open("eclipses.txt");
    for (auto &i: eclipse_dates) {
        myfile << i << std::endl;
    }
}
