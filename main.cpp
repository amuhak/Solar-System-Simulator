#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <vector>

#include "celestial_objects.h"

constexpr uint64_t starting_time_UNIX = 1704067200;
constexpr uint64_t ending_time_UNIX = 2020291200;
constexpr long double ticks_per_second = 1. / 1600;
constexpr uint64_t ticks = (ending_time_UNIX - starting_time_UNIX) * ticks_per_second;
std::vector<std::string> eclipse_dates;
long double current_time = starting_time_UNIX;
std::vector<std::pair<int, int> > celestial_object_pairs;


std::pair<vec3<long double>, vec3<long double> > calculate_acceleration(CelestialObject &body_one,
                                                                        CelestialObject &body_two) {
    vec3<long double> displacement_vector = body_two.coordinates - body_one.coordinates;
    long double distance = displacement_vector.norm();
    long double grav_force = 6.67430E-11 * body_one.mass * body_two.mass / (distance * distance);
    vec3<long double> force_vector = displacement_vector * (grav_force / distance);
    vec3<long double> body_one_acceleration = force_vector / body_one.mass;
    vec3<long double> body_two_acceleration = force_vector * -1 / body_two.mass;
    return {body_one_acceleration, body_two_acceleration};
}


std::tuple<int, vec3<long double>, int, vec3<long double> > calculate_acceleration_for_pair(int a, int b) {
    CelestialObject body_one = celestial_objects[a];
    CelestialObject body_two = celestial_objects[b];
    auto [body_one_acceleration, body_two_acceleration] = calculate_acceleration(body_one, body_two);
    return {a, body_one_acceleration, b, body_two_acceleration};
}

void iterate_simulation(long double tps) {
    vec3<long double> direction_vector = celestial_objects[3].coordinates - celestial_objects[0].coordinates;
    direction_vector = direction_vector / direction_vector.norm();
    int alignment_threshold = 1000;
    vec3<long double> sol_to_luna_vector = celestial_objects[4].coordinates - celestial_objects[0].coordinates;
    long double perpendicular_distance = sol_to_luna_vector.cross(direction_vector).norm();
    bool is_aligned = perpendicular_distance <= alignment_threshold + celestial_objects[4].radius +
                      celestial_objects[3].radius;

    if (is_aligned) {
        std::string eclipse_type = "Solar";
        if (distance_calc(celestial_objects[0].coordinates, celestial_objects[4].coordinates) >
            distance_calc(celestial_objects[0].coordinates, celestial_objects[3].coordinates)) {
            eclipse_type = "Lunar";
        }
        eclipse_dates.push_back(eclipse_type + " Eclipse on: " + std::to_string(current_time));
    }
    std::vector<vec3<long double> > body_cumulative_acceleration(celestial_objects.size());
    std::vector<std::tuple<int, vec3<long double>, int, vec3<long double> > > acceleration_results(
        celestial_object_pairs.size());
    for (int i = 0; i < celestial_object_pairs.size(); ++i) {
        acceleration_results[i] = calculate_acceleration_for_pair(celestial_object_pairs[i].first,
                                                                  celestial_object_pairs[i].second);
    }
    for (auto &[body_one_index, body_one_accel, body_two_index, body_two_acc]: acceleration_results) {
        body_cumulative_acceleration[body_one_index] = body_cumulative_acceleration[body_one_index] + body_one_accel;
        body_cumulative_acceleration[body_two_index] = body_cumulative_acceleration[body_two_index] + body_two_acc;
    }

    for (int i = 0; i < celestial_objects.size(); ++i) {
        celestial_objects[i].applyAcceleration(body_cumulative_acceleration[i], tps);
        celestial_objects[i].applyVelocity(tps);
    }

    long double total_mass = 0;
    vec3<long double> barycenter = {0, 0, 0};
    for (auto &obj: celestial_objects) {
        long double mass = obj.mass;
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
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < celestial_objects.size(); ++i) {
        for (int j = i + 1; j < celestial_objects.size(); ++j) {
            celestial_object_pairs.emplace_back(i, j);
        }
    }

    std::cout << std::setprecision(20);

    for (int count = 0; count < ticks; ++count) {
        iterate_simulation(ticks_per_second);
        long double percent_finished = (count + 1.) / ticks * 100;
        std::cout << "\rProgress: " << percent_finished << "%";
        std::cout << std::flush;
        current_time += 1. / ticks_per_second;
    }

    for (auto &i: eclipse_dates) {
        std::cout << i << std::endl;
    }

    std::cout << "Final UNIX timestamp: " << current_time << std::endl;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "Elapsed time: " << elapsed_time.count() << "ms" << std::endl;

    std::ranges::sort(eclipse_dates);
    std::ofstream myfile;
    myfile.open("eclipses.txt");
    for (auto &i: eclipse_dates) {
        myfile << i << std::endl;
    }
}
