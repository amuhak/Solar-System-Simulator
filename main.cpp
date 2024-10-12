#include <algorithm>
#include <cstdio>
#include <fstream>
#include <vector>

#include "celestial_objects.h"

constexpr u_int32_t starting_time_UNIX = 1704067200;
constexpr u_int32_t ending_time_UNIX = 2020291200;
constexpr double ticks_per_second = 1. / 1600;
constexpr u_int32_t ticks = (ending_time_UNIX - starting_time_UNIX) * ticks_per_second;
std::vector<std::string> eclipse_dates;

int main() {
    double current_time = starting_time_UNIX;
    std::cout << sol << std::endl;
    std::vector<std::pair<int, int> > celestial_object_pairs;
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
