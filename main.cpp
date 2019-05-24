#include <iostream>
#include <vector>
#include "poisson.hpp"
#include <algorithm>
#include <cmath>
#include <utility>
#include <functional>

int main(int argc, char** argv) {

    int ur_x = argc > 1 ? atoi(argv[1]) : 1000;
    int ur_y = argc > 2 ? atoi(argv[2]) : 1000;
    double radius = argc > 3 ? (double)atoi(argv[3]) : 5.0;
    int seed = argc > 4 ? atoi(argv[4]) : 0;
    int k_tries = argc > 5 ? atoi(argv[5]) : 30;

    std::vector<dv2> points = poisson_disk_sampling(dv2(0, 0), dv2(ur_x, ur_y), seed, radius, k_tries);
    std::vector<double> xs;
    std::vector<double> ys;


    for(size_t i = 0; i < points.size(); i++) {
        xs.push_back(points[i].x);
        ys.push_back(points[i].y);
    }

    double shift_x = *std::max_element(xs.begin(), xs.end())/2;
    double shift_y = *std::max_element(ys.begin(), ys.end())/2;

    for(size_t i = 0; i < points.size(); i++) {
        points[i].x -= shift_x;
        points[i].y -= shift_y;
    }

    const double alpha = 250;
    const double beta = 1;

    std::function<std::pair<double, double>(double, double)> gauss = [&alpha, &beta](double x, double y) -> std::pair<double, double> {
        double factor = 1.0 - exp(-((x * x) / (2.0 * alpha * alpha) + (y * y) / (2.0 * alpha * alpha)));
        factor = pow(factor, 1.0 / beta);
        return {x * factor, y * factor};
    };

    for(size_t i = 0; i < points.size(); i++) {
        auto [nx, ny] = gauss(points[i].x, points[i].y);
        points[i].x = nx;
        points[i].y = ny;
    }

    const double crop = 0.75;

    std::vector<dv2> final_points;

    for(size_t i = 0; i < points.size(); i++) {
        if(fabs(points[i].x) <= crop * shift_x && fabs(points[i].y) <= crop * shift_y) {
            final_points.push_back(points[i]);
        }
    }

    shift_x = std::min_element(final_points.begin(), final_points.end(), [](const dv2& a, const dv2& b) { return a.x < b.x; })->x;
    shift_y = std::min_element(final_points.begin(), final_points.end(), [](const dv2& a, const dv2& b) { return a.y < b.y; })->y;

    for(size_t i = 0; i < final_points.size(); i++) {
        final_points[i].x -= shift_x;
        final_points[i].y -= shift_y;
    }

    shift_x = std::max_element(final_points.begin(), final_points.end(), [](const dv2& a, const dv2& b) { return a.x < b.x; })->x;
    shift_y = std::max_element(final_points.begin(), final_points.end(), [](const dv2& a, const dv2& b) { return a.y < b.y; })->y;

    for(size_t i = 0; i < final_points.size(); i++) {
        final_points[i].x /= shift_x;
        final_points[i].y /= shift_y;
    }

    for(size_t i = 0; i < final_points.size(); i++) {
        std::cout << final_points[i].x << " " << final_points[i].y << std::endl;
    }
    return 0;
}
