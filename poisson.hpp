#pragma once
#include "Vec2.hpp"
#include <vector>
#include <queue>
#include <random>
#include <cstdint>

using dv2 = Vec2<double>;
using i64 = int64_t;
using u32 = uint32_t;

std::vector<dv2> poisson_disk_sampling(dv2 ll, dv2 ur, int seed = 0, double radius = 5.0, u32 k = 30) {
    std::vector<dv2> points;
    double edge = radius/std::sqrt(2);
    i64 width = (ur.x - ll.x) / edge + 1;
    i64 height = (ur.y - ll.y) / edge + 1;
    std::vector<std::vector<i64>> grid(height, std::vector<i64>(width, -1));
    std::mt19937_64 rand(seed);
    std::uniform_real_distribution<double> gen_x(0, width - 1);
    std::uniform_real_distribution<double> gen_y(0, height - 1);
    std::uniform_real_distribution<double> gen_radius(radius, 2.0 * radius);
    std::uniform_real_distribution<double> gen_angle(0.0, 2 * M_PI);
    auto gen_point = [&rand, &gen_x, &gen_y]() {
        return dv2(gen_x(rand), gen_y(rand));
    };

    auto point2grid = [&ll, &ur, &edge](dv2 p) -> std::pair<i64, i64>{
        // centering
        p.x -= ll.x;
        p.y -= ll.y;
        // converting
        p.x /= edge;
        p.y /= edge;

        return {p.y, p.x};
    };

    std::queue<dv2> active_list;
    auto x0 = gen_point();
    auto [sy0, sx0] = point2grid(x0);
    grid[sy0][sx0] = 0;
    points.push_back(x0);
    active_list.push(x0);

    while(!active_list.empty()) {
        auto p = active_list.front();
        active_list.pop();
        bool inserted = false;
        for(u32 i = 0; i < k; i++) {
            double radius = gen_radius(rand);
            double theta = gen_angle(rand);
            double candidate_x = radius * std::cos(theta);
            double candidate_y = radius * std::sin(theta);
            dv2 new_point(p.x + candidate_x, p.y + candidate_y);
            // if it's outside the domain, we are not interested
            if(new_point.x < ll.x || new_point.x >= ur.x || new_point.y < ll.y || new_point.y >= ur.y) continue;
            auto [grid_y, grid_x] = point2grid(new_point);
            // if the point it's outside the grid we are not interested
            if(grid_x < 0 || grid_x >= width || grid_y < 0 || grid_y >= height) continue;
            // it there is already a point in that location of the grid we are not interested
            if(grid[grid_y][grid_x] >= 0) continue;
            // at this point we can insert the candidate if there are no neighbours within the given radius
            // if a neighbour is encountered, we skip all the loops setting the flag to false
            bool flag = true;
            for(i64 y = std::max(0L, grid_y - 1); y <= std::min(grid_y + 1, height - 1) && flag; y++)
            for(i64 x = std::max(0L, grid_x - 1); x <= std::min(grid_x + 1, width - 1) && flag; x++) {
                // if there is a point, than the grid will contain a semipositive index instead of -1
                if(grid[y][x] >= 0) {
                    // if a cell contain a point, it's required to check the radius constraint
                    // in order to decide if we can insert the point at that location or another
                    // the check for the query cell [grid_x, grid_y] it's executed above
                    auto query_point = points[grid[y][x]];
                    if(std::hypot(query_point.x - new_point.x, query_point.y - new_point.y) < radius) {
                        flag = false;
                    }
                }
            }
            // if the query point is closer than radius to another point, then we skip all the loops
            // and we are not interested, so we jump to the next iteration (if there is) like we've done above
            if(!flag) continue;
            grid[grid_y][grid_x] = points.size();
            points.push_back(new_point);
            active_list.push(new_point);
            // if we have inserted a point, p has to be reinserted in the queue because there may be
            // some other space that can give a sample
            inserted = true;
        }
        // after k iterations, if at least one candidate is found, then we have to reinsert p
        if(inserted) active_list.push(p);
    }

    return points;
}
