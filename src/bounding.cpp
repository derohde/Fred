/*
Copyright 2023 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "bounding.hpp"

std::pair<Point, distance_t> bounding_sphere(const Points &points) {
    if (points.size() < 1)
            return std::make_pair(Point(0), std::numeric_limits<distance_t>::infinity());

    const auto &x = points[0];
    
    if (points.size() < 2)
        return std::make_pair(x, distance_t(0));
    
    distance_t max_dist = distance_t(0), tdist, radius;
    
    Point y(x.dimensions()), z(x.dimensions()), center(x.dimensions());
    
    for (const auto &point : points) {
        tdist = x.dist_sqr(point);
        if (tdist > max_dist) {
            y = point;
            max_dist = tdist;
        }
    }
    
    if (points.size() < 3) 
        return std::make_pair((x + y) / coordinate_t(2), x.dist(y) / distance_t(2));
    
    max_dist = distance_t(0);
    
    for (const auto &point : points) {
        tdist = y.dist_sqr(point);
        if (tdist > max_dist) {
            z = point;
            max_dist = tdist;
        }
    }
    
    center = (y + z) / coordinate_t(2);
    radius = y.dist(z) / distance_t(2);
    
    if (points.size() > 3) {
        for (const auto &point : points) {
            tdist = center.dist(point);
            if (tdist > radius) {
                radius = (radius + tdist) / distance_t(2);
                center = (center * radius + point * (tdist - radius)) / tdist;
            }
        }
    }
    
    return std::make_pair(center, radius);

}
