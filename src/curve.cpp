/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "curve.hpp"

Curve::Curve(const Points& points, dimensions_t dimensions) : 
    points(points), number_dimensions(dimensions) {
    if (points.empty()) { 
        std::cerr << "warning: constructed empty curve" << std::endl;
        return; 
    }

    #if DEBUG
    std::cout << "constructed curve of complexity " << points.size() << std::endl;
    #endif
}

Curve::Curve(const np::ndarray &in): number_dimensions(in.shape(1)) {
    auto dimensions = in.get_nd();
    if (dimensions != 2 or in.get_dtype() != np::dtype::get_builtin<coordinate_t>()){
        std::cerr << "A Polygonal_Curve requires an 2-dimensional numpy array of type double."<< std::endl;
        std::cerr << "Current dimensiont: " << dimensions << std::endl;
        return;
    }
    auto number_points = in.shape(0);
    auto point_size = in.shape(1);
    
    #if DEBUG
    std::cout << "constructing curve of size " << number_points << " and " << point_size << " dimensions" << std::endl;
    #endif
    
    auto strides0 = in.strides(0) / sizeof(coordinate_t);
    auto strides1 = in.strides(1) / sizeof(coordinate_t);
    
    auto data = reinterpret_cast<const coordinate_t*>(in.get_data());
    
    points = Points(number_points);
    
    for (index_type i = 0; i < number_points; ++i, data += strides0) {
        points[i] = Point(std::vector<coordinate_t> (point_size));
        
        auto coord_data = data;
        
        for(index_type j = 0; j < point_size; ++j, coord_data += strides1){
          points[i].getCoordinates()[j] = *coord_data;
        }
    }
    
    if (points.empty()) { 
        std::cerr << "warning: constructed empty curve" << std::endl;
    return; 
    }
}

void Curve::push_back(const Point &point) {
    points.push_back(point);
}

std::ostream& operator<<(std::ostream &out, const Curve &curve)
{
    out << "[";
    for (const auto &point: curve) {
        out << point << ", ";
    }
    out << "]" << std::endl;

    return out;
}

