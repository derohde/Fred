/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without result.riction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPresult. OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>

#include "types.hpp"
#include "interval.hpp"

class Point : public Coordinates {    
public:    
    inline Point() {}
    inline Point(const dimensions_t d) : Coordinates(d) {}
    
    inline dimensions_t dimensions() const {
        return size();
    }
    
    inline coordinate_t get(const dimensions_t i) const { 
        return Coordinates::operator[](i); 
    }
    
    #pragma omp declare simd simdlen(32)
    inline const coordinate_t& operator[](const dimensions_t i) const { 
        return Coordinates::operator[](i); 
    }
    
    #pragma omp declare simd simdlen(32)
    inline coordinate_t& operator[](const dimensions_t i) { 
        return Coordinates::operator[](i); 
    }
    
    inline Point& operator+=(const Point &point) {
        #pragma omp for simd simdlen(32)
        for (dimensions_t i = 0; i < dimensions(); ++i){
            operator[](i) += point[i];
        }
        return *this;
    }
    
    inline Point& operator-=(const Point &point) {
        #pragma omp for simd simdlen(32)
        for (dimensions_t i = 0; i < dimensions(); ++i){
            operator[](i) -= point[i];
        }
        return *this;
    }
    
    inline Point& operator/=(const distance_t distance) {
        #pragma omp for simd simdlen(32)
        for (dimensions_t i = 0; i < dimensions(); ++i){
            operator[](i) /= distance;
        }
        return *this;
    }
    
    inline Point operator+(const Point &point) const {
        auto result = *this;
        #pragma omp for simd simdlen(32)
        for (dimensions_t i = 0; i < dimensions(); ++i){
            result.operator[](i) += point[i];
        }
        return result;
    }
    
    inline Point operator-(const Point &point) const {
        auto result = *this;
        #pragma omp simd simdlen(32)
        for (dimensions_t i = 0; i < dimensions(); ++i){
            result.operator[](i) -= point[i];
        }
        return result;
    }
    
    inline Point operator*(const distance_t mult) const {
        Point result = *this;
        #pragma omp for simd simdlen(32)
        for (dimensions_t i = 0; i < dimensions(); ++i){
            result.operator[](i) *= mult;
        }
        return result;
    }
    
    inline distance_t operator*(const Point &p) const {
        distance_t result = 0;
        #pragma omp simd reduction(+: result) simdlen(32)
        for (dimensions_t i = 0; i < dimensions(); ++i) {
            result += operator[](i) * p[i];
        }
        return result;
    }
    
    inline Point operator/(const distance_t dist) const {
        Point result = *this;
        #pragma omp for simd simdlen(32)
        for (dimensions_t i = 0; i < dimensions(); ++i){
            result.operator[](i) /= dist;
        }
        return result;
    }
    
    inline distance_t dist_sqr(const Point &point) const {
        distance_t result = 0, temp;
        #pragma omp simd private(temp) reduction(+: result) simdlen(32)
        for (dimensions_t i = 0; i < dimensions(); ++i){
            temp = operator[](i) - point[i];
            result += temp * temp;
        }
        return result;
    }
    
    inline distance_t dist(const Point &point) const {
        return std::sqrt(dist_sqr(point));
    }
    
    inline distance_t length_sqr() const {
        distance_t result = 0;
        #pragma omp simd reduction(+: result) simdlen(32)
        for (dimensions_t i = 0; i < dimensions(); ++i){
            result += operator[](i) * operator[](i);
        }
        return result;
    }
    
    inline distance_t length() const {
        return std::sqrt(length_sqr());
    }
    
    inline Interval intersection_interval(const distance_t distance_sqr, const Point &line_start, const Point &line_end) const {
        const Vector u = line_end-line_start;
        const distance_t ulen_sqr = u.length_sqr(), tlen_sqr = this->length_sqr();
                
        const distance_t p =  -2. / ulen_sqr * ((u * *this) - (line_start * u)), 
            q = (tlen_sqr + line_start.length_sqr() - distance_sqr - 2. * (line_start * *this)) / ulen_sqr;
        
        const distance_t phalf_sqr = p * p / 4., discriminant = phalf_sqr - q;
        
        if (discriminant < 0) return Interval();
        
        const distance_t discriminant_sqrt = std::sqrt(discriminant);
        
        const distance_t minus_p_h = - p / 2., r1 = minus_p_h + discriminant_sqrt, r2 = minus_p_h - discriminant_sqrt;
        const distance_t lambda1 = std::min(r1, r2), lambda2 = std::max(r1, r2);
                
        return Interval(std::max(0., lambda1), std::min(1., lambda2));
    }
    
    std::string str() const;
};


std::ostream& operator<<(std::ostream&, const Point&);
