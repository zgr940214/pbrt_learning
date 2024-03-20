#include "core/geometry.hpp"
#include <stdexcept>

namespace eric
{
/*
** Vector2
*/
template<typename Ty>
Ty Vector2<Ty>::dot(const Vector2<Ty> &v) const {
    return x * v.x + y * v.y;
};

template<typename Ty>
Vector2<Ty>& Vector2<Ty>::abs() {
    x = std::abs(x);
    y = std::abs(y);
    return *this;
};

template<typename Ty>
Ty Vector2<Ty>::length() const {
    return std::sqrt(length_squared());
};

template<typename Ty>
Ty Vector2<Ty>::length_squared() const {
    return x*x + y*y;
};

template<typename Ty>
Vector2<Ty>& Vector2<Ty>::normalize() {
    *this /= length();
};

template<typename Ty>
Vector2<Ty> Vector2<Ty>::normalized() const {
    Vector2<Ty> result = *this;
    return result.normalize();
};

template<typename Ty>
bool Vector2<Ty>::has_nans() const{
    return std::isnan(x) || std::isnan(y);
}

template<typename Ty>
Ty Vector2<Ty>::operator [](int index) const {
    assert(index < 3 && index >= 0);
    if (index == 0) {
        return x;
    } else if (index == 1) {
        return y;
    } else {
        return z;
    }
}

template<typename Ty>
Vector2<Ty> Vector2<Ty>::operator+(const Vector2<Ty> &v) const {
    return Vector2<Ty> {x + v.x , y + v.y};
};

template<typename Ty>
Vector2<Ty> & Vector2<Ty>::operator+=(const Vector2<Ty> &v) {
    x += v.x;
    y += v.y;
    return *this;
};

template<typename Ty>
Vector2<Ty> Vector2<Ty>::operator-(const Vector2<Ty> &v) const {
    return Vector2<Ty> {x - v.x, y - v.y};
};

template<typename Ty>
Vector2<Ty> & Vector2<Ty>::operator-=(const Vector2<Ty> &v) {
    x -= v.x;
    y -= v.y;
    return *this;
};

template<typename Ty>
Vector2<Ty> Vector2<Ty>::operator*(Ty v) const {
    return Vector2<Ty> {x * v.x, y * v.y};
};

template<typename Ty>
Vector2<Ty>& Vector2<Ty>::operator*=(Ty v) {
    x *= v.x;
    y *= v.y;
    return *this;
}; 

template<typename Ty>
Vector2<Ty> Vector2<Ty>::operator/(Ty v) const {
    return Vector2<Ty> {x / v.x, y / v.y};
};

template<typename Ty>
Vector2<Ty>& Vector2<Ty>::operator/=(Ty v) {
    x /= v.x;
    y /= v.y;
    return *this;
};

template<typename Ty>
Vector2<Ty>& Vector2<Ty>::operator-() {// unary operator -, negated all val
    x = -x;
    y = -y;
    return *this;
};

/*
** Vector3
*/

template <typename Ty>
Vector3<Ty>::Vector3(const Normal3<Ty> &n):x(n.x), y(n.y), z(n.x){};

template<typename Ty>
Ty Vector3<Ty>::dot(const Vector3<Ty> &v) const {
    return x * v.x + y * v.y + z * v.z;
};

template<typename Ty>
Vector3<Ty> Vector3<Ty>::cross(const Vector3<Ty> &v) const {
    return Vector3<Ty>{
        y * v.z - z * v.y,
        z * v.x - x * v.z,
        x * v.y - y * v.x
    };
};

template<typename Ty>
Vector3<Ty>& Vector3<Ty>::abs() {
    x = std::abs(x);
    y = std::abs(y);
    z = std::abs(z);
};

template<typename Ty>
Ty Vector3<Ty>::length() const {
    return std::sqrt(length_squared());
};

template<typename Ty>
Ty Vector3<Ty>::length_squared() const {
    return x * x + y * y + z * z;
};

template<typename Ty>
Vector3<Ty> Vector3<Ty>::normalized() const {
    Vector3<Ty> result = *this;
    result /= length();
    return result;
};

template<typename Ty>
Vector3<Ty>& Vector3<Ty>::normalize() {
    *this /= length();
    return *this;
};

template<typename Ty>
bool Vector3<Ty>::has_nans() const{
    return std::isnan(x) || std::isnan(y) || std::isnan(z);
};

template<typename Ty>
Vector3<Ty>::operator Normal3<Ty>() const{
    return Normal3<Ty>(x, y, z);
};

template<typename Ty>
Ty Vector3<Ty>::operator [](int index) const {
    assert(index < 3 && index >= 0);
    if (index == 0) {
        return x;
    } else if (index == 1) {
        return y;
    } else if (index == 2){
        return z;
    } else {
        throw std::runtime_error("wrong index for retrieve vector component\n");
    }
}

template<typename Ty>
Vector3<Ty> Vector3<Ty>::operator+(const Vector3<Ty> &v) const {
    return Vector3<Ty> {x + v.x , y + v.y, z + v.z};
};

template<typename Ty>
Vector3<Ty> & Vector3<Ty>::operator+=(const Vector3<Ty> &v) {
    x += v.x;
    y += v.y;
    z += v.z
    return *this;
};

template<typename Ty>
Vector3<Ty> Vector3<Ty>::operator-(const Vector3<Ty> &v) const {
    return Vector3<Ty> {x - v.x, y - v.y, z - v.z};
};

template<typename Ty>
Vector3<Ty> & Vector3<Ty>::operator-=(const Vector3<Ty> &v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
};

template<typename Ty>
Vector3<Ty> Vector3<Ty>::operator*(Ty v) const {
    return Vector3<Ty> {x * v.x, y * v.y, z * v.z};
};

template<typename Ty>
Vector3<Ty>& Vector3<Ty>::operator*=(Ty v) {
    x *= v.x;
    y *= v.y;
    z *= v.z;
    return *this;
}; 

template<typename Ty>
Vector3<Ty> Vector3<Ty>::operator/(Ty v) const {
    return Vector3<Ty> {x / v.x, y / v.y, z / v.z};
};

template<typename Ty>
Vector3<Ty>& Vector3<Ty>::operator/=(Ty v) {
    x /= v.x;
    y /= v.y;
    z /= v.z;
    return *this;
};

template<typename Ty>
Vector3<Ty> Vector3<Ty>::operator-() const {// unary operator -, negated all val
    Vector3f v;
    v.x = -x;
    v.y = -y;
    v.z = -z;
    return v;
};


/*
** Point2
*/

template<typename Ty>
bool Point2<Ty>::has_nans() const{
    return std::isnan(x) || std::isnan(y);
};

template<typename Ty>
Ty Point2<Ty>::length() const {
    return std::sqrt(length_squared());
};

template<typename Ty>
Ty Point2<Ty>::length_squared() const {
    return x * x + y * y ;
};

template<typename Ty>
Point2<Ty>::operator Vector2<Ty>() {
    return Vector2<Ty>(x, y);
};

template<typename Ty>
Point2<Ty>::operator Vector3<Ty>() {
    return Vector3<Ty>(x, y, 0);
};

template<typename Ty>
Point2<Ty> Point2<Ty>::operator+(const Point2<Ty> &v) const {
    return Point2<Ty>(x + v.x, y + v.y);
};

template<typename Ty>
Point2<Ty> & Point2<Ty>::operator+=(const Point2<Ty> &v) {
    x += v.x;
    y += v.y;
    return *this;
};

template<typename Ty>
Point2<Ty> Point2<Ty>::operator-(const Vector2<Ty> &v) const {
    return Point2<Ty>(x - v.x, y - v.y);
};

template<typename Ty>
Vector2<Ty> Point2<Ty>::operator-(const Point2<Ty> &v) const {
    return Vector2<Ty>(x - v.x, y - v.y);
};

template<typename Ty>
Point2<Ty> & Point2<Ty>::operator-=(const Point2<Ty> &v) {
    x -= v.x;
    y -= v.y;
    return *this;
};

template<typename Ty>
Point2<Ty> Point2<Ty>::operator*(Ty v) const {
    return Point2<Ty>(x * v.x, y * v.y);
};

template<typename Ty>
Point2<Ty>& Point2<Ty>::operator*=(Ty v) {
    x *= v.x;
    y *= v.y;
    return *this;
}; 

template<typename Ty>
Point2<Ty> Point2<Ty>::operator/(Ty v) const {
    return Point2<Ty>(x / v.x, y / v.y);
};

template<typename Ty>
Point2<Ty>& Point2<Ty>::operator/=(Ty v) {
    x /= v.x;
    y /= v.y;
    return *this;
};

template<typename Ty>
Point2<Ty>& Point2<Ty>::operator-() {
    x = -x;
    y = -y;
    return *this;
}; // unary operator -, negated all val;


/*
** Point3
*/

template<typename Ty>
bool Point3<Ty>::has_nans() const{
    return std::isnan(x) || std::isnan(y) || std::isnan(z);
};

template<typename Ty>
Ty Point3<Ty>::length() const {
    return std::sqrt(length_squared());
};

template<typename Ty>
Ty Point3<Ty>::length_squared() const {
    return x * x + y * y + z * z;
};

template<typename Ty>
Point3<Ty>::operator Vector3<Ty>() {
    return Vector3<Ty>(x, y, z);
};

template<typename Ty>
Point3<Ty> Point3<Ty>::operator+(const Vector3<Ty> &v) const {
    return Point3<Ty>(x + v.x, y + v.y, z + v.z);
};

template<typename Ty>
Point3<Ty> Point3<Ty>::operator+(const Point3<Ty> &p) const {
    return Point3<Ty>(x + p.x, y + p.y, z + p.z);
};

template<typename Ty>
Point3<Ty> & Point3<Ty>::operator+=(const Point3<Ty> &p) {
    x += p.x;
    y += p.y;
    z += p.z;
    return *this;
};

template<typename Ty>
Point3<Ty> & Point3<Ty>::operator+=(const Vector3<Ty> &v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
};

template<typename Ty>
Point3<Ty> Point3<Ty>::operator-(const Vector3<Ty> &v) const {
    return Point3<Ty>(x - v.x, y - v.y, z - v.z);
};

template<typename Ty>
Vector3<Ty> Point3<Ty>::operator-(const Point3<Ty> &v) const {
    return Vector3<Ty>(x - v.x, y - v.y, z - v.z);
};

template<typename Ty>
Point3<Ty> & Point3<Ty>::operator-=(const Point3<Ty> &v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
};

template<typename Ty>
Point3<Ty> Point3<Ty>::operator*(Ty v) const {
    return Point3<Ty>(x * v.x, y * v.y, z * v.z);
};

template<typename Ty>
Point3<Ty>& Point3<Ty>::operator*=(Ty v) {
    x *= v.x;
    y *= v.y;
    z *= v.z;
    return *this;
}; 

template<typename Ty>
Point3<Ty> Point3<Ty>::operator/(Ty v) const {
    return Point3<Ty>(x / v.x, y / v.y, z / v.z);
};

template<typename Ty>
Point3<Ty>& Point3<Ty>::operator/=(Ty v) {
    x /= v.x;
    y /= v.y;
    z /= v.z;
    return *this;
};

template<typename Ty>
Point3<Ty>& Point3<Ty>::operator-() {
    x = -x;
    y = -y;
    z = -z;
    return *this;
}; // unary operator -, negated all val;

template<typename Ty>
Ty Point3<Ty>::operator[](int index) const {
    if (index == 0) {
        return x;
    } else if (index == 1) {
        return y;
    } else if (index == 2) {
        return z;
    } else {
        throw std::runtime_error("wrong index for retrieve point coordinates\n");
    }
};

/*
** Normal3
*/

template <typename Ty>
bool Normal3<Ty>::has_nans() const{
    return std::isnan(x) || std::isnan(y) || std::isnan(z);
}

template <typename Ty>
Ty Normal3<Ty>::dot(const Vector3<Ty> &v) const {
    return x * v.x + y * v.y + z * v.z;
};
        
template <typename Ty>        
Vector3<Ty> Normal3<Ty>::cross(const Vector3<Ty> &v) const {
    return Vector3<Ty> (
            y * v.z - z * v.y, 
            z * v.x - x * v.z;
            x * v.y - y * v.x;
            );
};

template <typename Ty>
Normal3<Ty>& Normal3<Ty>::scale(Ty s) {
    return this->operator*=(s); 
};

template <typename Ty>
Normal3<Ty> Normal3<Ty>::scaled(Ty s) const {
    Normal3<Ty> result = *this;
    result *= s;
    return result;
};

template <typename Ty>
Normal3<Ty>& Normal3<Ty>::normalize() {
    this->operator/=(length());
};

template <typename Ty>
Normal3<Ty> Normal3<Ty>::normalized() const {
    Normal3<Ty> result = *this;
    result /= length();
    return result;
};

template <typename Ty>
Ty Normal3<Ty>::length() const {
    return std::sqrt(length_squared());
};

template <typename Ty>
Ty Normal3<Ty>::length_squared() const {
    return x * x + y * y + z * z;
};


} // namespace eric