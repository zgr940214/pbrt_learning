#include "ray.hpp"

/*
-----------------------------------------------------------------------------------------
** Ray class implementation
----------------------------------------------------------------------------------------- 
*/  
namespace eric {
Point3f Ray::operator()(Float t) const {
    return (o + t*d);
};
}