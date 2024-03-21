#include "curves.hpp"

namespace eric {
Bounds3f Curve::ObjectBound() const {
    Bounds3f b = Union(Bounds3f(common->cpOjb[0], common->cpOjb[1]),
                    Bounds3f(common->cpOjb[2], common->cpOjb[3]));
    Float width[2] = {
        Lerp(uMin, common->width[0], common->width[1]),
        Lerp(uMax, common->width[0], common->width[1])
    };
    return Expand(b, std::max(width[0], width[1]) * 0.5f);
};

bool Curve::Intersect(const Ray &ray, Float *tHit, 
    SurfaceInteraction *isect, bool testAlphaTexture = true) const {
        // calculate segment control points 
        Point3f ctlPoint[4];
        ctlPoint[0] = BlossomBezier(common->cpOjb, uMin, uMin, uMin);
        ctlPoint[1] = BlossomBezier(common->cpOjb, uMin, uMin, uMax);
        ctlPoint[2] = BlossomBezier(common->cpOjb, uMin, uMax, uMax);
        ctlPoint[3] = BlossomBezier(common->cpOjb, uMax, uMax, uMax);

        Transform object2Ray;

        // transform segment to ray's coordinate system
        // decide the refined depth for recursive

    };

bool Curve::IntersectP(const Ray &ray,
    bool testAlphaTexture = true) const {

    };

Float Curve::SurfaceArea() const {
    
};

}