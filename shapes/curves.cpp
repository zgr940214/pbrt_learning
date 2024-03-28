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

bool Curve::Intersect(const Ray &r, Float *tHit, 
    SurfaceInteraction *isect, bool testAlphaTexture = true) const {
        // calculate segment control points 
        // transform segment to ray's coordinate system
        Vector3f oErr, dErr;
        Ray ray = (*WorldToObject)(r, &oErr, &dErr);
        Point3f ctlPoint[4];
        Vector3f v1,v2;
        CoordinateSystem(ray.d, &v1, &v2);
        Transform object2Ray(Matrix4f(
            v1.x, v2.x, ray.d.x, -ray.o.x,
            v1.y, v2.y, ray.d.y, -ray.o.y,
            v1.z, v2.z, ray.d.z, -ray.o.z,
            0.f, 0.f, 0.f, 1.f));

        ctlPoint[0] = object2Ray(BlossomBezier(common->cpOjb, uMin, uMin, uMin));
        ctlPoint[1] = object2Ray(BlossomBezier(common->cpOjb, uMin, uMin, uMax));
        ctlPoint[2] = object2Ray(BlossomBezier(common->cpOjb, uMin, uMax, uMax));
        ctlPoint[3] = object2Ray(BlossomBezier(common->cpOjb, uMax, uMax, uMax));

        Float maxWidth = std::max(
            Lerp(uMin, common->width[0], common->width[1]),
            Lerp(uMax, common->width[0], common->width[1])  
        );
        if (std::max(std::max(ctlPoint[0].y, ctlPoint[1].y), 
                std::max(ctlPoint[2].y, ctlPoint[3].y)) + 0.5 * maxWidth < 0.f
            || std::min(std::min(ctlPoint[0].y, ctlPoint[1].y), 
                std::min(ctlPoint[2].y, ctlPoint[3].y)) - 0.5 * maxWidth > 0.f) {
                    return false;
                }

        if (std::max(std::max(ctlPoint[0].x, ctlPoint[1].x), 
                std::max(ctlPoint[2].x, ctlPoint[3].x)) + 0.5 * maxWidth < 0.f
            || std::min(std::min(ctlPoint[0].x, ctlPoint[1].x), 
                std::min(ctlPoint[2].x, ctlPoint[3].x)) - 0.5 * maxWidth > 0.f) {
                    return false;
                }
        
        Float rayLength = ray.d.length() * ray.tMax;
        if (std::max(std::max(ctlPoint[0].z, ctlPoint[1].z), 
                std::max(ctlPoint[2].z, ctlPoint[3].z)) + 0.5 * maxWidth < 0.f
            || std::min(std::min(ctlPoint[0].z, ctlPoint[1].z), 
                std::min(ctlPoint[2].z, ctlPoint[3].z)) - 0.5 * maxWidth > rayLength) {
                    return false;
                }

        // decide the refined depth for recursive
        Float L0 = 0;
        for (int i = 0; i < 2; ++i)
            L0 = std::max(
                L0, std::max(
                        std::max(std::abs(ctlPoint[i].x - 2 * ctlPoint[i + 1].x + ctlPoint[i + 2].x),
                                std::abs(ctlPoint[i].y - 2 * ctlPoint[i + 1].y + ctlPoint[i + 2].y)),
                        std::abs(ctlPoint[i].z - 2 * ctlPoint[i + 1].z + ctlPoint[i + 2].z)));

        Float eps =
            std::max(common->width[0], common->width[1]) * .05f;  // width / 20
        auto Log2 = [](float v) -> int {
            if (v < 1) return 0;
            uint32_t bits = FloatToBits(v);
            // https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
            // (With an additional add so get round-to-nearest rather than
            // round down.)
            return (bits >> 23) - 127 + (bits & (1 << 22) ? 1 : 0);
        };
        // Compute log base 4 by dividing log2 in half.
        int r0 = Log2(1.41421356237f * 6.f * L0 / (8.f * eps)) / 2;
        int maxDepth = Clamp(r0, 0, 10);
        return recursiveIntersect(ray, tHit, isect, ctlPoint, Inverse(object2Ray), uMin,
                              uMax, maxDepth);
    };

 bool Curve::recursiveIntersect(const Ray &ray, Float *tHit, 
        SurfaceInteraction *isect, Point3f* cp, Transform ray2Object, 
        Float uMin, Float uMax, int maxDepth) const {
            

    };
        

bool Curve::IntersectP(const Ray &ray,
    bool testAlphaTexture = true) const {

    };

Float Curve::SurfaceArea() const {
    
};

}