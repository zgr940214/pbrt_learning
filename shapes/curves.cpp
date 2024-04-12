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
        // we would like to make start point always on the left side , 
        // that way we could determine the order of cross product for edge function
        v1 = common->cpOjb[3] - common->cpOjb[0]; 
        v1 = v1 - Dot(v1, ray.d) * ray.d;
        v1 = Normalize(v1);
        v2 = Normalize(Cross(v1, ray.d));
        //LookAt row major
        Transform object2Ray(Matrix4f(
            v1.x, v1.y, v1.z, -Dot(v1, Vector3f(ray.o)),
            v2.x, v2.y, v2.z, -Dot(v2, Vector3f(ray.o)),
            ray.d.x, ray.d.y, ray.d.z, -Dot(ray.d, Vector3f(ray.o)),
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
            Bounds3f box = Union(Bounds3f(cp[0], cp[1]), Bounds3f(cp[2], cp[3]));
            if (!box._intersect_p(ray, NULL, NULL)) {
                return false;
            }
            if(maxDepth > 0) {
                Point3f ctlPoint[7];
                Point3f *cps = ctlPoint;
                Float uu[3] = {
                    uMin,
                    (uMin + uMax) / 2,
                    uMax
                };
                Float *up = uu;
                SubdivideBezier(cp, ctlPoint);
                bool hit = false;
                for (int i = 0; i < 2 ; i++, cps+=3, up++) {
                    Float maxWidth = std::max(
                        Lerp(up[0], common->width[0], common->width[1]),
                        Lerp(up[1], common->width[0], common->width[1])  
                    );
                    if (std::max(std::max(cps[0].y, cps[1].y), 
                            std::max(cps[2].y, cps[3].y)) + 0.5 * maxWidth < 0.f
                        || std::min(std::min(cps[0].y, cps[1].y), 
                            std::min(cps[2].y, cps[3].y)) - 0.5 * maxWidth > 0.f) {
                               continue;
                            }

                    if (std::max(std::max(cps[0].x, cps[1].x), 
                            std::max(cps[2].x, cps[3].x)) + 0.5 * maxWidth < 0.f
                        || std::min(std::min(cps[0].x, cps[1].x), 
                            std::min(cps[2].x, cps[3].x)) - 0.5 * maxWidth > 0.f) {
                                continue;
                            }
                    
                    Float rayLength = ray.d.length() * ray.tMax;
                    if (std::max(std::max(cps[0].z, cps[1].z), 
                            std::max(cps[2].z, cps[3].z)) + 0.5 * maxWidth < 0.f
                        || std::min(std::min(cps[0].z, cps[1].z), 
                            std::min(cps[2].z, cps[3].z)) - 0.5 * maxWidth > rayLength) {
                                continue;
                            }
                    hit |= recursiveIntersect(ray, tHit, isect, cps, ray2Object, up[0], up[1], maxDepth - 1);
                    return hit;
                }
            } else {
                // hit point vector cross bitangent should be z positive (we adjust the ray coordinate system to make sure the start
                // point's x coordinate is smaller than end point's x coordinate system which make shoule the anti clock wise cross product produce positive z
                // if the hit point is on the right side of start point's bitangent)
                Float e1 = (cp[1].y - cp[0].y) * cp[0].x - (cp[1].x - cp[0].x) * cp[0].y; 
                if (e1 < 0) {
                    return false;
                }
                // same as e1 but we adjust the cross product order , since the hit point should be left side of end point's bitangent
                Float e2 = (cp[3].x - cp[2].x) * cp[3].y - (cp[3].y - cp[2].y) * cp[3].x;
                if (e2 < 0) {
                    return false;
                }
                
                // calculete the parameter u where the point is closest to ray
                Vector2f v1 = Vector2f(cp[3] - cp[0]);
                Vector2f v2 = Vector2f(-cp[0].x, cp[0].y);
                Float w = Dot(v1, v2) / v1.length_squared();

                Float u = Lerp(w, uMin, uMax);

                Float d = Lerp(u, common->width[0], common->width[1]); 

                Float sin0 = std::sin((1 - u) * common->normalAngle) * common->invSinNormalAngle;
                Float sin1 = std::sin(u * common->normalAngle) * common->invSinNormalAngle;
                Normal3f nHit = sin0 * common->n[0] + sin1 * common->n[1];
                //test curve width 
                if (common->type == CurveType::RIBBON) {
                    Float cosTheta = Dot(nHit, ray.d) / ray.d.length(); // Slerp keeps normal length unchanged
                    d *= cosTheta; 
                }

                Point3f pHit = BlossomBezier(common->cpOjb, u, u, u);
                if (Vector2f(pHit.x, pHit.y).length() > d / 2) {
                    return false;
                }

                Float t = pHit.z;
                if (t < 0 || t > ray.tMax) {
                    return false;
                }

                if (tHit) {
                    *tHit = t;
                }

                Float edge = v2.x * v1.y - v2.y * v1.x > 0 ? 1 : -1;
                Float v = (Vector2f(pHit.x, pHit.y).length() / (d));
                if (edge < 0) {
                    v = - 0.5 + v;
                }

                Vector3f dpdu, dpdv;
                EvalBezier(common->cpOjb, u, dpdu);
                if (common->type == CurveType::RIBBON) {
                    dpdv = Normalize(Cross(nHit, dpdu)) * d;
                } else if (common->type == CurveType::FLAT) {
                    dpdv = Normalize(Cross(ray.d, dpdu)) * d;
                } else { // cylinder
                    dpdv = Normalize(Cross(ray.d, dpdu)) * d;
                    // rotate the dpdv along dpdu over theta
                    Float cosTheta, sinTheta;
                    if (v < 0) {
                        cosTheta = (0.5 + v) / (d / 2);
                    } else {
                        cosTheta =  2 * v / d;
                    }
                    sinTheta = std::sqrt(1 - cosTheta * cosTheta);
                    //Rodrigues
                    Vector3f K = Normalize(dpdu);

                    dpdv = dpdv * cosTheta + Cross(K, dpdv) * sinTheta + Dot(K, dpdv) * (1 - cosTheta) * K;
                    dpdv = Normalize(dpdv);
                    dpdv /= sinTheta;
                }

                Vector3f pErr;
                //calculate error bound
                *isect = SurfaceInteraction(pHit,
                    pErr, Point2f(u, v), -ray.d, 
                    dpdu, dpdv, Normal3f(0.f, 0.f, 0.f),
                    Normal3f(0.f, 0.f, 0.f), ray.time, this);
                
                return true;

            }

    };
        
Float Curve::SurfaceArea() const {
    Point3f cp[4];
    cp[0] = BlossomBezier(common->cpOjb, uMin, uMin, uMin);
    cp[1] = BlossomBezier(common->cpOjb, uMin, uMin, uMax);
    cp[2] = BlossomBezier(common->cpOjb, uMin, uMax, uMax);
    cp[3] = BlossomBezier(common->cpOjb, uMax, uMax, uMax);

    Float w0 = Lerp(uMin, common->width[0], common->width[1]);
    Float w1 = Lerp(uMax, common->width[0], common->width[1]);

    Float averageWidth = (w0 + w1) / 2;
    Float eLength = 0.f;
    for (int i = 0; i < 3; i++) {
        eLength += Distance(cp[i], cp[i + 1]);
    }
    return eLength * averageWidth;
};

}