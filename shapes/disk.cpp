#include "disk.hpp"

namespace eric {
Bounds3f Disk::ObjectBound() const {
    return Bounds3f(Point3f(-radius, -radius, height), Point3f(radius, radius, height));
};

bool Disk::Intersect (const Ray &ray, Float *tHit, 
    SurfaceInteraction *isect, bool testAlphaTexture = true) const {
        Vector3f dErr, oErr;
        Ray r = (*WorldToObject)(ray);

        if (r.d.z == 0.f) {
            return false;
        }

        Float t = (height - r.o.z) / r.d.z;
        
        if (t < 0.f || t > r.tMax) {
            return false;
        }

        Point3f pHit = r(t);
        Float distance = pHit.x * pHit.x + pHit.y * pHit.y;
        if (distance < innerRadius * innerRadius || distance > radius * radius) {
            return false;
        }

        Float phi = atan2(pHit.y, pHit.x);
        if (phi < 0.f) {
            phi += 2 * M_PI;
        }

        if (phi > phiMax) {
            return false;
        }

        Float rd = sqrt(distance);
        Float u = phi / phiMax;
        Float oneMinusV = (rd - innerRadius) / (radius - innerRadius);
        Float v = 1 - oneMinusV;

        Vector3f dpdu = Vector3f(
            -1.f * phiMax * pHit.y,
            phiMax * pHit.x,
            0.f
        );

        Vector3f dpdv = Vector3f(
            pHit.x / rd * (innerRadius - radius),
            pHit.y / rd * (innerRadius - radius),
            0.f
        );

        Normal3f dndu = Normal3f(0.f, 0.f, 0.f);
        Normal3f dndv = Normal3f(0.f, 0.f, 0.f);

        Vector3f pError(0.f, 0.f, 0.f);

        *isect = (*ObjectToWorld)(
            SurfaceInteraction(pHit, pError, Point2f(u, v), -r.d, dpdu, dpdv,
            dndu, dndv, ray.time, this));
        *tHit = t;
        return true;
    };

bool Disk::IntersectP (const Ray &ray, bool testAlphaTexture = true) const {
    Vector3f dErr, oErr;
    Ray r = (*WorldToObject)(ray);

    if (r.d.z == 0.f) {
        return false;
    }

    Float t = (height - r.o.z) / r.d.z;
    
    if (t < 0.f || t > r.tMax) {
        return false;
    }

    Point3f pHit = r(t);
    Float distance = pHit.x * pHit.x + pHit.y * pHit.y;
    if (distance < innerRadius * innerRadius || distance > radius * radius) {
        return false;
    }

    Float phi = atan2(pHit.y, pHit.x);
    if (phi < 0.f) {
        phi += 2 * M_PI;
    }

    if (phi > phiMax) {
        return false;
    }
    return true;
};

Float Disk::SurfaceArea() const {
    return phiMax * 0.5f * (radius * radius - innerRadius * innerRadius);
};

}