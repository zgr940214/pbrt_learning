#include "cylinder.hpp"

namespace eric
{

Bounds3f Cylinder::ObjectBound() const {
    if (phiMax < 2 * M_PI && phiMax > 3 / 2 * M_PI) {
        return Bounds3f(Point3f(-radius, -radius, zMin), Point3f(radius, radius, zMax));
    } else if (phiMax > M_PI && phiMax < 3 / 2 * M_PI){
        return Bounds3f(Point3f(cos(1.5f * M_PI - phiMax) * (-radius), -radius, zMin), Point3f(radius, radius, zMax));
    } else if (phiMax > 0.5f * M_PI && phiMax < M_PI) {
        return Bounds3f(Point3f(0.f, cos(M_PI - phiMax) * (-radius), zMin), Point3f(radius, radius, zMax));
    } else if (phiMax > 0.0005f) {
        return Bounds3f(Point3f(0.f, 0.f, zMin), Point3f(radius, cos(0.5f * M_PI - phiMax) * radius, zMax));
    } else {
        return Bounds3f(Point3f(0.f, 0.f, 0.f), Point3f(0.f, 0.f, 0.f));
    }
}

bool Cylinder::Intersect(const Ray &ray, Float* tHit,  
    SurfaceInteraction *isect, bool testAlphaTexture = true) const {
        Vector3f oErr, dErr;
        Ray r = (*WorldToObject)(ray, &oErr, &dErr); // object space ray
        EFloat dx = EFloat(r.d.x, dErr.x), dy = EFloat(r.d.y, dErr.y), dz = EFloat(r.d.z, dErr.z);
        EFloat ox = EFloat(r.o.x, oErr.x), oy = EFloat(r.o.y, oErr.y), oz = EFloat(r.o.z, oErr.z);

        //compute quadratics parameters
        EFloat a = dx * dx + dy * dy;
        EFloat b = 2.f * (dx * ox + dy * oy);
        EFloat c = ox * ox + oy * oy - EFloat(radius) * EFloat(radius);

        EFloat t0, t1;
        if (!Quadratic(a, b, c, &t0, &t1)) {
            return false;
        }
        
        if (t1.UpperBound() < 0.f || t0.LowerBound() > r.tMax) {
            return false;
        } 
        EFloat tShapeHit = t0;
        if (tShapeHit.LowerBound() <= 0.f) {
            tShapeHit = t1;
        }
        if (tShapeHit.UpperBound() > r.tMax) {
            return false;
        }

        Point3f pHit = r(tShapeHit);
        Float x = pHit.x, y = pHit.y, z = pHit.z;
        Float phi = atan2(y, x);
        if (phi < 0.f) {
            phi += M_PI * 2.f;
        }
        //check intersection point in z phi range
        if (z < zMin || z > zMax || phi > phiMax) {
            if (tShapeHit == t1) {
                return false;
            }
            tShapeHit = t1;
            pHit = r(tShapeHit);
            x = pHit.x;
            y = pHit.y;
            z = pHit.z;
            if (z < zMin || z > zMax || phi > phiMax) {
                return false;
            }
        }

        Float v = (z - zMin) / (zMax - zMin);
        Float u = phi / phiMax;

        Vector3f dpdu = Vector3f(-1.f * phiMax * y, 
                        phiMax * x,
                        0.f);
        Vector3f dpdv = Vector3f( 0.f,
                        0.f,
                        (zMax - zMin));
        
        Vector3f d2Pduu = -phiMax * phiMax * Vector3f(pHit.x, pHit.y, 0);
        Vector3f d2Pduv(0, 0, 0), d2Pdvv(0, 0, 0);

        Float E = Dot(dpdu, dpdu);
        Float F = Dot(dpdu, dpdv);
        Float G = Dot(dpdv, dpdv);
        Vector3f N = Normalize(Cross(dpdu, dpdv));
        Float e = Dot(N, d2Pduu);
        Float f = Dot(N, d2Pduv);
        Float g = Dot(N, d2Pdvv);

        Float invEGF2 = 1 / (E * G - F * F);
        Normal3f dndu = Normal3f((f * F - e * G) * invEGF2 * dpdu +
        (e * F - f * E) * invEGF2 * dpdv);
        Normal3f dndv = Normal3f((g * F - f * G) * invEGF2 * dpdu +
        (f * F - g * E) * invEGF2 * dpdv);

        Vector3f pError = gamma(3) * (Vector3f)Abs(pHit); // refine error 
        *isect = (*ObjectToWorld)(
            SurfaceInteraction(pHit, pError, Point2f(u, v), -r.d, dpdu, dpdv,
            dndu, dndv, ray.time, this));
        *tHit = tShapeHit;
        return true;
    }

bool Cylinder::IntersectP(const Ray &ray, bool testAlphaTexture = true) const {
    Vector3f oErr, dErr;
        Ray r = (*WorldToObject)(ray, &oErr, &dErr); // object space ray
        EFloat dx = EFloat(r.d.x, dErr.x), dy = EFloat(r.d.y, dErr.y), dz = EFloat(r.d.z, dErr.z);
        EFloat ox = EFloat(r.o.x, oErr.x), oy = EFloat(r.o.y, oErr.y), oz = EFloat(r.o.z, oErr.z);

        //compute quadratics parameters
        EFloat a = dx * dx + dy * dy;
        EFloat b = 2.f * (dx * ox + dy * oy);
        EFloat c = ox * ox + oy * oy - EFloat(radius) * EFloat(radius);

        EFloat t0, t1;
        if (!Quadratic(a, b, c, &t0, &t1)) {
            return false;
        }
        
        if (t1.LowerBound() < 0.f || t0.UpperBound() > r.tMax) {
            return false;
        } 

        EFloat tShapeHit = t0;
        if (tShapeHit.LowerBound() <= 0.f) {
            tShapeHit = t1;
            if (tShapeHit.UpperBound() > r.tMax) {
                return false;
            }
        }

        Point3f pHit = r(tShapeHit);
        Float x = pHit.x, y = pHit.y, z = pHit.z;
        Float phi = atan2(y, x);
        if (phi < 0.f) {
            phi += M_PI * 2.f;
        }
        //check intersection point in z phi range
        if (z < zMin || z > zMax || phi > phiMax) {
            if (tShapeHit == t1) {
                return false;
            }
            tShapeHit = t1;
            pHit = r(tShapeHit);
            x = pHit.x;
            y = pHit.y;
            z = pHit.z;
            if (z < zMin || z > zMax || phi > phiMax) {
                return false;
            }
        }
        return true;
} 

Float Cylinder::SurfaceArea() const {
    return (zMax - zMin) * phiMax * radius;
}

} // namespace eric