#include "paraboloios.hpp"

namespace eric {
Bounds3f Paraboloios::ObjectBound() const {
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
};

bool Paraboloios::Intersect(const Ray &ray, Float *tHit, 
    SurfaceInteraction *isect, bool testAlphaTexture = true) const {
        Vector3f oErr, dErr;
        Ray r = (*WorldToObject)(ray, &oErr, &dErr);
        EFloat ox = EFloat(r.o.x), oy = EFloat(r.o.y), oz = EFloat(r.o.z);
        EFloat dx = EFloat(r.d.x), dy = EFloat(r.o.y), dz = EFloat(r.d.z);

        EFloat k = height / (radius * radius);
        EFloat a = k * (dx * dx + dy * dy);
        EFloat b = 2.f * k * (dx * ox + dy * oy) - dz;
        EFloat c = k * (ox * ox + oy * oy) - oz;

        EFloat t0, t1;
        if (!Quadratic(a, b, c, &t0, &t1)) {
            return false;
        }
        if (t0.UpperBound() > ray.tMax || t1.LowerBound() < 0.f) {
            return false;
        }

        Float tShapeHit = t0;
        if (tShapeHit < 0.f) {
            tShapeHit = t1;
            if (tShapeHit > ray.tMax) {
                return false;
            }
        }

        Point3f pHit = r(tShapeHit);
        if (pHit.x == 0.f && pHit.y == 0.f) {
            pHit.x = 1e-5 * radius;
        }
        Float phi = atan2(pHit.y, pHit.x);
        if (phi < 0.f) {
            phi += 2 * M_PI;
        }

        if (pHit.z < zMin || pHit.z > zMax || phi > phiMax) {
            if (tShapeHit == (Float)t0) {
                return false;
            }
            tShapeHit = t1;
            pHit = r(tShapeHit);
            if (pHit.x == 0.f && pHit.y == 0.f) {
                pHit.x = 1e-5 * radius;
            }
            Float phi = atan2(pHit.y, pHit.x);
            if (phi < 0.f) {
                phi += 2 * M_PI;
            }
            if (pHit.z < zMin || pHit.z > zMax || phi > phiMax) {
                return false;
            }
        }

        Float u = phi / phiMax;
        Float v = (pHit.z - zMin) / (zMax - zMin);


        Vector3f dpdu = Vector3f(
            -pHit.y * phiMax,
            pHit.x * phiMax,
            0.f
        );

        Vector3f dpdv = Vector3f(
            pHit.x / pHit.z * 0.5f,
            pHit.y / pHit.z * 0.5f,
            1.f
        ) * (zMax - zMin);

        Vector3f d2Pduu = Vector3f(
            pHit.x,
            pHit.y,
            0.f
        ) * -1.f * phiMax * phiMax;

        Vector3f d2Pdvv = Vector3f(
            pHit.x / (pHit.z * pHit.z * 4.f),
            pHit.y / (pHit.z * pHit.z * 4.f),
            0.f
        ) * -1.f * (zMax - zMin) * (zMax - zMin);

        Vector3f d2Pduv = Vector3f(
            -0.5f * pHit.y / pHit.z ,
            0.5f * pHit.y / pHit.z,
            0.f
        ) * phiMax * (zMax - zMin);

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

        Vector3f pError;
        *isect = (*ObjectToWorld)(SurfaceInteraction(
            pHit, pError, Point2f(u, v), -r.d, dpdu, dpdv, dndu, dndv, ray.time, this));
        *tHit = tShapeHit;
        return true;
    };

bool Paraboloios::IntersectP(const Ray &ray,
    bool testAlphaTexture = true) const {
         Vector3f oErr, dErr;
        Ray r = (*WorldToObject)(ray, &oErr, &dErr);
        EFloat ox = EFloat(r.o.x), oy = EFloat(r.o.y), oz = EFloat(r.o.z);
        EFloat dx = EFloat(r.d.x), dy = EFloat(r.o.y), dz = EFloat(r.d.z);

        EFloat k = height / (radius * radius);
        EFloat a = k * (dx * dx + dy * dy);
        EFloat b = 2.f * k * (dx * ox + dy * oy) - dz;
        EFloat c = k * (ox * ox + oy * oy) - oz;

        EFloat t0, t1;
        if (!Quadratic(a, b, c, &t0, &t1)) {
            return false;
        }
        if (t0.UpperBound() > ray.tMax || t1.LowerBound() < 0.f) {
            return false;
        }

        Float tShapeHit = t0;
        if (tShapeHit < 0.f) {
            tShapeHit = t1;
            if (tShapeHit > ray.tMax) {
                return false;
            }
        }

        Point3f pHit = r(tShapeHit);
        if (pHit.x == 0.f && pHit.y == 0.f) {
            pHit.x = 1e-5 * radius;
        }
        Float phi = atan2(pHit.y, pHit.x);
        if (phi < 0.f) {
            phi += 2 * M_PI;
        }

        if (pHit.z < zMin || pHit.z > zMax || phi > phiMax) {
            if (tShapeHit == (Float)t0) {
                return false;
            }
            tShapeHit = t1;
            pHit = r(tShapeHit);
            if (pHit.x == 0.f && pHit.y == 0.f) {
                pHit.x = 1e-5 * radius;
            }
            Float phi = atan2(pHit.y, pHit.x);
            if (phi < 0.f) {
                phi += 2 * M_PI;
            }
            if (pHit.z < zMin || pHit.z > zMax || phi > phiMax) {
                return false;
            }
        }
        return true;
    };

/*
 **  z = h / r0^2 * r^2; => Φmax*∫{zMin, zMax}sqrt(r0^2/h *z + r0^2/4h^2) dz
 ** =>  A = 2Φmax * r0 / (3 *sqrt(h)) *(z + r0^2 /4h)^3/2 |{zMin, zMax}
 ** => A = Φmax*r0^4 / (12 * Zmax) * (4*Zmax/r0^2 * z + 1)^3/2 | {zMin, zMax}   
*/

Float Paraboloios::SurfaceArea() const {
    Float radius2 = radius * radius;
    Float k = 4 * zMax / radius2;
    return (radius2 * radius2 * phiMax / (12 * zMax * zMax)) *
           (std::pow(k * zMax + 1, 1.5f) - std::pow(k * zMin + 1, 1.5f));
};

}