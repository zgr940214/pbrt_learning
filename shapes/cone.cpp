#include "cone.hpp"

namespace eric {
Bounds3f Cone::ObjectBound() const {
    if (phiMax < 2 * M_PI && phiMax > 3 / 2 * M_PI) {
        return Bounds3f(Point3f(-radius, -radius, 0), Point3f(radius, radius, height));
    } else if (phiMax > M_PI && phiMax < 3 / 2 * M_PI){
        return Bounds3f(Point3f(cos(1.5f * M_PI - phiMax) * (-radius), -radius, 0), Point3f(radius, radius, height));
    } else if (phiMax > 0.5f * M_PI && phiMax < M_PI) {
        return Bounds3f(Point3f(0.f, cos(M_PI - phiMax) * (-radius), 0), Point3f(radius, radius, height));
    } else if (phiMax > 0.0005f) {
        return Bounds3f(Point3f(0.f, 0.f, 0), Point3f(radius, cos(0.5f * M_PI - phiMax) * radius, height));
    } else {
        return Bounds3f(Point3f(0.f, 0.f, 0.f), Point3f(0.f, 0.f, 0.f));
    }
};

bool Cone::Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect, bool testAlphaTexture = true) const {
    Vector3f dErr, oErr;
    Ray r = (*WorldToObject)(ray);

    EFloat dx = EFloat(r.d.x), dy = EFloat(r.d.y), dz = EFloat(r.d.z);
    EFloat ox = EFloat(r.o.x), oy = EFloat(r.o.y), oz = EFloat(r.o.z);

    EFloat a = dx * dx + dy * dy + EFloat((radius / height) * (radius / height)) * dz * dz;
    EFloat b = 2.f * (ox * dx + oy * dy - EFloat(radius / height) * (oz - EFloat(height)) * dz);
    EFloat c = ox * ox + oy * oy - EFloat(radius / height) * EFloat(radius / height) * (oz - EFloat(height));

    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1)) {
        return false;
    }
    if (t1.LowerBound() < 0.f || t0.UpperBound() > r.tMax) {
        return false;
    }

    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() < 0.f) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > r.tMax) {
            return false;
        }
    }

    Point3f pHit = r(tShapeHit);
    
    if (pHit.x == 0.f && pHit.y == 0.f) {
        pHit.x = radius * 1e-5;
    }

    Float phi = atan2(pHit.y, pHit.x);
    if (phi < 0.f) {
        phi += 2 * M_PI;
    }
    if (phi > phiMax) {
        if (tShapeHit == t1) {
            return false;
        }
        tShapeHit = t0;
        pHit = r(tShapeHit);
        phi = atan2(pHit.y, pHit.x);
        if (phi < 0.f) {
            phi += 2 * M_PI;
        }
        if (phi > phiMax) {
            return false;
        }
    }

    Float u = phi / phiMax;
    Float v = pHit.z / height;

    Vector3f dpdu = Vector3f(
        -1.f * phiMax * pHit.y,
        phiMax * pHit.x,
        0.f
    );
    
    Vector3f dpdv = Vector3f(
        -1.f * (pHit.x) / (1.f - v),
        -1.f * (pHit.y) / (1.f - v),
        height
    );

    Vector3f d2Pduu = Vector3f(pHit.x, pHit.y, 0.f) * (phiMax * phiMax * -1.f);
    Vector3f d2Pduv = Vector3f(pHit.y, pHit.x, 0.f) * phiMax * (1.f - v);
    Vector3f d2Pdvv = Vector3f(0.f, 0.f, 0.f);

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

bool Cone::IntersectP(const Ray &ray, bool testAlphaTexture = true) const {
    Vector3f dErr, oErr;
    Ray r = (*WorldToObject)(ray);

    EFloat dx = EFloat(r.d.x), dy = EFloat(r.d.y), dz = EFloat(r.d.z);
    EFloat ox = EFloat(r.o.x), oy = EFloat(r.o.y), oz = EFloat(r.o.z);

    EFloat a = dx * dx + dy * dy + EFloat((radius / height) * (radius / height)) * dz * dz;
    EFloat b = 2.f * (ox * dx + oy * dy - EFloat(radius / height) * (oz - EFloat(height)) * dz);
    EFloat c = ox * ox + oy * oy - EFloat(radius / height) * EFloat(radius / height) * (oz - EFloat(height));

    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1)) {
        return false;
    }
    if (t1.LowerBound() < 0.f || t0.UpperBound() > r.tMax) {
        return false;
    }

    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() < 0.f) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > r.tMax) {
            return false;
        }
    }
    Point3f pHit = r(tShapeHit);
    if (pHit.x == 0.f && pHit.y == 0.f) {
        pHit.x = radius * 1e-5;
    }

    Float phi = atan2(pHit.y, pHit.x);
    if (phi < 0.f) {
        phi += 2 * M_PI;
    }
    if (phi > phiMax) {
        if (tShapeHit == t1) {
            return false;
        }
        tShapeHit = t0;
        pHit = r(tShapeHit);
        phi = atan2(pHit.y, pHit.x);
        if (phi < 0.f) {
            phi += 2 * M_PI;
        }
        if (phi > phiMax) {
            return false;
        }
    }
    return true;
};

Float Cone::SurfaceArea() const {
    return phiMax * sqrt((radius * radius) + (height * height)) * 0.5f * radius;
};

} // namespace eric