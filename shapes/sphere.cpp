#include "sphere.hpp"

namespace eric {

Bounds3f Sphere::ObjectBound() const {
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

bool Sphere::Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect, bool testAlphaTexture) const {
    Vector3f oErr, dErr;
    Ray r = (*WorldToObject)(ray, &oErr, &dErr);
    EFloat dx(r.d.x, dErr.x), dy(r.d.y, dErr.y), dz(r.d.z, dErr.z);
    EFloat ox(r.o.x, oErr.x), oy(r.o.y, oErr.y), oz(r.o.z, oErr.z);

    EFloat a = dx * dx + dy * dy + dz * dz;
    EFloat b = EFloat(2.f) * (dx * oy + dy * oz + dz * ox);
    EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1)) {
        return false;
    }
    // check the intersection point
    if (t0.UpperBound() > r.tMax || t1.LowerBound() < 0.f) {
        return false;
    }

    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() < 0.f) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > r.tMax) {
            return false;
        }
    }

    //compute hit position and phi for partial sphere
    auto pHit = r((Float)tShapeHit);
    //refine hit point
    if (pHit.x == 0.f && pHit.y == 0.f) { // if x== 0.f y == 0.f  atan2 may not result a stable value so refine x to a small number
        pHit.x = 1e-5f * radius;
    }

    auto phi = atan2(pHit.y, pHit.x);
    if (phi < 0.f) {
        phi += 2 * M_PI;
    }

    if ((zMax < radius && pHit.z > zMax) 
        || (zMax > -radius && pHit.z < zMin)
        || phi > phiMax) {
            if (tShapeHit == t1) {
                return false;
            }
            tShapeHit = t1;
            // reCalculate pHit and check the condition
            pHit = r(tShapeHit);
            if (pHit.x == 0.f && pHit.y == 0.f) { // if x== 0.f y == 0.f  atan2 may not result a stable value so refine x to a small number
                pHit.x = 1e-5f * radius;
            }

            phi = atan2(pHit.y, pHit.x);
            if (phi < 0.f) {
                phi += 2 * M_PI;
            }

            if ((zMax < radius && pHit.z > zMax) 
                || (zMax > -radius && pHit.z < zMin)
                || phi > phiMax) {
                    return false;
                }
    }

    Float u = phi / phiMax;
    Float theta = acos(Clamp((Float)pHit.z / radius, -1, 1));
    Float v = (theta - thetaMin) / (thetaMax - thetaMin);
    // calculate partial derivative dp/du dp/dv
    Vector3f dpdv, dpdu;

    dpdu.x = -1.f * phiMax * pHit.y;
    dpdu.y = phiMax * pHit.x;
    dpdu.z = 0.f;

    const Float sinPhi = sin(phi);
    const Float cosPhi = cos(phi);
    const Float sinTheta = sin(theta);
    const Float delta_theta = thetaMax - thetaMin;
    dpdv.x = pHit.z * sinPhi * delta_theta;
    dpdv.y = pHit.z * cosPhi * delta_theta;
    dpdv.z = -1.f * radius * sinTheta * delta_theta;

    //calculate partial derivative of surface normal dn/du dn/dv (Weingarten equations)

    Vector3f d2Pduu = -phiMax * phiMax * Vector3f(pHit.x, pHit.y, 0);
    Vector3f d2Pduv =
        (thetaMax - thetaMin) * pHit.z * phiMax * Vector3f(-sinPhi, cosPhi, 0.);
    Vector3f d2Pdvv = -(thetaMax - thetaMin) * (thetaMax - thetaMin) *
                      Vector3f(pHit.x, pHit.y, pHit.z);

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

    //Calculate SurfaceInteraction
    Vector3f pError = gamma(5) * (Vector3f)Abs(pHit); 
    *isect = (*ObjectToWorld)(
        SurfaceInteraction(pHit, pError, Point2f(u, v), -r.d, dpdu, dpdv,
        dndu, dndv, ray.time, this));
    *tHit = tShapeHit;
    return true;
};

bool Sphere::IntersectP(const Ray &ray, bool testAlphaTexture = true) const {
    Vector3f oErr, dErr;
    Ray r = (*WorldToObject)(ray, &oErr, &dErr);
    EFloat dx(r.d.x), dy(r.d.y), dz(r.d.z);
    EFloat ox(r.o.x), oy(r.o.y), oz(r.o.z);

    EFloat a = dx * dx + dy * dy + dz * dz;
    EFloat b = EFloat(2.f) * (dx * oy + dy * oz + dz * ox);
    EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1)) {
        return false;
    }
    // check the intersection point
    if (t0.LowerBound() > r.tMax || t1.UpperBound() < 0.f) {
        return false;
    }

    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() < 0.f) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > r.tMax) {
            return false;
        }
    }

    //compute hit position and phi for partial sphere
    auto pHit = r((Float)tShapeHit);
    //refine hit point
    if (pHit.x == 0.f && pHit.y == 0.f) { // if x== 0.f y == 0.f  atan2 may not result a stable value so refine x to a small number
        pHit.x = 1e-5f * radius;
    }

    auto phi = atan2(pHit.y, pHit.x);
    if (phi < 0.f) {
        phi += 2 * M_PI;
    }

    if ((zMax < radius && pHit.z > zMax) 
        || (zMax > -radius && pHit.z < zMin)
        || phi > phiMax) {
            if (tShapeHit == t1) {
                return false;
            }
            tShapeHit = t1;
            // reCalculate pHit and check the condition
            pHit = r(tShapeHit);
            if (pHit.x == 0.f && pHit.y == 0.f) { // if x== 0.f y == 0.f  atan2 may not result a stable value so refine x to a small number
                pHit.x = 1e-5f * radius;
            }

            phi = atan2(pHit.y, pHit.x);
            if (phi < 0.f) {
                phi += 2 * M_PI;
            }

            if ((zMax < radius && pHit.z > zMax) 
                || (zMax > -radius && pHit.z < zMin)
                || phi > phiMax) {
                    return false;
                }
    }
    return true;
}

Float Sphere::SurfaceArea() const {
    return phiMax * radius * (zMax - zMin);
};

};