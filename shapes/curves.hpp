#include "shape.hpp"
#include <memory>

namespace eric {

enum class CurveType {
    FLAT,
    CYLINDER,
    RIBBON,
};

struct CurveCommon {
    const CurveType type;
    const Point3f cpOjb[4];
    const Float width[2];
    Normal3f n[2];
    Float normalAngle, invSinNormalAngle;
};

class Curve : public Shape {
    private:
        std::shared_ptr<CurveCommon> common;
        Float uMin, uMax;
    public:
        Curve(const Transform *object2World, const Transform *world2Ojbect, bool reverseOrientation,
            std::shared_ptr<CurveCommon> com, Float min, Float max): 
                Shape(object2World, world2Ojbect, reverseOrientation), 
                common(com), uMin(min), uMax(max) {};
        ~Curve() {};

        Bounds3f ObjectBound() const override;

        bool Intersect(const Ray &ray, Float *tHit, 
            SurfaceInteraction *isect, bool testAlphaTexture = true) const;
        
        bool IntersectP(const Ray &ray,
            bool testAlphaTexture = true) const;
        
        Float SurfaceArea() const;

        static Point3f BlossomBezier(const Point3f p[4], Float u0, Float u1,
                Float u2) {
                    Point3f a[3] = { Lerp(u0, p[0], p[1]),
                    Lerp(u0, p[1], p[2]),
                    Lerp(u0, p[2], p[3]) };
                    Point3f b[2] = { Lerp(u1, a[0], a[1]), Lerp(u1, a[1], a[2]) };
                    return Lerp(u2, b[0], b[1]);
                }
        
        static void SubdivideBezier(const Point3f p[4], Point3f cpSplit[7]) {
            cpSplit[0] = p[0];
            cpSplit[1] = (p[0] + p[1]) / 2;
            cpSplit[2] = (p[0] + 2 * p[1] + p[2]) / 4;
            cpSplit[3] = (p[0] + 3 * p[1] + 3 * p[2] + p[3]) / 8;
            cpSplit[4] = (p[1] + 2 * p[2] + p[3]) / 4;
            cpSplit[5] = (p[2] + p[3]) / 2;
            cpSplit[6] = p[3];
        }
};
}