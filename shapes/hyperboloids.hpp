#include "core/math.hpp"
#include "core/shape.hpp"

namespace eric {

// Hyperboloid Declarations
class Hyperboloid : public Shape {
  public:
    // Hyperboloid Public Methods
    Hyperboloid(const Transform *o2w, const Transform *w2o, bool ro,
                const Point3f &point1, const Point3f &point2, Float tm);
    Bounds3f ObjectBound() const;
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture) const;
    bool IntersectP(const Ray &ray, bool testAlphaTexture) const;
    Float Area() const;
    Interaction Sample(const Point2f &u, Float *pdf) const;

  protected:
    // Hyperboloid Private Data
    Point3f p1, p2;
    Float zMin, zMax;
    Float phiMax;
    Float rMax;
    Float ah, ch;
};
}