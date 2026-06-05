#include "Object3D.h"

#include "VecUtils.h"

#include <cmath>

bool Sphere::intersect(const Ray& r, float tmin, Hit& h) const
{
    // BEGIN STARTER

    // We provide sphere intersection code for you.
    // You should model other intersection implementations after this one.

    // Locate intersection point ( 2 pts )
    const Vector3f& rayOrigin = r.getOrigin(); // Ray origin in the world coordinate
    const Vector3f& dir = r.getDirection();

    Vector3f origin = rayOrigin - _center; // Ray origin in the sphere coordinate

    float a = dir.absSquared();
    float b = 2 * Vector3f::dot(dir, origin);
    float c = origin.absSquared() - _radius * _radius;

    // no intersection
    if (b * b - 4 * a * c < 0) {
        return false;
    }

    float d = sqrt(b * b - 4 * a * c);

    float tplus = (-b + d) / (2.0f * a);
    float tminus = (-b - d) / (2.0f * a);

    // the two intersections are at the camera back
    if ((tplus < tmin) && (tminus < tmin)) {
        return false;
    }

    float t = 10000;
    // the two intersections are at the camera front
    if (tminus > tmin) {
        t = tminus;
    }

    // one intersection at the front. one at the back
    if ((tplus > tmin) && (tminus < tmin)) {
        t = tplus;
    }

    if (t < h.getT()) {
        Vector3f normal = r.pointAtParameter(t) - _center;
        normal = normal.normalized();
        h.set(t, this->material, normal);
        return true;
    }
    // END STARTER
    return false;
}

// Add object to group
void Group::addObject(Object3D* obj) { m_members.push_back(obj); }

// Return number of objects in group
int Group::getGroupSize() const { return (int)m_members.size(); }

bool Group::intersect(const Ray& r, float tmin, Hit& h) const
{
    // BEGIN STARTER
    // we implemented this for you
    bool hit = false;
    for (Object3D* o : m_members) {
        if (o->intersect(r, tmin, h)) {
            hit = true;
        }
    }
    return hit;
    // END STARTER
}

Plane::Plane(const Vector3f& normal, float d, Material* m)
    : Object3D(m)
{
    _normal = normal.normalized();
    _d = d;
}
bool Plane::intersect(const Ray& r, float tmin, Hit& h) const
{
    const Vector3f& dir = r.getDirection();
    float denom = Vector3f::dot(_normal, dir);
    if (std::abs(denom) < 1e-6f) {
        return false;
    }

    float t = (_d - Vector3f::dot(_normal, r.getOrigin())) / denom;
    if (t < tmin || t >= h.getT()) {
        return false;
    }

    h.set(t, this->material, _normal);
    return true;
}
bool Triangle::intersect(const Ray& r, float tmin, Hit& h) const
{
    const Vector3f& orig = r.getOrigin();
    const Vector3f& dir = r.getDirection();

    Vector3f edge1 = _v[1] - _v[0];
    Vector3f edge2 = _v[2] - _v[0];

    Vector3f pvec = Vector3f::cross(dir, edge2);
    float det = Vector3f::dot(edge1, pvec);
    if (std::abs(det) < 1e-6f) {
        return false;
    }

    float invDet = 1.0f / det;
    Vector3f tvec = orig - _v[0];
    float u = Vector3f::dot(tvec, pvec) * invDet;
    if (u < 0.0f || u > 1.0f) {
        return false;
    }

    Vector3f qvec = Vector3f::cross(tvec, edge1);
    float v = Vector3f::dot(dir, qvec) * invDet;
    if (v < 0.0f || (u + v) > 1.0f) {
        return false;
    }

    float t = Vector3f::dot(edge2, qvec) * invDet;
    if (t < tmin || t >= h.getT()) {
        return false;
    }

    Vector3f normal = (1.0f - u - v) * _normals[0] + u * _normals[1] + v * _normals[2];
    normal = normal.normalized();
    h.set(t, this->material, normal);
    return true;
}

Transform::Transform(const Matrix4f& m, Object3D* obj)
    : _object(obj)
{
    _matrix = m;
    _inverse = m.inverse();
}
bool Transform::intersect(const Ray& r, float tmin, Hit& h) const
{
    Vector3f localOrigin = VecUtils::transformPoint(_inverse, r.getOrigin());
    Vector3f localDir = VecUtils::transformDirection(_inverse, r.getDirection());
    Ray localRay(localOrigin, localDir);

    Hit localHit;
    if (!_object->intersect(localRay, tmin, localHit)) {
        return false;
    }

    Vector3f localPoint = localRay.pointAtParameter(localHit.getT());
    Vector3f worldPoint = VecUtils::transformPoint(_matrix, localPoint);
    Vector3f worldToPoint = worldPoint - r.getOrigin();
    float dirLen2 = r.getDirection().absSquared();
    float tWorld = dirLen2 > 0.0f ? Vector3f::dot(worldToPoint, r.getDirection()) / dirLen2
                                  : localHit.getT();

    if (tWorld < tmin || tWorld >= h.getT()) {
        return false;
    }

    Matrix4f inverseTranspose = _inverse.transposed();
    Vector3f worldNormal = VecUtils::transformDirection(inverseTranspose, localHit.getNormal());
    worldNormal = worldNormal.normalized();

    h.set(tWorld, localHit.getMaterial(), worldNormal);
    return true;
}