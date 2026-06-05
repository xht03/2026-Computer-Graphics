#include "Material.h"

#include <algorithm>
#include <cmath>
Vector3f Material::shade(
    const Ray& ray, const Hit& hit, const Vector3f& dirToLight, const Vector3f& lightIntensity)
{
    Vector3f normal = hit.getNormal().normalized();
    if (Vector3f::dot(normal, ray.getDirection()) > 0.0f) {
        normal = -normal;
    }

    Vector3f lightDir = dirToLight.normalized();
    Vector3f viewDir = (-ray.getDirection()).normalized();

    float nDotL = std::max(0.0f, Vector3f::dot(normal, lightDir));
    Vector3f diffuse = _diffuseColor * lightIntensity * nDotL;

    Vector3f specular = Vector3f::ZERO;
    if (_shininess > 0.0f) {
        Vector3f reflectDir = (2.0f * nDotL) * normal - lightDir;
        float rDotV = std::max(0.0f, Vector3f::dot(reflectDir.normalized(), viewDir));
        specular = _specularColor * lightIntensity * std::pow(rDotV, _shininess);
    }

    return diffuse + specular;
}
