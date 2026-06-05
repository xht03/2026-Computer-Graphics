#include "Renderer.h"

#include "ArgParser.h"
#include "Camera.h"
#include "Image.h"
#include "Ray.h"
#include "VecUtils.h"

#include <limits>

Renderer::Renderer(const ArgParser& args)
    : _args(args)
    , _scene(args.input_file)
{
}

void Renderer::Render()
{
    int w = _args.width;
    int h = _args.height;

    Image image(w, h);
    Image nimage(w, h);
    Image dimage(w, h);

    // loop through all the pixels in the image
    // generate all the samples

    // This look generates camera rays and callse traceRay.
    // It also write to the color, normal, and depth images.
    // You should understand what this code does.
    Camera* cam = _scene.getCamera();
    for (int y = 0; y < h; ++y) {
        float ndcy = 2 * (y / (h - 1.0f)) - 1.0f;
        for (int x = 0; x < w; ++x) {
            float ndcx = 2 * (x / (w - 1.0f)) - 1.0f;
            // Use PerspectiveCamera to generate a ray.
            // You should understand what generateRay() does.
            Ray r = cam->generateRay(Vector2f(ndcx, ndcy));

            Hit h;
            Vector3f color = traceRay(r, cam->getTMin(), _args.bounces, h);

            image.setPixel(x, y, color);
            nimage.setPixel(x, y, (h.getNormal() + 1.0f) / 2.0f);
            float range = (_args.depth_max - _args.depth_min);
            if (range) {
                dimage.setPixel(x, y, Vector3f((h.t - _args.depth_min) / range));
            }
        }
    }
    // END SOLN

    // save the files
    if (_args.output_file.size()) {
        image.savePNG(_args.output_file);
    }
    if (_args.depth_file.size()) {
        dimage.savePNG(_args.depth_file);
    }
    if (_args.normals_file.size()) {
        nimage.savePNG(_args.normals_file);
    }
}

Vector3f Renderer::traceRay(const Ray& r, float tmin, int bounces, Hit& h) const
{
    // The starter code only implements basic drawing of sphere primitives.
    // You will implement phong shading, recursive ray tracing, and shadow rays.

    if (!_scene.getGroup()->intersect(r, tmin, h)) {
        return _scene.getBackgroundColor(r.getDirection());
    }

    const float epsilon = 1e-4f;
    Vector3f hitPoint = r.pointAtParameter(h.getT());
    Vector3f normal = h.getNormal().normalized();
    if (Vector3f::dot(normal, r.getDirection()) > 0.0f) {
        normal = -normal;
    }

    Material* material = h.getMaterial();
    Vector3f color = _scene.getAmbientLight() * material->getDiffuseColor();

    int numLights = _scene.getNumLights();
    for (int i = 0; i < numLights; ++i) {
        Vector3f toLight;
        Vector3f intensity;
        float distToLight = 0.0f;
        _scene.getLight(i)->getIllumination(hitPoint, toLight, intensity, distToLight);

        if (_args.shadows) {
            Ray shadowRay(hitPoint + normal * epsilon, toLight);
            Hit shadowHit;
            if (_scene.getGroup()->intersect(shadowRay, epsilon, shadowHit)) {
                if (shadowHit.getT() < distToLight) {
                    continue;
                }
            }
        }

        color += material->shade(r, h, toLight, intensity);
    }

    if (bounces > 0) {
        Vector3f reflectDir = r.getDirection().normalized();
        reflectDir = reflectDir - 2.0f * Vector3f::dot(reflectDir, normal) * normal;
        Ray reflectRay(hitPoint + normal * epsilon, reflectDir);
        Hit reflectHit;
        Vector3f reflectedColor = traceRay(reflectRay, epsilon, bounces - 1, reflectHit);
        color += material->getSpecularColor() * reflectedColor;
    }

    return color;
}
