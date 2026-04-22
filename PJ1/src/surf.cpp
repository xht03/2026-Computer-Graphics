#include "surf.h"
#include "vertexrecorder.h"
#include <cmath>
using namespace std;

namespace {
const float c_pi = 3.14159265358979323846f;

inline bool approx(const Vector3f& lhs, const Vector3f& rhs)
{
    const float eps = 1e-8f;
    return (lhs - rhs).absSquared() < eps;
}

inline Vector3f safeNormalized(const Vector3f& v, const Vector3f& fallback)
{
    const float eps = 1e-12f;
    if (v.absSquared() < eps)
        return fallback;
    Vector3f out = v;
    out.normalize();
    return out;
}

// We're only implenting swept surfaces where the profile curve is
// flat on the xy-plane.  This is a check function.
static bool checkFlat(const Curve& profile)
{
    for (unsigned i = 0; i < profile.size(); i++)
        if (profile[i].V[2] != 0.0 || profile[i].T[2] != 0.0 || profile[i].N[2] != 0.0)
            return false;

    return true;
}
}

// DEBUG HELPER
Surface quad()
{
    Surface ret;
    ret.VV.push_back(Vector3f(-1, -1, 0));
    ret.VV.push_back(Vector3f(+1, -1, 0));
    ret.VV.push_back(Vector3f(+1, +1, 0));
    ret.VV.push_back(Vector3f(-1, +1, 0));

    ret.VN.push_back(Vector3f(0, 0, 1));
    ret.VN.push_back(Vector3f(0, 0, 1));
    ret.VN.push_back(Vector3f(0, 0, 1));
    ret.VN.push_back(Vector3f(0, 0, 1));

    ret.VF.push_back(Tup3u(0, 1, 2));
    ret.VF.push_back(Tup3u(0, 2, 3));
    return ret;
}

Surface makeSurfRev(const Curve& profile, unsigned steps)
{
    Surface surface;

    if (steps < 3)
        steps = 3;
    if (profile.size() < 2)
        return surface;
    if (!checkFlat(profile)) {
        cerr << "surfRev profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    const unsigned profCount = (unsigned)profile.size();
    const unsigned rings = steps + 1; // duplicate seam vertex ring

    surface.VV.reserve(rings * profCount);
    surface.VN.reserve(rings * profCount);

    for (unsigned i = 0; i < rings; ++i) {
        const float theta = 2.0f * c_pi * float(i) / float(steps);
        const float c = cos(theta);
        const float s = sin(theta);

        for (unsigned j = 0; j < profCount; ++j) {
            const Vector3f p = profile[j].V;
            const Vector3f n2 = profile[j].N;

            // Rotate around y-axis
            Vector3f V(p.x() * c + p.z() * s, p.y(), -p.x() * s + p.z() * c);
            Vector3f N(n2.x() * c + n2.z() * s, n2.y(), -n2.x() * s + n2.z() * c);
            N = safeNormalized(N, Vector3f(c, 0, -s));

            // Ensure outward normal (away from y-axis)
            Vector3f radial(V.x(), 0, V.z());
            if (radial.absSquared() > 1e-12f && Vector3f::dot(N, radial) < 0)
                N = -N;

            surface.VV.push_back(V);
            surface.VN.push_back(N);
        }
    }

    // Triangulate between rings; profile direction is not wrapped (open)
    for (unsigned i = 0; i < steps; ++i) {
        for (unsigned j = 0; j + 1 < profCount; ++j) {
            const unsigned a = i * profCount + j;
            const unsigned b = (i + 1) * profCount + j;
            const unsigned c0 = (i + 1) * profCount + (j + 1);
            const unsigned d = i * profCount + (j + 1);

            // Two triangles (a,b,c0) and (a,c0,d)
            surface.VF.push_back(Tup3u(a, b, c0));
            surface.VF.push_back(Tup3u(a, c0, d));
        }
    }

    return surface;
}

Surface makeGenCyl(const Curve& profile, const Curve& sweep)
{
    Surface surface;
    if (profile.size() < 2 || sweep.size() < 2)
        return surface;
    if (!checkFlat(profile)) {
        cerr << "genCyl profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    bool profileClosed = approx(profile.front().V, profile.back().V);
    bool sweepClosed = approx(sweep.front().V, sweep.back().V);

    unsigned profCount = (unsigned)profile.size();
    unsigned sweepCount = (unsigned)sweep.size();
    if (profileClosed && profCount > 1)
        profCount -= 1;
    if (sweepClosed && sweepCount > 1)
        sweepCount -= 1;

    if (profCount < 2 || sweepCount < 2)
        return surface;

    surface.VV.reserve(profCount * sweepCount);
    surface.VN.reserve(profCount * sweepCount);

    auto idx = [profCount](unsigned i, unsigned j) { return i * profCount + j; };

    for (unsigned i = 0; i < sweepCount; ++i) {
        const CurvePoint& sw = sweep[i];

        // Approximate sweep derivatives with finite differences (scale irrelevant)
        unsigned ip = (i == 0) ? (sweepClosed ? (sweepCount - 1) : 0) : (i - 1);
        unsigned in = (i + 1 == sweepCount) ? (sweepClosed ? 0 : (sweepCount - 1)) : (i + 1);

        const Vector3f dV = sweep[in].V - sweep[ip].V;
        const Vector3f dN = sweep[in].N - sweep[ip].N;
        const Vector3f dB = sweep[in].B - sweep[ip].B;

        for (unsigned j = 0; j < profCount; ++j) {
            const CurvePoint& pr = profile[j];
            const float x = pr.V.x();
            const float y = pr.V.y();

            const Vector3f offset = x * sw.N + y * sw.B;
            const Vector3f V = sw.V + offset;

            // Surface tangents
            const Vector3f Su = pr.T.x() * sw.N + pr.T.y() * sw.B;
            const Vector3f Sv = dV + x * dN + y * dB;

            Vector3f N = Vector3f::cross(Sv, Su);
            N = safeNormalized(N, sw.N);

            // Ensure outward (away from the sweep curve centerline)
            if (offset.absSquared() > 1e-12f && Vector3f::dot(N, offset) < 0)
                N = -N;

            surface.VV.push_back(V);
            surface.VN.push_back(N);
        }
    }

    const unsigned sweepSegs = sweepClosed ? sweepCount : (sweepCount - 1);
    const unsigned profSegs = profileClosed ? profCount : (profCount - 1);

    for (unsigned i = 0; i < sweepSegs; ++i) {
        const unsigned inext = (i + 1) % sweepCount;
        for (unsigned j = 0; j < profSegs; ++j) {
            const unsigned jnext = (j + 1) % profCount;

            const unsigned a = idx(i, j);
            const unsigned b = idx(inext, j);
            const unsigned c0 = idx(inext, jnext);
            const unsigned d = idx(i, jnext);

            surface.VF.push_back(Tup3u(a, b, c0));
            surface.VF.push_back(Tup3u(a, c0, d));
        }
    }

    return surface;
}

void recordSurface(const Surface& surface, VertexRecorder* recorder)
{
    const Vector3f WIRECOLOR(0.4f, 0.4f, 0.4f);
    for (int i = 0; i < (int)surface.VF.size(); i++) {
        recorder->record(surface.VV[surface.VF[i][0]], surface.VN[surface.VF[i][0]], WIRECOLOR);
        recorder->record(surface.VV[surface.VF[i][1]], surface.VN[surface.VF[i][1]], WIRECOLOR);
        recorder->record(surface.VV[surface.VF[i][2]], surface.VN[surface.VF[i][2]], WIRECOLOR);
    }
}

void recordNormals(const Surface& surface, VertexRecorder* recorder, float len)
{
    const Vector3f NORMALCOLOR(0, 1, 1);
    for (int i = 0; i < (int)surface.VV.size(); i++) {
        recorder->record_poscolor(surface.VV[i], NORMALCOLOR);
        recorder->record_poscolor(surface.VV[i] + surface.VN[i] * len, NORMALCOLOR);
    }
}

void outputObjFile(ostream& out, const Surface& surface)
{

    for (int i = 0; i < (int)surface.VV.size(); i++)
        out << "v  " << surface.VV[i][0] << " " << surface.VV[i][1] << " " << surface.VV[i][2]
            << endl;

    for (int i = 0; i < (int)surface.VN.size(); i++)
        out << "vn " << surface.VN[i][0] << " " << surface.VN[i][1] << " " << surface.VN[i][2]
            << endl;

    out << "vt  0 0 0" << endl;

    for (int i = 0; i < (int)surface.VF.size(); i++) {
        out << "f  ";
        for (unsigned j = 0; j < 3; j++) {
            unsigned a = surface.VF[i][j] + 1;
            out << a << "/" << "1" << "/" << a << " ";
        }
        out << endl;
    }
}
