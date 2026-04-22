#include "curve.h"
#include "vertexrecorder.h"
#include <cmath>
using namespace std;

const float c_pi = 3.14159265358979323846f;

namespace {
// Approximately equal to.  We don't want to use == because of
// precision issues with floating point.
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

inline void orthonormalizeTNB(Vector3f& T, Vector3f& N, Vector3f& B)
{
    T = safeNormalized(T, Vector3f(1, 0, 0));
    // Make N orthogonal to T
    N = N - Vector3f::dot(N, T) * T;
    N = safeNormalized(N, Vector3f(0, 1, 0));
    B = Vector3f::cross(T, N);
    B = safeNormalized(B, Vector3f(0, 0, 1));
    // Recompute N to remove accumulated error and ensure right-handed basis
    N = Vector3f::cross(B, T);
    N = safeNormalized(N, Vector3f(0, 1, 0));
}

inline void computeFrameFromDerivatives(const Vector3f& d1, const Vector3f& d2, bool havePrev,
    const CurvePoint& prev, Vector3f& outT, Vector3f& outN, Vector3f& outB)
{
    const float eps = 1e-12f;
    // Tangent from first derivative
    Vector3f T = d1;
    if (T.absSquared() < eps)
        T = havePrev ? prev.T : Vector3f(1, 0, 0);
    T.normalize();

    // Special-case: planar profile curves on the xy-plane.
    // Swept surface code assumes V/T/N all have z==0 for 2D profiles.
    const bool planarXY = (fabs(d1.z()) < 1e-8f) && (fabs(d2.z()) < 1e-8f) && (fabs(T.z()) < 1e-8f);
    if (planarXY) {
        Vector3f B = havePrev ? prev.B : Vector3f(0, 0, 1);
        if (B.absSquared() < eps || fabs(B.z()) < 0.5f)
            B = Vector3f(0, 0, 1);
        // Use a planar normal so that N.z == 0
        Vector3f N = Vector3f::cross(B, T);
        orthonormalizeTNB(T, N, B);
        if (havePrev && Vector3f::dot(N, prev.N) < 0) {
            N = -N;
            B = -B;
        }
        outT = T;
        outN = N;
        outB = B;
        return;
    }

    // Preferred: Frenet normal from curvature (second derivative projected off tangent)
    Vector3f N = d2 - Vector3f::dot(d2, T) * T;
    Vector3f B;

    if (N.absSquared() < eps) {
        // Near-straight segment or inflection: transport previous frame if possible
        if (havePrev && prev.B.absSquared() > eps) {
            B = prev.B;
            N = Vector3f::cross(B, T);
        } else {
            // Choose an arbitrary axis not parallel to T
            Vector3f a = (fabs(T.z()) < 0.9f) ? Vector3f(0, 0, 1) : Vector3f(0, 1, 0);
            B = Vector3f::cross(T, a);
            N = Vector3f::cross(B, T);
        }
    }

    orthonormalizeTNB(T, N, B);

    // Keep frame direction continuous (avoid random flips)
    if (havePrev && Vector3f::dot(N, prev.N) < 0) {
        N = -N;
        B = -B;
    }

    outT = T;
    outN = N;
    outB = B;
}

}

Curve evalBezier(const vector<Vector3f>& P, unsigned steps)
{
    // Check
    if (P.size() < 4 || P.size() % 3 != 1) {
        cerr << "evalBezier must be called with 3n+1 control points." << endl;
        return Curve();
    }
    if (steps == 0)
        steps = 1;

    // TODO:
    // You should implement this function so that it returns a Curve
    // (e.g., a vector< CurvePoint >).  The variable "steps" tells you
    // the number of points to generate on each piece of the spline.
    // At least, that's how the sample solution is implemented and how
    // the SWP files are written.  But you are free to interpret this
    // variable however you want, so long as you can control the
    // "resolution" of the discretized spline curve with it.

    // Make sure that this function computes all the appropriate
    // Vector3fs for each CurvePoint: V,T,N,B.
    // [NBT] should be unit and orthogonal.

    // Also note that you may assume that all Bezier curves that you
    // receive have G1 continuity.  Otherwise, the TNB will not be
    // be defined at points where this does not hold.

    const unsigned nPieces = (unsigned)((P.size() - 1) / 3);
    Curve result;
    result.reserve(nPieces * steps + 1);

    for (unsigned piece = 0; piece < nPieces; ++piece) {
        const Vector3f p0 = P[3 * piece + 0];
        const Vector3f p1 = P[3 * piece + 1];
        const Vector3f p2 = P[3 * piece + 2];
        const Vector3f p3 = P[3 * piece + 3];

        for (unsigned s = 0; s <= steps; ++s) {
            if (piece > 0 && s == 0)
                continue; // avoid duplicate endpoint

            const float t = float(s) / float(steps);
            const float u = 1.0f - t;

            // Position
            const Vector3f V
                = u * u * u * p0 + 3.0f * u * u * t * p1 + 3.0f * u * t * t * p2 + t * t * t * p3;

            // First derivative
            const Vector3f d1
                = 3.0f * u * u * (p1 - p0) + 6.0f * u * t * (p2 - p1) + 3.0f * t * t * (p3 - p2);

            // Second derivative
            const Vector3f d2 = 6.0f * u * (p2 - 2.0f * p1 + p0) + 6.0f * t * (p3 - 2.0f * p2 + p1);

            CurvePoint cp;
            cp.V = V;
            computeFrameFromDerivatives(
                d1, d2, !result.empty(), result.empty() ? cp : result.back(), cp.T, cp.N, cp.B);
            result.push_back(cp);
        }
    }

    return result;
}

Curve evalBspline(const vector<Vector3f>& P, unsigned steps)
{
    // Check
    if (P.size() < 4) {
        cerr << "evalBspline must be called with 4 or more control points." << endl;
        return Curve();
    }
    if (steps == 0)
        steps = 1;

    // TODO:
    // It is suggested that you implement this function by changing
    // basis from B-spline to Bezier.  That way, you can just call
    // your evalBezier function.

    const unsigned nPieces = (unsigned)(P.size() - 3);
    Curve result;
    result.reserve(nPieces * steps + 1);

    for (unsigned piece = 0; piece < nPieces; ++piece) {
        const Vector3f p0 = P[piece + 0];
        const Vector3f p1 = P[piece + 1];
        const Vector3f p2 = P[piece + 2];
        const Vector3f p3 = P[piece + 3];

        for (unsigned s = 0; s <= steps; ++s) {
            if (piece > 0 && s == 0)
                continue; // avoid duplicate endpoint

            const float u = float(s) / float(steps);
            const float um = 1.0f - u;

            // Uniform cubic B-spline basis on [0,1]
            const float N0 = (um * um * um) / 6.0f;
            const float N1 = (3.0f * u * u * u - 6.0f * u * u + 4.0f) / 6.0f;
            const float N2 = (-3.0f * u * u * u + 3.0f * u * u + 3.0f * u + 1.0f) / 6.0f;
            const float N3 = (u * u * u) / 6.0f;

            const Vector3f V = N0 * p0 + N1 * p1 + N2 * p2 + N3 * p3;

            // First derivative of basis
            const float dN0 = -(um * um) / 2.0f;
            const float dN1 = (3.0f * u * u - 4.0f * u) / 2.0f;
            const float dN2 = (-3.0f * u * u + 2.0f * u + 1.0f) / 2.0f;
            const float dN3 = (u * u) / 2.0f;
            const Vector3f d1 = dN0 * p0 + dN1 * p1 + dN2 * p2 + dN3 * p3;

            // Second derivative of basis
            const float ddN0 = um;
            const float ddN1 = 3.0f * u - 2.0f;
            const float ddN2 = -3.0f * u + 1.0f;
            const float ddN3 = u;
            const Vector3f d2 = ddN0 * p0 + ddN1 * p1 + ddN2 * p2 + ddN3 * p3;

            CurvePoint cp;
            cp.V = V;
            computeFrameFromDerivatives(
                d1, d2, !result.empty(), result.empty() ? cp : result.back(), cp.T, cp.N, cp.B);
            result.push_back(cp);
        }
    }

    return result;
}

Curve evalCircle(float radius, unsigned steps)
{
    // This is a sample function on how to properly initialize a Curve
    // (which is a vector< CurvePoint >).

    // Preallocate a curve with steps+1 CurvePoints
    Curve R(steps + 1);

    // Fill it in counterclockwise
    for (unsigned i = 0; i <= steps; ++i) {
        // step from 0 to 2pi
        float t = 2.0f * c_pi * float(i) / steps;

        // Initialize position
        // We're pivoting counterclockwise around the y-axis
        R[i].V = radius * Vector3f(cos(t), sin(t), 0);

        // Tangent vector is first derivative
        R[i].T = Vector3f(-sin(t), cos(t), 0);

        // Normal vector is second derivative
        R[i].N = Vector3f(-cos(t), -sin(t), 0);

        // Finally, binormal is facing up.
        R[i].B = Vector3f(0, 0, 1);
    }

    return R;
}

void recordCurve(const Curve& curve, VertexRecorder* recorder)
{
    const Vector3f WHITE(1, 1, 1);
    for (int i = 0; i < (int)curve.size() - 1; ++i) {
        recorder->record_poscolor(curve[i].V, WHITE);
        recorder->record_poscolor(curve[i + 1].V, WHITE);
    }
}
void recordCurveFrames(const Curve& curve, VertexRecorder* recorder, float framesize)
{
    Matrix4f T;
    const Vector3f RED(1, 0, 0);
    const Vector3f GREEN(0, 1, 0);
    const Vector3f BLUE(0, 0, 1);

    const Vector4f ORGN(0, 0, 0, 1);
    const Vector4f AXISX(framesize, 0, 0, 1);
    const Vector4f AXISY(0, framesize, 0, 1);
    const Vector4f AXISZ(0, 0, framesize, 1);

    for (int i = 0; i < (int)curve.size(); ++i) {
        T.setCol(0, Vector4f(curve[i].N, 0));
        T.setCol(1, Vector4f(curve[i].B, 0));
        T.setCol(2, Vector4f(curve[i].T, 0));
        T.setCol(3, Vector4f(curve[i].V, 1));

        // Transform orthogonal frames into model space
        Vector4f MORGN = T * ORGN;
        Vector4f MAXISX = T * AXISX;
        Vector4f MAXISY = T * AXISY;
        Vector4f MAXISZ = T * AXISZ;

        // Record in model space
        recorder->record_poscolor(MORGN.xyz(), RED);
        recorder->record_poscolor(MAXISX.xyz(), RED);

        recorder->record_poscolor(MORGN.xyz(), GREEN);
        recorder->record_poscolor(MAXISY.xyz(), GREEN);

        recorder->record_poscolor(MORGN.xyz(), BLUE);
        recorder->record_poscolor(MAXISZ.xyz(), BLUE);
    }
}
