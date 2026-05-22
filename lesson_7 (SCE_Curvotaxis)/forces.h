#ifndef FORCES_H
#define FORCES_H

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <stdexcept>

// -----------------------------------------------------------------------------
// forces.h
//
// Force utilities for the SIAM Workshop cell motility code.
//
// This header adapts the MATLAB force structure:
//   1. membrane-membrane Morse interaction, all-to-all
//   2. membrane-membrane spring interaction, nearest neighbors
//   3. membrane bending force using bending_force_fun logic
//   4. membrane-interior Morse interaction
//   5. interior-interior Morse interaction
//
// Assumptions:
// - Positions are stored in separate x, y, z vectors.
// - Forces are stored in separate Fx, Fy, Fz vectors.
// - Membrane nodes use periodic indexing.
// - Morse force sign follows the MATLAB convention:
//       f += z * (-dx) / r
// -----------------------------------------------------------------------------

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct MorseParameters {
    double U0 = 1.0;
    double V0 = 1.0;
    double k1 = 1.0;
    double k2 = 1.0;
    double cutoff = 1.0;
};

struct SpringParameters {
    double k = 1.0;
    double restLength = 0.1;
};

struct BendingParameters {
    double k = 1.0;
    double angleEq = 3.141592653589793;
};

struct AdhesionParams {
    int maxSubSitePerNode = 20;       // Maximum adhesion sites per membrane node

    double charUnbindDist = 0.1;      // Baseline characteristic distance for unbinding
    double kadh = 2.0;                // Adhesion spring constant

    // Unbinding asymmetry knobs.
    // Negative curvature stabilizes adhesions; positive curvature destabilizes them.
    double betaConc = 60.0;
    double betaConv = 0.02;

    // Curvature-dependent binding probability knobs.
    double pNeg = 0.0;
    double p0 = 0.5;
    double pPos = 1.0;
    double kappaPar = 0.5;
    double epsK = 1.0e-6;

    // Curvature-dependent lognormal adhesion length knobs.
    double meanAdhesionLength = 0.1;
    double gammaLen = 0.2;
    double lengthLognormalSigma = 0.25;
    double minMeanAdhesionLength = 1.0e-6;
};


// -----------------------------------------------------------------------------
// Default mechanical parameters from your MATLAB version.
// -----------------------------------------------------------------------------

constexpr double PI_FORCE = 3.141592653589793;

const SpringParameters DEFAULT_MEMBRANE_SPRING{
    300.0,     // k_spring
    0.03125    // L_spring
};

const BendingParameters DEFAULT_MEMBRANE_BENDING{
    1.0,       // k_angular
    PI_FORCE   // angle_eq; change this if your MATLAB code uses a different value
};

const MorseParameters DEFAULT_MI_MORSE{
    0.3125,    // U0_MI
    0.05,      // V0_MI
    0.05125,   // k1_MI
    0.625,     // k2_MI
    0.35625    // cut_off_MI
};

const MorseParameters DEFAULT_II_MORSE{
    0.188,     // U0_II
    0.146484,  // V0_II
    0.125,     // k1_II
    1.5625,    // k2_II
    0.3        // cut_off_II
};

const MorseParameters DEFAULT_MM_MORSE{
    0.0,       // U0_MM
    0.0,       // V0_MM; note: U0=V0=0 means no MM Morse force
    0.124,     // k1_MM
    0.625,     // k2_MM
    0.3635     // cut_off_MM
};

// -----------------------------------------------------------------------------
// Basic vector helpers
// -----------------------------------------------------------------------------

int periodicIndex(int i, int N) {
    return (i % N + N) % N;
}

Vec3 makeVec3(double x, double y, double z) {
    return Vec3{x, y, z};
}

Vec3 addVec3(const Vec3& a, const Vec3& b) {
    return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 subtractVec3(const Vec3& a, const Vec3& b) {
    return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 scaleVec3(const Vec3& v, double c) {
    return Vec3{c * v.x, c * v.y, c * v.z};
}

double dotVec3(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 crossVec3(const Vec3& a, const Vec3& b) {
    return Vec3{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

double normVec3(const Vec3& v) {
    return std::sqrt(dotVec3(v, v));
}

double clampValue(double value, double low, double high) {
    return std::max(low, std::min(value, high));
}

double distance3D(double xi, double yi, double zi,
                         double xj, double yj, double zj) {
    const double dx = xj - xi;
    const double dy = yj - yi;
    const double dz = zj - zi;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

inline double dist3D(const Vec3& a, const Vec3& b) {
    return distance3D(a.x, a.y, a.z, b.x, b.y, b.z);
}

inline double curvatureFromSurfaceDerivatives(double dzdx, double d2zdx2) {
    // MATLAB sign convention:
    //     kappa = -z''(x) / (1 + z'(x)^2)^(3/2)
    // For z = a*cos(x/b), peaks have positive curvature and valleys negative.
    return -d2zdx2 / std::pow(1.0 + dzdx * dzdx, 1.5);
}

template <typename SurfaceDzDxFunction, typename SurfaceD2zDx2Function>
inline double surfaceCurvatureFromCallbacks(
    double x,
    double y,
    SurfaceDzDxFunction surfaceDzDxFunc,
    SurfaceD2zDx2Function surfaceD2zDx2Func
) {
    const double dzdx = surfaceDzDxFunc(x, y);
    const double d2zdx2 = surfaceD2zDx2Func(x, y);
    return curvatureFromSurfaceDerivatives(dzdx, d2zdx2);
}

inline double medianAbsValue(std::vector<double> values) {
    if (values.empty()) {
        return 0.0;
    }

    for (double& v : values) {
        v = std::abs(v);
    }

    std::sort(values.begin(), values.end());

    const std::size_t n = values.size();
    if (n % 2 == 1) {
        return values[n / 2];
    }

    return 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

inline double curvatureModifiedUnbindDistance(
    double kappa,
    double kappaScale,
    const AdhesionParams& adh
) {
    kappaScale = std::max(kappaScale, adh.epsK);

    const double r = std::abs(kappa) / (kappaScale + std::abs(kappa));
    double mult = 1.0;

    if (kappa < 0.0) {
        mult = 1.0 + adh.betaConc * r;   // concave/valley: stronger, longer-lived
    } else {
        mult = 1.0 - adh.betaConv * r;   // convex/peak: weaker, shorter-lived
    }

    return std::max(adh.minMeanAdhesionLength, adh.charUnbindDist * mult);
}

inline double curvatureModifiedBindThreshold(
    double kappa,
    double kappaScale,
    const AdhesionParams& adh
) {
    const double k0 = std::max(adh.kappaPar * kappaScale, adh.epsK);

    if (std::abs(kappa) < adh.epsK) {
        kappa = 0.0;
    }

    const double r = std::abs(kappa) / (k0 + std::abs(kappa));

    double p = adh.p0;
    if (kappa >= 0.0) {
        p = adh.p0 + (adh.pPos - adh.p0) * r;
    } else {
        p = adh.p0 - (adh.p0 - adh.pNeg) * r;
    }

    return clampValue(p, 0.0, 1.0);
}

template <typename RNG>
inline double curvatureModifiedAdhesionLength(
    double kappa,
    double kappaScale,
    const AdhesionParams& adh,
    RNG& rng
) {
    kappaScale = std::max(kappaScale, adh.epsK);
    const double s = kappa / (kappaScale + std::abs(kappa)); // signed saturating factor in [-1, 1]

    double meanLength = adh.meanAdhesionLength * (1.0 + adh.gammaLen * s);
    meanLength = std::max(meanLength, adh.minMeanAdhesionLength);

    const double sigma = adh.lengthLognormalSigma;
    const double mu = std::log(meanLength) - 0.5 * sigma * sigma;

    std::normal_distribution<double> normal01(0.0, 1.0);
    return std::exp(mu + sigma * normal01(rng));
}

inline double outwardMembraneBindAngle(
    int i,
    const std::vector<double>& membraneX,
    const std::vector<double>& membraneY
) {
    double centerX = 0.0;
    double centerY = 0.0;

    const int N = static_cast<int>(membraneX.size());
    for (int k = 0; k < N; ++k) {
        centerX += membraneX[k];
        centerY += membraneY[k];
    }

    centerX /= static_cast<double>(N);
    centerY /= static_cast<double>(N);

    const double vx = membraneX[i] - centerX;
    const double vy = membraneY[i] - centerY;

    if (vx == 0.0 && vy == 0.0) {
        return 0.0;
    }

    return std::atan2(vy, vx);
}

// -----------------------------------------------------------------------------
// Morse forces
// -----------------------------------------------------------------------------

double morseForceMagnitude(double r, const MorseParameters& params) {
    // MATLAB:
    // z = U0/k1 * exp(-r/k1) - V0/k2 * exp(-r/k2)
    return (params.U0 / params.k1) * std::exp(-r / params.k1)
         - (params.V0 / params.k2) * std::exp(-r / params.k2);
}

void addMorseForceBetweenNodes(
    int i,
    int j,
    const std::vector<double>& xA,
    const std::vector<double>& yA,
    const std::vector<double>& zA,
    const std::vector<double>& xB,
    const std::vector<double>& yB,
    const std::vector<double>& zB,
    std::vector<double>& FxA,
    std::vector<double>& FyA,
    std::vector<double>& FzA,
    const MorseParameters& params
) {
    const double dx = xB[j] - xA[i];
    const double dy = yB[j] - yA[i];
    const double dz = zB[j] - zA[i];

    // Your MATLAB force block used r2 = dx*dx + dy*dy.
    // Since this C++ code stores z, this uses the full 3D distance.
    const double r2 = dx * dx + dy * dy + dz * dz;
    const double cutoff2 = params.cutoff * params.cutoff;

    if (r2 > 0.0 && r2 < cutoff2) {
        const double r = std::sqrt(r2);
        const double F = morseForceMagnitude(r, params);
        const double invr = 1.0 / r;

        FxA[i] += F * (-dx) * invr;
        FyA[i] += F * (-dy) * invr;
        FzA[i] += F * (-dz) * invr;
    }
}

void addMembraneMembraneMorseForces(
    const std::vector<double>& membraneX,
    const std::vector<double>& membraneY,
    const std::vector<double>& membraneZ,
    std::vector<double>& membraneFx,
    std::vector<double>& membraneFy,
    std::vector<double>& membraneFz,
    const MorseParameters& params = DEFAULT_MM_MORSE
) {
    const int N = static_cast<int>(membraneX.size());

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (j == i) {
                continue;
            }

            addMorseForceBetweenNodes(
                i, j,
                membraneX, membraneY, membraneZ,
                membraneX, membraneY, membraneZ,
                membraneFx, membraneFy, membraneFz,
                params
            );
        }
    }
}

void addMembraneInteriorMorseForcesOnMembrane(
    const std::vector<double>& membraneX,
    const std::vector<double>& membraneY,
    const std::vector<double>& membraneZ,
    const std::vector<double>& interiorX,
    const std::vector<double>& interiorY,
    const std::vector<double>& interiorZ,
    std::vector<double>& membraneFx,
    std::vector<double>& membraneFy,
    std::vector<double>& membraneFz,
    const MorseParameters& params = DEFAULT_MI_MORSE
) {
    const int Nmembrane = static_cast<int>(membraneX.size());
    const int Ninterior = static_cast<int>(interiorX.size());

    for (int i = 0; i < Nmembrane; i++) {
        for (int j = 0; j < Ninterior; j++) {
            addMorseForceBetweenNodes(
                i, j,
                membraneX, membraneY, membraneZ,
                interiorX, interiorY, interiorZ,
                membraneFx, membraneFy, membraneFz,
                params
            );
        }
    }
}

// Backward-compatible shorter name.
void addMembraneInteriorMorseForces(
    const std::vector<double>& membraneX,
    const std::vector<double>& membraneY,
    const std::vector<double>& membraneZ,
    const std::vector<double>& interiorX,
    const std::vector<double>& interiorY,
    const std::vector<double>& interiorZ,
    std::vector<double>& membraneFx,
    std::vector<double>& membraneFy,
    std::vector<double>& membraneFz,
    const MorseParameters& params = DEFAULT_MI_MORSE
) {
    addMembraneInteriorMorseForcesOnMembrane(
        membraneX, membraneY, membraneZ,
        interiorX, interiorY, interiorZ,
        membraneFx, membraneFy, membraneFz,
        params
    );
}

void addMembraneInteriorMorseForcesOnInterior(
    const std::vector<double>& membraneX,
    const std::vector<double>& membraneY,
    const std::vector<double>& membraneZ,
    const std::vector<double>& interiorX,
    const std::vector<double>& interiorY,
    const std::vector<double>& interiorZ,
    std::vector<double>& interiorFx,
    std::vector<double>& interiorFy,
    std::vector<double>& interiorFz,
    const MorseParameters& params = DEFAULT_MI_MORSE
) {
    const int Ninterior = static_cast<int>(interiorX.size());
    const int Nmembrane = static_cast<int>(membraneX.size());

    for (int i = 0; i < Ninterior; i++) {
        for (int j = 0; j < Nmembrane; j++) {
            addMorseForceBetweenNodes(
                i, j,
                interiorX, interiorY, interiorZ,
                membraneX, membraneY, membraneZ,
                interiorFx, interiorFy, interiorFz,
                params
            );
        }
    }
}

void addInteriorInteriorMorseForces(
    const std::vector<double>& interiorX,
    const std::vector<double>& interiorY,
    const std::vector<double>& interiorZ,
    std::vector<double>& interiorFx,
    std::vector<double>& interiorFy,
    std::vector<double>& interiorFz,
    const MorseParameters& params = DEFAULT_II_MORSE
) {
    const int N = static_cast<int>(interiorX.size());

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (j == i) {
                continue;
            }

            addMorseForceBetweenNodes(
                i, j,
                interiorX, interiorY, interiorZ,
                interiorX, interiorY, interiorZ,
                interiorFx, interiorFy, interiorFz,
                params
            );
        }
    }
}

// -----------------------------------------------------------------------------
// Membrane spring forces
// -----------------------------------------------------------------------------

void addSpringForceBetweenMembraneNeighbors(
    int i,
    int j,
    const std::vector<double>& membraneX,
    const std::vector<double>& membraneY,
    const std::vector<double>& membraneZ,
    std::vector<double>& membraneFx,
    std::vector<double>& membraneFy,
    std::vector<double>& membraneFz,
    const SpringParameters& params
) {
    const double dx = membraneX[j] - membraneX[i];
    const double dy = membraneY[j] - membraneY[i];
    const double dz = membraneZ[j] - membraneZ[i];
    const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (distance > 0.0) {
        // MATLAB:
        // f_i += -k * (distance - L) * (x_i - x_j) / distance
        // This is equivalent to:
        // f_i +=  k * (distance - L) * (x_j - x_i) / distance
        const double F = params.k * (distance - params.restLength);
        const double invDistance = 1.0 / distance;

        membraneFx[i] += F * dx * invDistance;
        membraneFy[i] += F * dy * invDistance;
        membraneFz[i] += F * dz * invDistance;
    }
}

void addMembraneSpringForces(
    const std::vector<double>& membraneX,
    const std::vector<double>& membraneY,
    const std::vector<double>& membraneZ,
    std::vector<double>& membraneFx,
    std::vector<double>& membraneFy,
    std::vector<double>& membraneFz,
    const SpringParameters& params = DEFAULT_MEMBRANE_SPRING
) {
    const int N = static_cast<int>(membraneX.size());

    for (int i = 0; i < N; i++) {
        const int leftNeighbor = periodicIndex(i - 1, N);
        const int rightNeighbor = periodicIndex(i + 1, N);

        addSpringForceBetweenMembraneNeighbors(
            i, leftNeighbor,
            membraneX, membraneY, membraneZ,
            membraneFx, membraneFy, membraneFz,
            params
        );

        addSpringForceBetweenMembraneNeighbors(
            i, rightNeighbor,
            membraneX, membraneY, membraneZ,
            membraneFx, membraneFy, membraneFz,
            params
        );
    }
}

// -----------------------------------------------------------------------------
// Bending forces
// -----------------------------------------------------------------------------

double orientedAngleFromTriple(const Vec3& left,
                                      const Vec3& center,
                                      const Vec3& right) {
    const Vec3 leftvec = subtractVec3(left, center);
    const Vec3 rightvec = subtractVec3(right, center);

    const double leftNorm = normVec3(leftvec);
    const double rightNorm = normVec3(rightvec);

    if (leftNorm <= 0.0 || rightNorm <= 0.0) {
        return 0.0;
    }

    const double cosine = clampValue(dotVec3(leftvec, rightvec) / (leftNorm * rightNorm), -1.0, 1.0);
    double angle = std::acos(cosine);

    // MATLAB uses kvec = [0; 0; 1] and checks dot(cross(leftvec, rightvec), kvec).
    const double crossdotN = crossVec3(leftvec, rightvec).z;
    if (crossdotN < 0.0) {
        angle = 2.0 * PI_FORCE - angle;
    }

    return angle;
}

Vec3 bendingForceFun(const std::string& flag,
                            double kBend,
                            double myAngle,
                            double angleEq,
                            const Vec3& p1,
                            const Vec3& p2,
                            const Vec3& p3) {
    // Direct C++ translation of your MATLAB bending_force_fun.
    const double sinAngle = std::sin(myAngle);
    const double eps = 1e-12;

    if (std::abs(sinAngle) < eps) {
        return Vec3{0.0, 0.0, 0.0};
    }

    const double coef = kBend * (myAngle - angleEq) / std::abs(sinAngle);
    const double cosAngle = std::cos(myAngle);

    if (flag == "center") {
        // p1 = left, p2 = center, p3 = right
        const Vec3 CLvec = subtractVec3(p1, p2);
        const Vec3 CRvec = subtractVec3(p3, p2);
        const double CLvecLen = normVec3(CLvec);
        const double CRvecLen = normVec3(CRvec);

        if (CLvecLen <= eps || CRvecLen <= eps) {
            return Vec3{0.0, 0.0, 0.0};
        }

        const Vec3 termL1 = scaleVec3(CLvec, -1.0 / (CLvecLen * CRvecLen));
        const Vec3 termL2 = scaleVec3(CLvec, cosAngle / (CLvecLen * CLvecLen));
        const Vec3 termR1 = scaleVec3(CRvec, -1.0 / (CLvecLen * CRvecLen));
        const Vec3 termR2 = scaleVec3(CRvec, cosAngle / (CRvecLen * CRvecLen));

        return scaleVec3(addVec3(addVec3(termL1, termL2), addVec3(termR1, termR2)), coef);
    }

    if (flag == "left") {
        // p1 = leftleft, p2 = left, p3 = center
        const Vec3 LLvec = subtractVec3(p1, p2);
        const Vec3 CLvec = subtractVec3(p2, p3);
        const double LLvecLen = normVec3(LLvec);
        const double CLvecLen = normVec3(CLvec);

        if (LLvecLen <= eps || CLvecLen <= eps) {
            return Vec3{0.0, 0.0, 0.0};
        }

        const Vec3 term1 = scaleVec3(LLvec, 1.0 / (LLvecLen * CLvecLen));
        const Vec3 term2 = scaleVec3(CLvec, cosAngle / (CLvecLen * CLvecLen));

        return scaleVec3(addVec3(term1, term2), coef);
    }

    if (flag == "right") {
        // p1 = center, p2 = right, p3 = rightright
        const Vec3 CRvec = subtractVec3(p2, p1);
        const Vec3 RRvec = subtractVec3(p3, p2);
        const double CRvecLen = normVec3(CRvec);
        const double RRvecLen = normVec3(RRvec);

        if (CRvecLen <= eps || RRvecLen <= eps) {
            return Vec3{0.0, 0.0, 0.0};
        }

        const Vec3 term1 = scaleVec3(RRvec, 1.0 / (RRvecLen * CRvecLen));
        const Vec3 term2 = scaleVec3(CRvec, cosAngle / (CRvecLen * CRvecLen));

        return scaleVec3(addVec3(term1, term2), coef);
    }

    return Vec3{0.0, 0.0, 0.0};
}

void addMembraneBendingForces(
    const std::vector<double>& membraneX,
    const std::vector<double>& membraneY,
    const std::vector<double>& membraneZ,
    std::vector<double>& membraneFx,
    std::vector<double>& membraneFy,
    std::vector<double>& membraneFz,
    const BendingParameters& params = DEFAULT_MEMBRANE_BENDING
) {
    const int N = static_cast<int>(membraneX.size());

    for (int i = 0; i < N; i++) {
        const int im1 = periodicIndex(i - 1, N);
        const int im2 = periodicIndex(i - 2, N);
        const int ip1 = periodicIndex(i + 1, N);
        const int ip2 = periodicIndex(i + 2, N);

        const Vec3 p_im2 = makeVec3(membraneX[im2], membraneY[im2], membraneZ[im2]);
        const Vec3 p_im1 = makeVec3(membraneX[im1], membraneY[im1], membraneZ[im1]);
        const Vec3 p_i   = makeVec3(membraneX[i],   membraneY[i],   membraneZ[i]);
        const Vec3 p_ip1 = makeVec3(membraneX[ip1], membraneY[ip1], membraneZ[ip1]);
        const Vec3 p_ip2 = makeVec3(membraneX[ip2], membraneY[ip2], membraneZ[ip2]);

        const double angleCenter = orientedAngleFromTriple(p_im1, p_i,   p_ip1);
        const double angleLeft   = orientedAngleFromTriple(p_im2, p_im1, p_i);
        const double angleRight  = orientedAngleFromTriple(p_i,   p_ip1, p_ip2);

        const Vec3 bendingCenter = bendingForceFun("center", params.k, angleCenter, params.angleEq,
                                                   p_im1, p_i, p_ip1);
        const Vec3 bendingLeft = bendingForceFun("left", params.k, angleLeft, params.angleEq,
                                                 p_im2, p_im1, p_i);
        const Vec3 bendingRight = bendingForceFun("right", params.k, angleRight, params.angleEq,
                                                  p_i, p_ip1, p_ip2);

        Vec3 total = addVec3(addVec3(bendingCenter, bendingLeft), bendingRight);

        // Preserve the MATLAB sign logic as closely as possible:
        // the original code flips the final bending force if the last computed
        // crossdotN value, from the right angle, is negative.
        const Vec3 rightLeftVec = subtractVec3(p_i, p_ip1);
        const Vec3 rightRightVec = subtractVec3(p_ip2, p_ip1);
        const double lastCrossdotN = crossVec3(rightLeftVec, rightRightVec).z;

        if (lastCrossdotN < 0.0) {
            total = scaleVec3(total, -1.0);
        }

        membraneFx[i] += total.x;
        membraneFy[i] += total.y;
        membraneFz[i] += total.z;
    }
}

// -----------------------------------------------------------------------------
// General helpers
// -----------------------------------------------------------------------------

void resetForces(std::vector<double>& Fx,
                        std::vector<double>& Fy,
                        std::vector<double>& Fz) {
    std::fill(Fx.begin(), Fx.end(), 0.0);
    std::fill(Fy.begin(), Fy.end(), 0.0);
    std::fill(Fz.begin(), Fz.end(), 0.0);
}


// ------------------------------------------------------------
// 2D distance helper
// ------------------------------------------------------------
inline double dist2D(double x0, double y0, double x1, double y1) {
    double dx = x1 - x0;
    double dy = y1 - y0;
    return std::sqrt(dx * dx + dy * dy);
}

// ------------------------------------------------------------
// Count how many adhesion sites are currently attached
// for membrane node i
// ------------------------------------------------------------
inline int countAttachedSites(
    int i,
    const std::vector<bool>& isAttached,
    const AdhesionParams& adh
) {
    int count = 0;

    for (int j = 0; j < adh.maxSubSitePerNode; ++j) {
        int siteIndex = i * adh.maxSubSitePerNode + j;

        if (isAttached[siteIndex]) {
            count++;
        }
    }

    return count;
}

// ------------------------------------------------------------
// Initialize adhesion vectors
//
// Call this in your main code after you know the number of
// membrane nodes.
// ------------------------------------------------------------
inline void initializeAdhesionSites(
    int numMembraneNodes,
    std::vector<double>& bindX,
    std::vector<double>& bindY,
    std::vector<double>& bindZ,
    std::vector<bool>& isAttached,
    const AdhesionParams& adh
) {
    int totalSites = numMembraneNodes * adh.maxSubSitePerNode;

    bindX.assign(totalSites, 0.0);
    bindY.assign(totalSites, 0.0);
    bindZ.assign(totalSites, 0.0);
    isAttached.assign(totalSites, false);
}


// ------------------------------------------------------------
// Add adhesion forces to membrane nodes
//
// This function does three things:
// 1. Removes stretched adhesion sites with a distance-dependent
//    unbinding probability.
// 2. Creates new adhesion sites if a node has open sites.
// 3. Adds spring-like adhesion forces from bound sites.
//
// For a curved substrate, new adhesion sites are placed at
//
//     z = surfaceZFunc(x, y)
//
// Distances are computed using ordinary 3D Euclidean distance, not
// surface arclength. This is faster and is usually good enough for
// this workshop model.
//
// The force from an attached site is:
//
//     F = kadh * (bindPosition - nodePosition)
//
// so the adhesion site pulls the node toward the substrate
// anchor point.
// ------------------------------------------------------------
template <typename RNG, typename SurfaceZFunction, typename CurvatureFunction>
inline void addAdhesionForcesWithCurvature(
    const std::vector<double>& membraneX,
    const std::vector<double>& membraneY,
    const std::vector<double>& membraneZ,

    std::vector<double>& bindX,
    std::vector<double>& bindY,
    std::vector<double>& bindZ,
    std::vector<bool>& isAttached,

    std::vector<double>& forceX,
    std::vector<double>& forceY,
    std::vector<double>& forceZ,

    const AdhesionParams& adh,
    RNG& rng,
    SurfaceZFunction surfaceZFunc,
    CurvatureFunction curvatureFunc
) {

    int numMembraneNodes = static_cast<int>(membraneX.size());

    if (
        membraneY.size() != membraneX.size() ||
        membraneZ.size() != membraneX.size() ||
        forceX.size() != membraneX.size() ||
        forceY.size() != membraneX.size() ||
        forceZ.size() != membraneX.size()
    ) {
        throw std::runtime_error("Membrane position and force vectors must have the same size.");
    }

    const std::size_t expectedSites = static_cast<std::size_t>(numMembraneNodes) * static_cast<std::size_t>(adh.maxSubSitePerNode);

    if (
        bindX.size() != expectedSites ||
        bindY.size() != expectedSites ||
        bindZ.size() != expectedSites ||
        isAttached.size() != expectedSites
    ) {
        throw std::runtime_error("Adhesion vectors have the wrong size. Call initializeAdhesionSites first.");
    }

    if (adh.charUnbindDist <= 0.0) {
        throw std::runtime_error("charUnbindDist must be positive.");
    }

    std::uniform_real_distribution<double> uniform01(0.0, 1.0);

    std::vector<double> curvature(numMembraneNodes, 0.0);
    for (int i = 0; i < numMembraneNodes; ++i) {
        curvature[i] = curvatureFunc(membraneX[i], membraneY[i]);
    }

    const double kappaScale = std::max(medianAbsValue(curvature), adh.epsK);

    for (int i = 0; i < numMembraneNodes; ++i) {

        const double kappa_i = curvature[i];

        // ----------------------------------------------------
        // 1. Unbinding logic
        // ----------------------------------------------------
        const double charUnbindDist_i =
            curvatureModifiedUnbindDistance(kappa_i, kappaScale, adh);

        for (int j = 0; j < adh.maxSubSitePerNode; ++j) {
            int siteIndex = i * adh.maxSubSitePerNode + j;

            if (isAttached[siteIndex]) {
                double distNodeSite = distance3D(
                    membraneX[i],
                    membraneY[i],
                    membraneZ[i],
                    bindX[siteIndex],
                    bindY[siteIndex],
                    bindZ[siteIndex]
                );

                double unbindProb = 1.0 - std::exp(-distNodeSite / charUnbindDist_i);

                if (uniform01(rng) < unbindProb) {
                    isAttached[siteIndex] = false;
                }
            }
        }

        // Recount after unbinding
        int bindSiteCount = countAttachedSites(i, isAttached, adh);

        // ----------------------------------------------------
        // 2. Binding logic
        // ----------------------------------------------------
        const double siteBindThreshold_i =
            curvatureModifiedBindThreshold(kappa_i, kappaScale, adh);

        if (bindSiteCount < adh.maxSubSitePerNode) {
            int openSites = adh.maxSubSitePerNode - bindSiteCount;

            for (int attempt = 0; attempt < openSites; ++attempt) {
                if (uniform01(rng) < siteBindThreshold_i) {

                    const double randAngle = outwardMembraneBindAngle(i, membraneX, membraneY);
                    const double randLen =
                        curvatureModifiedAdhesionLength(kappa_i, kappaScale, adh, rng);

                    const double x0 = membraneX[i];
                    const double y0 = membraneY[i];

                    // Performance change from MATLAB:
                    // use ordinary Euclidean xy offsets instead of converting
                    // through an arclength coordinate.
                    const double nodeSiteX = x0 + randLen * std::cos(randAngle);
                    const double nodeSiteY = y0 + randLen * std::sin(randAngle);
                    const double nodeSiteZ = surfaceZFunc(nodeSiteX, nodeSiteY);

                    // Store the new adhesion site in the first open slot
                    for (int j = 0; j < adh.maxSubSitePerNode; ++j) {
                        int siteIndex = i * adh.maxSubSitePerNode + j;

                        if (!isAttached[siteIndex]) {
                            bindX[siteIndex] = nodeSiteX;
                            bindY[siteIndex] = nodeSiteY;
                            bindZ[siteIndex] = nodeSiteZ;

                            isAttached[siteIndex] = true;
                            break;
                        }
                    }
                }
            }
        }

        // ----------------------------------------------------
        // 3. Adhesion force logic
        // ----------------------------------------------------
        for (int j = 0; j < adh.maxSubSitePerNode; ++j) {
            int siteIndex = i * adh.maxSubSitePerNode + j;

            if (isAttached[siteIndex]) {
                double distNodeSite = distance3D(
                    membraneX[i],
                    membraneY[i],
                    membraneZ[i],
                    bindX[siteIndex],
                    bindY[siteIndex],
                    bindZ[siteIndex]
                );

                if (distNodeSite > 0.0) {
                    forceX[i] += adh.kadh * (bindX[siteIndex] - membraneX[i]);
                    forceY[i] += adh.kadh * (bindY[siteIndex] - membraneY[i]);
                    forceZ[i] += adh.kadh * (bindZ[siteIndex] - membraneZ[i]);
                }
            }
        }
    }
}

// Preferred curved-substrate overload when analytic derivatives are available.
// Pass surfaceZFunc(x,y), dz/dx, and d2z/dx2 from your .cpp file.
template <typename RNG, typename SurfaceZFunction, typename SurfaceDzDxFunction, typename SurfaceD2zDx2Function>
inline void addAdhesionForces(
    const std::vector<double>& membraneX,
    const std::vector<double>& membraneY,
    const std::vector<double>& membraneZ,

    std::vector<double>& bindX,
    std::vector<double>& bindY,
    std::vector<double>& bindZ,
    std::vector<bool>& isAttached,

    std::vector<double>& forceX,
    std::vector<double>& forceY,
    std::vector<double>& forceZ,

    const AdhesionParams& adh,
    RNG& rng,
    SurfaceZFunction surfaceZFunc,
    SurfaceDzDxFunction surfaceDzDxFunc,
    SurfaceD2zDx2Function surfaceD2zDx2Func
) {
    auto curvatureFunc = [&](double x, double y) {
        return surfaceCurvatureFromCallbacks(x, y, surfaceDzDxFunc, surfaceD2zDx2Func);
    };

    addAdhesionForcesWithCurvature(
        membraneX, membraneY, membraneZ,
        bindX, bindY, bindZ, isAttached,
        forceX, forceY, forceZ,
        adh, rng,
        surfaceZFunc, curvatureFunc
    );
}


#endif
