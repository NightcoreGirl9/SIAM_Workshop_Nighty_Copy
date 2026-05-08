#ifndef FORCES_H
#define FORCES_H

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

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

#endif
