#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <chrono>

#include "forces.h"

// Curved substrate parameters
const double surfacePar_a = 0.1;
const double surfacePar_b = 3.141592653589793 / 6.0;

// Curved substrate function
double surfaceZ(double x, double y) {
    (void)y; // What does this do? Our func doesnt use y parameter so we jus get rid of it
    return surfacePar_a * std::cos(x / surfacePar_b);
}

// First derivative dz/dx for z = a*cos(x/b).
double surfaceDzDx(double x, double y) {
    (void)y;
    return -(surfacePar_a / surfacePar_b) * std::sin(x / surfacePar_b);
}

// Second derivative d2z/dx2 for z = a*cos(x/b).
double surfaceD2zDx2(double x, double y) {
    (void)y;
    return -(surfacePar_a / (surfacePar_b * surfacePar_b)) * std::cos(x / surfacePar_b);
}

int main() {
    // Time Keeping Code
    auto startTime = std::chrono::high_resolution_clock::now();

    const int N_membrane = 60;
    const int N_interior = 60;
    const double PI = 3.141592653589793;

    const double R_membrane = .5;
    const double R_interior_max = 0.4;

    // Euler motion parameters
    const int N_steps = 100000;
    const int output_every = 1000;

    // dt controls how strongly forces move the nodes each time step.
    const double dt = 0.0005;

    // Node Position Vectors
    std::vector<double> membraneX, membraneY, membraneZ;
    std::vector<double> interiorX, interiorY, interiorZ;

    //Bind Position Vectors
    std::vector<double> bindX, bindY, bindZ;
    std::vector<bool> isAttached;     

    // Force vectors
    std::vector<double> membraneFx(N_membrane, 0.0);
    std::vector<double> membraneFy(N_membrane, 0.0);
    std::vector<double> membraneFz(N_membrane, 0.0);

    std::vector<double> interiorFx(N_interior, 0.0);
    std::vector<double> interiorFy(N_interior, 0.0);
    std::vector<double> interiorFz(N_interior, 0.0);

    std::ofstream membraneFile("membrane.dat");
    std::ofstream interiorFile("interior.dat");

    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> angleDist(0.0, 2.0 * PI);
    std::uniform_real_distribution<double> radiusDist(0.0, 1.0);

    // Cell center: choose this off-center to break symmetry
    const double x_c = PI*PI / 12.0;   // change as needed
    const double y_c = 0.0;   // change as needed

    // Create cell membrane around (x_c, y_c)
    for (int i = 0; i < N_membrane; i++) {
        double theta = 2.0 * PI * i / N_membrane;

        double x = x_c + R_membrane * std::cos(theta);
        double y = y_c + R_membrane * std::sin(theta);
        double z = surfaceZ(x, y);

        membraneX.push_back(x);
        membraneY.push_back(y);
        membraneZ.push_back(z);
    }

    // Create interior nodes around (x_c, y_c)
    for (int i = 0; i < N_interior; i++) {
        double theta = angleDist(gen);

        // Generate radius strictly inside the membrane
        double r = R_interior_max * std::sqrt(radiusDist(gen));

        double x = x_c + r * std::cos(theta);
        double y = y_c + r * std::sin(theta);
        double z = surfaceZ(x, y);

        interiorX.push_back(x);
        interiorY.push_back(y);
        interiorZ.push_back(z);
    }
    // Helper function for writing one animation frame
    auto writeFrame = [&]() {
        for (int i = 0; i < N_membrane; i++) {
            membraneFile << membraneX[i] << " "
                         << membraneY[i] << " "
                         << membraneZ[i] << "\n";
        }
        membraneFile << "\n\n";

        for (int i = 0; i < N_interior; i++) {
            interiorFile << interiorX[i] << " "
                         << interiorY[i] << " "
                         << interiorZ[i] << "\n";
        }
        interiorFile << "\n\n";
    };

    // Write initial cell position
    writeFrame();

    //Initialize bind site parameters
    AdhesionParams adh;

    adh.maxSubSitePerNode = 20;
    adh.charUnbindDist = 0.1;
    adh.meanAdhesionLength = 0.1;
    adh.kadh = 2.0;

    // Curvature-dependent unbinding parameters
    adh.betaConc = 60.0;
    adh.betaConv = 0.02;

    // Curvature-dependent binding probability parameters
    // pNeg: target binding probability for strong negative curvature
    // p0:   binding probability near flat regions
    // pPos: target binding probability for strong positive curvature
    adh.pNeg = 0.0;
    adh.p0   = 0.5;
    adh.pPos = 1.0;

    // Curvature scale tuning
    // Larger kappaPar means curvature effects saturate more slowly
    adh.kappaPar = .5;

    // Curvature-dependent adhesion length parameters
    // gammaLen controls how strongly curvature modifies adhesion length
    adh.gammaLen = 0.2;


    initializeAdhesionSites(
        N_membrane,
        bindX,
        bindY,
        bindZ,
        isAttached,
        adh
    );
    // Time stepping loop
    for (int step = 1; step <= N_steps; step++) {

        // Reset all forces to zero at the start of each time step
        resetForces(membraneFx, membraneFy, membraneFz);
        resetForces(interiorFx, interiorFy, interiorFz);

        // Add membrane forces
        addMembraneMembraneMorseForces(
            membraneX, membraneY, membraneZ,
            membraneFx, membraneFy, membraneFz
        );

        addMembraneSpringForces(
            membraneX, membraneY, membraneZ,
            membraneFx, membraneFy, membraneFz
        );

        addMembraneBendingForces(
            membraneX, membraneY, membraneZ,
            membraneFx, membraneFy, membraneFz
        );

        // Add membrane-interior forces on the membrane nodes
        addMembraneInteriorMorseForces(
            membraneX, membraneY, membraneZ,
            interiorX, interiorY, interiorZ,
            membraneFx, membraneFy, membraneFz
        );

        // Add membrane-interior forces on the interior nodes
        addMembraneInteriorMorseForcesOnInterior(
            membraneX, membraneY, membraneZ,
            interiorX, interiorY, interiorZ,
            interiorFx, interiorFy, interiorFz
        );

        // Add interior-interior forces
        addInteriorInteriorMorseForces(
            interiorX, interiorY, interiorZ,
            interiorFx, interiorFy, interiorFz
        );

        // Add Substrate Adhesion Forces
        addAdhesionForces(
            membraneX,
            membraneY,
            membraneZ,

            bindX,
            bindY,
            bindZ,
            isAttached,

            membraneFx,
            membraneFy,
            membraneFz,

            adh,
            gen,
            surfaceZ,
            surfaceDzDx,
            surfaceD2zDx2
        );

        // Euler update for membrane nodes
        for (int i = 0; i < N_membrane; i++) {
            membraneX[i] += dt * membraneFx[i];
            membraneY[i] += dt * membraneFy[i];
            membraneZ[i] = surfaceZ(membraneX[i], membraneY[i]);
        }

        // Euler update for interior nodes
        for (int i = 0; i < N_interior; i++) {
            interiorX[i] += dt * interiorFx[i];
            interiorY[i] += dt * interiorFy[i];
            interiorZ[i] = surfaceZ(interiorX[i], interiorY[i]);
        }

        if (step % output_every == 0) {
            writeFrame();
        }
    }
    membraneFile.close();
    interiorFile.close();

    std::cout << "Created membrane.dat and interior.dat\n";

    auto endTime = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsedSeconds = endTime - startTime;

    std::cout << "Runtime: " << elapsedSeconds.count() / 60.0 << " minutes\n";
    return 0;
}