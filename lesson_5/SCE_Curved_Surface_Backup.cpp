#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <vector>

#include "forces.h"

// Curved Substrate Parameters
const double surfacePar_a = 0.20;
const double surfacePar_b = 1.00;

double surfaceZ(double x, double y) {
    return surfacePar_a * std::cos(x / surfacePar_b);
}

int main() {
    const int N_membrane = 120;
    const int N_interior = 120;
    const double PI = 3.141592653589793;

    const double R_membrane = 1.0;
    const double R_interior_max = 0.95;

    // Euler motion parameters
    const int N_steps = 10000;
    const int output_every = 1000;

    // dt controls how strongly forces move the nodes each time step.
    const double dt = 0.0001;

    // Random motion parameters
    // Each node moves in a random direction with random length between 0 and 0.001.
    const double random_step_max = 0.00005;

    std::vector<double> membraneX, membraneY, membraneZ;     // Membrane node positions
    std::vector<double> interiorX, interiorY, interiorZ;     // Interior node positions

    std::ofstream membraneFile("membrane.dat");
    std::ofstream interiorFile("interior.dat");

    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> angleDist(0.0, 2.0 * PI);
    std::uniform_real_distribution<double> radiusDist(0.0, 1.0);
    std::uniform_real_distribution<double> randomLengthDist(0.0, random_step_max);

    // Create cell membrane
    for (int i = 0; i < N_membrane; i++) {
        double theta = 2.0 * PI * i / N_membrane;

        double x = R_membrane * std::cos(theta);
        double y = R_membrane * std::sin(theta);
        double z = surfaceZ(x, y);

        membraneX.push_back(x);
        membraneY.push_back(y);
        membraneZ.push_back(z);
    }

    // Create interior nodes
    for (int i = 0; i < N_interior; i++) {
        double theta = angleDist(gen);

        // Generate radius strictly inside the membrane
        double r = R_interior_max * std::sqrt(radiusDist(gen));

        double x = r * std::cos(theta);
        double y = r * std::sin(theta);
        double z = surfaceZ(x, y);

        interiorX.push_back(x);
        interiorY.push_back(y);
        interiorZ.push_back(z);
    }

    // Force vectors
    std::vector<double> membraneFx(N_membrane, 0.0);
    std::vector<double> membraneFy(N_membrane, 0.0);
    std::vector<double> membraneFz(N_membrane, 0.0);

    std::vector<double> interiorFx(N_interior, 0.0);
    std::vector<double> interiorFy(N_interior, 0.0);
    std::vector<double> interiorFz(N_interior, 0.0);

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

        // Euler update for membrane nodes
        for (int i = 0; i < N_membrane; i++) {
            double theta_random = -PI;
            double random_length = randomLengthDist(gen);

            double random_dx = random_length * std::cos(theta_random);
            double random_dy = random_length * std::sin(theta_random);
            double random_dz = 0.0;

            membraneX[i] += dt * membraneFx[i] + random_dx;
            membraneY[i] += dt * membraneFy[i] + random_dy;
            membraneZ[i] = surfaceZ(membraneX[i], membraneY[i]);
        }

        // Euler update for interior nodes
        for (int i = 0; i < N_interior; i++) {
            double theta_random = PI; // This actually doesnt have any effect on the code because random_length is 0!!!
            double random_length = 0;

            double random_dx = random_length * std::cos(theta_random);
            double random_dy = random_length * std::sin(theta_random);
            double random_dz = 0.0;

            interiorX[i] += dt * interiorFx[i] + random_dx;
            interiorY[i] += dt * interiorFy[i] + random_dy;
            interiorZ[i] = surfaceZ(interiorX[i], interiorY[i]);
        }

        if (step % output_every == 0) {
            writeFrame();
        }
    }
    membraneFile.close();
    interiorFile.close();

    std::cout << "Created membrane.dat and interior.dat\n";

    return 0;
}