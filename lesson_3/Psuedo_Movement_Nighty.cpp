#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <vector>

int main() {

    // Membrane Vars
    const int N_membrane = 120;
    const int N_interior = 120;
    const double PI = 3.141592653589793;

    const double R_membrane = 1.0;
    const double R_interior_max = 0.95;

    // Time Step Vars; Can change vars if wanted
    const int N_steps = 1000;
    const int output_every = 100;
    const double dx = 0.001;
    const double dy = 0.001;
    const double dz = 0.001;

    // Vec Vars
    std::vector<double> membraneX;
    std::vector<double> membraneY;
    std::vector<double> membraneZ;

    std::vector<double> interiorX;
    std::vector<double> interiorY;
    std::vector<double> interiorZ;

    // Output File Stream Vars + Checking

    std::ofstream membraneFile("membrane.dat");
    std::ofstream interiorFile("interior.dat");

    if ((!membraneFile.is_open()) || (!interiorFile.is_open())) {
        if (!membraneFile.is_open()) {
            std::cout << "membrane.dat has failed to open" << std::endl;
        }

        if (!interiorFile.is_open()) {
            std::cout << "interior.dat has failed to open" << std::endl;
        }

        exit(EXIT_FAILURE);
    }

    // Random number generator for interior nodes
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> angleDist(0.0, 2.0 * PI);
    std::uniform_real_distribution<double> radiusDist(0.0, 1.0);

    // Create cell membrane
    for (int i = 0; i < N_membrane; i++) {
        double theta = 2.0 * PI * i / N_membrane;

        double x = R_membrane * std::cos(theta);
        double y = R_membrane * std::sin(theta);
        double z = 0.0;

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
        double z = 0.0;

        interiorX.push_back(x);
        interiorY.push_back(y);
        interiorZ.push_back(z);
    }

    // Time Step
    for (int step = 0; step < N_steps; step++) {

        // Updates membrane node
        for (int i = 0; i < N_membrane; i++) {
            membraneX.at(i) += dx;
            membraneY.at(i) += dy;
            membraneZ.at(i) += dz;
        }

        // Updates interior node
        for (int i = 0; i < N_interior; i++) {
            interiorX.at(i) += dx;
            interiorY.at(i) += dy;
            interiorZ.at(i) += dz;
        }

        // Writes to file
        if (step % output_every == 0) {
            for (int i = 0; i < N_membrane; i++) {
                membraneFile << membraneX.at(i) << " " 
                << membraneY.at(i) << " " << membraneZ.at(i) << "\n";
            }
            membraneFile << "\n\n";

            for (int i = 0; i < N_interior; i++) {
                interiorFile << interiorX.at(i) << " " 
                << interiorY.at(i) << " " << interiorZ.at(i) << "\n";
            }
            interiorFile << "\n\n";
        }
    }

    // Close ofstream Files

    membraneFile.close();
    interiorFile.close();

    std::cout << "Created membrane.dat and interior.dat\n";

    return 0;
}