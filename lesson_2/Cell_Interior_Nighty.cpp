#include <fstream>
#include <iostream>
#include <cmath>
#include <random>

int main() {
    const int N_membrane = 120;
    const int N_interior = 120;

    const double PI = 3.141592653589793;
    const double R_membrane = 1.0;
    const double R_interior_max = 0.95;

    // std::ofstream file("circle.dat");
    std::ofstream membrane_file;
    std::ofstream interior_file;

    membrane_file.open("membrane.dat");
    interior_file.open("interior.dat");

    if ((!membrane_file.is_open()) || (!interior_file.is_open())) {
        if (!membrane_file.is_open()) {
            std::cout << "membrane.dat has failed to open" << std::endl;
        }

        if (!interior_file.is_open()) {
            std::cout << "interior.dat has failed to open" << std::endl;
        }

        exit(EXIT_FAILURE);
    }

    // Random number generator for interior nodes
    std::random_device rd; // This creates a source of randomness used to generate a seed
    std::mt19937 gen(rd()); // This creates the actual random number generator. The generator is called `gen`, and it is initialized using the seed from `rd()`
    std::uniform_real_distribution<double> angleDist(0.0, 2.0 * PI); // This creates a distribution for the angle `theta`, so each angle between `0` and `2π` is equally likely
    std::uniform_real_distribution<double> radiusDist(0.0, 1.0); // This creates a distribution for a number between `0` and `1`, which we will use to generate the radius.


    // Create membrane Nodes
    for (int i = 0; i < N_membrane; i++) {
        double theta = 2.0 * PI * i / N_membrane;

        //Function to Change//
        // double R = 1;
        // const double R_membrane = 1.0;
        //-----------------//

        double x = R_membrane * std::cos(theta);
        double y = R_membrane * std::sin(theta);
        double z = 0;

        membrane_file << x << " " << y << " " << z << std::endl;
    }

    // Create interior nodes
    for (int i = 0; i < N_interior; i++) {
        double theta = angleDist(gen);

        // Generate radius strictly inside the membrane
        double r_interior_pt = R_interior_max * std::sqrt(radiusDist(gen));

        double x = r_interior_pt * std::cos(theta);
        double y = r_interior_pt * std::sin(theta);
        double z = 0.0;

        interior_file << x << " " << y << " " << z << std::endl;
    }

    membrane_file.close();
    interior_file.close();

    std::cout << "membrane.dat and interior.dat have been created!" << std::endl;

    return 0;
}