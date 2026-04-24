#include <fstream>
#include <iostream>
#include <cmath>

int main() {
    const int N = 120;
    const double PI = 3.141592653589793;

    std::ofstream given_file("circle.dat");

    /*if (given_file.is_open() == false) {
        std::cout << given_file << " does not exist or did not successfully open!" << std::endl;
        exit(EXIT_FAILURE);
    }*/

    for (int i = 0; i < N; i++) {
        double theta = 2.0 * PI * i / N;

        //Function to Change//
        double R = 0.5; // FIXME: Add data variable?
        //-----------------//

        double x = R * std::cos(theta);
        double y = R * std::sin(theta);
        double z = 0;

        given_file << x << " " << y << " " << z << "\n";
    }

    given_file.close();

    return 0;
}