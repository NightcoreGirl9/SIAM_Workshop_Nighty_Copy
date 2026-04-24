#include <fstream>
#include <iostream>
#include <cmath>

int main() {
    const int N = 120;
    const double PI = 3.141592653589793;

    std::ofstream file("circle.dat");

    for (int i = 0; i < N; i++) {
        double theta = 2.0 * PI * i / N;

        //Function to Change//
        double R =1;
        //-----------------//

        double x = R * std::cos(theta);
        double y = R * std::sin(theta);
        double z = 0;

        file << x << " " << y << " " << z << "\n";
    }

    file.close();

    return 0;
}