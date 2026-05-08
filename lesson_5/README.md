# Lesson 5: Moving the Cell on a Curved Substrate

In this lesson, we move our pseudo-cell model from a flat surface to a **curved substrate**.

Up to this point, all node positions had `z = 0`, so the cell lived on a flat plane. Now we will place the cell on the curved surface

```math
z = a\cos\left(\frac{x}{b}\right)
```

and study how the cell behaves on a non-flat geometry.

---

## Main Goals

In this lesson, we will:

1. Reduce randomness in the model so the motion is easier to understand.
2. Replace the flat surface `z = 0` with a curved surface.
3. Update each node so that after moving in `x` and `y`, its `z`-value is recomputed from the surface.
4. Visualize the motion in **3D** using **gnuplot**.
5. Test other surfaces by changing the surface function in both C++ and gnuplot.

---

## Mathematical Idea

We define the substrate by

```math
z = a\cos\left(\frac{x}{b}\right).
```

### Meaning of the parameters

- `a` controls the **height** of the hills and valleys.
- `b` controls the **horizontal spacing** of the hills.

In our code, we use:

```cpp
const double surfacePar_a = 0.20;
const double surfacePar_b = 1.00;
```

so the substrate is

```math
z = 0.2\cos(x).
```

This surface varies only in the `x`-direction, so it looks like ridges extending along the `y`-direction.

---

## Simplification for the Workshop

On a truly curved surface, distances along the substrate should be computed using **arc length**.

However, to keep the code simple and fast for this workshop, we will **not** compute arc length. Instead, all forces and distances are still computed using the usual **Euclidean distance** in 3D.

So this lesson is mainly about:

- placing the cell on a curved surface,
- updating the node heights correctly,
- and visualizing the result.

---

## Main Changes from the Previous Lesson

### 1. Add curved substrate parameters

Near the top of the code, we define:

```cpp
const double surfacePar_a = 0.20;
const double surfacePar_b = 1.00;
```

These are the parameters in

```math
z = a\cos\left(\frac{x}{b}\right).
```

---

### 2. Add a curved surface function

We define a function that returns the substrate height:

```cpp
double surfaceZ(double x, double y) {
    return surfacePar_a * std::cos(x / surfacePar_b);
}
```

Even though this function takes both `x` and `y`, the current surface depends only on `x`.

We include `y` so it is easy to replace this later with a more general surface:

```math
z = f(x,y).
```

---

### 3. Initialize nodes on the curved surface

Previously, each node was created with:

```cpp
double z = 0.0;
```

Now we use:

```cpp
double z = surfaceZ(x, y);
```

so every membrane node and interior node starts on the curved substrate.

---

### 4. Project nodes back onto the surface after each update

After updating `x` and `y`, we do **not** update `z` independently. Instead, we recompute it from the surface:

```cpp
membraneZ[i] = surfaceZ(membraneX[i], membraneY[i]);
```

and similarly for interior nodes:

```cpp
interiorZ[i] = surfaceZ(interiorX[i], interiorY[i]);
```

This keeps the cell constrained to the curved substrate.

---

## Reducing Randomness to Test the Curved Surface

Before adding more complicated motion on a curved substrate, we take a step back and simplify the movement.

The goal is to make the behavior easier to interpret.

### Membrane nodes

For membrane nodes, we use a very small random step length:

```cpp
const double random_step_max = 0.0001;
```

Then, instead of choosing a random direction, we force the membrane nodes to move in a fixed direction:

```cpp
double theta_random = -PI;
```

Since

```math
\cos(-\pi) = -1
```

and

```math
\sin(-\pi) = 0,
```

this makes the membrane's random step point in the negative `x`-direction.

So the random movement term becomes approximately:

```math
(random\_dx, random\_dy) = (-random\_length, 0).
```

This gives the cell a small consistent push in one direction.

---

### Interior nodes

For interior nodes, we set:

```cpp
double random_length = 0.0;
```

This means the interior nodes receive **no direct random movement**.

However, the interior nodes can still move because they interact with the rest of the cell through the force terms. In particular, the code still computes:

- membrane-interior forces acting on interior nodes,
- interior-interior forces acting on interior nodes.

This is an important modeling point:

> Even if the interior nodes are not directly pushed by random motion, the whole cell can still move because the membrane pulls the interior along through cell-cell force interactions.


---

## Files in This Lesson

This lesson uses:

- the C++ simulation code,
- `forces.h` for force calculations,
- and a gnuplot script to animate the results.

The simulation writes two output files:

- `membrane.dat`
- `interior.dat`

These files are then plotted in gnuplot.

---

## Compiling and Running

Compile the C++ file as usual.

For example, with `g++`:

```bash
g++ -std=c++11 SCE_Curved_Surface.cpp -o lesson5_cell
./lesson5_cell
```

After running, you should get:

- `membrane.dat`
- `interior.dat`

and the program should print:

```text
Created membrane.dat and interior.dat
```
---

## Gnuplot Script Explanation

This lesson uses `splot` instead of `plot` because we are now plotting in 3D.

### Main changes from the 2D script

#### Add a `z` label and `z` range

```gnuplot
set zlabel "z"
set zrange [-0.4:0.4]
```

Since the cell is no longer on the plane `z = 0`, we need to show the vertical direction.

---

#### Place the xy-plane at `z = 0`

```gnuplot
set xyplane at 0
```

This makes the reference plane appear at height `z = 0`.

---

#### Set the viewing angle

```gnuplot
set view 70, 225
```

The first number controls the elevation angle. The second number controls the rotation around the vertical axis.

If the cell is hidden behind the surface, changing this line is usually the easiest fix.

---

#### Define the same surface in gnuplot

```gnuplot
surface(x,y) = 0.2*cos(x)
```

This must match the surface in the C++ code.

In C++, the surface is:

```cpp
return surfacePar_a * std::cos(x / surfacePar_b);
```

Since `surfacePar_a = 0.20` and `surfacePar_b = 1.00`, the matching gnuplot surface is:

```gnuplot
surface(x,y) = 0.2*cos(x)
```

---

#### Use `splot` to plot the surface and the cell

```gnuplot
splot \
    surface(x,y) with pm3d, \
    'membrane.dat' index i using 1:2:($3 + 0.02) with points pt 7 ps 1.4 lc rgb "black", \
    'interior.dat' index i using 1:2:($3 + 0.02) with points pt 7 ps 1.0 lc rgb "red"
```

The surface is plotted with `pm3d`, while the membrane and interior nodes are plotted as points.

---

## Why Do We Plot `($3 + 0.02)`?

In the gnuplot script, we use:

```gnuplot
using 1:2:($3 + 0.02)
```

instead of just:

```gnuplot
using 1:2:3
```

This slightly lifts the plotted cell points above the surface so they are easier to see.

It is only a **visualization trick**. It does **not** change the simulation itself.

The actual simulation data still lies on:

```math
z = 0.2\cos(x).
```
---

## Full C++ Code

```cpp
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <vector>

#include "forces.h"

// Curved substrate parameters
const double surfacePar_a = 0.20;
const double surfacePar_b = 1.00;

// Curved substrate function
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
    const double random_step_max = 0.0001;

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
            // Move membrane nodes in a fixed direction.
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
            // Interior nodes do not receive direct random motion.
            // They move because they are connected to the rest of the cell through forces.
            double theta_random = PI;
            double random_length = 0.0;

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
```

---

## Full Gnuplot Script

```gnuplot
set terminal gif animate delay 50 size 800,800
set output "pseudo_movement.gif"

set title "Pseudo Cell Movement on Curved Surface"
set xlabel "x"
set ylabel "y"
set zlabel "z"

set xrange [-4:4]
set yrange [-4:4]
set zrange [-0.4:0.4]

set grid
set border 4095
unset key

set xyplane at 0
set view 70, 225

set samples 80
set isosamples 80

surface(x,y) = 0.2*cos(x)

# Gradient surface
set pm3d at s
set palette rgbformulae 22,13,-31
set style fill transparent solid 0.35

Nsteps = 10

do for [i=0:Nsteps-1] {
    set title sprintf("Pseudo Cell Movement on Curved Surface - Frame %d", i)

    splot \
        surface(x,y) with pm3d, \
        'membrane.dat' index i using 1:2:($3 + 0.02) with points pt 7 ps 1.4 lc rgb "black", \
        'interior.dat' index i using 1:2:($3 + 0.02) with points pt 7 ps 1.0 lc rgb "red"
}

set output
```

---

## Exercises: Try Different Surfaces

For each exercise, change the surface function in **both places**:

1. In the C++ code:

```cpp
double surfaceZ(double x, double y) {
    return /* new surface */;
}
```

2. In the gnuplot script:

```gnuplot
surface(x,y) = # matching gnuplot version
```

The C++ and gnuplot surfaces must match. Otherwise, your data and your plotted surface will not agree.

---

### Exercise 1: Try different parameter values for `a` and `b`

Current surface:

```math
z = a\cos\left(\frac{x}{b}\right)
```

Try changing:

```cpp
const double surfacePar_a = 0.10;
const double surfacePar_b = 1.00;
```

or:

```cpp
const double surfacePar_a = 0.30;
const double surfacePar_b = 2.00;
```

Matching gnuplot examples:

```gnuplot
surface(x,y) = 0.1*cos(x)
```

or:

```gnuplot
surface(x,y) = 0.3*cos(x/2.0)
```

Things to test:

- What happens when `a` is larger?
- What happens when `b` is larger?
- Does the cell look like it is moving over taller hills or wider hills?

---

### Exercise 2: Ridges in the `y`-direction

Instead of changing height with `x`, make the height change with `y`.

```math
z = 0.2\cos(y)
```

C++:

```cpp
double surfaceZ(double x, double y) {
    return 0.20 * std::cos(y);
}
```

Gnuplot:

```gnuplot
surface(x,y) = 0.2*cos(y)
```

This creates ridges running parallel to the `x`-axis.

---

### Exercise 3: Diagonal ridges

Make the ridges run diagonally across the domain.

```math
z = 0.2\cos(x+y)
```

C++:

```cpp
double surfaceZ(double x, double y) {
    return 0.20 * std::cos(x + y);
}
```

Gnuplot:

```gnuplot
surface(x,y) = 0.2*cos(x + y)
```

You can also try:

```math
z = 0.2\cos(x-y).
```

C++:

```cpp
double surfaceZ(double x, double y) {
    return 0.20 * std::cos(x - y);
}
```

Gnuplot:

```gnuplot
surface(x,y) = 0.2*cos(x - y)
```

---

### Exercise 4: Eggcrate surface

This surface bends in both the `x` and `y` directions.

```math
z = 0.2\cos(x)\cos(y)
```

C++:

```cpp
double surfaceZ(double x, double y) {
    return 0.20 * std::cos(x) * std::cos(y);
}
```

Gnuplot:

```gnuplot
surface(x,y) = 0.2*cos(x)*cos(y)
```

This creates alternating hills and valleys like an egg carton.

---

### Exercise 5: Gaussian hill and Gaussian valley

A Gaussian hill gives one smooth bump in the center.

```math
z = 0.3e^{-(x^2+y^2)}
```

C++:

```cpp
double surfaceZ(double x, double y) {
    return 0.30 * std::exp(-(x*x + y*y));
}
```

Gnuplot:

```gnuplot
surface(x,y) = 0.3*exp(-(x**2 + y**2))
```

A Gaussian valley gives one smooth pit in the center.

```math
z = -0.3e^{-(x^2+y^2)}
```

C++:

```cpp
double surfaceZ(double x, double y) {
    return -0.30 * std::exp(-(x*x + y*y));
}
```

Gnuplot:

```gnuplot
surface(x,y) = -0.3*exp(-(x**2 + y**2))
```

You may need to adjust:

```gnuplot
set zrange [-0.4:0.4]
```

depending on the amplitude you choose.

---

### Exercise 6: Parabolic bowl and parabolic hill

A parabolic bowl curves upward away from the origin.

```math
z = 0.03(x^2+y^2)
```

C++:

```cpp
double surfaceZ(double x, double y) {
    return 0.03 * (x*x + y*y);
}
```

Gnuplot:

```gnuplot
surface(x,y) = 0.03*(x**2 + y**2)
```

For this one, update the gnuplot z-range:

```gnuplot
set zrange [0:1.0]
```

A parabolic hill curves downward away from the origin.

```math
z = -0.03(x^2+y^2)
```

C++:

```cpp
double surfaceZ(double x, double y) {
    return -0.03 * (x*x + y*y);
}
```

Gnuplot:

```gnuplot
surface(x,y) = -0.03*(x**2 + y**2)
```

For this one, try:

```gnuplot
set zrange [-1.0:0]
```

---

### Exercise 7: Saddle surface

A saddle curves upward in one direction and downward in the other.

```math
z = 0.03(x^2-y^2)
```

C++:

```cpp
double surfaceZ(double x, double y) {
    return 0.03 * (x*x - y*y);
}
```

Gnuplot:

```gnuplot
surface(x,y) = 0.03*(x**2 - y**2)
```

Try:

```gnuplot
set zrange [-0.6:0.6]
```

This is a useful surface because it has both hill-like and valley-like behavior.

---

### Exercise 8: Rippled surface

This creates circular waves centered at the origin.

```math
z = 0.2\cos\left(\sqrt{x^2+y^2}\right)
```

C++:

```cpp
double surfaceZ(double x, double y) {
    return 0.20 * std::cos(std::sqrt(x*x + y*y));
}
```

Gnuplot:

```gnuplot
surface(x,y) = 0.2*cos(sqrt(x**2 + y**2))
```

A more textured version is a decaying ripple:

```math
z = 0.25e^{-(x^2+y^2)/4}\cos\left(3\sqrt{x^2+y^2}\right)
```

C++:

```cpp
double surfaceZ(double x, double y) {
    double r = std::sqrt(x*x + y*y);
    return 0.25 * std::exp(-(r*r)/4.0) * std::cos(3.0*r);
}
```

Gnuplot:

```gnuplot
surface(x,y) = 0.25*exp(-(x**2 + y**2)/4.0)*cos(3.0*sqrt(x**2 + y**2))
```

---

## Final Reminder

Whenever you change the surface in the C++ code, also change the surface in the gnuplot script.

For example, if C++ uses:

```cpp
return 0.20 * std::cos(x) * std::cos(y);
```

then gnuplot should use:

```gnuplot
surface(x,y) = 0.2*cos(x)*cos(y)
```

The simulation and plot need to agree.
