# Week 3: Adding Pseudo-Movement

Welcome to Week 3.

Last week we initialized the membrane and interior nodes of a simple cell. This week we add **pseudo-movement**: the whole cell will translate in the positive `x` direction over time. The point of this lesson is not yet to build a biologically realistic motility model. Instead, the goal is to introduce three core simulation ideas:

1. **Store node locations in memory** using vectors.
2. **Update those locations repeatedly** with a simple Euler-style step.
3. **Write snapshots to file** so the motion can be visualized in GNUplot.

In this version, the cell shape stays the same and only its position changes. That makes it a good first example of a time-stepping simulation.

---

# What this lesson adds

Compared to the previous week, the code now does four new things:

- adds the C++ `vector` library,
- stores membrane and interior node coordinates in six vectors,
- updates all node positions over many time steps,
- saves the data in a format GNUplot can animate.

The final program writes two files:

- `membrane.dat`
- `interior.dat`

Each saved frame is separated by **two blank lines**, so GNUplot can read each frame as a separate data block.

---

# Add the vector library

In the earlier version of the program, each node was created and written directly to file. That is enough for a single static plot, but not enough for motion. To move the nodes later, we first need to store their coordinates in memory.

Add this header near the top of the file:

```cpp
#include <vector>
```

---

# What is a vector?

You can think of a vector as a dynamic list. It stores many values of the same type and can grow as needed.

For this lesson, we store the `x`, `y`, and `z` coordinates separately for:

- membrane nodes,
- interior nodes.

That means we need six vectors total:

```cpp
std::vector<double> membraneX, membraneY, membraneZ;
std::vector<double> interiorX, interiorY, interiorZ;
```

A membrane node with index `i` has location

```cpp
(membraneX[i], membraneY[i], membraneZ[i])
```

and an interior node with index `i` has location

```cpp
(interiorX[i], interiorY[i], interiorZ[i])
```

Remember that C++ indexing starts at `0`.

For example, the 21st membrane node is stored at index `20`:

```cpp
(membraneX[20], membraneY[20], membraneZ[20])
```

To add a value to the end of a vector, use `push_back`:

```cpp
double x = 1.0;
membraneX.push_back(x);
```

---

# Storing the initialized nodes

When building the membrane and interior nodes, do **not** write them directly to file inside the initialization loops. Instead, store them in the vectors.

For the membrane nodes, use:

```cpp
membraneX.push_back(x);
membraneY.push_back(y);
membraneZ.push_back(z);
```

For the interior nodes, use:

```cpp
interiorX.push_back(x);
interiorY.push_back(y);
interiorZ.push_back(z);
```

This means the program first builds the full cell in memory. After that, the time-stepping loop can update the coordinates.

---

# Euler-style motion update

Now we add motion.

The update rule is based on a simple Euler step. The idea is:

```text
new value = old value + change
```

For this lesson, we will only move in the `x`-direction:

- move only in the `x` direction,
- keep `y` fixed,
- keep `z` fixed.

So each update looks like

```text
x_new = x_old + dx
```

with `dy = 0` and `dz = 0`.

This gives a rigid translation of the whole cell rather than a deformation.

Use these parameters:

```cpp
const int N_steps = 1000;
const int output_every = 100;
const double dx = 0.001;
const double dy = 0;
const double dz = 0;
```

## What these values mean

- `N_steps = 1000` means the simulation runs for 1,000 update steps.
- `output_every = 100` means the code saves one frame every 100 steps instead of every single step.
- `dx = 0.001` means every node moves by `0.001` in the positive `x` direction each time step.
- `dy = 0` and `dz = 0` mean there is no motion in the `y` or `z` directions.

So the total displacement in `x` after all 1,000 steps is

```text
1000 × 0.001 = 1
```

The cell starts centered near the origin and ends about 1 units to the right.

---

# Why we only save every 100 steps

If we wrote to file at every time step, the output files would be much larger and the animation would contain many nearly identical frames.

Instead, we use

```cpp
if (step % output_every == 0)
```

This means:

- save when `step = 0`,
- save when `step = 100`,
- save when `step = 200`,
- and so on.

Because `1000 / 100 = 100`, the simulation saves about 10 frames.

---

# Important GNUplot formatting detail

To make GNUplot read each saved frame as a separate `index`, the data blocks must be separated by **two blank lines**.

That is why the code uses:

```cpp
membraneFile << "\n\n";
interiorFile << "\n\n";
```

instead of a single blank line.

---

# Final code: 

```cpp
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <vector>

int main() {
    const int N_membrane = 120;
    const int N_interior = 120;
    const double PI = 3.141592653589793;

    const double R_membrane = 1.0;
    const double R_interior_max = 0.95;

    // Euler motion parameters
    const int N_steps = 1000;
    const int output_every = 100;
    const double dx = .001;
    const double dy = 0;
    const double dz = 0;

    std::vector<double> membraneX, membraneY, membraneZ;     // Vectors for membrane node positions
    std::vector<double> interiorX, interiorY, interiorZ;     // Vectors for interior node positions

    std::ofstream membraneFile("membrane.dat");
    std::ofstream interiorFile("interior.dat");

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
    membraneFile << "\n\n";

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
    interiorFile << "\n\n";

    // Time stepping loop
    for (int step = 0; step < N_steps; step++) {
        for (int i = 0; i < N_membrane; i++) {
            membraneX[i] += dx;
            membraneY[i] += dy;
            membraneZ[i] += dz;
        }

        for (int i = 0; i < N_interior; i++) {
            interiorX[i] += dx;
            interiorY[i] += dy;
            interiorZ[i] += dz;
        }

        if (step % output_every == 0) {
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
        }
    }
    membraneFile.close();
    interiorFile.close();

    std::cout << "Created membrane.dat and interior.dat\n";

    return 0;
}
```

---

# Visualizing the motion with GNUplot

This lesson also includes a GNUplot script called `plot_animation.gnu`.

Instead of typing the plotting commands by hand, open GNUplot and run:

```gnuplot
load 'plot_animation.gnu'
```

This script reads the saved blocks from `membrane.dat` and `interior.dat` and plays them as a 3D animation.

##GNUplot script

```gnuplot
set terminal wxt size 800,800
set title "Pseudo Cell Movement in 3D"
set xlabel "x"
set ylabel "y"
set zlabel "z"

set xrange [-4:4]
set yrange [-4:4]
set zrange [-1:1]

set xyplane at 0
set view 60,30
set view equal xyz

set size 1,1
set lmargin 0
set rmargin 0
set bmargin 0
set tmargin 0

set grid
set border 4095
unset key

Nsteps = 10

do for [i=0:Nsteps-1] {
    splot \
        'membrane.dat' index i using 1:2:3 with points pt 7 ps 1.0 lc rgb "black", \
        'interior.dat' index i using 1:2:3 with points pt 7 ps 0.6 lc rgb "red"

    pause .5
}
```

---

# Notes on the GNUplot script

A few lines are especially important:

- `set view equal xy` keeps the `x`, `y`, scales proportional.
- `set xyplane at 0` places the base plane at `z = 0`.
- the `index i` syntax tells GNUplot to read frame `i` from the `.dat` files.
- `pause .5` displays each frame for half a second.

## One practical note

For now, you need to manually match the number of frames from your simulation with `GNUplot`

```gnuplot
Nsteps = 10
```

If you want to show more motion, you will need to make changes to the code and the script..

---
## Try-It-Yourself Ideas

Now that the code is working, try changing some of the motion parameters and see what happens. This is a good way to build intuition for how the simulation loop and the plotting script work together.

### Suggestions to test
- **Change the total number of time steps**
  - Try making `N_steps` smaller, such as `500`, or larger, such as `20000`.
  - How does this change the total distance traveled by the cell?

- **Change the x-direction motion**
  - Try
    ```cpp
    const double dx = -0.001;
    ```
  - This should make the cell move in the negative `x`-direction instead of the positive `x`-direction.

- **Add motion in the y-direction**
  - Try
    ```cpp
    const double dy = 0.001;
    ```
  - Now the cell should move diagonally in the `xy`-plane.

- **Add motion in the z-direction**
  - Try
    ```cpp
    const double dz = 0.001;
    ```
  - This will lift the cell out of the `xy`-plane and make the motion fully 3D.

- **Change how often data is saved**
  - Try changing
    ```cpp
    const int output_every = 100;
    ```
    to something like `50` or `200`.
  - Smaller values save more frames. Larger values save fewer frames.

- **Change the number of membrane or interior nodes**
  - Try changing `N_membrane` or `N_interior`.
  - How does this affect the appearance of the cell?

### Important note
If you change the simulation parameters, you may also need to edit `plot_animation.gnu`.

For example, depending on the changes you make, you may need to adjust:
- the number of plotted frames,
- the `x`, `y`, or `z` axis ranges,
- the viewing window size,
- the pause time between frames.

In particular:
- if the cell moves farther, you may need a larger `xrange`, `yrange`, or `zrange`,
- if you save more or fewer frames, you may need to update the frame count in the script,
- if you add motion in the `z` direction, you may want to expand the `zrange`.

A good habit is to make one small change at a time, rerun the code, and then update the GNUplot script only as needed.
---

# Summary

At this stage, the program does the following:

1. builds a circular membrane,
2. fills the interior with random nodes,
3. stores all node coordinates in vectors,
4. updates those coordinates through time,
5. writes snapshots to `.dat` files,
6. animates the saved frames in GNUplot.

This is the first complete version of a simple time-stepping simulation pipeline.

The motion is still very basic, but the structure is now in place for future upgrades such as:

- different membrane and interior velocities,
- random perturbations,
- deformation instead of rigid translation,
- curved substrates,
- simple force laws.
