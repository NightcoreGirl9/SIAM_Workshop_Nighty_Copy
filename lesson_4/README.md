# Week 4: Adding Forces to the Cell

Welcome to Week 4.

Last week we added pseudo-movement to the cell. The membrane nodes and interior nodes were stored in vectors, and the program used an Euler-style update to move the nodes through time.

This week we make the movement more realistic. Instead of letting every node move freely with only random motion, we add **node-to-node forces**. These forces help the cell behave like a connected object instead of a cloud of unrelated points.

The main idea is:

```text
Random motion alone makes the cell fall apart.
Forces help hold the cell together.
```

In this lesson, the cell still does not represent a fully realistic biological cell. However, this is an important step toward a force-based cell motility model.

---

# What this lesson adds

Compared to `Psuedo_movement.cpp`, the new file `Psuedo_movement_with_forces.cpp` adds several major changes:

- includes a custom header file called `forces.h`,
- creates force vectors for each node,
- computes membrane forces and interior forces at each time step,
- keeps the random movement from the previous version,
- updates each node using both forces and randomness,
- uses a helper function called `writeFrame()` to avoid repeated output code.

The final program still writes two files:

- `membrane.dat`
- `interior.dat`

These files are still used by GNUplot to animate the cell.

---

# Why random motion is not enough

In the previous version, each node moved in a random direction with a random step length. Conceptually, each node did something like this:

```text
x_new = x_old + random_dx
y_new = y_old + random_dy
z_new = z_old + random_dz
```

This creates motion, but it has a problem: every node is acting independently.

That means the membrane nodes do not know they are supposed to form a connected membrane. The interior nodes do not know they are supposed to stay inside or near the cell. After enough time steps, the cell can spread out and fall apart.

To fix this, we add forces.

The new update has the form:

```text
new position = old position + force-based movement + random movement
```

For example:

```text
x_new = x_old + dt * Fx + random_dx
y_new = y_old + dt * Fy + random_dy
z_new = z_old + dt * Fz + random_dz
```

Here:

- `Fx`, `Fy`, and `Fz` are the total force components on the node,
- `dt` controls how strongly the force affects the node position,
- `random_dx`, `random_dy`, and `random_dz` keep the pseudo-random movement.

---

# Including the force header file

The new code uses a custom header file:

```cpp
#include "forces.h"
```

This file contains the force functions. Keeping the force formulas in a separate file helps keep the main simulation code readable.

The main file focuses on the simulation flow:

```text
initialize nodes
create force vectors
for each time step:
    reset forces
    compute forces
    update node positions
    write output frames
```

The header file focuses on the force calculations.

Make sure `forces.h` is in the same folder as your `.cpp` file.

Your folder should look something like this:

```text
lesson_4/
├── forces.h
├── plot_animation.gnu
├── Psuedo_movement_random.cpp
└── Psuedo_movement_with_forces.cpp
```

---

# Basic force idea

A force tells a node which direction to move and how strongly to move in that direction.

In 3D, a force has three components:

```text
F = (Fx, Fy, Fz)
```

For our code, each node needs its own force in the `x`, `y`, and `z` directions.

So for membrane nodes, we create:

```cpp
std::vector<double> membraneFx(N_membrane, 0.0);
std::vector<double> membraneFy(N_membrane, 0.0);
std::vector<double> membraneFz(N_membrane, 0.0);
```

For interior nodes, we create:

```cpp
std::vector<double> interiorFx(N_interior, 0.0);
std::vector<double> interiorFy(N_interior, 0.0);
std::vector<double> interiorFz(N_interior, 0.0);
```

A membrane node with index `i` has position:

```cpp
(membraneX[i], membraneY[i], membraneZ[i])
```

and force:

```cpp
(membraneFx[i], membraneFy[i], membraneFz[i])
```

The same idea applies to interior nodes.

---

# Resetting the force vectors

At the start of every time step, all forces must be reset to zero:

```cpp
resetForces(membraneFx, membraneFy, membraneFz);
resetForces(interiorFx, interiorFy, interiorFz);
```

This is important.

The force vectors store the force for the **current** time step only. If we forget to reset them, forces from old time steps would keep accumulating forever, which would make the motion incorrect.

The pattern is:

```text
start time step
set all forces to zero
add spring forces
add bending forces
add Morse forces
use total forces to update positions
end time step
```

---

# The forces used in this lesson

This lesson uses three main types of forces:

1. spring forces,
2. Morse forces,
3. bending forces.

Each force has a different purpose.

---

# Spring forces

The spring forces connect neighboring membrane nodes.

The membrane is a loop of points, so each membrane node has a left neighbor and a right neighbor:

```text
left neighbor -- current node -- right neighbor
```

The spring force helps neighboring membrane nodes stay close to a preferred distance.

A simple spring force has the form:

```text
F_spring = k * (distance - restLength)
```

where:

- `k` is the spring stiffness,
- `distance` is the current distance between two neighboring nodes,
- `restLength` is the preferred distance between neighboring nodes.

If the distance is larger than the rest length, the spring pulls the nodes together.

If the distance is smaller than the rest length, the spring pushes the nodes apart.

In this project, the spring force helps keep the membrane connected.

The function used directly in the code is:

```cpp
addMembraneSpringForces(
    membraneX, membraneY, membraneZ,
    membraneFx, membraneFy, membraneFz
);
```

This adds spring force contributions to the membrane force vectors.

---

# Morse forces

The Morse force is used for node-to-node interactions. It can represent attraction and repulsion depending on the distance between nodes and the parameter values.

The force magnitude has the form:

```text
F_Morse(r) = (U0 / k1) * exp(-r / k1) - (V0 / k2) * exp(-r / k2)
```

where:

- `r` is the distance between two nodes,
- `U0` controls the attraction force,
- `V0` controls replusion force,
- `k1` and `k2` control how quickly the interaction decays with distance.

The exponential terms make the force strongest when nodes are close and weaker when nodes are farther apart.

The code also uses a cutoff distance. If two nodes are farther apart than the cutoff, the Morse force is not computed.

This keeps the simulation more efficient and prevents faraway nodes from interacting unnecessarily.

The general idea is:

```text
if distance < cutoff:
    compute Morse force
else:
    force = 0
```

---

# Membrane-membrane Morse force

The membrane-membrane Morse force allows membrane nodes to interact with other membrane nodes.

The function used in the code is:

```cpp
addMembraneMembraneMorseForces(
    membraneX, membraneY, membraneZ,
    membraneFx, membraneFy, membraneFz
);
```

This adds Morse force contributions to the membrane nodes.

In your current parameters, the membrane-membrane Morse constants are set to zero, so this term is effectively turned off for now. The function is still included because it gives us a place to turn on this force later.

---

# Membrane-interior Morse forces

The membrane and interior nodes also interact.

There are two directions to think about:

```text
interior nodes act on membrane nodes
membrane nodes act on interior nodes
```

The first function adds the effect of interior nodes on the membrane:

```cpp
addMembraneInteriorMorseForces(
    membraneX, membraneY, membraneZ,
    interiorX, interiorY, interiorZ,
    membraneFx, membraneFy, membraneFz
);
```

This updates:

```cpp
membraneFx
membraneFy
membraneFz
```

The second function adds the effect of membrane nodes on the interior nodes:

```cpp
addMembraneInteriorMorseForcesOnInterior(
    membraneX, membraneY, membraneZ,
    interiorX, interiorY, interiorZ,
    interiorFx, interiorFy, interiorFz
);
```

This updates:

```cpp
interiorFx
interiorFy
interiorFz
```

Together, these forces help couple the membrane and interior nodes.

---

# Interior-interior Morse force

Interior nodes also interact with other interior nodes.

The function used in the code is:

```cpp
addInteriorInteriorMorseForces(
    interiorX, interiorY, interiorZ,
    interiorFx, interiorFy, interiorFz
);
```

This adds forces between interior nodes.

This force helps the interior nodes behave more like a connected material instead of completely independent random particles.

---

# Bending force

The bending force is a structural force on the membrane.

The spring force keeps neighboring membrane nodes close to the correct distance. The bending force helps keep the membrane smooth.

Without bending, the membrane may become jagged or form sharp corners.

The bending force compares local angles in the membrane to a preferred angle. Conceptually:

```text
F_bending depends on (current angle - preferred angle)
```

A simplified way to think about it is:

```text
if the membrane bends too sharply:
    bending force pushes it toward a smoother shape
```

The function used in the code is:

```cpp
addMembraneBendingForces(
    membraneX, membraneY, membraneZ,
    membraneFx, membraneFy, membraneFz
);
```

This updates the membrane force vectors.

Bending is important because the membrane is not just a collection of separate springs. It should also resist sharp local changes in shape.

---

# The helper `writeFrame` function

The previous code wrote the membrane and interior data directly inside the time-stepping loop.

That works, but it creates repeated code. We need to write a frame at the beginning and then again every `output_every` time steps.

Instead of copying and pasting the same output code multiple times, the new code defines a helper function inside `main()`:

```cpp
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
```

This is a C++ lambda function.

For our purposes, you can think of it as a small local function called `writeFrame`.

The line:

```cpp
auto writeFrame = [&]() {
```

means:

- `auto` lets C++ figure out the function type,
- `[&]` lets the function use nearby variables like `membraneFile`, `membraneX`, and `interiorX`,
- `()` means the function takes no input arguments,
- `{ ... }` contains the code that runs when we call it.

After defining it, we can write one frame using:

```cpp
writeFrame();
```

The code writes the initial cell position before the time loop:

```cpp
writeFrame();
```

Then it writes later frames inside the time loop:

```cpp
if (step % output_every == 0) {
    writeFrame();
}
```

This keeps the program cleaner.

---

# Warning: Large simulations can take time

This version of the code is more expensive to run than the previous pseudo-movement code. The reason is that the nodes are no longer just moving randomly. At every time step, the code now computes several **node-to-node interactions**.

For this lesson, we have:

```cpp
const int N_membrane = 120;
const int N_interior = 120;
```

So the cell has:

```text
120 membrane nodes + 120 interior nodes = 240 total nodes
```

Some force calculations only use nearby membrane nodes, such as the spring and bending forces. These are relatively cheap. However, the Morse forces compare many nodes against many other nodes.

For example:

```text
membrane-membrane interactions: 120 × 120 = 14,400 checks per time step
membrane-interior interactions: 120 × 120 = 14,400 checks per time step
interior-interior interactions: 120 × 120 = 14,400 checks per time step
```

So the code performs roughly:

```text
14,400 + 14,400 + 14,400 = 43,200 node-to-node checks per time step
```

That means the cost grows quickly as we increase the number of time steps.

## Example timing

On one test run, the code took approximately:

```text
1,000 steps  → 1.52 seconds
10,000 steps → 15 seconds
```

Using this estimate:

| Number of time steps | Estimated runtime |
|---:|---:|
| 1,000 | about 1.5 seconds |
| 100,000 | about 150 seconds, or 2.5 minutes |
| 1,000,000 | about 1,500 seconds, or 25 minutes |

These times will depend on your computer, compiler, and how many other programs are running.

## Important note about saved frames

The number of frames saved to file is controlled by:

```cpp
const int output_every = 100000;
```

For example, if

```cpp
const int N_steps = 1000000;
const int output_every = 100000;
```

then the code only saves about 10 frames.

This keeps the output files smaller, but it does **not** make the force calculations much cheaper. The forces are still computed at every time step.

## Recommendation

When testing changes, start with a smaller simulation first:

```cpp
const int N_steps = 1000;
const int output_every = 100;
```

or

```cpp
const int N_steps = 10000;
const int output_every = 1000;
```

Once the code is working and the cell behaves reasonably, then try larger values such as:

```cpp
const int N_steps = 100000;
const int output_every = 10000;
```

A good habit is to test small first, check the animation, and only then run a larger simulation.

---

# The new time-stepping loop

The old time-stepping loop was based mostly on direct motion.

The new time-stepping loop has a force-based structure:

```text
for each time step:
    reset all forces to zero
    compute membrane forces
    compute membrane-interior forces
    compute interior forces
    update membrane node positions
    update interior node positions
    write frame if needed
```

In code, the main force-computation part looks like this:

```cpp
resetForces(membraneFx, membraneFy, membraneFz);
resetForces(interiorFx, interiorFy, interiorFz);

addMembraneMembraneMorseForces(...);
addMembraneSpringForces(...);
addMembraneBendingForces(...);
addMembraneInteriorMorseForces(...);
addMembraneInteriorMorseForcesOnInterior(...);
addInteriorInteriorMorseForces(...);
```

The force functions add to the force vectors. After all of them have been called, each node has a total force.

---

# The Euler update with forces and randomness

The old version moved nodes using a fixed change, or in the random version, a random displacement.

The new version uses both force and randomness.

For membrane nodes:

```cpp
membraneX[i] += dt * membraneFx[i] + random_dx;
membraneY[i] += dt * membraneFy[i] + random_dy;
membraneZ[i] += dt * membraneFz[i] + random_dz;
```

For interior nodes:

```cpp
interiorX[i] += dt * interiorFx[i] + random_dx;
interiorY[i] += dt * interiorFy[i] + random_dy;
interiorZ[i] += dt * interiorFz[i] + random_dz;
```

This has the form:

```text
new position = old position + force motion + random motion
```

The force motion is:

```text
dt * force
```

The random motion is generated by choosing:

- a random angle,
- a random length between `0` and `0.01`.

The code uses:

```cpp
double theta_random = angleDist(gen);
double random_length = randomLengthDist(gen);

double random_dx = random_length * std::cos(theta_random);
double random_dy = random_length * std::sin(theta_random);
double random_dz = 0.0;
```

This makes each node take a small random step in the `xy`-plane.

---

# Why `dt` is small

The code uses:

```cpp
const double dt = 0.00001;
```

This value is small because the force constants can be large. For example, the membrane spring stiffness is large enough that using a large `dt` may make the simulation unstable.

If `dt` is too large, the cell may explode, fly apart, or produce strange numerical behavior.

If that happens, try reducing `dt`:

```cpp
const double dt = 0.000001;
```

You can also reduce the random step size:

```cpp
const double random_step_max = 0.001;
```

---

# Important note about stability

This code is now a force-based simulation. Force-based simulations are more sensitive than simple translation.

Small parameter changes can have large effects.

If the animation looks unstable, try changing only one thing at a time:

- decrease `dt`,
- decrease `random_step_max`,
- decrease the number of steps,
- increase `output_every` while testing,
- temporarily turn off one force at a time by commenting out one force function call.

For example, to test only membrane spring forces, you could temporarily comment out:

```cpp
// addMembraneBendingForces(...);
// addMembraneInteriorMorseForces(...);
// addInteriorInteriorMorseForces(...);
```

Then add the forces back one at a time.

---

# Full code

```cpp
```cpp
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <vector>

#include "forces.h"

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
    const double random_step_max = 0.001;

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
            double theta_random = angleDist(gen);
            double random_length = randomLengthDist(gen);

            double random_dx = random_length * std::cos(theta_random);
            double random_dy = random_length * std::sin(theta_random);
            double random_dz = 0.0;

            membraneX[i] += dt * membraneFx[i] + random_dx;
            membraneY[i] += dt * membraneFy[i] + random_dy;
            membraneZ[i] += dt * membraneFz[i] + random_dz;
        }

        // Euler update for interior nodes
        for (int i = 0; i < N_interior; i++) {
            double theta_random = angleDist(gen);
            double random_length = randomLengthDist(gen);

            double random_dx = random_length * std::cos(theta_random);
            double random_dy = random_length * std::sin(theta_random);
            double random_dz = 0.0;

            interiorX[i] += dt * interiorFx[i] + random_dx;
            interiorY[i] += dt * interiorFy[i] + random_dy;
            interiorZ[i] += dt * interiorFz[i] + random_dz;
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

# Visualizing with GNUplot
The steps are the same as in previous week. This week however, the `plot_animation.gnu` is in 2D to make viewing the simulation easier. It will also save the file as a `gif` which you can view in VScode. It should work for Windows, Mac and Linux.
To view the simulation, open `GNUplot` in the folder for `lesson_4` and run
```gnuplot
load 'plot_animation.gnu'
```
Then click on `psuedo_movement.gif` to view the animation.
---

# Try-It-Yourself Ideas

Once the code runs, try changing one feature at a time.

## Change the random step size

Try:

```cpp
const double random_step_max = 0.0005;
```

or:

```cpp
const double random_step_max = 0.00001;
```

Question: does the cell stay together better when the random motion is smaller?

---

## Change the time step

Try:

```cpp
const double dt = 0.000001;
```

Question: does the motion become more stable?

---

## Turn off random motion

Try:

```cpp
const double random_step_max = 0.0;
```

Question: what do the forces do by themselves?

---

## Turn off one force at a time

Comment out one force function call and rerun the code.

For example:

```cpp
// addMembraneBendingForces(...);
```

Question: what happens to the membrane when the bending force is removed?

---

## Change the number of nodes

Try changing:

```cpp
const int N_membrane = 120;
const int N_interior = 120;
```

Question: does the simulation look smoother with more nodes? Does it run slower?

---

