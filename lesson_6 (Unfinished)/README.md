# Lesson 6: Adding Substrate Adhesion

## Review

So far, we have built up our cell model step by step.

First, we created the **cell membrane** by placing membrane nodes around a circle. Then we added **interior nodes** randomly distributed inside the cell. After that, we added the **Euler method** so that we could update the position of each node over time and create motion.

However, when every node moves independently, the cell falls apart. To fix this, we added **intercellular forces** that keep the cell together. These forces include attraction, repulsion, spring forces, and bending forces. Together, these forces help the membrane keep its shape and help the interior nodes stay connected to the rest of the cell.

We also added a **substrate**, which represents the surface the cell moves on.

This week, we will add **adhesion mechanics** so that the cell can interact with the substrate.

The main idea is:

> The cell does not just float freely. It grabs onto the substrate, pulls on those attachment points, and then updates its position according to the forces acting on it.

---

## Substrate Adhesion

Most of the new code has already been written for you and is included in the `forces.h` header file. Here is the main idea behind the adhesion mechanics.

Each membrane node is allowed to form a certain number of possible adhesion sites. These adhesion sites represent places where the cell attaches to the substrate.

For each membrane node, the code does three things:

1. **Check whether existing adhesion sites detach.**
2. **Create new adhesion sites if there is room.**
3. **Add adhesion forces from attached sites.**

### 1. Adhesion Sites Can Detach

If a membrane node is already attached to the substrate, the code checks whether the adhesion site should detach.

The farther the membrane node moves from its binding site, the more likely it is to detach. This is controlled by the parameter

```cpp
adh.charUnbindDist
```

A small value means adhesions detach more easily. A larger value means adhesions last longer.

The probability of detachment is based on the distance between the membrane node and the adhesion site.

---

### 2. New Adhesion Sites Can Form

If a membrane node has open adhesion slots, it may create a new adhesion site.

The probability of forming a new site is controlled by

```cpp
adh.siteBindThreshold
```

For example,

```cpp
adh.siteBindThreshold = 0.1;
```

means that each available adhesion attempt has a 10% chance of forming a new adhesion site.

For now, we will hard code the direction in which new adhesion sites form. Later, this can be changed to depend on curvature, randomness, or other biological rules.

The adhesion length is randomly chosen using an exponential distribution. This means most adhesion sites form close to the node, but occasionally a longer adhesion site can form.

This parameter controls the typical adhesion length:

```cpp
adh.lambda
```

Larger values of `lambda` create shorter average adhesion lengths.

---

### 3. Adhesion Sites Add Forces

Once an adhesion site is attached, it acts like a spring pulling the membrane node toward the binding location.

The force has the form

```cpp
F = kadh * (bindPosition - nodePosition)
```

In component form, this becomes

```cpp
Fx = kadh * (bindX - membraneX)
Fy = kadh * (bindY - membraneY)
Fz = kadh * (bindZ - membraneZ)
```

The parameter

```cpp
adh.kadh
```

controls the strength of the adhesion force.

A larger value makes adhesions pull harder. A smaller value makes adhesions weaker.

---

## Main Code Changes

To use these new forces, make the following changes to the main code.

---

## 1. Hard Code the Adhesion Direction

For now, we will hard code the direction of new adhesion sites.

Add this near your constant variable section:

```cpp
const double direction = PI;
```

This means new adhesion sites will tend to form in the direction

```cpp
theta = PI
```

which points in the negative \(x\)-direction.

Later, you can try changing this value.

For example:

```cpp
const double direction = 0.0;        // positive x-direction
const double direction = PI;         // negative x-direction
const double direction = PI / 2.0;   // positive y-direction
const double direction = -PI / 2.0;  // negative y-direction
```

Only keep one of these lines in your actual code.

---

## 2. Add Bind Position Vectors

Place these near your other position vectors:

```cpp
// Bind Position Vectors
std::vector<double> bindX, bindY, bindZ;
std::vector<bool> isAttached;
```

These vectors store the adhesion information.

The vectors

```cpp
bindX, bindY, bindZ
```

store the position of each adhesion binding site.

This is different from the membrane node positions. The membrane node moves, but the binding site represents the location where the cell attached to the substrate.

So if a membrane node attaches to the substrate, we store the binding location in these vectors.

The vector

```cpp
std::vector<bool> isAttached;
```

stores whether each possible adhesion site is currently active.

A `bool` is a variable that can only be one of two values:

```cpp
true
false
```

In this case:

```cpp
isAttached[siteIndex] = true;
```

means that adhesion site is currently attached to the substrate.

Meanwhile,

```cpp
isAttached[siteIndex] = false;
```

means that adhesion site is currently empty or detached.

Each membrane node can have several possible adhesion sites. For example, if

```cpp
adh.maxSubSitePerNode = 5;
```

then each membrane node has 5 possible adhesion slots.

So if we have `N_membrane` membrane nodes, then the total number of possible adhesion sites is

```cpp
N_membrane * adh.maxSubSitePerNode
```

The indexing is handled inside the header file, but conceptually the adhesion sites are organized like this:

```cpp
siteIndex = nodeIndex * maxSubSitePerNode + localSiteIndex;
```

For example, if each node has 5 possible adhesion sites, then membrane node 0 uses sites 0 through 4, membrane node 1 uses sites 5 through 9, membrane node 2 uses sites 10 through 14, and so on.

---

## 3. Initialize the Adhesion Site Parameters

After we initialize the membrane and interior nodes, we need to initialize the adhesion bind vectors.

Add this before the main time loop:

```cpp
// Initialize bind site parameters
AdhesionParams adh;

adh.maxSubSitePerNode = 5;
adh.siteBindThreshold = 0.1;
adh.charUnbindDist = 0.2;
adh.lambda = 10.0;
adh.kadh = 1.0;
adh.substrateZ = 0.0;

initializeAdhesionSites(
    N_membrane,
    bindX,
    bindY,
    bindZ,
    isAttached,
    adh
);
```

Here is what each parameter means.

```cpp
adh.maxSubSitePerNode = 5;
```

Each membrane node can have up to 5 adhesion sites.

```cpp
adh.siteBindThreshold = 0.1;
```

This controls the probability that a new adhesion site forms.

```cpp
adh.charUnbindDist = 0.2;
```

This controls how easily adhesion sites detach. If the node moves far away from the binding site, detachment becomes more likely.

```cpp
adh.lambda = 10.0;
```

This controls the length scale for newly created adhesion sites. Larger values usually create shorter adhesion lengths.

```cpp
adh.kadh = 1.0;
```

This controls the strength of the adhesion force.

```cpp
adh.substrateZ = 0.0;
```

This places the substrate at \(z = 0\).

The function

```cpp
initializeAdhesionSites(...)
```

creates enough storage space for all possible adhesion sites and initially sets them to detached.

---

## 4. Add Adhesion Forces Inside the Time Step Loop

Inside the time step loop, but before the Euler update, add the new adhesion force function.

```cpp
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
    direction,
    gen;
    surfaceZ
);
```

This function updates the adhesion sites and adds the adhesion forces to the membrane force vectors.

The function uses:

```cpp
membraneX, membraneY, membraneZ
```

to know where the membrane nodes currently are.

It uses:

```cpp
bindX, bindY, bindZ
```

to store where the substrate adhesion sites are located.

It uses:

```cpp
isAttached
```

to check whether each adhesion site is currently attached.

It adds forces into:

```cpp
membraneFx, membraneFy, membraneFz
```

These force vectors are then used in the Euler update.

---

## 5. Euler Update

After all forces have been added, the Euler update moves the membrane nodes.

The idea is still the same:

```cpp
new position = old position + dt * force
```

For example:

```cpp
membraneX[i] += dt * membraneFx[i];
membraneY[i] += dt * membraneFy[i];
membraneZ[i] += dt * membraneFz[i];
```

The difference is that now the total force includes the adhesion force from the substrate.

So the cell is no longer only being held together by internal forces. It is also interacting with the surface underneath it.

---

## Full Code
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

    const double direction = PI;

    // Euler motion parameters
    const int N_steps = 200000;
    const int output_every = 2000;

    // dt controls how strongly forces move the nodes each time step.
    const double dt = 0.0001;

    // Random motion parameters
    // Each node moves in a random direction with random length between 0 and 0.001.
    const double random_step_max = 0.0001;

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

    adh.maxSubSitePerNode = 5;
    adh.siteBindThreshold = 0.5;
    adh.charUnbindDist = 0.2;
    adh.lambda = 10.0;
    adh.kadh = 1.0;
    adh.substrateZ = 0.0;

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
            direction,
            gen,
            surfaceZ
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

## Things to Try

Once your code runs, try changing the adhesion parameters.

### 1. Change the Direction

Try:

```cpp
const double direction = 0.0;
```

or

```cpp
const double direction = PI / 2.0;
```

What direction does the cell move?

---

### 2. Change the Number of Adhesion Sites

Try:

```cpp
adh.maxSubSitePerNode = 1;
```

then try:

```cpp
adh.maxSubSitePerNode = 10;
```

What changes when each membrane node can have more adhesion sites?

---

### 3. Change the Binding Probability

Try:

```cpp
adh.siteBindThreshold = 0.01;
```

then try:

```cpp
adh.siteBindThreshold = 0.5;
```

Does the cell attach more strongly to the substrate?

---

### 4. Change the Adhesion Strength

Try:

```cpp
adh.kadh = 0.1;
```

then try:

```cpp
adh.kadh = 5.0;
```

What happens when adhesion forces are weak? What happens when they are strong?

---

### 5. Change the Unbinding Distance

Try:

```cpp
adh.charUnbindDist = 0.05;
```

then try:

```cpp
adh.charUnbindDist = 1.0;
```

Do the adhesions detach quickly, or do they stay attached for a long time?

---

## Summary

In this lesson, we added a new biological mechanism to the model: **substrate adhesion**.

The cell can now attach to the substrate, pull on those attachment points, detach from old sites, and form new sites.

This is important because real cells do not move by simply sliding around freely. They move by creating attachments, generating forces, and releasing attachments over time.

This adhesion mechanism gives us a foundation for more interesting movement rules, including movement on curved surfaces and eventually curvotaxis.
