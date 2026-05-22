# Lesson 7: Curvotaxis on a Curved Surface

## Goal

The goal of this lesson is to modify our cell motility model so that it can reproduce the main qualitative behavior from the curvotaxis paper by Pieuchot et al.: cells can sense curvature and tend to migrate toward valleys on a curved surface.

In earlier lessons, we built a cell from membrane and interior nodes, added force laws to help the cell keep its shape, and placed the cell on a curved substrate. This week, we add curvature-dependent adhesion rules. These rules allow the cell to respond differently depending on whether a membrane node is located near a peak, a valley, or a nearly flat part of the surface.

The code this week is more complex and may take longer to compile and run, so the complete code has already been uploaded for you. Your main task is to understand the logic, run the simulation, and experiment with the tuning parameters.

---
## Biological background

Pieuchot et al. introduced **curvotaxis**, the tendency of cells to respond to cell-scale curvature and preferentially migrate toward **concave valleys** rather than remaining on **convex peaks**.

A key biological observation from the paper is that focal adhesions behave differently depending on local surface curvature:

- **On convex peaks/ridges:** cells form **more focal adhesions**, but those adhesions are **shorter lived** and less stable.
- **In concave valleys:** cells form **fewer focal adhesions**, but those adhesions are **more stable** and longer lived.

In our model, we represent this idea by making three adhesion quantities depend on curvature:

1. the probability that an existing adhesion unbinds,
2. the probability that a new adhesion forms,
3. the length/strength scale of a newly formed adhesion.

The goal is not to exactly reproduce every biological detail, but to encode the main qualitative idea:

> Curvature changes adhesion dynamics, and asymmetric adhesion dynamics can bias cell motion.

---

## Curvature sign convention

Our surface is

$$
z = f(x).
$$

The signed curvature used in the code is

$$
\kappa = \frac{-f''(x)}{\left(1 + (f'(x))^2\right)^{3/2}}.
$$

For the surface

$$
z = a\cos(x/b),
$$

this gives the convention:

$$
\kappa > 0 \quad \text{on convex peaks/ridges},
$$

and

$$
\kappa < 0 \quad \text{in concave valleys}.
$$

So in the code:

```cpp
kappa > 0.0   // convex peak/ridge
kappa < 0.0   // concave valley
```

---

## 1. Curvature-dependent unbinding

The unbinding probability for an existing adhesion is modeled as

$$ 
P_{\text{unbind}}=1 - \exp\left(-\frac{d}{L(\kappa)}\right),
$$S

where:

- $d$ is the distance between the node and the adhesion site,
- $L(\kappa)$ is the curvature-dependent characteristic unbinding distance.

The curvature-modified unbinding distance is

$$
L(\kappa)=L_0 M_{\text{unbind}}(\kappa),
$$

where $L_0$ is the base value `charUnbindDist`.

The multiplier is

$$
M_{\text{unbind}}(\kappa)=\begin{cases} 1 + \beta_{\text{conc}} r, & \kappa < 0, \\
1 - \beta_{\text{conv}} r, & \kappa \ge 0,\end{cases}
$$

with

$$
r = \frac{|\kappa|}{\kappa_{\text{scale}} + |\kappa|}.
$$

### Interpretation

For concave valleys, $\kappa < 0$, so

$$
L(\kappa) = L_0(1 + \beta_{\text{conc}}r).
$$

This makes $L(\kappa)$ larger. Since

$$
P_{\text{unbind}} = 1 - e^{-d/L},
$$

a larger $L$ gives a smaller unbinding probability. Therefore:

> Increasing `betaConc` makes adhesions in valleys more stable.

For convex peaks, $\kappa > 0$, so

$$
L(\kappa) = L_0(1 - \beta_{\text{conv}}r).
$$

This makes $L(\kappa)$ smaller, which increases the unbinding probability. Therefore:

> Increasing `betaConv` makes adhesions on peaks more likely to unbind.

### Parameters

```cpp
adh.charUnbindDist
adh.betaConc
adh.betaConv
```

- Increase `charUnbindDist`: adhesions generally live longer everywhere.
- Decrease `charUnbindDist`: adhesions generally unbind more easily everywhere.
- Increase `betaConc`: valleys become more adhesive/stable.
- Decrease `betaConc`: valleys become less stabilizing.
- Increase `betaConv`: peaks become more slippery because adhesions unbind more easily.
- Decrease `betaConv`: adhesions on peaks become less likely to unbind.

This is the main model component that supports Pieuchot-style migration toward valleys.

---

## 2. Curvature-dependent new binding probability

The probability of forming a new adhesion site is also curvature-dependent.

We define

$$
k_0 = \kappa_{\text{par}}\kappa_{\text{scale}},
$$

and

$$
r = \frac{|\kappa|}{k_0 + |\kappa|}.
$$

The binding probability is

$$
p(\kappa) = \begin{cases} p_0 + (p_{\text{pos}} - p_0)r, & \kappa \ge 0, \\ 
p_0 - (p_0 - p_{\text{neg}})r, & \kappa < 0. \end{cases}
$$

### Interpretation

For positive curvature, $\kappa > 0$, the probability moves from $p_0$ toward $p_{\text{pos}}$. For negative curvature, $\kappa < 0$, the probability moves from $p_0$ toward $p_{\text{neg}}$.

With the common default values

$$
p_{\text{neg}} = 0, \qquad p_0 = 0.5, \qquad p_{\text{pos}} = 1,
$$

we get:

- peaks/ridges have high binding probability,
- flat regions have intermediate binding probability,
- valleys have low binding probability.

This matches the biological idea:

> More adhesions form on peaks, while fewer adhesions form in valleys.

### Parameters

```cpp
adh.pNeg
adh.p0
adh.pPos
adh.kappaPar
```

- Increase `pPos`: more new adhesions form on peaks.
- Decrease `pPos`: weakens the peak-binding preference.
- Increase `pNeg`: more new adhesions form in valleys.
- Decrease `pNeg`: fewer new adhesions form in valleys.
- Increase `p0`: raises the baseline binding probability everywhere.
- Decrease `p0`: lowers the baseline binding probability everywhere.
- Increase `kappaPar`: curvature effects saturate more slowly, so curvature has a weaker effect unless curvature is large.
- Decrease `kappaPar`: curvature effects saturate faster, so curvature has a stronger effect even for moderate curvature.

This mechanism can favor peaks because more adhesions form there. In the biological interpretation, those peak adhesions should be shorter lived, so this mechanism should be balanced against curvature-dependent unbinding.

---

## 3. Curvature-dependent adhesion length / strength scale

When a new adhesion forms, its initial length is sampled from a lognormal distribution. The mean of this distribution depends on curvature.

First define the signed curvature saturation factor

$$
s = \frac{\kappa}{\kappa_{\text{scale}} + |\kappa|}.
$$

Unlike $r$, this quantity keeps the sign of curvature:

$$
s > 0 \quad \text{on peaks},
$$

and

$$
s < 0 \quad \text{in valleys}.
$$

The curvature-dependent mean adhesion length is

$$
m(\kappa)
=
m_0(1 + \gamma_{\text{len}}s),
$$

where $m_0$ is the base value `meanAdhesionLength`.

The actual adhesion length is sampled as

$$
\ell = \exp(\mu + \sigma Z),
$$

where

$$
Z \sim \mathcal{N}(0,1),
$$

and

$$
\mu = \log(m(\kappa)) - \frac{\sigma^2}{2}.
$$

This choice of $\mu$ makes the expected adhesion length equal to $m(\kappa)$.

### Interpretation

If $\gamma_{\text{len}} > 0$, then:

- peaks/ridges have larger mean adhesion length,
- valleys have smaller mean adhesion length.

This can make peak adhesions mechanically stronger or longer at formation, depending on how the length interacts with the force law.

If this effect is too strong, the cell may climb toward peaks instead of stabilizing in valleys.

### Parameters

```cpp
adh.meanAdhesionLength
adh.gammaLen
adh.lognormalSigma
adh.minMeanAdhesionLength
```

- Increase `meanAdhesionLength`: adhesions are generally longer when they form.
- Decrease `meanAdhesionLength`: adhesions are generally shorter when they form.
- Increase `gammaLen`: curvature has a stronger effect on adhesion length.
- Set `gammaLen = 0.0`: turns off curvature-dependent adhesion length.
- Increase `lognormalSigma`: adhesion lengths become more variable.
- Decrease `lognormalSigma`: adhesion lengths become more consistent.
- Increase `minMeanAdhesionLength`: prevents adhesion lengths from becoming too small.

This mechanism should be tuned carefully because it can compete with the valley-stabilizing unbinding mechanism.

---

## Summary of the three mechanisms

| Model component | Equation | Biological meaning | Likely bias |
|---|---|---|---|
| Curvature-dependent unbinding | $P_{\text{unbind}} = 1 - e^{-d/L(\kappa)}$ | Adhesions in valleys are more stable; adhesions on peaks are shorter lived | Valley-favoring |
| Curvature-dependent binding | $p(\kappa)$ moves toward $p_{\text{pos}}$ on peaks and $p_{\text{neg}}$ in valleys | More adhesions form on peaks; fewer form in valleys | Peak-favoring |
| Curvature-dependent adhesion length | $m(\kappa)=m_0(1+\gamma_{\text{len}}s)$ | Curvature changes the size/strength scale of new adhesions | Depends on tuning |

The important modeling point is that these effects compete. The biological observation is not simply “more adhesions means the cell goes there.” Pieuchot et al. found that adhesions can be more numerous on convex regions but more stable in concave regions. In the model, cell behavior depends on the balance between adhesion formation and adhesion lifetime.

---

## Suggested tuning exercises

### Exercise 1: Isolate curvature-dependent unbinding

```cpp
adh.pNeg = 0.5;
adh.p0   = 0.5;
adh.pPos = 0.5;
adh.gammaLen = 0.0;

adh.betaConc = 60.0;
adh.betaConv = 0.02;
```

This turns off curvature-dependent new binding and curvature-dependent adhesion length. Only adhesion lifetime depends on curvature.

Expected behavior: the cell should be more likely to stabilize near valleys.

---

### Exercise 2: Strengthen peak unbinding

```cpp
adh.betaConv = 0.1;
```

Increasing `betaConv` makes adhesions on peaks shorter lived.

Expected behavior: peaks become more slippery, so the cell should be less likely to stay on peaks.

---

### Exercise 3: Strengthen valley stabilization

```cpp
adh.betaConc = 100.0;
```

Increasing `betaConc` makes adhesions in valleys longer lived.

Expected behavior: once the cell reaches a valley, it should remain there more strongly.

---

### Exercise 4: Weaken peak binding

```cpp
adh.pNeg = 0.25;
adh.p0   = 0.5;
adh.pPos = 0.75;
```

This reduces the difference between valley and peak binding probabilities.

Expected behavior: the cell should be less biased toward forming many adhesions on peaks.

---

### Exercise 5: Turn off curvature-dependent adhesion length

```cpp
adh.gammaLen = 0.0;
```

This makes new adhesion lengths independent of curvature.

Expected behavior: the result will depend mainly on the competition between new binding probability and unbinding lifetime.

---

## Files for This Lesson

This lesson uses two main files:

```text
main.cpp        // Main simulation file
forces.h        // Force and adhesion helper functions
```

The main file sets up the cell, surface, time-stepping loop, output files, and simulation parameters.

The header file contains the force functions, adhesion data structures, curvature calculations, and curvature-dependent adhesion logic.

---

# Header File: `forces.h`

The header file contains most of the reusable mechanics of the model. The major parts are described below.

---

## `struct AdhesionParams`

`AdhesionParams` stores the parameters used by the adhesion model.

```cpp
struct AdhesionParams {
    int maxSubSitePerNode = 20;

    double charUnbindDist = 0.1;
    double kadh = 2.0;

    double betaConc = 60.0;
    double betaConv = 0.02;

    double pNeg = 0.0;
    double p0 = 0.5;
    double pPos = 1.0;
    double kappaPar = 0.5;
    double epsK = 1.0e-6;

    double meanAdhesionLength = 0.1;
    double gammaLen = 0.2;
    double lengthLognormalSigma = 0.25;
    double minMeanAdhesionLength = 1.0e-6;
};
```

### Basic adhesion parameters

| Parameter | Meaning |
|---|---|
| `maxSubSitePerNode` | Maximum number of possible adhesion sites per membrane node. |
| `charUnbindDist` | Baseline characteristic distance controlling unbinding probability. |
| `kadh` | Adhesion spring constant. Larger values make adhesion sites pull more strongly. |

### Curvature-dependent unbinding parameters

| Parameter | Meaning |
|---|---|
| `betaConc` | Strength of adhesion stabilization in concave regions. Larger values make valley adhesions longer-lived. |
| `betaConv` | Strength of adhesion destabilization in convex regions. Larger values make peak adhesions shorter-lived. |

The unbinding rule is one of the most important mechanisms for getting valley-directed migration. If `betaConc` is large, adhesions in valleys persist longer.

### Curvature-dependent binding parameters

| Parameter | Meaning |
|---|---|
| `pNeg` | Target binding probability for strong negative curvature. |
| `p0` | Binding probability near flat regions. |
| `pPos` | Target binding probability for strong positive curvature. |
| `kappaPar` | Controls how quickly curvature effects saturate. |
| `epsK` | Small numerical threshold used to avoid division by zero and ignore tiny curvature values. |

These parameters control how likely a membrane node is to form new adhesion sites.

### Curvature-dependent adhesion length parameters

| Parameter | Meaning |
|---|---|
| `meanAdhesionLength` | Baseline mean adhesion length. |
| `gammaLen` | Strength of curvature dependence in the adhesion length. |
| `lengthLognormalSigma` | Spread of the lognormal adhesion length distribution. |
| `minMeanAdhesionLength` | Safety lower bound so the mean length does not become zero or negative. |

New adhesion lengths are sampled from a lognormal distribution. The mean of that distribution is modified by curvature.

---

## `curvatureFromSurfaceDerivatives`

```cpp
inline double curvatureFromSurfaceDerivatives(double dzdx, double d2zdx2) {
    return -d2zdx2 / std::pow(1.0 + dzdx * dzdx, 1.5);
}
```

This function computes the signed curvature of a surface of the form

$$
z = f(x).
$$

The formula is

$$
\kappa = \frac{-z''(x)}{\left(1 + (z'(x))^2\right)^{3/2}}.
$$

The negative sign is important because we are using the same sign convention as the MATLAB version:

- `kappa > 0`: convex peak or ridge.
- `kappa < 0`: concave valley.
- `kappa ≈ 0`: nearly flat region.

For the surface

$$
z = a\cos(x/b),
$$

peaks have positive curvature and valleys have negative curvature.

---

## `surfaceCurvatureFromCallbacks`

```cpp
template <typename SurfaceDzDxFunction, typename SurfaceD2zDx2Function>
inline double surfaceCurvatureFromCallbacks(
    double x,
    double y,
    SurfaceDzDxFunction surfaceDzDxFunc,
    SurfaceD2zDx2Function surfaceD2zDx2Func
) {
    const double dzdx = surfaceDzDxFunc(x, y);
    const double d2zdx2 = surfaceD2zDx2Func(x, y);
    return curvatureFromSurfaceDerivatives(dzdx, d2zdx2);
}
```

This function connects the surface functions from the main file to the curvature formula in the header file.

The header does not know what surface you chose. Instead, the main file provides:

```cpp
double surfaceZ(double x, double y);
double surfaceDzDx(double x, double y);
double surfaceD2zDx2(double x, double y);
```

Then `surfaceCurvatureFromCallbacks` evaluates the first and second derivatives at the current node location and sends those values to `curvatureFromSurfaceDerivatives`.

This design is important because if you change the surface, you must also change the derivative functions in the main file.

---

## `curvatureModifiedUnbindDistance`

```cpp
inline double curvatureModifiedUnbindDistance(
    double kappa,
    double kappaScale,
    const AdhesionParams& adh
) {
    kappaScale = std::max(kappaScale, adh.epsK);

    const double r = std::abs(kappa) / (kappaScale + std::abs(kappa));
    double mult = 1.0;

    if (kappa < 0.0) {
        mult = 1.0 + adh.betaConc * r;
    } else {
        mult = 1.0 - adh.betaConv * r;
    }

    return std::max(adh.minMeanAdhesionLength, adh.charUnbindDist * mult);
}
```

This function computes the curvature-modified characteristic unbinding distance.

The adhesion unbinding probability has the form

$$
P_{unbind} = 1 - e^{-d/L},
$$

where:

- `d` is the distance between the membrane node and the adhesion site,
- `L` is the local unbinding distance.

If `L` is larger, then the adhesion is less likely to unbind. If `L` is smaller, then the adhesion is more likely to unbind.

The logic is:

- In valleys, `kappa < 0`, so the code increases `charUnbindDist`.
- Near peaks, `kappa > 0`, so the code decreases `charUnbindDist`.

This makes adhesions longer-lived in valleys and shorter-lived near peaks.

---

## `curvatureModifiedBindThreshold`

```cpp
inline double curvatureModifiedBindThreshold(
    double kappa,
    double kappaScale,
    const AdhesionParams& adh
) {
    const double k0 = std::max(adh.kappaPar * kappaScale, adh.epsK);

    if (std::abs(kappa) < adh.epsK) {
        kappa = 0.0;
    }

    const double r = std::abs(kappa) / (k0 + std::abs(kappa));

    double p = adh.p0;
    if (kappa >= 0.0) {
        p = adh.p0 + (adh.pPos - adh.p0) * r;
    } else {
        p = adh.p0 - (adh.p0 - adh.pNeg) * r;
    }

    return clampValue(p, 0.0, 1.0);
}
```

This function computes the probability threshold for forming a new adhesion site.

The logic is:

- Near flat regions, the binding probability is close to `p0`.
- For strong positive curvature, the binding probability moves toward `pPos`.
- For strong negative curvature, the binding probability moves toward `pNeg`.

With the default tuning

```cpp
adh.pNeg = 0.0;
adh.p0   = 0.5;
adh.pPos = 1.0;
```

positive curvature creates more new adhesions, while negative curvature creates fewer new adhesions.

This mechanism alone may bias motion toward convex regions. The full curvotaxis behavior depends on the balance between binding, unbinding, and adhesion length.

---

## `curvatureModifiedAdhesionLength`

```cpp
template <typename RNG>
inline double curvatureModifiedAdhesionLength(
    double kappa,
    double kappaScale,
    const AdhesionParams& adh,
    RNG& rng
) {
    kappaScale = std::max(kappaScale, adh.epsK);
    const double s = kappa / (kappaScale + std::abs(kappa));

    double meanLength = adh.meanAdhesionLength * (1.0 + adh.gammaLen * s);
    meanLength = std::max(meanLength, adh.minMeanAdhesionLength);

    const double sigma = adh.lengthLognormalSigma;
    const double mu = std::log(meanLength) - 0.5 * sigma * sigma;

    std::normal_distribution<double> normal01(0.0, 1.0);
    return std::exp(mu + sigma * normal01(rng));
}
```

This function samples a new adhesion length using a lognormal distribution.

The baseline mean length is `meanAdhesionLength`, but curvature modifies it through `gammaLen`.

The signed factor

$$
s = \frac{\kappa}{\kappaScale + |\kappa|}
$$

lies approximately between `-1` and `1`.

Then the mean length is modified by

$$
meanLength = meanAdhesionLength(1 + gammaLen \cdot s).
$$

Increasing `gammaLen` makes adhesion placement more sensitive to curvature.

---

## `outwardMembraneBindAngle`

```cpp
inline double outwardMembraneBindAngle(
    int i,
    const std::vector<double>& membraneX,
    const std::vector<double>& membraneY
) {
    double centerX = 0.0;
    double centerY = 0.0;

    const int N = static_cast<int>(membraneX.size());
    for (int k = 0; k < N; ++k) {
        centerX += membraneX[k];
        centerY += membraneY[k];
    }

    centerX /= static_cast<double>(N);
    centerY /= static_cast<double>(N);

    const double vx = membraneX[i] - centerX;
    const double vy = membraneY[i] - centerY;

    if (vx == 0.0 && vy == 0.0) {
        return 0.0;
    }

    return std::atan2(vy, vx);
}
```

This function makes new membrane adhesion sites point outward from the cell center.

The old version used a direction variable. That direction has now been removed because adhesion direction is determined by geometry, not by a manually chosen direction.

For membrane node `i`, the function computes the cell center by averaging all membrane node positions. Then it points from the center of the cell to the current membrane node.

This gives an outward adhesion direction.

---

## `addAdhesionForces`

The adhesion force function now takes the surface function and the first and second derivative functions.

The new version should be called like this:

```cpp
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
```

Notice that the old `direction` argument has been removed.

Inside the function, the adhesion update has three steps:

1. **Unbinding**
   - Compute curvature at each membrane node.
   - Modify the unbinding distance using `curvatureModifiedUnbindDistance`.
   - Randomly remove stretched adhesions.

2. **Binding**
   - Compute the curvature-dependent binding threshold using `curvatureModifiedBindThreshold`.
   - If a node has open adhesion slots, try to create new adhesion sites.
   - New sites are placed outward from the membrane node.

3. **Force accumulation**
   - Each attached adhesion site acts like a spring.
   - The force has the form

$$
F = k_{adh}(x_{bind} - x_{node}).
$$

The code uses ordinary 3D Euclidean distance for adhesion stretch. This is faster than the MATLAB arclength calculation and is better for this workshop version.

---

# Main File Changes

---

## Added `<chrono>` for Runtime Tracking

The main file now includes:

```cpp
#include <chrono>
```

This lets us measure how long the simulation takes.

At the start of `main`, we can write:

```cpp
auto startTime = std::chrono::high_resolution_clock::now();
```

At the end of the simulation, we can write:

```cpp
auto endTime = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsedSeconds = endTime - startTime;

std::cout << "Simulation runtime: "
          << elapsedSeconds.count() << " seconds\n";

std::cout << "Simulation runtime: "
          << elapsedSeconds.count() / 60.0 << " minutes\n";
```

Since the simulation is now more complex, it is good practice to begin with a small simulation before scaling up.

For reference:

- A simulation with `20,000` steps took about `20–30` seconds.
- A simulation with `200,000` steps took about `3` minutes.

This was for a single cell with a total of about `120` nodes. If you increase the number of membrane nodes, interior nodes, adhesion sites, or cells, the simulation time will increase accordingly.

---

## Surface Derivatives Added

The main file now explicitly defines the surface, its first derivative, and its second derivative.

For the surface

$$
z = a\cos(x/b),
$$

we use:

```cpp
// Curved substrate function
double surfaceZ(double x, double y) {
    (void)y;
    return surfacePar_a * std::cos(x / surfacePar_b);
}

// First derivative dz/dx
double surfaceDzDx(double x, double y) {
    (void)y;
    return -(surfacePar_a / surfacePar_b) * std::sin(x / surfacePar_b);
}

// Second derivative d2z/dx2
double surfaceD2zDx2(double x, double y) {
    (void)y;
    return -(surfacePar_a / (surfacePar_b * surfacePar_b)) * std::cos(x / surfacePar_b);
}
```

If you change the surface function, you must also change these derivative functions. The curvature calculation depends on them.

For example, if you replace the surface with a Gaussian valley, you need to write the correct first and second derivatives for that Gaussian surface.

---

## Cell Center Initialization Added

The cell is now initialized around a chosen center:

```cpp
const double x_c = PI * PI / 12.0;
const double y_c = 0.0;
```

The membrane nodes are created around this center:

```cpp
double x = x_c + R_membrane * std::cos(theta);
double y = y_c + R_membrane * std::sin(theta);
double z = surfaceZ(x, y);
```

The interior nodes are also created around this center:

```cpp
double x = x_c + r * std::cos(theta);
double y = y_c + r * std::sin(theta);
double z = surfaceZ(x, y);
```

This is important because if the cell starts exactly on a peak or in a perfectly symmetric position, the forces may balance and the cell may not move. Placing the cell slightly off-center breaks symmetry and allows the cell to start migrating.

---

## Adhesion Tuning Parameters Added in the Main File

The main file now includes the tuning parameters in one place so they are easy to modify:

```cpp
// Initialize bind site parameters
AdhesionParams adh;

adh.maxSubSitePerNode = 20;
adh.charUnbindDist = 0.1;
adh.meanAdhesionLength = 0.1;
adh.kadh = 2.0;

// Curvature-dependent unbinding parameters
adh.betaConc = 60.0;
adh.betaConv = 0.02;

// Curvature-dependent binding probability parameters
adh.pNeg = 0.0;
adh.p0   = 0.5;
adh.pPos = 1.0;

// Curvature scale tuning
adh.kappaPar = 0.5;

// Curvature-dependent adhesion length parameter
adh.gammaLen = 0.2;
```

These are the main parameters to modify when tuning the simulation.

---
# Runtime Scaling Notes

The expensive parts of the simulation are the pairwise force loops and adhesion updates.

The approximate cost per time step is

$$
O(N_m^2 + N_i^2 + N_mN_i + N_mS),
$$

where:

- `N_m` is the number of membrane nodes,
- `N_i` is the number of interior nodes,
- `S` is the number of adhesion sites per membrane node.

This means doubling the number of nodes can make the force calculations much more expensive. Increasing `N_steps` increases runtime almost linearly.

For debugging, start with something like:

```cpp
const int N_steps = 20000;
const int output_every = 200;
```

Then scale up to:

```cpp
const int N_steps = 200000;
const int output_every = 2000;
```

Both settings output about 100 frames, but the second simulation takes about 10 times more time steps.

---
