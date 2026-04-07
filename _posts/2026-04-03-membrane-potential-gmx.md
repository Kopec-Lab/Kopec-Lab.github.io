---
layout: post
title: Calculating membrane potential from GROMACS trajectories
date: 2026-04-05 12:00:00+0000
description: How the -correct flag works, its mathematical basis, and what the tool actually computes in simulations with an applied electric field.
tags: [membrane potential, GROMACS, electrostatics, molecular dynamics, computational electrophysiology]
categories: methods
related_posts: false
toc:
  sidebar: left

---

### Intro

If you have ever run an MD simulation of a lipid membrane, one of the analyses likely included calculating the membrane potential. In GROMACS, this is conveniently done using 'gmx potential'. However, what you likely observed was an asymmetric potential profile, varying wildly depending on the number of slices used for integration and the simulation length (as shown below; figures are made using our simple [gplot tool](https://github.com/Kopec-Lab/gplot)). 

{% include figure.liquid 
   path="assets/img/blog/potentials-raw-gmx.png" 
   caption="Figure 1: Membrane potentials calculated by gmx potential with a different number of slices (-sl)" 
   width="80%" %}

You may then have tried the -correct flag in the same tool, which appears to fix these issues, although it is still debated whether, and under what circumstances, it should be used. Below, I present some thoughts on this, together with a new method for calculating the membrane potential.

{% include figure.liquid 
   path="assets/img/blog/potential-correct-gmx.png" 
   caption="Figure 2: Membrane potentials calculated by gmx potential with a different number of slices (-sl) and with the -correct flag" 
   width="80%" %}


### The problem

When computing the electrostatic potential from an MD simulation, the standard recipe is:

1. Divide the simulation box into thin slices along one axis (typically Z)
2. Compute the charge density rho(z) in each slice
3. Integrate rho(z) twice to obtain the potential psi(z) via the Poisson equation

In the classical form (Eq. 2 in Gurtovenko 2009), the boundary conditions are psi(0) = 0 and E(0) = 0, and the integration proceeds from z = 0 to z = L (the box length):

```
E(z)   = (1/epsilon_0) * integral_0^z rho(z') dz'
psi(z) = -integral_0^z E(z') dz'
```

For a symmetric POPC bilayer (centered so the membrane is in the middle of the box and water is at the edges), this approach produces potential profiles that are **strongly dependent on the number of slices**, with large asymmetries that vary erratically. For example, in my hands, applying gmx potential without the -correct flag to a trajectory of a symmetric POPC bilayer produced:

```
Slices   Peak (V)   Asymmetry psi(L)-psi(0) (V)
  50      0.55       -0.08
 100      0.81       +0.47
 200      0.69       +0.17
 300      0.75       +0.25
 500      0.69       +0.13
 800      0.61       -0.00
1000      0.65       +0.06
```

For a system that is perfectly symmetric, psi(L) should equal psi(0). Instead, the potential drop across the box fluctuates between -0.08 and +0.47 V depending on how many slices are used. The peak potential varies between 0.55 and 0.81 V.

### Root cause: molecule splitting at the periodic boundary

After some investigation, it seems that the artifact originates from how atoms are assigned to slices when the system is periodic.

After centering (shifting so the bilayer center of mass is at z = L/2), each atom's z-coordinate is wrapped into [0, L) using modular arithmetic (`z % L`). This ensures all atoms fall within the box. However, molecules that straddle the wrapping boundary at z = 0 / z = L get split: some atoms of the molecule end up near z = 0 and others near z = L.

In my test simulation, approximately 400 water molecules are split this way in each frame. For a TIP3P water molecule at the boundary:

- Before wrapping: O at z = -0.3 A, H1 at z = 0.5 A, H2 at z = 0.7 A (molecule is intact)
- After wrapping: O at z = L - 0.3, H1 at z = 0.5, H2 at z = 0.7 (molecule is torn apart)

The oxygen (charge -0.834e) is now at the right edge of the box, while the hydrogens (+0.417e each) remain at the left edge. Each split molecule creates an artificial dipole spanning nearly the entire box. While the total charge of each molecule is zero, its charges are now distributed across opposite ends of the integration path.

When the Poisson equation is integrated from z = 0 to z = L, this artificial dipole contributes a linear drift to the potential. The magnitude of this drift depends on:

1. **The number of split molecules** -- varies between frames as molecules diffuse
2. **The exact positions of split atoms within their slices** -- changes with slice width, hence with the number of slices

This is why the artifact varies erratically with the slice count: different slice widths cause the split charges to land in different bins, changing the effective artificial dipole moment and thus the potential drift.

### What doesn't fix it

I tested several intuitive approaches and found that they were ineffective:

**Making molecules whole before wrapping.** Even when using a preprocessed trajectory where molecules are intact (e.g., via `gmx trjconv -pbc whole`), the centering step of gmx potential followed by `z % L` wrapping re-splits molecules at the new boundary. The artifact is identical to using a raw trajectory.

**Wrapping by molecule center of mass.** Instead of wrapping individual atoms, one can shift entire molecules based on their center of mass. This keeps molecules intact but atoms of molecules near the boundary extend beyond [0, L). These out-of-range atoms must be handled somehow:

- **Clipping** (placing them in the edge slice): Concentrates charge at the boundary, making the artifact much worse. At 1000 slices, the potential reached 47 V.
- **Modulo index wrapping** (`index % n_slices`): Effectively re-splits the molecule -- the out-of-range atoms land on the opposite side. Results are identical to simple atom wrapping.

**Scaling coordinates to a fixed box size.** Gurtovenko & Vattulainen recommend scaling all coordinates to the initial box size to avoid effects from box size fluctuations under NPT pressure coupling. Testing showed this has negligible effect: the artifact is dominated by molecule splitting, not box fluctuations. The asymmetry pattern was essentially unchanged with or without scaling.

**Placing the wrapping boundary through the membrane.** Moving the boundary from the water phase (where many molecules cross it) to the membrane center (where fewer cross). This made the artifact much worse, because lipid headgroups carry much larger partial charges than water, and splitting a lipid headgroup creates a correspondingly larger artificial dipole.

### Fix 1: Sachs et al. correction

Sachs et al. (J. Chem. Phys. 121, 10847, 2004) proposed an alternative form of the Poisson equation where the potential is constrained to be equal on both sides of the box: psi(0) = psi(L). This amounts to subtracting a linear correction term:

```
psi_Sachs(z) = psi_classical(z) - (z/L) * psi_classical(L)
```

Gurtovenko & Vattulainen (Eq. 4-11) showed that the difference between the two forms is proportional to the total dipole moment of the system. For electroneutral systems, the correction is:

```
Delta(z) = z / (epsilon_0 * V) * P_total
```

where P_total is the system dipole moment and V is the box volume.

The Sachs correction eliminates the linear drift and produces stable results:

```
Slices   Peak (V)   Asymmetry (V)
  50      0.56       -0.001
 100      0.58       +0.005
 200      0.61       +0.001
 300      0.63       +0.001
 500      0.62       +0.000
 800      0.61       -0.000
1000      0.62       +0.000
```

However, it has a fundamental limitation: **it removes ALL linear drift, including any real transmembrane potential difference**. For asymmetric bilayers (e.g., POPC/POPE), the transmembrane potential is a real physical quantity that should not be subtracted. The Sachs form is therefore only appropriate for symmetric systems.

### Fix 2: Dipole correction (subtracting the artifact only)

An alternative approach is to compute the artificial dipole moment caused by molecule splitting and subtract only that contribution:

1. For each frame, compute the dipole moment of the system with intact (unwrapped) molecules: P_whole = sum(q_i * z_i_whole)
2. Compute the dipole moment after wrapping: P_wrapped = sum(q_i * z_i_wrapped)
3. The artifact dipole is: P_artifact = P_wrapped - P_whole
4. Subtract the corresponding linear potential: Delta(z) = z / (epsilon_0 * V) * P_artifact

In principle, this should remove only the wrapping artifact while preserving the real transmembrane potential. In practice, testing showed this **overcorrects**: the "whole" and "wrapped" dipole moments are not cleanly separable because the centering operation itself changes the dipole moment relative to the box origin. The corrected asymmetry was worse than uncorrected for many slice counts.

### Fix 3: Fourier-space integration - our new tool

One idea was to solve the Poisson equation in Fourier (reciprocal) space, where periodicity is built in by construction. I tried the implementation here: [Fourier potential tool](https://github.com/Kopec-Lab/membrane-potential-tool)

The 1D Poisson equation in real space:

```
d^2 psi(z) / dz^2 = -rho(z) / epsilon_0
```

becomes, after Fourier transform:

```
-k^2 * psi_hat(k) = -rho_hat(k) / epsilon_0
```

which gives:

```
psi_hat(k) = rho_hat(k) / (k^2 * epsilon_0)    for k != 0
```

The k = 0 component corresponds to an arbitrary overall offset of the potential and is set to zero. The electric field is obtained similarly:

```
E_hat(k) = -i*k * psi_hat(k)
```

The inverse Fourier transform gives the real-space potential and field.


{% include figure.liquid 
   path="assets/img/blog/potential-fourier.png" 
   caption="Figure 3: Membrane potential calculated by our tool with the Fourier-space integration" 
   width="80%" %}


This approach has several nice advantages:

1. **Inherently periodic.** The Fourier transform assumes the input signal is periodic. The potential automatically satisfies psi(0) = psi(L), with no need for corrections or special boundary conditions. There is no "starting point" for integration that could accumulate errors.

2. **Immune to the boundary artifact.** Split molecules create charge density features at both edges of the box. In the classical method, the double integration treats these as separated by the full box length, amplifying their effect. In Fourier space, periodicity means the edges are the same point: the charge distribution from a split molecule is treated as if the molecule were whole, because the Fourier basis functions wrap around naturally.

3. **Insensitive to the number of slices.** Because there is no error accumulation from sequential integration, the result converges rapidly and varies minimally with slice count:

```
Slices   Peak (V)   Asymmetry (V)
  50      0.65       +0.003
 100      0.63       +0.003
 200      0.64       +0.001
 300      0.64       +0.001
 500      0.63       +0.000
 800      0.63       +0.000
1000      0.63       +0.000
```

The peak potential is consistent at ~0.63 V regardless of slice count, and the asymmetry is essentially zero (< 0.003 V). Compare this with the classical method's range of 0.55-0.81 V for the peak and up to 0.47 V of asymmetry.

### How the tool works

### Step 1: Charge density

The box is divided into `N` slices along the chosen axis. For each trajectory frame, every atom is assigned to a slice based on its coordinate. The partial charge (from the `.tpr` file) is added to that slice. The charge density is the total charge in each slice divided by the slice volume, averaged over all frames. Conveniently, the tool automatically detects the optimal number of slices so there is no need of fiddling with these anymore.

### Step 2: Solve the Poisson equation

**Fourier method (default):** The charge density is Fourier-transformed, divided by `k^2 * epsilon_0` in reciprocal space (skipping k=0), and inverse-transformed to obtain the potential. The electric field is computed analogously using `-i*k * psi_hat(k)`.

**Classical method (`-classical`):** The charge density is numerically integrated twice using the cumulative trapezoidal rule, with boundary conditions psi(0) = 0 and E(0) = 0. With `-sachs`, a linear correction term is subtracted to enforce psi(0) = psi(L).

Both methods use the `gmx potential` sign convention, where the dipole potential inside the membrane is positive relative to water.

### Known issues and considerations

- **Centering matters**: Without `-center`, bilayer center-of-mass fluctuations smear the charge density and produce asymmetric profiles. Always use `-center` for membrane systems.
- **Convergence**: Potential profiles require long trajectories (100+ ns for bilayers) to converge, particularly in the water region near lipid headgroups.
- **Dielectric constant**: The calculation uses vacuum permittivity (epsilon_r = 1). The output represents the bare electrostatic potential from the explicit charge distribution, not a dielectrically screened quantity.
- **Undulations**: For large membranes with significant undulations, the 1D slab decomposition may not accurately capture the local potential. Using small, planar bilayer patches avoids this issue.
- **Preprocessing**: It is recommended (but not strictly necessary) to preprocess the trajectory with `gmx trjconv -pbc whole` to make molecules whole. While the Fourier method handles split molecules much better than the classical method, intact molecules produce marginally cleaner charge density profiles.







### What gmx potential `-correct` actually does

Inspection of the GROMACS source code reveals that `-correct` performs **two mean subtractions** during the classical double integration:

1. **Before the first integration**: subtract the mean charge density from all slices that have nonzero charge density. This makes the effective total charge exactly zero, removing a constant drift term from the electric field.

2. **After the first integration**: subtract the mean electric field (again, only from slices with nonzero charge density). This removes a constant offset from E(z), which would otherwise produce a linear drift in the potential.

I believe this is *not* the same as "zeroing net charge per slab" as sometimes described. It is a global correction applied to the entire charge density profile and then the entire electric field profile.

The key insight is that Step 2 (mean E subtraction) is doing essentially the same thing as the Sachs correction and the Fourier method's k=0 exclusion, but from a different angle.

In a periodic system, the electric field must also be periodic: E(0) = E(L). The classical integration starts with E(0) = 0, so it requires E(L) = 0 as well. If the total integrated charge density is not exactly zero (due to molecule splitting at the boundary, finite precision, or slight charge imbalance from the slicing), E(L) != 0, and this mismatch propagates as a linear drift in the potential.

Subtracting the mean electric field shifts E(z) so that its average is zero. For a periodic function, this is equivalent to enforcing that the integral of E over the full box vanishes, which is exactly the condition for the potential to be periodic (psi(0) = psi(L)).

Step 1 (mean rho subtraction) has a negligible effect in practice because the total system charge is already very close to zero. Its main effect is to eliminate tiny numerical imbalances from the binning.

### Test results

Testing on a symmetric POPC bilayer with the `-correct` flag produces results nearly identical to the Fourier method:

```
Classical + correct:
Slices   Peak (V)   Asymmetry (V)
  50      0.56       +0.002
 100      0.61       -0.003
 200      0.63       +0.001
 300      0.63       +0.001
 500      0.63       +0.000
 800      0.62       +0.000
1000      0.63       +0.000
```

Compare with Fourier (peak ~0.63 V, asymmetry < 0.003 V) and uncorrected classical (peak 0.55-0.81 V, asymmetry up to 0.47 V). The `-correct` flag effectively eliminates the slice-count sensitivity.

The `-correct` procedure is mathematically equivalent to enforcing periodic boundary conditions on the electric field, which is physically correct for a periodic system. In this sense, I believe it is a valid fix.

---

### Simulations with applied electric field

A common approach to simulate a transmembrane voltage in MD is to apply a constant external electric field E perpendicular to the membrane (GROMACS `electric-field-z` mdp parameter). The resulting voltage is V = E * L_z, where L_z is the box length in the field direction. After spending many years trying to fully understand the physics here, this section explains what the tool computes in this scenario, how to interpret the results, and what the `-efield` flag in our new tool does.

### Background: the constant electric field method

In GROMACS, a constant electric field E_z applies a force F = q_i * E_z to every charged particle. The applied voltage across the simulation box is:

```
V = E_z * L_z
```

This is equivalent to the influence of two infinite baths connected by an electromotive force (EMF), as shown by Roux (Biophys. J. 95, 4205, 2008). The constant field method has been validated extensively (Gumbart et al., BBA 1818, 294, 2012; Jensen et al., J. Gen. Physiol. 141, 619, 2013) and is widely used for studying ion channels, electroporation, and voltage-gated conformational changes.

### Decomposition of the total potential

The total electrostatic potential in a simulation with an applied field has two components (Jensen et al. Eq. 2):

```
Phi_total(z) = Phi_reaction(z) - E_applied * z
```

where:

- **Phi_reaction(z)** is the reaction potential arising from the explicit charges in the system (ions, water dipoles, lipid headgroups, protein residues). This is computed by the particle-mesh Ewald (PME) method during the simulation and is what our tool calculates from the charge density by solving the Poisson equation. **This potential is periodic**: Phi_reaction(0) = Phi_reaction(L).

- **-E_applied * z** is the linear ramp from the applied external field. This component has no source charges: it is imposed externally. **This breaks periodicity**: the total potential differs by V = E * L_z between the two ends of the box.

The charge density rho(z) that the tool computes from atom positions reflects the system's *response* to the applied field: the rearrangement of ions, water dipoles, and other charges. It does not contain the external field itself. Therefore, the Poisson equation solved by the tool yields only the reaction potential, which is periodic. The Fourier method remains the correct solver for this component as shown below.

### What happens in the bulk water

Bulk aqueous electrolyte behaves approximately as a conductor: it self-organizes to screen the applied field. In the water region:

```
E_total = E_reaction + E_applied ≈ 0
```

Therefore the reaction electric field in water is approximately:

```
E_reaction ≈ -E_applied
```

and the reaction potential develops a slope in the water region:

```
dPhi_reaction/dz ≈ E_applied    (in water)
```

The total potential, being the sum of the sloped reaction potential and the opposite-sloped applied ramp, is approximately flat in bulk water, as expected for a conductor. The entire voltage drop concentrates across the low-dielectric membrane region, which is the biologically meaningful quantity.

### Recovering the voltage from the charge density

One question I wanted to test was if the applied voltage V can be recovered purely from the charge density profile, without knowing E_applied a priori. There should be two ways to do that:

**From the electric field profile**: The average reaction electric field in the bulk water region should be approximately -E_applied. Reading this value and multiplying by L_z gives an estimate of V.

**From the potential slope**: The reaction potential has a slope ≈ E_applied in the water region. Fitting a line to the potential in the water phase and multiplying the slope by L_z gives V.

However, these are **approximate** estimates. I tested it on a gramicidin A channel system (E_applied = 0.032 V/nm, expected V = 293 mV) showed.


```
Reaction potential slope in water (per region):
  Region 1 (left water):  slope = 0.0411 V/nm
  Region 2 (right water): slope = 0.0137 V/nm
  Average slope:           0.0274 V/nm
  -> V estimate:           251 mV (85.7% of expected 293 mV)

Average reaction E field in water: -0.0278 V/nm
Expected (-E_applied):             -0.0320 V/nm
Average total E field in water:     0.0042 V/nm (ideal: 0)
```

The ~85% recovery is consistent and reproducible. The underestimate reflects the fact that water is not a perfect conductor:

1. **Finite water phase**: The bulk water regions (~2 nm on each side) may not be thick enough for complete screening.
2. **Ionic current**: In a conducting channel (like gramicidin), steady-state ionic current requires a residual driving field in the water, so the screening is inherently incomplete.
3. **Finite system size**: Gumbart et al. (2012) emphasize that V = E * L_z is exact by construction, but extracting V from the potential profile introduces finite-size errors.

**Therefore, the exact voltage is always V = E * L_z.** The slope-based estimate is useful as a consistency check but should not be used as the primary voltage determination.

### What the `-efield` flag does in our tool

When `-efield VALUE` is specified (where VALUE is E_z in V/nm, matching the GROMACS `electric-field-z` mdp parameter):

1. **Reports the exact voltage**: V = E_z * L_z, using the average box length from the trajectory.

2. **Writes the total potential**: Phi_total(z) = Phi_reaction(z) - E * z is written to the file specified by `-ot` (default: `potential_total.xvg`). This non-periodic profile shows the voltage drop across the membrane.

3. **Detects bulk water regions**: Water molecules are automatically identified (by residue name: SOL, TIP3, HOH, WAT, SPC) and their positions are binned during trajectory processing. Contiguous regions with high water density are identified as bulk water.

4. **Measures the reaction potential slope**: A linear fit is performed in each bulk water region independently (since different regions are at different absolute potential levels in the periodic reaction potential), and the average slope is reported.

5. **Reports voltage estimates**: Both the exact V = E * L_z and the slope-based estimate are printed, along with the recovery percentage. The average reaction electric field in water is also reported and compared to the expected -E_applied.

{% include figure.liquid 
   path="assets/img/blog/potential-fourier-efield.png" 
   caption="Figure 4: Membrane potential calculated by our tool with the Fourier-space integration in a system with applied electric field" 
   width="80%" %}


### Classical integration with `-correct` also works

Perhaps surprisingly, `-classical -correct` (and gmx potential with -correct) produces the same water slope and voltage estimate as the Fourier method. Testing on the gramicidin system shows identical slopes (0.0270 V/nm) for both methods.

The reason: `-correct` subtracts the global mean E field across all slices. This mean (-0.051 V/nm in the test system) reflects contributions from both the water regions (where E_reaction ≈ -0.029 V/nm) and the membrane region (where E_reaction ≈ +0.017 V/nm, reflecting the large headgroup dipole fields). Because these contributions partially cancel, the global mean is dominated by the integration drift artifact, not by the physical field from the applied voltage. Subtracting it removes only the artifact, preserving the spatial structure — exactly as the Fourier method does by setting k=0 to zero.

```
                           Uncorrected    After -correct    Fourier
Mean E, all slices:        -0.0512         0.0000           0.0000  V/nm
Mean E, left water:        -0.0826        -0.0314          -0.0314  V/nm
Mean E, right water:       -0.0775        -0.0263          -0.0263  V/nm
Mean E, membrane:          -0.0338        +0.0174          +0.0174  V/nm
```

Both corrected methods give identical E field profiles. The `-sachs` correction, by contrast, enforces psi(0) = psi(L) by subtracting a linear term, which would remove the physical slope and should not be used with applied fields.

### Practical considerations

- **Always use `-center`** when analyzing systems with applied fields. Without centering, membrane drift smears the charge density and corrupts the slope measurement.

- **The total potential (`-ot` output) is not periodic.** It shows a linear ramp across the box with the voltage drop concentrated at the membrane. This is the physically meaningful potential profile.

- **The reaction potential (`-o` output) is periodic.** It shows the intrinsic membrane dipole potential plus the system's response to the applied field. The slope in the water region reflects the screening of the applied field.

- **Both Fourier and `-classical -correct` work correctly** with applied fields. Both methods remove only the unphysical DC offset (k=0 mode) from the electric field while preserving all spatial structure, including the slope in water that carries the voltage information. The `-correct` procedure subtracts the global mean E field across all slices — this mean (-0.051 V/nm in the gramicidin test) is NOT the same as the E field in water (-0.029 V/nm), because the membrane region has a very different E field that pulls the average away. Therefore, the physical slope is preserved. The `-sachs` correction should still be avoided with applied fields, as it explicitly enforces psi(0) = psi(L), which would remove the linear component of the reaction potential that encodes the voltage.

### Why the tool computes the reaction potential, not the total potential

A potentially confusing aspect of electrostatic potential calculations from MD simulations — both in `gmx potential` and in this tool — is that what gets computed is **the reaction potential**, not the total electrostatic potential. It has confused me for a long time.

**What the tool does:** It bins atomic charges into slices, constructs the charge density rho(z), and solves the Poisson equation:

```
d²Phi/dz² = -rho(z) / epsilon_0
```

The charge density rho(z) contains only the explicit charges in the system: ions, water partial charges, lipid headgroup charges, protein residues, etc. The resulting potential Phi(z) is the **reaction potential** — the electrostatic potential generated by the system's own charges in response to whatever conditions are imposed.

**Why not the total potential?** When an external electric field is applied in GROMACS (`electric-field-z`), the field acts as a force F = q * E_z on every particle, but it has **no source charges**. It is an external boundary condition, not a charge distribution. Since the applied field contributes no charges to rho(z), it does not appear in the Poisson equation, and the tool cannot "see" it. The total potential is:

```
Phi_total(z) = Phi_reaction(z) - E_applied * z
```

The linear ramp from the applied field must be added manually (which is what the `-efield` flag does).

**Why this hasn't caused problems historically:** In the vast majority of membrane simulations — symmetric bilayers, asymmetric bilayers, systems with ion gradients — there is no applied external field. In these cases, Phi_reaction(z) *is* the total potential. The distinction between reaction and total potential simply doesn't arise, and the output of `gmx potential` (with the -correct flag) is the complete answer.

**Where confusion can arise:** When we apply an external field (e.g., to study voltage-gated channels or electroporation) and then run `gmx potential` on the trajectory, we get the reaction potential, which is periodic and shows no net voltage drop across the box. This can be mistaken for a calculation error or a sign that the applied voltage had no effect. In fact, the tool is working correctly; the voltage drop is carried entirely by the applied field ramp, which is not part of the charge density. The reaction potential shows the system's *response* to the applied voltage: screening in the water phase (visible as a slope) and the intrinsic membrane dipole potential.

**To obtain the total potential** in simulations with an applied field, use `-efield VALUE` in our tool where VALUE matches the `electric-field-z` parameter from the GROMACS mdp file. The tool will then write the total potential (reaction + applied ramp) to `potential_total.xvg`, which shows the non-periodic voltage profile with the full voltage drop across the membrane.

### Why is the reaction potential periodic?

At first glance, it seems paradoxical that the reaction potential Phi_reaction(z) is periodic even when a directional electric field is applied. One might expect charges to accumulate on one side of the membrane — as in a real capacitor — breaking the symmetry and making the charge density (and hence the potential) non-periodic. But this does not happen, for a fundamental reason: **periodic boundary conditions prevent true charge accumulation**.

**There is only one water phase.** In a membrane simulation with PBC, the water region above the membrane and the water region below it are not separate compartments — they are the same continuous water phase, connected through the periodic boundary. A water molecule (or ion) that exits the box at z = L_z immediately re-enters at z = 0. There are no edges, walls, or capacitor plates where charge can build up.

**Charge density inherits the periodicity of the box.** Because every atom's position is defined modulo L_z, the ensemble-averaged charge density satisfies rho(z + L_z) = rho(z) by construction. No matter how strong the applied field, the charge density profile repeats with period L_z. The system responds to the field through *local* polarization: water dipoles align preferentially, ions redistribute slightly within the periodic cell — but these rearrangements are periodic.

**Poisson's equation preserves periodicity.** The reaction potential is obtained by solving the Poisson equation with the periodic charge density as input. If rho(z) is periodic, Phi_reaction(z) is also periodic (the Fourier method enforces this explicitly; the classical method with `-correct` achieves it by removing the DC drift). The reaction potential develops a slope in the water regions (reflecting the screening of the applied field) but this slope is compensated by the opposite response in the membrane region, so that the potential returns to its starting value over one full period.

**How does the system sustain a voltage without charge accumulation?** This is the key insight from Roux (2008). In PBC, the applied voltage cannot be represented by accumulated surface charge (there are no surfaces). Instead, it enters as an external parameter — the EMF in Roux's formulation — that adds a force F = q * E_z to every charged particle. The system's response is captured by the *displacement charge* Q_d = sum(q_i * z_i) / L_z, which measures how much charge has been displaced along z relative to the periodic cell, without ever breaking periodicity. The displacement charge is analogous to the charge on a capacitor plate, but it arises from the collective shift of all charges within the periodic cell, not from accumulation at a boundary.

**In summary:** PBC means the simulation box has no boundaries, only a periodic unit cell. Charges cannot accumulate at edges that don't exist. The charge density and reaction potential are therefore necessarily periodic. The entire voltage drop is encoded in the external field ramp (-E * Lz), which must be added manually to obtain the non-periodic total potential.

### Special case: pure water under applied electric field

Testing the tool on a pure water box with and without an applied electric field reveals a counterintuitive but I believe a correct result: **the reaction potential is identical in both cases**, and the total potential under the applied field is a perfectly linear ramp.

**Why the charge density doesn't change:** The applied field does orient water dipoles: each molecule tilts slightly so its dipole aligns with E_z. However, this polarization is *spatially uniform*: every slice of the water box has the same degree of alignment. Uniform polarization produces no charge density variation:

```
rho_bound = -dP/dz = 0    (when P is constant)
```

Microscopically: if all hydrogens shift slightly "up", each slice gains H charge from molecules below but loses the same amount to the slice above. For a uniform density of molecular centers, these contributions cancel exactly in every slice. The net dipole moment is real, but it is invisible in rho(z).

**Why there is no screening:** In a finite dielectric slab (or at a membrane interface), the edges of the polarized region carry bound surface charges (sigma_b = P dot n) that create a depolarization field screening E_applied by a factor of 1/epsilon_r. In PBC, there are no surfaces — the "+sigma" at z = L_z wraps around and cancels the "-sigma" at z = 0. Without bound surface charges, there is no screening:

| System | E_reaction in water | E_total in water | Screening? |
|--------|-------------------|------------------|------------|
| Pure water (PBC) | 0 | E_applied | No |
| Water + membrane | ~ -E_applied | ~ 0 | Yes |

This means that in pure water with PBC, the total electric field equals the full applied field everywhere. The dielectric constant of water is irrelevant — epsilon_r only matters when interfaces are present for bound charges to accumulate on. The voltage V = E * L_z is exact regardless of epsilon_r, which is a fundamental property of the constant-field method in PBC (Roux, 2008).

**What this means in practice:** The perfectly linear total potential ramp in pure water is NOT because water is "a good conductor that screens the field" (leaving only the ramp). It is the opposite — water *cannot* screen the field at all in a homogeneous system with PBC, so the reaction potential is flat and the total potential is the bare applied ramp. Screening requires interfaces (like a membrane) where the polarization changes spatially, creating bound charges that generate the reaction field.

**Connection to uncompensated dipole artifacts:** This physics also explains the PBC-induced water ordering artifact described by us recently (Kopec & Gapsys in revision): when a membrane has an uncompensated electric dipole (from asymmetric charged lipids or an embedded protein), PBC constrains the voltage across the box to zero. The system compensates by creating a counter-field — water dipoles order and ions redistribute. This ordering is an artifact of the periodic boundary, not a physical response. Adding salt (mobile ions) alleviates it by providing charge carriers that can redistribute to screen the membrane dipole without forcing bulk water to order.

### References

- B. Roux, "The membrane potential and its representation by a constant electric field in computer simulations," *Biophys. J.* **95**, 4205-4216 (2008).
- J. Gumbart, F. Khalili-Araghi, M. Sotomayor, B. Roux, "Constant electric field simulations of the membrane potential illustrated with simple systems," *BBA - Biomembranes* **1818**, 294-302 (2012).
- M.O. Jensen, V. Jogini, M.P. Eastwood, D.E. Shaw, "Atomic-level simulation of current-voltage relationships in single-file ion channels," *J. Gen. Physiol.* **141**, 619-632 (2013).
- W. Kopec, V. Gapsys, "Periodic boundaries in molecular dynamics simulations: why do we need salt?," (2026).
