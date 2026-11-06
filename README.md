# From Data to Control: Companion Code

Companion code for the paper:

> **From Data to Control: A Formal Compositional Framework for Large-Scale Interconnected Networks**

From a **single noisy input–state trajectory** of one representative subsystem with an **unknown model**, the synthesis scripts construct a control storage certificate (CSC) `S_i(x_i) = x_i' P_i x_i` and a local polynomial controller via SOS programming, verify the dissipativity-based compositionality condition (LMI) over the full network, and compute the initial/unsafe level sets `(eta_i, mu_i)`. The remaining scripts reproduce the paper's simulations, complexity figure, and the comparison against a small-gain compositionality approach.

## Contents

- `Spacecraft_Fully_Network.m`
- `Spacecraft_Fully_Simulation.m`
- `Lorenz_Fully_Network.m`
- `Lorenz_Fully_Simulation.m`
- `Complexity_Comparison.m`
- `Signed_Lorenz_Comparison.m`

## Computational environment

Our simulations were executed in **MATLAB R2023b** on a **MacBook M2** with **32 GB RAM**.
The **RT** and **MU** values reported in the paper refer only to the synthesis scripts (`*_Network.m`). In contrast, the network simulation scripts (`*_Simulation.m`) for generating the corresponding figures may take **considerably longer** to finish, depending on the number of subsystems and the available computational resources.

> **Note** The reported **running times (RT)** and **memory usage (MU)** correspond to the time and memory needed to solve the SOS problems for each subsystem, as detailed in **Algorithm 1** of the paper; they are measured before any network-level matrix is formed.
>
> **Note on memory units**
> MATLAB's `whos` reports usage in **bytes**. For clarity, the scripts convert these values to **megabytes (MB)** before printing them.

## Prerequisites & Installation

### Supported operating systems
Any relatively recent version of **MATLAB** on any OS should work, but for consistent results, we recommend the version used in our tests: **MATLAB R2023b**.

### Required optimization toolboxes & solvers
Add the following packages to your MATLAB path:

- **SOSTOOLS** 4.03
- **MOSEK** 10.2.3

(`Spacecraft_Fully_Simulation.m`, `Lorenz_Fully_Simulation.m`, and `Complexity_Comparison.m` need only base MATLAB.)

### Setup steps
1. Download
   - **SOSTOOLS** <https://github.com/oxfordcontrol/SOSTOOLS/archive/refs/tags/v4.03.zip>
   - **MOSEK**    <https://www.mosek.com/downloads/>
2. Extract each `.zip` archive.
3. In MATLAB, open **Set Path… → Add with Subfolders…** and add the extracted folders.
4. Save the updated path and restart MATLAB (optional but recommended).

> **Note** Both packages **must** be on the MATLAB path before running `Spacecraft_Fully_Network.m`, `Lorenz_Fully_Network.m`, or `Signed_Lorenz_Comparison.m`.

### MATLAB Path Instructions

Follow the steps below **for each toolbox** (e.g., **SOSTOOLS 4.03**):

1. Open **MATLAB**.
2. Go to the **Home** tab in the MATLAB toolbar.
3. Click **Set Path** in the **Environment** section.
4. In the dialog that appears, select **Add with Subfolders…**.
5. Locate and select the extracted toolbox folder (after unzipping, if applicable).
6. Click **Save** to confirm the changes.

### Alternative method

1. Place the toolbox folder in the same location as your other `.m` files.
2. In MATLAB's **Current Folder** browser, right-click each folder and choose
   **Add to Path → Selected Folders and Subfolders**.

### MOSEK setup

1. **Get a license**
   - Request a free academic license → [MOSEK Academic Licenses](https://www.mosek.com/products/academic-licenses/)
   - Or request a time-limited trial license → [MOSEK Trial Licenses](https://www.mosek.com/products/trial/)
   - Watch for an email (check your spam folder) containing the `mosek.lic` file.

2. **Save the license file**
   - Download the attached `mosek.lic` from the email.

3. **Place `mosek.lic` in the correct folder**

   | OS       | Destination path (create the folder if it doesn't exist) |
   |----------|-----------------------------------------------------------|
   | Linux    | `/home/<YOUR_USER_NAME>/mosek/mosek.lic`                   |
   | macOS    | `/Users/<YOUR_USER_NAME>/mosek/mosek.lic`                  |
   | Windows  | `C:\Users\<YOUR_USER_NAME>\mosek\mosek.lic`                |

4. **Install MOSEK binaries**
   Follow the OS-specific installation guide in the MOSEK documentation: <https://docs.mosek.com/latest/install/installation.html>.

### How to run the code

To regenerate the results presented in the paper:

1. **Download** the relevant script.
2. **Open MATLAB**.
3. **Add SOSTOOLS and MOSEK** (and their subfolders) to your MATLAB path.
4. **Run** the script by name in the Command Window, e.g.:

   ```matlab
   Spacecraft_Fully_Network
   ```

## Reproducibility

- `rng(1)` fixes the collected data and initial conditions in every script. Certificate coefficients may differ in the last digits across solver versions and platforms.
- The sampled ranges of the internal-input data `w` are immaterial for the soundness of the certificates: `w` enters the data equations exactly.
