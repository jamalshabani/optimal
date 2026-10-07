# FireTop — Topology Optimization with Firedrake and PETSc/TAO

FireTop is Jamal Shabani's research code for the computational design of compliant morphing structures. It contains finite-element and phase-field optimization experiments involving structural and responsive materials, stimulus, and prescribed motion.

The implementation is organized into individual Python solvers and experiment launchers. Each solver defines its own model, boundary conditions, objective, and optimization settings.

## Research scope

The repository explores:

- Two-material and three-material topology optimization.
- Phase-field regularization of material distributions.
- Responsive materials and stimulus-driven deformation.
- Compliant motion and blocking-load experiments.
- Numerical optimization using PETSc/TAO.
- Adjoint-based formulations in selected experiments.

## Related repositories

| Repository | Purpose |
| --- | --- |
| [optimal](https://github.com/jamalshabani/optimal) | Main FireTop research code and numerical experiments. |
| [first_paper](https://github.com/jamalshabani/first_paper) | Computations used for the first paper and doctoral thesis. |
| [trajectory](https://github.com/jamalshabani/trajectory) | Ongoing work on trajectory topology optimization. |

## Repository structure

| Files | Description |
| --- | --- |
| `2_mat_prob.py`, `2_mat_res_prob.py` | Two-material experiments. |
| `2_mat_prob_adjoint.py` | Adjoint variant using `firedrake_adjoint` and `pyadjoint`. |
| `3_mat_prob.py`, `3_mat_prob_explicit.py`, `3_mat_res_prob.py` | Three-material formulations and variants. |
| `2_mat_motion_stimulus_prob.py`, `3_mat_motion_*.py` | Motion and stimulus experiments. |
| `blocking_load_2_mat.py` | Blocking-load experiment. |
| `motion_problem.py`, `motion_nof_problem.py` | Additional motion formulations. |
| `run_*.py` | Experiment launchers containing parameter combinations. |
| `*.msh` | Input meshes, including beam aspect-ratio variants. |

## Dependencies

The main numerical dependencies are:

- Python 3.
- Firedrake.
- PETSc with TAO, accessed through `petsc4py`.
- NumPy.
- A VTK-compatible viewer, such as ParaView, for inspecting results.

The adjoint variant additionally requires compatible versions of `firedrake_adjoint` and `pyadjoint`.

Run the solvers inside a compatible Firedrake environment. Check the core imports with:

```bash
python3 -c "import firedrake, numpy; from petsc4py import PETSc; print(PETSc.Sys.getVersion())"
```

The scripts include uses of Firedrake's `File` API. Compatibility depends on the installed Firedrake version; a fully pinned computational environment is not specified here.

## Getting started

Run commands from the repository root so relative mesh paths resolve correctly.

Inspect the available options:

```bash
python3 2_mat_res_prob.py --help
```

Run a two-material example:

```bash
mkdir -p runs/two-material-example

python3 2_mat_res_prob.py \
  -m 1_to_1_mesh.msh \
  -o runs/two-material-example \
  -tao_monitor \
  -tao_max_it 20 \
  -er 0.1 \
  -es 1.0 \
  -l -0.005 \
  -k 1.0e-5 \
  -e 4.0e-3 \
  -p 1.0
```

This example uses a parameter set from `run_2_mat_res_prob.py` with a reduced iteration limit. It is intended as an initial execution example, rather than a converged reproduction of a published result.

### Selected parameters

The following options apply to `2_mat_res_prob.py`:

| Option | Meaning |
| --- | --- |
| `-m`, `--mesh` | Input mesh file. |
| `-o`, `--output` | Output directory. |
| `-er`, `-es` | Responsive and structural elastic moduli. |
| `-l` | Lagrange multiplier. |
| `-k` | Modica–Mortola regularization weight. |
| `-e` | Phase-field regularization parameter. |
| `-v` | Responsive-material volume fraction parameter. |
| `-p` | Elasticity interpolation exponent. |
| `-tao_max_it` | Maximum optimization iterations. |
| `-tao_monitor` | Display optimizer progress. |

Options and defaults vary between scripts. Inspect the relevant argument parser before adapting a command to another formulation.

## Output and visualization

The example writes results under the specified output directory, including:

- `rho_initial.pvd`.
- Intermediate density fields.
- `final-rho.pvd`.
- `displacement.pvd`.

Open the PVD collections in ParaView or another VTK-compatible viewer. Keep the associated data files together with the PVD files.

Use a separate output directory for every experiment.

## Reproducibility

For each computation, record:

1. The repository commit.
2. The exact command and parameter values.
3. The input mesh.
4. Python, Firedrake, PETSc, and NumPy versions.
5. Solver logs and convergence information.

Retain the supplied mesh boundary labels. Replacing a mesh requires checking the boundary conditions in the source.

Several solvers create TAO on `PETSc.COMM_SELF`; review the implementation before attempting distributed execution.

The repository includes substantial simulation output, so a complete download may require significant storage. For computations associated with the paper and thesis, see the companion `first_paper` repository.

## Publications

- J. Shabani, K. Bhattacharya, and B. Bourdin, *Systematic Design of Compliant Morphing Structures: A Phase-Field Approach*. [Preprint](https://arxiv.org/abs/2411.06289).
- Jamal Shabani, *Systematic design of compliant morphing structures with stimulus as design and state variable*, doctoral thesis, McMaster University. [Thesis Publication Link](https://macsphere.mcmaster.ca/items/8108f6e2-92a0-4529-a892-baa3c2526f2d).

If you use this code in research, cite the associated publications and identify the repository commit used.

## Author

**Jamal Shabani**

[GitHub](https://github.com/jamalshabani)
