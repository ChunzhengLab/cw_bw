# cw_bw

- **Some issues are known to exist.**
- **So although it is under the MIT license, please contact the author before use.**

`cw_bw` is a lightweight simulation and fitting framework for investigating local baryon conservation (LBC) effects in heavy-ion collisions. It features a blast-wave event generator, with a standalone blast wave fitter on spectrum and elliptic flow.

## Project Structure

```
cw_bw-main/
├── CMakeLists.txt
├── include/                    # Core C++ headers
├── src/                        # Core implementations
├── configs/                    # YAML configuration files
├── refdata/                    # Experimental data (ROOT/CSV)
├── scripts/
│   ├── BlastWaveFitter/        # Blast-wave spectrum fitting module (C++ ROOT)
│   ├── FittedData/             # Model vs data comparison (Python)
│   └── ScanPars/               # Parameter scan & plotting (Python)
└── LICENSE
```

## Build Instructions

Requires C++17, ROOT, and CMake.

```bash
mkdir build && cd build
cmake ..
make -j
```

## Usage

### 1. Generate Events

```bash
./bwgen -C 55 -n 100000 -c ../configs/default.yaml -o results.root
```

Arguments:

- `-C`: Centrality class (e.g., 15, 25, 35, 45, 55)
- `-n`: Number of events
- `-c`: YAML config file
- `-o`: Output ROOT file

### 2. Compare with Data

```bash
python3 scripts/FittedData/compare_model_alice_pre.py
```

### 3. Parameter Scan and Plot

```bash
python3 scripts/ScanPars/analysis.py
```

### 4. Blast-Wave Spectra Fit

The `scripts/BlastWaveFitter/` directory contains a standalone ROOT macro for blast-wave fits:

```bash
root -l -b -q run_batch.C
```

This module fits identified particle spectra using the blast-wave model:

$$
\\frac{1}{2\\pi p_T}\\frac{d^2N}{dp_T\\,dy} \\propto \\int_0^R rdr \\int_0^{2\\pi} d\\phi \\, m_T I_0\\left(\\frac{p_T\\sinh\\rho}{T}\\right) K_1\\left(\\frac{m_T\\cosh\\rho}{T}\\right)
$$

You can adjust flow parameters like `T_kin`, `beta_T`, and `n` inside `BlastWaveFitter.C`.

## Reference Data

- ROOT files from [HEPData](https://www.hepdata.net/)
- CSV files from ALICE preliminary results

##  Config Examples

Available in `configs/`:

- `default.yaml`
- `no_v2.yaml`
- `default_badboost.yaml`

## License

This project is licensed under the MIT License. See [LICENSE](./LICENSE) for details.
