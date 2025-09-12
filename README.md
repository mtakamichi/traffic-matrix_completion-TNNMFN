# traffic-matrix_completion-TNNMFN
Official implementation of "Completion of Traffic Matrix by Tensor Nuclear Norm Minus Frobenius Norm Minimization and Time Slicing" IEEE NOMS Workshop 2024.


> **TL;DR**: This repo completes missing entries in Internet traffic matrices by combining **time slicing** with a **tensor nuclear norm minus Frobenius norm (TNMF)** regularizer, solved via a specialized ADMM routine. It achieves accuracy comparable to (or better than) recent methods **without any training data**.

## Paper
Miyata, T. “Completion of Traffic Matrix by Tensor Nuclear Norm Minus Frobenius Norm Minimization and Time Slicing,” **IEEE NOMS Workshop 2024**.  
DOI: 10.1109/NOMS59830.2024.10575433  
IEEE Xplore: https://ieeexplore.ieee.org/document/10575433

---

## Method at a Glance
- **Goal**: Recover missing entries in a spatio-temporal **traffic matrix (TM)**.
- **Time Slicing (TS)**: Split the TM along time into a 3-way tensor to better capture temporal periodicity.
- **Regularizer**: **Tensor Nuclear Norm − α · Frobenius Norm (TNMF)** over mode-unfoldings, a simple non-convex surrogate that avoids manual singular-value weighting.
- **Solver**: A specialized **ADMM** with a closed-form proximal operator for the NNMF term on each unfolding.
- **Why it works**: Encodes low-rank spatio-temporal structure *without* pretraining or tuning many weights; strong performance even under high missing rates. (See Figs. 1–3 in the paper.)

---

## Requirements
- **MATLAB** (R2020b or newer recommended)
- Standard MATLAB toolboxes (no special dependencies expected)

---
## Getting Started

### 1) Download Abilene TM data
This project uses the **Abilene** dataset.

Download the dataset from Matthew Roughan’s website and extract it **inside this folder**:
```bash
# Direct link
curl -LO https://roughan.info/data/Abilene.tar.gz

# Extract
tar -xvf Abilene.tar.gz
```
For background and documentation on the dataset, see: https://roughan.info/project/traffic_matrix/


### 2) Preprocess

Run the MATLAB preprocessing script to convert the raw Abilene data into a .mat file used by the experiments:

```
% In MATLAB
Preprocess_Abilene_data
```

This will create `abilene_tm_2016.mat` in the repository directory.

###  3) Run TM Completion

Execute the main script to perform traffic matrix completion with TNMF + Time Slicing:
```
% In MATLAB
TNNMFN_Abilene_ANNET
```
If everything is set correctly, you should obtain the same Figure 1 as reported in the paper (NMAE vs. missing probability).

## Citation
If you use this code, please cite:

```bibtex
@INPROCEEDINGS{Miyata24_TNNMFN,
  author    = {Miyata, Takamichi},
  booktitle = {IEEE Network Operations and Management Symposium Workshop},
  title     = {Completion of Traffic Matrix by Tensor Nuclear Norm Minus Frobenius Norm Minimization and Time Slicing},
  year      = {2024},
  pages     = {1-5},
  doi       = {10.1109/NOMS59830.2024.10575433}
}
```

---

## License
This project is released under the Apache 2.0 License.  
See the [LICENSE](./LICENSE) file for details.

---

## Contact
For questions, feedback, or collaboration inquiries:  
📧 **takamichi.miyata@it-chiba.ac.jp**
