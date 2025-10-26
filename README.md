# Analysis Pipelines for *Evolutionarily Conserved Transcriptional Regulators Control Monoaminergic Neuron Development*  
**Lewis et al., 2025**

**Contact:** Clifton Lewis — [lewis.clifton1597@gmail.com](mailto:lewis.clifton1597@gmail.com)

---

This repository contains all analysis and image-processing scripts developed for the study  
*“Evolutionarily Conserved Transcriptional Regulators Control Monoaminergic Neuron Development”* (Lewis et al., 2025).  

It accompanies the processed datasets hosted on Zenodo:  
- **scRNA-seq datasets** – [Zenodo DOI](https://doi.org/XXXXXXX)  [tbc]
- **Confocal imaging data** – [Zenodo DOI](https://doi.org/XXXXXXX)   [tbc]
- **FACS datasets** – [Zenodo DOI](https://doi.org/XXXXXXX)   [tbc]

---

## Repository Overview

### Confocal Image Processing Scripts  
FIJI (ImageJ) macros and workflows for confocal embryo imaging, alignment, and orientation correction.

### HCR Image Processing Scripts  
Pipelines for hybridisation chain reaction (HCR) image processing, including channel separation, and figure generation.

### scRNA-seq Analysis  
R and Python scripts for:  
- Pre-processing and integration of time-resolved *Drosophila* scRNA-seq data (Seurat ≥ v5.2)  
- RNA velocity and trajectory inference (scVelo, CellRank v2)  
- Gene ontology and transcription factor dynamics analyses  
- Orthogroup mapping across bilaterian species (eggNOG v5.0)

### Targeted Metabolomics  
R scripts for processing and visualising targeted LC–MS/MS metabolomics data, corresponding to the analyses presented in **Figure 5** of the manuscript.

---

## Environment and Dependencies

- **R ≥ 4.3.0**  
- **Python ≥ 3.10**  
- **Seurat ≥ 5.2**  
- **scVelo ≥ 0.3**  
- **CellRank ≥ 2.0**  
- **FIJI ≥ 2.9** (for ImageJ macros)  

---

## License

All scripts and code are distributed under the **MIT License**.  
Please cite this repository and the associated publication when reusing or adapting any materials.

---

## Citation

**Lewis, C.**, Goulty, M., Wroblewska, A., Croxall, N., Onion, D., Robinson, S. W., Zinzen, R. P., Solana, J., Kyriacou, C. P., Rosato, E., & Feuda, R. (2025).  
*Evolutionarily Conserved Transcriptional Regulators Control Monoaminergic Neuron Development.*  
University of Leicester

---
#
