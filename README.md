# LTS-BiodiversityAnalysisCAO
Code repository for the support project P_Snoeijs-Leijonmalm_2205, specifically 
related to identifying and mapping COI reads.

## Installation

1. Clone the repository
```bash
git clone git@github.com:NBISweden/LTS-BiodiversityAnalysisCAO.git
```

2. Install the software environment
```
conda env create -f environment.yml
```

3. Activate the environment
```
conda activate biodivcao
```

## Testrun

```bash
snakemake --profile test -c 1
```