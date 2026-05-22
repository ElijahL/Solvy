# Description

Source code for the project. All utility functions as well as backend and frontend servers are defined here

## Data Processing
The entry point for data preprocessing is `/preprocess.py` which clears the data (standardize any SMILES that are present here) and splits the given dataset into train / test / val according to the given split method. 

It is a cli application, for more detailed overview please run `uv run -m solvy.preprocess -h` from the root of the project

### Data preparation
We're trying to standardize the SMILES that we get from the dataset. 

Standardization is implemented inside `solvy.standardize`

Currently the following preprocessing of molecules is done:
1. **Cleanup**. Fixes common issues (rdMolStandardize.Cleanup)
2. **Fragment parent**. Strips any redundant parts (i.e. salts, e.g. CC(=O)O.[Na+] -> CC(=O)O)
3. **Uncharge**. Neutralizes the molecule (adds or remove H^+ to remove any ionization)
4. **Canonical tautomer**. Picks canonical tautomer

The order is quite important as almost every step potentially does some modification to the structure of a molecule. For instance, if uncharge is done before parent fragmentation, we may end up with charged molecule.

Take `CC(=O)[O-].[Na+]` (Sodium Acetate) for example. 
- Uncharge -> fragment parent
1. uncharge: `CC(=O)[O-].[Na+]` -> `CC(=O)[O-].[Na+]` - Nothing happens as the molecule is neutral
2. fragment_parent: `CC(=O)[O-].[Na+]` -> `CC(=O)[O-]` - The base molecule was found (which is charged!)

- Fragment parent -> uncharge
1. fragment_parent: `CC(=O)[O-].[Na+]` -> `CC(=O)[O-]` - The base molecule was found (which is charged!)
2. uncharge: `CC(=O)[O-]` -> `CC(=O)O` - Charged molecule was neutralized

As a result we end up with 2 different structures after performing steps in different order. A code for experimentation can be found in `notebooks/02_scaffolds.ipynb`

For each molecule a `StandardizationResult` dataclass is created, which stores every steps which caused the molecule structure to change.

### Data splits
When a `preprocess` is run - a folder with the name of the split method inside `data/processed` directory is created and each split (i.e. train / test / val ) is saved inside that directory.

Split methods available:
- scaffold (MurckoScaffoldSmiles)

For each group of splits a `manifest.json` file is created, describing how the split was performed: number of entries per split, split method used, source file hash sum, etc.

An example `manifest.json` file:
```json
{
    "dataset_name": "esol",
    "source_file": "data/raw/esol.csv",
    "source_sha256": "768a05b6ae300d5dc0750dce68e97b697b43a0d69de2bbb83d5915b3e09d0064",
    "split_method": "scaffold",
    "split_fractions": {
        "train": 0.8,
        "val": 0.1,
        "test": 0.1
    },
    "n_train": 902,
    "n_val": 112,
    "n_test": 114,
    "created_at": "2026-05-21T14:02:01.030708+00:00",
    "rdkit_version": "2026.03.1"
}
```
