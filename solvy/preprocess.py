"""
Preprocesses the dataset by checking the checksum, clearing data and splitting it into train, val and test files
"""
import sys
import os
import json
import pandas as pd
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import rdkit
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
from collections import defaultdict
from solvy.utils import load_config


def scaffold_splits(
    smiles_list: list[str], 
    frac_train: float = 0.8, 
    frac_val: float = 0.1, 
    frac_test: float = 0.1
) -> tuple[list[int], list[int], list[int]]:
    """
    Uses scaffolds to split the data
    
    Parameters
    ----------
    smiles_list: list of strings
        List of SMILES
    frac_train: float, default = 0.8
        Fraction of train dataset
    frac_val: float, default = 0.1
        Fraction of val dataset
    frac_test: float, default = 0.1
        Fraction of test dataset

    Returns
    -------
    tuple of indices (list[int)
        Returns three lists of indices, i.e. results of splits
    """
    assert frac_train + frac_test + frac_val == 1, "Train / Val / Test fractions do not sum up to 1!"

    scaffold_indices = defaultdict(list)
    for idx, smiles in enumerate(smiles_list):
        scaffold_indices[MurckoScaffoldSmiles(smiles=smiles)].append(idx)

    train_idx, val_idx, test_idx = [], [], []
    for indices in sorted(scaffold_indices.values(), key=len, reverse=True):
        if len(train_idx) + len(indices) <= frac_train * len(smiles_list):
            train_idx += indices
        elif len(val_idx) + len(indices) <= frac_val * len(smiles_list):
            val_idx += indices
        else:
            test_idx += indices

    return train_idx, val_idx, test_idx


def preprocess(
    name: str,
    source_file: str,
    source_sha256: str,
    smiles_col: str = "SMILES",
    split_method: str = "scaffold",
    frac_train: float = 0.8, 
    frac_val: float = 0.1, 
    frac_test: float = 0.1
) -> None:
    """
    Runs the whole preprocess
    1. Verify checksums
    2. Clear data
    3. Split into train, val, test

    Parameters
    ----------
    source_file: str,
        A filepath to the source datatset
    name: str,
        Dataset name
    smiles_col: str = "SMILES",
        Column name of SMILES inside the dataset
    split_method: str = "scaffold",
        Split method to use
    frac_train: float, default = 0.8
        Fraction of train dataset
    frac_val: float, default = 0.1
        Fraction of val dataset
    frac_test: float, default = 0.1
        Fraction of test dataset

    Returns
    -------
    None
    """
    df = pd.read_csv(fr"data/raw/esol.csv")
    available_split_methods = {"scaffold"}
    assert split_method in available_split_methods, (
        f"Split method {split_method} doesn't exist!"
        f" Please choose one of the following: {available_split_methods}"
    )

    # Identify split methods
    if split_method == "scaffold":
        if not smiles_col:
            raise ValueError("Smiles column is not provided for scaffold method!")

        smiles_list = df[smiles_col]
        train_idx, val_idx, test_idx = scaffold_splits(
            smiles_list=smiles_list, frac_train=frac_train, frac_test=frac_test, frac_val=frac_val
        )

    file_name_fmt = "{name}_{split_type}.csv"
    folder_name = Path("data/processed") / split_method
    folder_name.mkdir(parents=True, exist_ok=True)

    # Create manifest
    manifest_params = {
        "dataset_name": name,
        "source_file": str(source_file),
        "source_sha256": source_sha256,
        "split_method": split_method,
        "split_fractions": {
            "train": frac_train,
            "val": frac_val,
            "test": frac_test,
        },
        "n_train": len(train_idx),
        "n_val": len(val_idx),
        "n_test": len(test_idx),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "rdkit_version": rdkit.__version__,
    }

    # Write manifest
    with open(folder_name / "manifest.json", "w") as f:
        json.dump(manifest_params, f, indent=4)

    # Write splits
    splits = [
        ("train", train_idx),
        ("val", val_idx),
        ("test", test_idx),
    ]
    for split_type, split_indices in splits:
        split_filepath = folder_name / file_name_fmt.format(name=name, split_type=split_type)
        df.iloc[split_indices].to_csv(split_filepath, compression="gzip")


if __name__ == "__main__":
    cfg = load_config("datasets")
    dataset_name = "esol"
    dataset_cfg = cfg[dataset_name]

    data_dir = Path("data/raw")
    source_file = data_dir / dataset_cfg["filename"]
    preprocess(
        name=dataset_name,
        source_file=source_file,
        source_sha256=dataset_cfg["sha256"],
        smiles_col="SMILES",
        split_method="scaffold",
        frac_train=0.8, 
        frac_val=0.1, 
        frac_test=0.1,
    )
