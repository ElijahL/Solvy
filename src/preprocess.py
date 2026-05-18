import sys
import os
import pandas as pd
from pathlib import Path

from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
from collections import defaultdict

def scaffold_splits(
    smiles_list: list[str], 
    frac_train: float = 0.8, 
    frac_val: float = 0.1, 
    frac_test: float = 0.1
) -> tuple[list[int], list[int], list[int]]:
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
    df: pd.DataFrame, 
    name: str,
    smiles_col: str = "SMILES",
    split_method: str = "scaffold"
):
    smiles_list = df[smiles_col]
    available_split_methods = {"scaffold"}
    assert split_method in available_split_methods, (
        f"Split method {split_method} doesn't exist!"
        f" Please choose one of the following: {available_split_methods}"
    )

    if split_method == "scaffold":
        train_idx, val_idx, test_idx = scaffold_splits(smiles_list=smiles_list)

    file_name_fmt = "{name}_{split_type}.csv"
    folder_name = Path("data/processed") / split_method
    if not folder_name.exists():
        os.makedirs(str(folder_name))

    df.iloc[train_idx].to_csv(folder_name / file_name_fmt.format(name=name, split_type="train"), compression="gzip")
    df.iloc[val_idx].to_csv(folder_name / file_name_fmt.format(name=name, split_type="val"), compression="gzip")
    df.iloc[test_idx].to_csv(folder_name / file_name_fmt.format(name=name, split_type="test"), compression="gzip")


if __name__ == "__main__":
    df = pd.read_csv(fr"data/raw/esol.csv")
    preprocess(df, "esol")