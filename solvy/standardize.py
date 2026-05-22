from typing import Optional, TypeAlias, Callable
from dataclasses import dataclass, field
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


class StandardizationException(Exception): ...

StandardizationStep: TypeAlias = dict[str, str, str] # keys: step, input_smiles, output_smiles

@dataclass 
class StandardizationResult:
    input_smiles: str
    output_smiles: str
    steps: list[StandardizationStep] = field(default_factory=list)
    is_dropped: bool = False
    reason: Optional[str] = None


def run_standardization_pipeline(mol: Chem.Mol, pipeline: list[tuple[str, Callable]]) -> list[StandardizationStep]:
    """
    Return standardization pipeline, raises StandardizationException if unprocesssable

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        rdkit molecule object. *Mutable*
    pipeline: list of tuples of step name and executable
        Defines the steps to standardize a molecule

    Return
    ------
    list[dict]
        A list of steps under which the molecule was modified, and what exactly was modified
    """
    steps = []
    for name, fn in pipeline:
        before = Chem.MolToSmiles(mol)
        try:
            mol = fn(mol)
        except Exception as e:
            raise StandardizationException(f"{name}_failed: {e}")

        after = Chem.MolToSmiles(mol)
        if before != after:
            steps.append({"step": name, "before": before, "after": after})

    return steps


def standardize(smiles: str) -> StandardizationResult:
    """
    Standardizes smiles by doing the following:
    - cleanup
        Fix common issues (rdMolStandardize.Cleanup)
    - fragment parent
        Strips any redundant parts (i.e. salts, e.g. CC(=O)O.[Na+] -> CC(=O)O)
    - uncharge
        Neutralizes the molecule (adds or remove H^+ to remove any ionization)
    - canonical tautomer
        Pick canonical tautomer

    Parameters
    ----------
    smiles: str
        Input smiles string

    Return
    ------
    StandardizationResult
        Result class that stores input_smiles, output_smiles 
        as well as all the standardization methods and steps
    """
    res = StandardizationResult(input_smiles=smiles, output_smiles=None)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        res.dropped = True
        res.reason = "unparseable_smiles"
        return res

    # Define standardization pipeline
    pipeline = [
        ("cleanup",             rdMolStandardize.Cleanup),
        ("fragment_parent",     rdMolStandardize.FragmentParent),
        ("uncharge",            rdMolStandardize.Uncharger().uncharge),
        ("canonical_tautomer",  rdMolStandardize.TautomerEnumerator().Canonicalize),
    ]

    # Standardize molecules recording any failes
    try:
        res.steps += run_standardization_pipeline(mol, pipeline)
    except StandardizationException as e:
        res.dropped = True
        res.reason = f"{name}_failed: {e}"
        return res

    res.output_smiles = Chem.MolToSmiles(mol)
    return res


