import yaml
from pathlib import Path


def load_config(config_name: str) -> dict:
    """
    Loads config file in yaml format which is located in configs/

    Parameters
    ----------
    config_name: str
        Name of the config, i.e. config file name (w/o extension)

    Return
    ------
    dict
        Loaded config in dict format
    """
    filename = f"{config_name}.yaml"
    cfg_path = Path("configs") / filename
    if not cfg_path.exists():
        raise ValueError(f"{str(cfg_path)} doesn't exist!")

    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    return cfg
