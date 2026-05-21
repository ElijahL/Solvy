Source code for the project. All utility functions as well as backend and frontend servers are defined here

# Description
## Data Processing
The entry point for data preprocessing is `/predict.py` which clears the data (standardize any SMILES that are present here) and splits the given dataset into train / test / val according to the given split method. 
It is a cli application, for more detailed overview please run `uv run -m solvy.preprocess -h` from the root of the project

