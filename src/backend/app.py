from fastapi import FastAPI

from src import predict


app = FastAPI()


@app.get("/")
async def root():
    return {"message": "Hello, World!"}


@app.get("/healthcheck")
async def healthcheck():
    return True


@app.get("/predict_solubility")
async def predict_solubility(smiles: str) -> float:
    return predict.predict_solubility(smiles=smiles)

