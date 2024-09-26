# backend.py
from fastapi import FastAPI, File, UploadFile, Form, HTTPException
from pydantic import BaseModel
from io import BytesIO
from fastapi.responses import StreamingResponse
import base64
from chem_converter import (
    # getStereoChoice,
    getCanonicalIsomerSmiles,
    CanonicalToSmiles,
    InChIKeyToSmiles,
    InchiToSmiles,
    sdfToSmiles,
    molToSmiles,
    iupacToSmiles,
    smartsToSmiles,
    peptideToSmiles,
    singleHandleChoice
)
from standardization import stnSmilesRdkit, stnSmilesMolvs, stnSmilesChembl
from make_iamge import rdkitSmilesToImage, pubchemSmilesToPng, chemblSmilesToPng


app = FastAPI()

# 요청에 사용할 데이터 모델 정의
class ChemicalConversionRequest(BaseModel):
    input_type: str
    input_value: str = None
    
class ChemicalStandardizationRequest(BaseModel):
    convert_smiles : str
    select_stn: str
    
class ChemicalCanonicalIsomerRequest(BaseModel):
    convert_smiles : str
    select_isomer: bool

class ChemicalImageRequest(BaseModel):
    convert_smiles : str

@app.post("/convert/")
def convertChemical(data: ChemicalConversionRequest):
    try:
        if data.input_type == "smiles":
            convert = CanonicalToSmiles(data.input_value)
        elif data.input_type == "inchi_key":
            convert = InChIKeyToSmiles(data.input_value)
        elif data.input_type == "inchi":
            convert = InchiToSmiles(data.input_value)
        elif data.input_type == "iupac":
            convert = iupacToSmiles(data.input_value)
        elif data.input_type == "peptide":
            convert = peptideToSmiles(data.input_value)
        elif data.input_type == "smarts":
            convert = smartsToSmiles(data.input_value)
        else:
            raise HTTPException(status_code=400, detail="Invalid input type")
        return {"convert": convert}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/convert_file/")
async def convertChemicalFile(
    input_type: str = Form(...),  # 폼 데이터로 입력 받음
    file: UploadFile = File(None)  # SDF 또는 MOL2 파일 업로드
):
    try:
        # SDF 또는 MOL2 파일 처리
        if input_type in ["sdf", "mol2"] and file:
            file_content = await file.read()
            print(file_content[:1000])

            if input_type == "sdf":
                try:
                    convert_file = sdfToSmiles(file_content)
                except ValueError as e:
                    raise HTTPException(status_code=400, detail=str(e))
            elif input_type == "mol2":
                try:
                    convert_file = molToSmiles(file_content)
                except ValueError as e:
                    raise HTTPException(status_code=400, detail=str(e))
            else:
                raise HTTPException(status_code=400, detail="Invalid input type for file")

        return {"convert": convert_file}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    
@app.post("/stn/")
def standardization(data: ChemicalStandardizationRequest):
    try:
        if data.select_stn == "rdkit":
            stn = stnSmilesRdkit(data.convert_smiles)
        elif data.select_stn == "molvs":
            stn = stnSmilesMolvs(data.convert_smiles)
        elif data.select_stn == "chembl_pipeline":
            stn = stnSmilesChembl(data.convert_smiles)
        else:
            raise HTTPException(status_code=400, detail="Invalid input type")
        return {"stn": stn}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/isomer/")
def canonicalIsomerSmiles(date: ChemicalCanonicalIsomerRequest):
    print(date.select_isomer)
    try:
        if date.select_isomer is not None  :
            canonical_isomer_smiles = getCanonicalIsomerSmiles(date.convert_smiles, date.select_isomer)
        else:
            raise HTTPException(status_code=400, detail="Invalid input type")
        return {"isomer": canonical_isomer_smiles}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    

@app.post("/images/")
def chemical_images(data: ChemicalImageRequest):
    smiles = data.convert_smiles
    if not smiles:
        raise HTTPException(status_code=400, detail="SMILES 문자열이 없습니다.")

    try:
        # RDKit 이미지 생성
        rdkit_image = rdkitSmilesToImage(smiles)
        if isinstance(rdkit_image, BytesIO):
            rdkit_image_base64 = base64.b64encode(rdkit_image.getvalue()).decode('utf-8')
        else:
            rdkit_image_base64 = None

        # PubChem 이미지 생성
        pubchem_image = pubchemSmilesToPng(smiles)
        if isinstance(pubchem_image, BytesIO):
            pubchem_image_base64 = base64.b64encode(pubchem_image.getvalue()).decode('utf-8')
        else:
            pubchem_image_base64 = None

        # ChEMBL 이미지 생성
        chembl_image = chemblSmilesToPng(smiles)
        if isinstance(chembl_image, BytesIO):
            chembl_image_base64 = base64.b64encode(chembl_image.getvalue()).decode('utf-8')
        else:
            chembl_image_base64 = None

        # 결과 반환
        return {
            "rdkit": rdkit_image_base64,
            "pubchem": pubchem_image_base64,
            "chembl": chembl_image_base64,
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))