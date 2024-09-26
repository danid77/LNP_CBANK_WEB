import argparse
from rdkit import Chem
from molvs import Standardizer as MolvsStandardizer, validate_smiles as molvsValidateSmiles
from chembl_structure_pipeline import standardizer as chemblStandardizer, checker

# RDKit만 사용하는 SMILES 처리 함수
def stnSmilesRdkit(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Error: Invalid SMILES string."
        stnSmiles = Chem.MolToSmiles(mol, canonical=False)
        print(f"standardization Canonical SMILES (RDKit): {stnSmiles}")
        
        return stnSmiles
    except Exception as e:
        return f"Error: {e}"

# molvs를 사용하는 SMILES 처리 함수
def stnSmilesMolvs(smiles):
    if not smiles or smiles.strip() == '':
        return 'Invalid SMILES', 'Invalid SMILES'
    
    validationMessages = molvsValidateSmiles(smiles)
    
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        
        standardizer = MolvsStandardizer()
        standardizedMol = standardizer.standardize(mol)
        
        stnSmiles = Chem.MolToSmiles(standardizedMol, canonical=False)
        print(f"standardization Canonical SMILES (Molvs): {stnSmiles}, validationMessages : {validationMessages}")
        
        return stnSmiles
    except Exception as e:
        return validationMessages, str(e)

# chembl_structure_pipeline을 사용하는 SMILES 처리 함수
def stnSmilesChembl(smiles):
    try:
        # SMILES에서 Mol 객체 생성
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            raise ValueError(f"SMILES 변환 실패: {smiles}")
        
        # Mol 객체를 MolBlock으로 변환
        molBlock = Chem.MolToMolBlock(mol)
        
        # 표준화 및 부모 구조 생성
        stdMolBlock = chemblStandardizer.standardize_molblock(molBlock)
        parentMolBlock, _ = chemblStandardizer.get_parent_molblock(stdMolBlock)
        
        # 문제 확인 (경고 발생 시 무시하고 계속 진행)
        issues = checker.check_molblock(parentMolBlock)
        if issues:
            print(f"Issues found during checking for : {issues}")
        
        # 부모 MolBlock을 다시 Mol 객체로 변환
        parentMol = Chem.MolFromMolBlock(parentMolBlock, sanitize=False)
        if parentMol is None:
            raise ValueError(f"MOL 블록 변환 실패: {parentMolBlock}")

        # Canonical SMILES 생성 (False로 지정하여 canonical 처리하지 않음)
        stnSmiles = Chem.MolToSmiles(parentMol, canonical=False)
        print(f"Standardization Canonical SMILES (chembl_structure_pipeline): {stnSmiles}")
        
        return stnSmiles
    
    # 모든 예외 처리 (오류가 발생해도 결과 반환)
    except Exception as e:
        return f"Error occurred: {e}"

# 메인 함수
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a SMILES string using RDKit, molvs, or chembl_structure_pipeline.")
    
    parser.add_argument("--smiles", type=str, required=True, help="The SMILES string to process")
    parser.add_argument(
        "--mode", type=str, choices=["rdkit", "molvs", "chembl"], default="rdkit",
        help="Choose which library to use: 'rdkit' for RDKit only, 'molvs' for MolVS validation and standardization, 'chembl' for chembl_structure_pipeline"
    )

    # 명령줄 인자 파싱
    args = parser.parse_args()

    # SMILES 처리
    if args.mode == "rdkit":
        result = stnSmilesRdkit(args.smiles)
    elif args.mode == "molvs":
        result = stnSmilesMolvs(args.smiles)
    elif args.mode == "chembl":
        result = stnSmilesChembl(args.smiles)

    # 최종 결과 출력
    print(result)
