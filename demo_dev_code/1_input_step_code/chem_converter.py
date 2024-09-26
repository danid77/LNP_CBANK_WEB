import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import SDWriter
from chembl_webresource_client.new_client import new_client
from openbabel import openbabel
import requests
import pandas as pd
import os
import re

# ChEMBL 클라이언트 설정
molecule = new_client.molecule

def getStereoChoice(stereo_value):
    while True:
        stereo_choice = stereo_value
        if stereo_choice in ['y', 'n']:
            return True if stereo_choice == 'y' else False
        else:
            print("잘못된 입력입니다. y 또는 n을 입력해주세요.")

def CanonicalToSmiles(smiles, consider_stereochemistry):
    mol = Chem.MolFromSmiles(smiles)
    try:
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=consider_stereochemistry)
    except Exception as e:
        print(f"Warning: stereochemistry handling failed. Error: {str(e)}")
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    return canonical_smiles

def InChIKeyToSmiles(inchi_key, consider_stereochemistry):
    compound = pcp.get_compounds(inchi_key, 'inchikey')
    if len(compound) == 0:
        result = molecule.filter(molecule_structures__standard_inchi_key=inchi_key)
        if result:
            return result[0]['molecule_structures']['canonical_smiles']
        else:
            return "No compound found in PubChem or ChEMBL."
    inchi = compound[0].inchi
    mol = Chem.MolFromInchi(inchi)
    try:
        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=consider_stereochemistry)
    except Exception as e:
        print(f"Warning: stereochemistry handling failed. Error: {str(e)}")
        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    return smiles

def InchiToSmiles(inchi, consider_stereochemistry):
    mol = Chem.MolFromInchi(inchi)
    try:
        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=consider_stereochemistry)
    except Exception as e:
        print(f"Warning: stereochemistry handling failed. Error: {str(e)}")
        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    return smiles

def sdfToSmiles(sdf_path, consider_stereochemistry):
    supplier = Chem.SDMolSupplier(sdf_path)
    sdf_smiles = ''
    for mol in supplier:
        if mol is not None:
            try:
                sdf_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=consider_stereochemistry)
            except Exception as e:
                print(f"Warning: stereochemistry handling failed for SDF molecule. Error: {str(e)}")
                sdf_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
            print('sdf : ', sdf_smiles)
    return sdf_smiles

def sdfDirToSmiles(sdf_directory, output_csv_path, consider_stereochemistry=True):
    data = []
    
    for filename in os.listdir(sdf_directory):
        if filename.endswith('.sdf'):
            file_path = os.path.join(sdf_directory, filename)
            supplier = Chem.SDMolSupplier(file_path)
            for mol in supplier:
                if mol is not None:
                    try:
                        sdf_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=consider_stereochemistry)
                    except Exception as e:
                        print(f"Warning: stereochemistry handling failed for SDF molecule in {filename}. Error: {str(e)}")
                        sdf_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
                    data.append({'filename': filename, 'smiles': sdf_smiles})
    
    # DataFrame으로 변환 후 지정된 경로에 CSV 파일로 저장
    df = pd.DataFrame(data)
    df.to_csv(output_csv_path, index=False)
    
    # 절대 경로 반환
    absolute_path = os.path.abspath(output_csv_path)
    print(f"SDF SMILES CSV 파일이 저장되었습니다: {absolute_path}")
    
    return absolute_path


def molToSmiles(mol_path, consider_stereochemistry):
    mol = Chem.MolFromMolFile(mol_path)
    try:
        mol_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=consider_stereochemistry)
    except Exception as e:
        print(f"Warning: stereochemistry handling failed. Error: {str(e)}")
        mol_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    return mol_smiles

def molDirToSmiles(mol_directory, output_csv_path, consider_stereochemistry):
    data = []
    
    for filename in os.listdir(mol_directory):
        if filename.endswith('.mol'):
            file_path = os.path.join(mol_directory, filename)
            mol = Chem.MolFromMolFile(file_path)
            if mol is not None:
                try:
                    mol_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=consider_stereochemistry)
                except Exception as e:
                    print(f"Warning: stereochemistry handling failed for MOL file {filename}. Error: {str(e)}")
                    mol_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
                data.append({'filename': filename, 'smiles': mol_smiles})
    
    # DataFrame으로 변환 후 지정된 경로에 CSV 파일로 저장
    df = pd.DataFrame(data)
    df.to_csv(output_csv_path, index=False)
    
    # 절대 경로 반환
    absolute_path = os.path.abspath(output_csv_path)
    print(f"MOL SMILES CSV 파일이 저장되었습니다: {absolute_path}")
    
    return absolute_path


def iupacToSmiles(chemical_IUPAC, consider_stereochemistry):
    try:
        # Step 1: Try PubChem
        compounds = pcp.get_compounds(chemical_IUPAC, 'name')
        if len(compounds) == 0:
            # Step 2: Try ChEMBL
            result = molecule.filter(pref_name__iexact=chemical_IUPAC)
            if result:
                return result[0]['molecule_structures']['canonical_smiles']
            else:
                # Step 3: If both fail, use OPSIN API to parse the IUPAC name
                print("No compound found in PubChem or ChEMBL. Attempting to parse IUPAC name using OPSIN API...")
                try:
                    response = requests.get(f'https://opsin.ch.cam.ac.uk/opsin/{chemical_IUPAC}.json')
                    if response.status_code == 200:
                        data = response.json()
                        smiles = data.get('smiles', None)
                        if smiles:
                            mol = Chem.MolFromSmiles(smiles)
                            smiles_with_stereochemistry = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=consider_stereochemistry)
                            return smiles_with_stereochemistry
                        else:
                            return "OPSIN failed to parse the IUPAC name."
                    else:
                        return f"OPSIN API request failed with status code {response.status_code}."
                except Exception as e:
                    return f"Error during OPSIN IUPAC parsing: {str(e)}"
        else:
            inchi = compounds[0].inchi6
            mol = Chem.MolFromInchi(inchi)
            try:
                smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=consider_stereochemistry)
            except Exception as e:
                print(f"Warning: stereochemistry handling failed. Error: {str(e)}")
                smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
            return smiles
    except Exception as e:
        return f"An error occurred: {str(e)}"

def smartsToSmiles(smarts, consider_stereochemistry):
    molecule = Chem.MolFromSmarts(smarts)
    
    if molecule is None:
        return "Invalid SMARTS pattern"
    
    try:
        smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=consider_stereochemistry)
    except Exception as e:
        print(f"Warning: stereochemistry handling failed. Error: {str(e)}")
        smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
    
    return smiles

def peptideToSmiles(sequence):
    smiles_list = []
    amino_acids_smiles = {
        'A': 'CC(C(=O)O)N',   # Alanine
        'C': 'C(C(=O)O)N',    # Cysteine
        'D': 'CC(C(=O)O)C(=O)O', # Aspartic Acid
        'E': 'CCC(C(=O)O)C(=O)O', # Glutamic Acid
        'F': 'C1=CC=C(C=C1)CC(C(=O)O)N', # Phenylalanine
        'G': 'C(C(=O)O)N',   # Glycine
        'H': 'C1C=NC=N1CC(C(=O)O)N', # Histidine
        'I': 'CC(C)C(C(=O)O)N', # Isoleucine
        'K': 'CCCCN(C(=O)O)N', # Lysine
        'L': 'CC(C)CC(C(=O)O)N', # Leucine
        'M': 'CSCC(C(=O)O)N', # Methionine
        'N': 'CC(C(=O)O)C(=O)N', # Asparagine
        'P': 'C1CCNC1C(=O)O', # Proline
        'Q': 'CCC(C(=O)O)C(=O)N', # Glutamine
        'R': 'C(N=C(N)N)C(=O)O', # Arginine
        'S': 'C(C(C(=O)O)N', # Serine
        'T': 'C(C(C)C(=O)O)N', # Threonine
        'V': 'CC(C)C(C(=O)O)N', # Valine
        'W': 'C1=CC=C(C=C1)C2=CC=CC2C(C(=O)O)N', # Tryptophan
        'Y': 'C1=CC=C(C=C1)C(C(=O)O)O', # Tyrosine
    }
    for aa in sequence:
        if aa in amino_acids_smiles:
            smiles_list.append(amino_acids_smiles[aa])
        else:
            raise ValueError(f"Unknown amino acid code: {aa}")
    peptide_smiles = ".".join(smiles_list)
    return peptide_smiles


def smilesToSDF(smiles, output_filename='output.sdf'):
    try:
        # 분자를 SMILES로부터 생성
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("입력된 값이 유효한 SMILES 형식이 아닙니다.")
        
        # 파일을 output_filename 이름으로 생성
        writer = SDWriter(output_filename)
        writer.write(mol)
        writer.close()
        
        # 절대 경로를 얻기 위해 os.path.abspath 사용
        sdf_file_path = os.path.abspath(output_filename)
        
        # 생성된 SDF 파일의 경로를 반환
        return sdf_file_path
    except Exception as e:
        raise ValueError(f"파일 생성 중 오류가 발생했습니다: {str(e)}")
    
def smilesToMol2(smiles: str, output_file: str):
    """
    SMILES 문자열을 .mol2 파일로 변환하고, 저장된 파일의 절대 경로를 반환하는 함수.

    :param smiles: 변환할 SMILES 문자열
    :param output_file: 저장할 .mol2 파일 경로
    :return: 저장된 .mol2 파일의 절대 경로
    """
    # SMILES를 Mol 객체로 변환
    mol = Chem.MolFromSmiles(smiles)

    # Mol 객체를 MolBlock 형식으로 변환 (중간 단계)
    mol_block = Chem.MolToMolBlock(mol)

    # OpenBabel을 사용하여 MolBlock을 .mol2로 변환
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("mol", "mol2")

    ob_mol = openbabel.OBMol()
    ob_conversion.ReadString(ob_mol, mol_block)

    # 변환된 .mol2 파일 저장
    ob_conversion.WriteFile(ob_mol, output_file)

    # 저장된 파일의 절대 경로 반환
    absolute_path = os.path.abspath(output_file)
    print(f"{absolute_path} 파일로 저장되었습니다.")
    return absolute_path

def getSmilesInput() -> str:
    """
    SMILES 입력을 받아 유효성을 검사하고 반환하는 함수.
    
    :return: 유효한 SMILES 문자열
    """
    smiles_valid = False
    smiles = ""

    while not smiles_valid:
        try:
            smiles = input("SMILES 입력: ")
            
            # 유효하지 않은 SMILES 입력 시 예외 발생
            if not smiles:
                raise ValueError("SMILES 문자열을 입력해주세요.")
            
            # SMILES 유효성 검사
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("입력된 값이 유효한 SMILES 형식이 아닙니다.")
            
            smiles_valid = True  # 유효한 SMILES 입력 시 루프 종료

        except ValueError as ve:
            print(f"오류: {str(ve)}. 다시 시도해주세요.")
        except Exception as e:
            print(f"예상치 못한 오류가 발생했습니다: {str(e)}. 다시 시도해주세요.")
    
    return smiles
    
def getOutputFilename(default_filename: str) -> str:
    """
    파일 이름을 입력받고, 유효성을 검사하여 반환하는 함수.
    
    :param default_filename: 기본값으로 사용할 파일 이름
    :return: 유효한 파일 이름
    """
    path_valid = False
    output_filename = ""

    while not path_valid:
        try:
            output_filename = input(f"파일 이름 (기본값: {default_filename}): ") or default_filename
            
            # 파일 경로의 유효성 검사
            if not output_filename.strip():
                raise ValueError("유효한 파일 경로를 입력해주세요.")
            
            path_valid = True

        except ValueError as ve:
            print(f"오류: {str(ve)}. 다시 시도해주세요.")
        except Exception as e:
            print(f"예상치 못한 오류가 발생했습니다: {str(e)}. 다시 시도해주세요.")
    
    return output_filename  

def singleHandleChoice(choice: str):
    """
    주어진 선택에 따라 SMILES를 mol2 또는 sdf 파일로 변환하여 저장.
    
    :param choice: 사용자가 선택한 옵션 ("9" 또는 "10")
    """
    smiles = getSmilesInput()

    if choice == "9":
        # SDF 파일 처리
        output_filename = getOutputFilename("output.sdf")
        sdf_path = smilesToSDF(smiles, output_filename)
        print(f"생성된 SDF 파일 경로: {sdf_path}")

    elif choice == "10":
        # mol2 파일 처리
        output_filename = getOutputFilename("output.mol2")
        mol2_path = smilesToMol2(smiles, output_filename)
        print(f"생성된 mol2 파일 경로: {mol2_path}")

def sanitize_filename(filename):
    return re.sub(r'[\\/*?:"<>|]', "", filename)
       
def csvHandleChoice(choice):
    input_csv_path = input("CSV 파일 경로 입력: ")
    
    if choice == "9":
        output_dir = input("파일이 저장될 폴더 경로 입력: ") or "SMILES_to_sdf"
    elif choice == "10":
        output_dir = input("파일이 저장될 폴더 경로 입력: ") or "SMILES_to_mol2"

    # CSV 파일 읽기
    df = pd.read_csv(input_csv_path)

    # SMILES 컬럼 확인
    if 'SMILES' not in df.columns:
        raise ValueError("CSV 파일에 'SMILES' 컬럼이 존재하지 않습니다.")

    # 파일명에 사용할 컬럼 확인 (Name, Title, origin_name 중 첫 번째로 발견된 컬럼 사용)
    filename_column = None
    for possible_col in ['Name', 'Title', 'origin_name', 'cid']:
        if possible_col in df.columns:
            filename_column = possible_col
            break

    if filename_column is None:
        raise ValueError("CSV 파일에 파일명을 생성할 'Name', 'Title', 'origin_name' 컬럼이 존재하지 않습니다.")

    # 폴더가 없으면 생성
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 각 행에 대해 SMILES 값을 처리
    for idx, row in df.iterrows():
        smiles = row['SMILES']
        filename = sanitize_filename(str(row[filename_column]))

        if choice == "9":
            sdf_filename = os.path.join(output_dir, f"{filename}.sdf")
            sdf_path = smilesToSDF(smiles, sdf_filename)
            if sdf_path:
                print(f"SDF 파일 생성: {sdf_path}")

        elif choice == "10":
            mol2_filename = os.path.join(output_dir, f"{filename}.mol2")
            mol2_path = smilesToMol2(smiles, mol2_filename)
            if mol2_path:
                print(f"MOL2 파일 생성: {mol2_path}")      
            
