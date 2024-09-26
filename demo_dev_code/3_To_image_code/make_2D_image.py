import os
import time
import requests
import pandas as pd
from urllib.parse import quote
from requests.exceptions import HTTPError, ConnectionError, Timeout, RequestException
from rdkit import Chem
from rdkit.Chem import Draw
import cairosvg

# 폴더가 없으면 생성하는 함수
def create_directory_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"폴더 {directory}가 생성되었습니다.")
    else:
        print(f"폴더 {directory}는 이미 존재합니다.")

# RDKit SMILES -> 이미지 변환
def rdkitSmilesToImage(smiles, filename, title, sanitize=True):
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=sanitize)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        img = Draw.MolToImage(mol, size=(300, 300))
        img.save(filename)
        return True
    except Exception as e:
        print(f"rdkit SMILES 변환에 실패했습니다: {smiles}, 에러: {str(e)}")
        return False

# PubChem SMILES -> 이미지 변환
def pubchemSmilesToPng(smiles, filename, retries=3, delay=2):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    if not smiles or len(smiles.strip()) == 0:
        print("SMILES 문자열이 비어 있습니다. 변환을 건너뜁니다.")
        return False

    encoded_smiles = quote(smiles)
    cid_url = f"{base_url}/smiles/{encoded_smiles}/cids/JSON"
    
    for attempt in range(retries):
        try:
            cid_response = requests.get(cid_url)
            cid_response.raise_for_status()
            
            cid_data = cid_response.json()
            if "IdentifierList" in cid_data and "CID" in cid_data["IdentifierList"]:
                cid = cid_data["IdentifierList"]["CID"][0]
                png_url = f"{base_url}/cid/{cid}/PNG"
                png_response = requests.get(png_url)
                
                if png_response.status_code == 200:
                    with open(filename, 'wb') as file:
                        file.write(png_response.content)
                    return True
                else:
                    print(f"pubchem PNG 가져오는 중 오류 발생: {png_response.status_code}")
                    return False
            else:
                print("pubchem CID를 찾을 수 없습니다.")
                return False
        except (HTTPError, ConnectionError, Timeout, RequestException) as err:
            print(f"pubchem 오류 발생: {err}")
            if cid_response.status_code in [429, 503]:
                time.sleep(delay)
                delay *= 2
            else:
                return False
    return False

# ChEMBL SMILES -> SVG -> PNG 변환
def chemblSmilesToPng(smiles, svg_filename, png_filename):
    os.system(f"curl -X POST -d \"{smiles}\" https://www.ebi.ac.uk/chembl/api/utils/smiles2svg > {svg_filename}")
    cairosvg.svg2png(url=svg_filename, write_to=png_filename)
    if os.path.exists(svg_filename):
        os.remove(svg_filename)

# CSV 데이터를 읽고 여러 형식의 이미지를 생성
def convertSmilesToImages(csv_file, output_dir):
    df = pd.read_csv(csv_file)
    
    # 실패한 SMILES 저장을 위한 딕셔너리
    failed_smiles_dict = {
        'pubchem': pd.DataFrame(columns=df.columns),
        'rdkit': pd.DataFrame(columns=df.columns),
        'chembl': pd.DataFrame(columns=df.columns)
    }
    
    # 출력 경로 설정
    pubchem_out_dir = f"{output_dir}/pubchem"
    rdkit_out_dir = f"{output_dir}/rdkit"
    chembl_out_dir = f"{output_dir}/chembl"

    # 폴더 생성
    for dir_path in [pubchem_out_dir, rdkit_out_dir, chembl_out_dir]:
        create_directory_if_not_exists(dir_path)
    
    # 파일과 SMILES 컬럼 찾기
    filename_column = next((col for col in ['Name', 'Title', 'origin_name', 'cid'] if col in df.columns), None)
    smiles_column = next((col for col in ['Canonical_SMILES', 'Canonical SMILES'] if col in df.columns), None)
    
    if not filename_column or not smiles_column:
        print("필요한 컬럼을 찾을 수 없습니다.")
        return
    
    # SMILES -> 이미지 변환
    for i in range(len(df)):
        title = df.loc[i, filename_column]
        smiles = df.loc[i, smiles_column]
        
        # 파일 경로 설정
        pubchem_png_path = f"{pubchem_out_dir}/{title}.png"
        rdkit_png_path = f"{rdkit_out_dir}/{title}.png"
        
        chembl_svg_path = f"{chembl_out_dir}/{title}.svg"
        chembl_png_path = f"{chembl_out_dir}/{title}.png"
        
        # PubChem 이미지 변환
        if not pubchemSmilesToPng(smiles, pubchem_png_path):
            failed_smiles_dict['pubchem'] = pd.concat([failed_smiles_dict['pubchem'], df.iloc[[i]]])
        
        # RDKit 이미지 변환
        if not rdkitSmilesToImage(smiles, rdkit_png_path, title):
            failed_smiles_dict['rdkit'] = pd.concat([failed_smiles_dict['rdkit'], df.iloc[[i]]])
        
        # ChEMBL 이미지 변환 (SVG -> PNG)
        try:
            chemblSmilesToPng(smiles, chembl_svg_path, chembl_png_path)
        except Exception as e:
            print(f"ChEMBL 이미지 변환 실패: {e}")
            failed_smiles_dict['chembl'] = pd.concat([failed_smiles_dict['chembl'], df.iloc[[i]]])
        
        time.sleep(1)  # API 제한 방지
    
    # 각 폴더별로 실패한 SMILES를 저장
    for key, failed_smiles_df in failed_smiles_dict.items():
        if not failed_smiles_df.empty:
            failed_csv_filename = f"{output_dir}/{key}_failed_smiles.csv"
            failed_smiles_df.to_csv(failed_csv_filename, index=False)
            print(f"{key} 이미지 생성이 실패한 데이터가 {failed_csv_filename}로 저장되었습니다.")
        else:
            print(f"{key}에 실패한 SMILES가 없습니다. 모든 이미지가 성공적으로 생성되었습니다.")
    
    return [pubchem_out_dir, rdkit_out_dir, chembl_out_dir]

def countFilesInDirectory(directory_path):
    try:
        # 경로가 존재하는지 확인
        if not os.path.exists(directory_path):
            print(f"경로 {directory_path} 가 존재하지 않습니다.")
            return 0
        
        # 디렉토리 내 파일 갯수 세기
        file_count = len([f for f in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, f))])
        
        return print(f"디렉토리 '{directory_path}'에 있는 파일 개수: {file_count}")
    
    except Exception as e:
        print(f"오류 발생: {e}")
        return 0

         
if __name__ == '__main__':
    # csv_list = ["/group/CBANK_01/0_data/240913/2_standardization_complete_data/molvs/240919_test_data_repair_molvs_sdtOut.csv", "/group/CBANK_01/0_data/240913/2_standardization_complete_data/rdkit/240919_test_data_repair_rdkit_sdtOut.csv"]
    # output_list = ["/group/CBANK_01/3_To_image_code_result/results/240919_origin_240913_testdata/molvs_std", "/group/CBANK_01/3_To_image_code_result/results/240919_origin_240913_testdata/rdkit_std"]
    
    # for i in range(len(csv_list)): 
    #     print("이 실행파일은 PubChem(NCBI) API, RDKit, ChEMBL API을 활용하여 이미지를 생성하는 코드 입니다.")
    #     csv_file = csv_list[i]
        
    #     print("저장할 경로를 입력하면 해당 폴더에 각각의 폴더(3개)가 생성됩니다. 그 안에 이미지가 저장됩니다. ")
    #     output_dir = output_list[i]
        
    #     file_path = convertSmilesToImages(csv_file, output_dir)
    #     countFilesInDirectory(file_path[0])
    #     countFilesInDirectory(file_path[1])
    #     countFilesInDirectory(file_path[2])
    
    
    print("이 실행파일은 PubChem(NCBI) API, RDKit, ChEMBL API을 활용하여 이미지를 생성하는 코드 입니다.")
    csv_file = input("csv 경로를 입력해주세요 : ")
    print("저장할 경로를 입력하면 해당 폴더에 각각의 폴더(3개)가 생성됩니다. 그 안에 이미지가 저장됩니다. ")
    output_dir = input("저장될 경로를 입력해주세요 : ")
    convertSmilesToImages(csv_file, output_dir)
        