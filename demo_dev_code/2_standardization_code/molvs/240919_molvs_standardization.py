import os
import sys
import pandas as pd
from rdkit import Chem
from molvs import Standardizer, validate_smiles
import glob
import time  # time 모듈 추가

# SMILES 검증, 표준화 및 캐노니컬 SMILES 변환 함수
def process_smiles(smiles, standardizer):
    # NaN 또는 빈 값을 건너뛰기
    if pd.isna(smiles) or smiles.strip() == '':
        return 'Invalid SMILES', 'Invalid SMILES', False
    
    # SMILES 검증
    print(f"Validating SMILES: {smiles}")
    validation_messages = validate_smiles(smiles)
    print(f"Validation result for {smiles}: {validation_messages}")
    
    # SMILES 표준화 및 캐노니컬 변환
    try:
        print(f"Standardizing SMILES: {smiles}")
        # sanitize=False를 사용하여 Mol 객체 생성
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        
        # RDKit의 SanitizeMol로 분자를 표준화
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        
        # 표준화된 분자로 변환
        standardized_mol = standardizer.standardize(mol)
        
        # 캐노니컬 SMILES로 변환
        canonical_smiles = Chem.MolToSmiles(standardized_mol, canonical=True)
        print(f"Canonical SMILES for {smiles}: {canonical_smiles}")
        
        return validation_messages, canonical_smiles, True  # 성공 시 True 반환
    except Exception as e:
        canonical_smiles = str(e)
        print(f"Error processing SMILES: {smiles}, Error: {e}")
        return validation_messages, canonical_smiles, False  # 실패 시 False 반환

# 디렉토리에서 CSV 파일을 읽고 처리하는 함수
def process_csv_files(directory):
    # Standardizer 생성
    standardizer = Standardizer()
    
    # 디렉토리 내 모든 CSV 파일 찾기
    csv_files = glob.glob(os.path.join(directory, "*.csv"))
    
    for file_path in csv_files:
        print(f"Processing file: {file_path}")
        
        # 시작 시간 기록
        start_time = time.time()
        
        # 성공 및 실패 카운트 변수 초기화
        success_count = 0
        failure_count = 0
        
        # CSV 파일 읽기 (Title, origin_name, SMILES 열 포함)
        df = pd.read_csv(file_path, delimiter=',')
        
        # 검증, 표준화 및 캐노니컬 SMILES 변환 적용
        results = df['SMILES'].apply(lambda x: process_smiles(x, standardizer))
        df['Validation'] = results.apply(lambda x: x[0])
        df['Canonical_SMILES'] = results.apply(lambda x: x[1])
        
        # 성공 및 실패 카운트 갱신
        success_count = results.apply(lambda x: x[2]).sum()  # True인 경우만 카운트
        failure_count = len(results) - success_count  # 전체에서 성공을 빼면 실패 수
        
        # 출력 파일 경로 설정 (원본 파일에 '_out' 추가)
        output_file_path = file_path.replace('.csv', '_out.csv')
        
        # 변환된 결과를 원본 파일과 동일한 형식으로 저장
        df.to_csv(output_file_path, index=False)
        
        # 종료 시간 기록 및 경과 시간 출력
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Processed file saved as: {output_file_path}, Time taken: {elapsed_time:.2f} seconds")
        print(f"Success: {success_count}, Failure: {failure_count} for SMILES in file {file_path}")

# 메인 함수
if __name__ == "__main__":
    if len(sys.argv) == 1:
        # 명령줄 인자가 없으면 현재 디렉토리에서 작업
        current_directory = os.getcwd()
        print(f"Processing files in current directory: {current_directory}")
        process_csv_files(current_directory)
    elif len(sys.argv) == 2:
        # 디렉토리 경로가 주어진 경우 해당 디렉토리에서 작업
        directory = sys.argv[1]
        if os.path.isdir(directory):
            print(f"Processing files in directory: {directory}")
            process_csv_files(directory)
        else:
            print(f"Error: {directory} is not a valid directory.")
    else:
        print("Usage:")
        print("  python process_smiles.py           # Process CSV files in current directory")
        print("  python process_smiles.py /path/to/directory  # Process CSV files in the specified directory")

