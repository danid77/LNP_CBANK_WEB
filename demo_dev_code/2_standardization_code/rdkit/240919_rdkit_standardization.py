import os
import sys
import pandas as pd
from rdkit import Chem
import glob
import time  # 시간 측정을 위한 모듈

# SMILES를 캐노니컬 SMILES로 표준화하는 함수
def rdkit_standardization(smiles):
    try:
        # SMILES를 RDKit 분자로 변환 후 다시 캐노니컬 SMILES로 변환
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "X", False  # 변환에 실패하면 'X'로 표기하고 False 반환
        return Chem.MolToSmiles(mol, canonical=True), True  # 성공하면 True 반환
    except:
        return "X", False  # 오류가 발생하면 'X'로 표기하고 False 반환

# 디렉토리에서 CSV 파일을 읽고 처리하는 함수
def process_csv_files(directory):
    # 디렉토리 내 모든 CSV 파일 찾기
    csv_files = glob.glob(os.path.join(directory, "*.csv"))
    
    for file_path in csv_files:
        print(f"Processing file: {file_path}")
        
        # 시작 시간 기록
        start_time = time.time()
        
        # CSV 파일 읽기 (Title, origin_name, SMILES 열 포함)
        df = pd.read_csv(file_path, delimiter=',')
        
        # 성공 및 실패 카운트 초기화
        success_count = 0
        failure_count = 0
        
        # SMILES 표준화 후 캐노니컬 SMILES로 변환
        results = df['SMILES'].apply(rdkit_standardization)
        
        # Canonical SMILES 및 성공 여부 분리
        df['Canonical_SMILES'] = results.apply(lambda x: x[0])
        success_count = results.apply(lambda x: x[1]).sum()
        failure_count = len(results) - success_count
        
        # 출력 파일 경로 설정 (원본 파일에 '_out' 추가)
        output_file_path = file_path.replace('.csv', '_out.csv')
        
        # 변환된 결과를 CSV 파일로 저장 (Title, origin_name, SMILES, Canonical SMILES 순서)
        df.to_csv(output_file_path, index=False)
        print(f"Processed file saved as: {output_file_path}")
        
        # 종료 시간 기록 및 경과 시간 출력
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Time taken to process {file_path}: {elapsed_time:.2f} seconds")
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

