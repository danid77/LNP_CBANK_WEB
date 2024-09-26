import csv
import os
import sys
import glob
import time
from rdkit import Chem
from chembl_structure_pipeline import standardizer, checker

# CSV 파일 처리 함수 정의
def process_csv(input_csv, output_csv, success_count, failure_count):
    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        # CSV 파일 구분자가 콤마임을 명시
        reader = csv.DictReader(infile, delimiter=',')
        fieldnames = reader.fieldnames + ['Canonical SMILES']  # 기존 필드에 'Canonical SMILES' 추가
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        
        # 첫 번째 row에 헤더 추가
        writer.writeheader()
        
        # 각 row 처리
        for row in reader:
            smiles = row['smiles']
            
            try:
                # SMILES를 Mol 객체로 변환 (sanitize=False로 RDKit 검증 생략)
                mol = Chem.MolFromSmiles(smiles, sanitize=False)
                if mol is None:
                    raise ValueError(f"SMILES 변환 실패: {smiles}")

                # Mol 객체를 MOL 블록으로 변환
                mol_block = Chem.MolToMolBlock(mol)
                
                # 표준화된 MOL 블록 생성
                std_mol_block = standardizer.standardize_molblock(mol_block)
                
                # 부모 화합물 추출
                parent_mol_block, _ = standardizer.get_parent_molblock(std_mol_block)
                
                # 체크 과정
                issues = checker.check_molblock(parent_mol_block)
                if issues:
                    print(f"Issues found during checking for {row['Title']}: {issues}")
                
                # 부모 화합물 MOL 블록을 RDKit Mol 객체로 변환
                parent_mol = Chem.MolFromMolBlock(parent_mol_block, sanitize=False)
                if parent_mol is None:
                    raise ValueError(f"MOL 블록 변환 실패: {parent_mol_block}")
                
                # 표준화된 canonical SMILES 생성
                canonical_smiles = Chem.MolToSmiles(parent_mol, canonical=False)
                
                # 성공 카운트 증가
                success_count += 1
            
            except Exception as e:
                # 변환 실패 시 canonical_smiles에 "X" 표기
                canonical_smiles = "X"
                print(f"Error occurred while processing {row['Title']}: {e}")
                # 실패 카운트 증가
                failure_count += 1
            
            # 결과를 CSV에 기록
            row['Canonical SMILES'] = canonical_smiles
            writer.writerow(row)
    
    return success_count, failure_count

# 디렉토리의 모든 CSV 파일을 처리하는 함수 정의
def process_directory(directory):
    start_time = time.time()  # 시작 시간 기록
    # 디렉토리 내 모든 CSV 파일 찾기
    csv_files = glob.glob(os.path.join(directory, '*.csv'))
    
    if not csv_files:
        print(f"No CSV files found in directory: {directory}")
        return
    
    total_success = 0
    total_failure = 0
    
    for csv_file in csv_files:
        # 각 CSV 파일에 대해 "_out.csv" 파일명 생성
        base_name = os.path.basename(csv_file)
        file_name, _ = os.path.splitext(base_name)
        output_file = os.path.join(directory, f"{file_name}_out.csv")
        
        print(f"Processing {csv_file} -> {output_file}")
        
        # 파일별 성공/실패 카운트
        success_count = 0
        failure_count = 0
        
        # CSV 파일 처리
        success_count, failure_count = process_csv(csv_file, output_file, success_count, failure_count)
        
        # 총 성공/실패 카운트에 합산
        total_success += success_count
        total_failure += failure_count
    
    end_time = time.time()  # 종료 시간 기록
    elapsed_time = end_time - start_time  # 걸린 시간 계산
    
    print(f"Directory processing completed in {elapsed_time:.2f} seconds.")
    print(f"Total Success: {total_success}, Total Failure: {total_failure}")

# 스크립트 실행
if __name__ == "__main__":
    # 명령줄 인자 처리
    if len(sys.argv) == 1:
        # 명령줄 인자가 없을 경우 현재 디렉토리에서 실행
        current_dir = os.getcwd()
        print(f"현재 디렉토리: {current_dir}")
        process_directory(current_dir)
    
    elif len(sys.argv) == 2:
        # 특정 디렉토리가 주어졌을 경우 해당 디렉토리에서 실행
        target_dir = sys.argv[1]
        
        # 디렉토리가 실제로 존재하는지 확인
        if os.path.isdir(target_dir):
            print(f"지정된 디렉토리: {target_dir}")
            process_directory(target_dir)
        else:
            print(f"오류: 지정된 디렉토리가 존재하지 않습니다 - {target_dir}")
    else:
        print("사용법: python standard.py [directory_path (선택사항)]")

