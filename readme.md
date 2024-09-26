# Chemical Converter with Images

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Chemical 입력 문자열 타입
- SMILES
- Inchi_key
- Inchi
- IUPAC
- Peptide (e.g. AQRYT)
- SMARTS
- SDF FILE
- MOL2 FILE

## 표준화 선택
- Chembl pipeline
- Rdkit
- Molvs

## 이미지 생성 
- Chembl API
- Pumchem API (안 되는경우도 있음)
- Rdkit

## web server 실행

1. **Clone the repository**:
    ```bash
    git clone https://github.com/danid77/LNP_CBANK_WEB.git
    ```

2. **Navigate to the project directory**:
    ```bash
    cd LNP_CBANK_WEB
    ```

3. **Create and activate a virtual environment** (recommended):
    ```bash
    conda env create -f enviroment.yml
    conda activate to_smiles
    ```

4. **run**
- [server 창 2개를 활용해야함]
- [backend]
    ```bash
    uvicorn backend:app --reload --host 0.0.0.0 --port 8000
    ```

- [frontend]
    ```bash
    streamlit run frontend.py --server.port 8501
    ```
