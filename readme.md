# Chemical Converter with Images

## Chemical 입력 문자열 타입
- [SMILES]
- [Inchi_key]
- [Inchi]
- [IUPAC]
- [Peptide](e.g. AQRYT)
- [SMARTS]
- [SDF FILE]
- [MOL2 FILE]

## 표준화 선택
- [Chembl pipeline]
- [Rdkit]
- [Molvs]

## 이미지 생성 
- [Chembl API]
- [Pumchem API](안 되는경우도 있음)
- [Rdkit]

## web server 실행
- [server 창 2개를 활용해야함]
- [backend](uvicorn backend:app --reload --host 0.0.0.0 --port 8000)
- [frontend](streamlit run frontend.py --server.port 8501) 
