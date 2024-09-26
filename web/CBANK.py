import streamlit as st
import requests
from io import BytesIO
from PIL import Image
import base64

# FastAPI 백엔드 URL
BACKEND_URL_CONVERT = "http://localhost:8503/convert/"
BACKEND_URL_CONVERT_FILE = "http://localhost:8503/convert_file/"
BACKEND_URL_STN = "http://localhost:8503/stn/"
BACKEND_URL_ISOMER = "http://localhost:8503/isomer/"
BACKEND_URL_IMAGES = "http://localhost:8503/images/"  # 이미지 생성용 백엔드 URL

# 페이지 제목
st.title("Chemical Converter with Multiple Images")

# 사용자에게 표시할 값과 실제로 사용할 값을 매핑
type_display = ["SMILES", "Inchi_key", "Inchi", "IUPAC", "Peptide", "SMARTS", "SDF FILE", "MOL2 FILE"]
type_value = ["smiles", "inchi_key", "inchi", "iupac", "peptide", "smarts", "sdf", "mol2"]

# 사용자가 선택한 값을 실제로 사용할 값으로 매핑
input_type_display = st.selectbox("Chemical 입력 문자열 타입", type_display)
input_type = type_value[type_display.index(input_type_display)]

# # 입력 형식 선택
# input_type = st.selectbox("Input 형식을 선택하세요", 
#                           ["smiles", "inchi_key", "inchi", "iupac", "peptide", "smarts", "sdf", "mol2"])

# 파일 업로드 및 입력 값 입력
input_file = None
input_value = None

# 입력 값 입력
if input_type in ["sdf", "mol2"]:
    input_file = st.file_uploader(f"{input_type.upper()} 파일을 업로드하세요", type=[input_type])
else:
    input_value = st.text_input(f"{input_type.upper()} 입력").strip()

# 표준화 처리
# 사용자에게 표시할 값과 실제로 사용할 값을 매핑
stn_display = ["Chembl pipeline", "Rdkit", "Molvs"]
stn_value = ["chembl_pipeline", "rdkit", "molvs"]

# 사용자가 선택한 값을 실제로 사용할 값으로 매핑
select_display = st.selectbox("표준화 모드 선택", stn_display)
select_stn = stn_value[stn_display.index(select_display)]

# Isomer 옵션 선택
select_isomer = st.radio("Isomer 유지 여부", ("Yes", "No"))

# API 요청 함수
def postRequest(url, payload, files=None):
    try:
        if files:
            response = requests.post(url, files=files, data=payload)
        else:
            response = requests.post(url, json=payload)
            
        if response.status_code == 200:
            return response.json(), None
        else:
            return None, f"오류 발생: {response.json().get('detail', '상세 정보 없음')}"
    except Exception as e:
        return None, f"요청 중 오류 발생: {str(e)}"

# 변환 요청 함수 (문자열 처리)
def requestConversion(input_type, input_value):
    payload_convert = {
        "input_type": input_type, 
        "input_value": input_value
    }
    return postRequest(BACKEND_URL_CONVERT, payload_convert)

# 변환 요청 함수 (파일 업로드 처리)
def requestConversionFile(input_type, input_file=None):
    files = {"file": (input_file.name, input_file, input_file.type)}
    form_data = {"input_type": input_type}
    return postRequest(BACKEND_URL_CONVERT_FILE, form_data, files=files)

# def requestConversionFile(input_type, input_value, input_file=None):
#     if input_type in ["sdf", "mol2"] and input_file:
#         files = {"file": (input_file.name, input_file, input_file.type)}
#         form_data = {"input_type": input_type}
#         return postRequest(BACKEND_URL_CONVERT, form_data, files=files)
#     else:
#         payload_convert = {
#             "input_type": input_type, 
#             "input_value": input_value
#         }
#         return postRequest(BACKEND_URL_CONVERT, payload_convert)
    
# 표준화 요청 함수
def requestStandardization(convert_smiles, select_stn):
    payload_stn = {
        "convert_smiles": convert_smiles, 
        "select_stn": select_stn
    }
    return postRequest(BACKEND_URL_STN, payload_stn)

# isomer 선택여부
def requestIsomer(convert_smiles, select_isomer):
    # print(select_isomer, type(select_isomer))
    payload_isomer = {
        "convert_smiles": convert_smiles, 
        "select_isomer": select_isomer == "Yes"
    }
    return postRequest(BACKEND_URL_ISOMER, payload_isomer)

# 표준화 및 Isomer 결과 표시
def displayResults(stn_result, isomer_result):
    # Section Header
    st.markdown("### 변환 및 표준화 결과")
    
    # Displaying the standardization result in a clean format
    with st.container():
        st.write("**Standardized SMILES (using selected method):**")
        st.code(stn_result, language="markdown")
    
    # Section for Isomer result
    st.markdown("### Isomer 선택 및 Canonical 변환 결과")
    
    with st.container():
        st.write("**Isomer-selected and Canonical SMILES (Rdkit 기준):**")
        st.code(isomer_result, language="markdown")

# Example function to display results
def displayResults():
    if st.session_state["stn_result"]:
        st.markdown("### 변환 및 표준화 결과")
        st.code(st.session_state["stn_result"], language="markdown")

    if st.session_state["isomer_result"]:
        st.markdown("### Isomer 선택 및 Canonical 변환 결과")
        st.code(st.session_state["isomer_result"], language="markdown")

# # 이미지 생성 요청 함수
# def requestImages(convert_smiles):
#     payload_images = {"convert_smiles": convert_smiles}
#     response_images = requests.post(BACKEND_URL_IMAGES, json=payload_images)

#     if response_images.status_code == 200:
#         return response_images.json(), None
#     else:
#         return None, response_images.json().get('detail')

# RDKit, PubChem, ChEMBL 이미지 표시 함수
# def displayImages(result_images):
#     if result_images.get("rdkit"):
#         rdkit_image_data = base64.b64decode(result_images["rdkit"])
#         rdkit_image = Image.open(BytesIO(rdkit_image_data))
#         st.image(rdkit_image, caption="RDKit Molecule Image")

#     if result_images.get("pubchem"):
#         pubchem_image_data = base64.b64decode(result_images["pubchem"])
#         pubchem_image = Image.open(BytesIO(pubchem_image_data))
#         st.image(pubchem_image, caption="PubChem Molecule Image")

#     if result_images.get("chembl"):
#         chembl_image_data = base64.b64decode(result_images["chembl"])
#         chembl_image = Image.open(BytesIO(chembl_image_data))
#         st.image(chembl_image, caption="ChEMBL Molecule Image")

# 이미지 생성 요청 함수
def requestImages(convert_smiles):
    payload_images = {"convert_smiles": convert_smiles}
    response_images = requests.post(BACKEND_URL_IMAGES, json=payload_images)

    if response_images.status_code == 200:
        return response_images.json(), None
    else:
        return None, response_images.json().get('detail')

# 이미지 및 다운로드 버튼 표시 함수
def displayImagesWithDownload(result_images):
    if result_images.get("rdkit"):
        try:
            rdkit_image_data = base64.b64decode(result_images["rdkit"])
            rdkit_image = Image.open(BytesIO(rdkit_image_data))

            # Display the image
            st.image(rdkit_image, caption="RDKit Molecule Image")

            # Create a download button for the RDKit image
            buffered = BytesIO()
            rdkit_image.save(buffered, format="PNG")
            st.download_button(
                label="Download RDKit Image",
                data=buffered.getvalue(),
                file_name="rdkit_molecule.png",
                mime="image/png"
            )
        except Exception as e:
            st.error("RDKit 이미지를 불러올 수 없습니다.")
    else:
        st.warning("RDKit 이미지 변환이 불가능 합니다.")

    if result_images.get("pubchem"):
        try:
            pubchem_image_data = base64.b64decode(result_images["pubchem"])
            pubchem_image = Image.open(BytesIO(pubchem_image_data))

            # Display the image
            st.image(pubchem_image, caption="PubChem Molecule Image")

            # Create a download button for the PubChem image
            buffered = BytesIO()
            pubchem_image.save(buffered, format="PNG")
            st.download_button(
                label="Download PubChem Image",
                data=buffered.getvalue(),
                file_name="pubchem_molecule.png",
                mime="image/png"
            )
        except Exception as e:
            st.error("PubChem 이미지를 불러올 수 없습니다.")
    else:
        st.warning("PubChem API의 변환이 불가능 합니다. (PubChem에서 해당 SMILES를 찾을 수 없음)")

    if result_images.get("chembl"):
        try:
            chembl_image_data = base64.b64decode(result_images["chembl"])
            chembl_image = Image.open(BytesIO(chembl_image_data))

            # Display the image
            st.image(chembl_image, caption="ChEMBL Molecule Image")

            # Create a download button for the ChEMBL image
            buffered = BytesIO()
            chembl_image.save(buffered, format="PNG")
            st.download_button(
                label="Download ChEMBL Image",
                data=buffered.getvalue(),
                file_name="chembl_molecule.png",
                mime="image/png"
            )
        except Exception as e:
            st.error("ChEMBL 이미지를 불러올 수 없습니다.")
    else:
        st.warning("ChEMBL API의 변환이 불가능 합니다.")

def resetSession():
    st.session_state["stn_result"] = None
    st.session_state["isomer_result"] = None
    st.session_state["result_images"] = None

# Use session state to store results and images
if "stn_result" not in st.session_state:
    st.session_state["stn_result"] = None
if "isomer_result" not in st.session_state:
    st.session_state["isomer_result"] = None
if "result_images" not in st.session_state:
    st.session_state["result_images"] = None

# Main function
def main():
    if st.button("변환 실행"):
        if input_type in ["sdf", "mol2"] and not input_file:
            st.error(f"{input_type.upper()} 파일을 업로드하세요")
            resetSession()
            return
         
        # 첫 번째 변환 요청
        if input_type in ["sdf", "mol2"]:
            result_convert, error_convert = requestConversionFile(input_type, input_file)
        else:
            result_convert, error_convert = requestConversion(input_type, input_value)

        if error_convert:
            st.error(error_convert)
            resetSession()
        else:
            convert_smiles = result_convert.get("convert")

            if convert_smiles:
                # 두 번째 표준화 요청
                result_stn, error_stn = requestStandardization(convert_smiles, select_stn)

                if error_stn:
                    st.error(error_stn)
                    resetSession()
                else:
                    st.session_state["stn_result"] = result_stn.get("stn")
                    
                    # if "InChI: Omitted undefined stereo" in st.session_state["stn_result"]:
                    #     st.warning("변경할 수 없는 SMILES입니다.")
                    if st.session_state["stn_result"]:
                        # 세번째 Isomer 선택 요청 추가
                        result_isomer, error_isomer = requestIsomer(st.session_state["stn_result"], select_isomer)

                        if error_isomer:
                            st.error(f"Isomer 선택 오류: {error_isomer}")
                            resetSession()
                        else:
                            st.session_state["isomer_result"] = result_isomer.get("isomer")
                            if st.session_state["isomer_result"]:
                                # 이미지 생성 요청
                                # result_images, error_images = requestImages(convert_smiles)
                                result_images, error_images = requestImages(st.session_state["isomer_result"])

                                if error_images:
                                    st.error(f"이미지 생성 오류: {error_images}")
                                    resetSession()
                                else:
                                    st.session_state["result_images"] = result_images

    # Display the results and images
    displayResults()

    # Display images with download buttons if available
    if st.session_state["result_images"]:
        displayImagesWithDownload(st.session_state["result_images"])

# 실행
if __name__ == "__main__":
    main()
