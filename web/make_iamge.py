import os
import time
from io import BytesIO
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageOps
import requests
import cairosvg
from urllib.parse import quote
from requests.exceptions import HTTPError, ConnectionError, Timeout, RequestException



# Function to convert RDKit SMILES to image with larger size and white background
def rdkitSmilesToImage(smiles: str, sanitize: bool = True, size=(600, 600)):
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=sanitize)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        # Generate image with larger size
        img = Draw.MolToImage(mol, size=size)
        
        # Convert image to RGBA and add a white background
        img = img.convert("RGBA")
        background = Image.new("RGBA", img.size, (255, 255, 255))  # White background
        img_with_bg = Image.alpha_composite(background, img).convert("RGB")
        
        # Save the image to a BytesIO object
        img_byte_arr = BytesIO()
        img_with_bg.save(img_byte_arr, format='PNG')
        img_byte_arr.seek(0)
        
        return img_byte_arr  # Return as BytesIO object
    
    except Exception as e:
        print(f"rdkit SMILES 변환에 실패했습니다: {smiles}, 에러: {str(e)}")
        return False

# Function to get PubChem PNG with larger size and white background
def pubchemSmilesToPng(smiles: str, retries: int = 3, delay: int = 2, size=(600, 600)):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    if not smiles or len(smiles.strip()) == 0:
        print("SMILES 문자열이 비어 있습니다. 변환을 건너뜁니다.")
        return False
    
    encoded_smiles = quote(smiles)
    cid_url = f"{base_url}/smiles/{encoded_smiles}/cids/JSON"
    
    for attempt in range(retries):
        try:
            # Get CID for the SMILES
            cid_response = requests.get(cid_url)
            cid_response.raise_for_status()
            
            cid_data = cid_response.json()
            if "IdentifierList" in cid_data and "CID" in cid_data["IdentifierList"]:
                cid = cid_data["IdentifierList"]["CID"][0]
                
                # Get image using CID
                png_url = f"{base_url}/cid/{cid}/PNG"
                png_response = requests.get(png_url)
                
                if png_response.status_code == 200:
                    img_byte_arr = BytesIO(png_response.content)
                    img = Image.open(img_byte_arr)
                    
                    # Resize the image
                    img = img.resize(size, Image.Resampling.LANCZOS)
                    
                    # Add a white background if the image has transparency
                    if img.mode == "RGBA":
                        background = Image.new("RGBA", img.size, (255, 255, 255))
                        img_with_bg = Image.alpha_composite(background, img).convert("RGB")
                    else:
                        img_with_bg = img
                    
                    # Save the image to a BytesIO object
                    final_img_byte_arr = BytesIO()
                    img_with_bg.save(final_img_byte_arr, format="PNG")
                    final_img_byte_arr.seek(0)
                    
                    return final_img_byte_arr  # Return as BytesIO object
                
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

# Function to get ChEMBL image and add white background
def chemblSmilesToPng(smiles: str, size=(600, 600)):
    svg_url = "https://www.ebi.ac.uk/chembl/api/utils/smiles2svg"
    
    try:
        # Post the SMILES to the ChEMBL API to get SVG data
        response = requests.post(svg_url, data=smiles)
        response.raise_for_status()
        
        svg_data = response.content
        
        # Convert the SVG to PNG and increase size
        png_byte_arr = BytesIO()
        cairosvg.svg2png(bytestring=svg_data, write_to=png_byte_arr, output_width=size[0], output_height=size[1])
        png_byte_arr.seek(0)
        
        # Open the PNG image and add a white background if needed
        img = Image.open(png_byte_arr)
        if img.mode == "RGBA":
            background = Image.new("RGBA", img.size, (255, 255, 255))  # White background
            img_with_bg = Image.alpha_composite(background, img).convert("RGB")
        else:
            img_with_bg = img
        
        # Save the final image to BytesIO and return
        final_img_byte_arr = BytesIO()
        img_with_bg.save(final_img_byte_arr, format="PNG")
        final_img_byte_arr.seek(0)
        
        return final_img_byte_arr  # Return as BytesIO object
    
    except Exception as e:
        print(f"chembl SMILES 변환에 실패했습니다: {smiles}, 에러: {str(e)}")
        return False
    
# def rdkitSmilesToImage(smiles: str, sanitize: bool = True):
#     try:
#         mol = Chem.MolFromSmiles(smiles, sanitize=sanitize)
#         if mol is None:
#             raise ValueError(f"Invalid SMILES: {smiles}")
        
#         # 이미지를 메모리로 생성
#         img = Draw.MolToImage(mol, size=(300, 300))
#         img_byte_arr = BytesIO()
#         img.save(img_byte_arr, format='PNG')
#         img_byte_arr.seek(0)
        
#         # 이미지를 바이트 스트림으로 반환
#         return img_byte_arr  # StreamingResponse 대신 BytesIO 객체를 반환
    
#     except Exception as e:
#         print(f"rdkit SMILES 변환에 실패했습니다: {smiles}, 에러: {str(e)}")
#         return False

    
# def pubchemSmilesToPng(smiles: str, retries: int = 3, delay: int = 2):
#     base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
#     if not smiles or len(smiles.strip()) == 0:
#         print("SMILES 문자열이 비어 있습니다. 변환을 건너뜁니다.")
#         return False
    
#     encoded_smiles = quote(smiles)
#     cid_url = f"{base_url}/smiles/{encoded_smiles}/cids/JSON"
    
#     for attempt in range(retries):
#         try:
#             # CID 가져오기
#             cid_response = requests.get(cid_url)
#             cid_response.raise_for_status()
            
#             cid_data = cid_response.json()
#             if "IdentifierList" in cid_data and "CID" in cid_data["IdentifierList"]:
#                 cid = cid_data["IdentifierList"]["CID"][0]
                
#                 # CID로 이미지 가져오기
#                 png_url = f"{base_url}/cid/{cid}/PNG"
#                 png_response = requests.get(png_url)
                
#                 if png_response.status_code == 200:
#                     img_byte_arr = BytesIO(png_response.content)
#                     img_byte_arr.seek(0)
                    
#                     # 이미지를 바이트 스트림으로 반환
#                     return img_byte_arr  # StreamingResponse 대신 BytesIO 객체를 반환
#                 else:
#                     print(f"pubchem PNG 가져오는 중 오류 발생: {png_response.status_code}")
#                     return False
#             else:
#                 print("pubchem CID를 찾을 수 없습니다.")
#                 return False
#         except (HTTPError, ConnectionError, Timeout, RequestException) as err:
#             print(f"pubchem 오류 발생: {err}")
#             if cid_response.status_code in [429, 503]:
#                 time.sleep(delay)
#                 delay *= 2
#             else:
#                 return False
#     return False

# def chemblSmilesToPng(smiles: str):
#     svg_url = "https://www.ebi.ac.uk/chembl/api/utils/smiles2svg"
    
#     try:
#         # ChEMBL API로 SMILES를 보내서 SVG 데이터 가져오기
#         response = requests.post(svg_url, data=smiles)
#         response.raise_for_status()
        
#         svg_data = response.content
        
#         # SVG 데이터를 PNG로 변환
#         png_byte_arr = BytesIO()
#         cairosvg.svg2png(bytestring=svg_data, write_to=png_byte_arr)
#         png_byte_arr.seek(0)
        
#         # PNG 이미지를 바이트 스트림으로 반환
#         return png_byte_arr  # Response 대신 BytesIO 객체를 반환
    
#     except Exception as e:
#         print(f"chembl SMILES 변환에 실패했습니다: {smiles}, 에러: {str(e)}")
#         return False