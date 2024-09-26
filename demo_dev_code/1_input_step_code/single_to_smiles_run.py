from chem_converter import (
    CanonicalToSmiles,
    InChIKeyToSmiles,
    InchiToSmiles,
    sdfToSmiles,
    molToSmiles,
    iupacToSmiles,
    smartsToSmiles,
    peptideToSmiles,
    getStereoChoice,
    singleHandleChoice
)

def main():
    while True:
        print("\nInput 형식을 선택하세요:")
        print("1. Just SMILES or isoSMILES")
        print("2. Inchi_key")
        print("3. Inchi")
        print("4. IUPAC")
        print("5. peptide sequence")
        print("6. SMARTS")
        print("7. sdf_path")
        print("8. mol_path")
        print("9. SMILES to sdf(smiles와 sdf 파일 경로와 이름 입력)")
        print("10. SMILES to mol2(smiles와 mol2 파일 경로와 이름 입력)")
        print("0. 종료")

        choice = input("선택: ").strip()
        
        if choice == "9" or choice == "10":
            singleHandleChoice(choice)
            # get_stereo_choice를 건너뛰기 위해 루프의 다음 단계로 이동
            continue
        elif choice == "0":
            pass
        elif choice not in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']:
            print("유효한 선택이 아닙니다. 다시 시도해주세요.")

        if choice != "0":
            stereo_value = input("Isomer 적용 여부 (y, n): ").strip()
            consider_stereochemistry = getStereoChoice(stereo_value)

        if choice == "1":
            isosmile = input("SMILES 입력: ").strip()
            print("CanonicalToSmiles: ", CanonicalToSmiles(isosmile, consider_stereochemistry))

        elif choice == "2":
            inchi_key = input("InChIKey 입력: ").strip()
            print("InChIKeyToSmiles: ", InChIKeyToSmiles(inchi_key, consider_stereochemistry))

        elif choice == "3":
            inchi = input("InChI 입력: ").strip()
            print("InchiToSmiles: ", InchiToSmiles(inchi, consider_stereochemistry))

        elif choice == "4":
            chemical_IUPAC = input("IUPAC 이름 입력: ").strip()
            print("IUPACToSmiles: ", iupacToSmiles(chemical_IUPAC, consider_stereochemistry))

        elif choice == "5":
            sequence = input("Peptide sequence 입력 (1-letter code): ").strip()
            print("Peptide to SMILES: ", peptideToSmiles(sequence))
        
        elif choice == "6":
            smarts = input("smarts 입력 : ").strip()
            print("Peptide to SMILES: ", smartsToSmiles(smarts, consider_stereochemistry))
                
        elif choice == "7":
            sdf_path = input("SDF 파일 경로 입력: ").strip()
            print("sdfToSmiles: ", sdfToSmiles(sdf_path, consider_stereochemistry))

        elif choice == "8":
            mol_path = input("Mol 파일 경로 입력: ").strip()
            print("molToSmiles: ", molToSmiles(mol_path, consider_stereochemistry))
            
        # elif choice == "9":
        #     smiles = input("smiles 입력 : ")
        #     sdf_path = input("sdf 파일 이름")
        #     print("생성된 SDF 파일 경로: ", smilesToSDF(smiles, sdf_path))

        elif choice == "0":
            print("프로그램을 종료합니다.")
            break

        else:
            print("잘못된 선택입니다. 다시 선택해주세요.")

if __name__ == '__main__':
    main()