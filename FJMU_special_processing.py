import os
import shutil
from dcmrtstruct2nii import dcmrtstruct2nii


def extract_and_save_dicoms(source_dir, target_dir):
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    for patient_case in os.listdir(source_dir):
        patient_case_path = os.path.join(source_dir, patient_case)
        if os.path.isdir(patient_case_path):
            target_patient_dir = os.path.join(target_dir, patient_case)
            if not os.path.exists(target_patient_dir):
                os.makedirs(target_patient_dir)

            for file in os.listdir(patient_case_path):
                if file.endswith('.dcm'):
                    source_file_path = os.path.join(patient_case_path, file)
                    target_file_path = os.path.join(target_patient_dir, file)
                    shutil.copy2(source_file_path, target_file_path)

    print("DICOM files have been copied to the target directory.")


def move_dicom_files(target_dir):
    for patient_case in os.listdir(target_dir):
        patient_case_path = os.path.join(target_dir, patient_case)
        if os.path.isdir(patient_case_path):
            image_folder = os.path.join(patient_case_path, 'img')
            mask_folder = os.path.join(patient_case_path, 'RT_Struct')

            if not os.path.exists(image_folder):
                os.makedirs(image_folder)
            if not os.path.exists(mask_folder):
                os.makedirs(mask_folder)

            for file in os.listdir(patient_case_path):
                if file.startswith('MR'):
                    source_file_path = os.path.join(patient_case_path, file)
                    target_file_path = os.path.join(image_folder, file)
                    shutil.move(source_file_path, target_file_path)
                elif file.startswith('RS'):
                    source_file_path = os.path.join(patient_case_path, file)
                    target_file_path = os.path.join(mask_folder, file)
                    shutil.move(source_file_path, target_file_path)

            if not os.listdir(image_folder) or not os.listdir(mask_folder):
                print(f"Warning: Missing files in image or mask folders for case {patient_case}")
            else:
                print(f"Case {patient_case}: Image and mask folders are populated.")


def process_all_cases(base_dir):
    for case in os.listdir(base_dir):
        case_path = os.path.join(base_dir, case)
        if os.path.isdir(case_path):
            image_folder = os.path.join(case_path, 'img')
            mask_folder = os.path.join(case_path, 'RT_Struct')
            merge_folder = os.path.join(case_path, 'merge')

            if not os.path.exists(merge_folder):
                os.makedirs(merge_folder)

            rtstruct_file = None
            for file in os.listdir(mask_folder):
                if file.startswith('RS') and file.endswith('.dcm'):
                    rtstruct_file = os.path.join(mask_folder, file)
                    break

            if rtstruct_file:
                dcmrtstruct2nii(rtstruct_file, image_folder, merge_folder)
                print(f"Processed case: {case}")
            else:
                print(f"No RTSTRUCT file found for case: {case}")


# Main execution
if __name__ == "__main__":
    source_directory = "D:\\projects\\ITHscore\\02_FJMU_cases\\FJMU_GBM_152"
    target_directory = "D:\\projects\\ITHscore\\02_FJMU_cases\\FJMU_GBM_152_beta"

    # 提取并保存DICOM文件
    extract_and_save_dicoms(source_directory, target_directory)

    # 移动DICOM文件
    move_dicom_files(target_directory)

    # 处理所有病例
    process_all_cases(target_directory)


####check point####
def check_mask_files_count(target_dir):
    """
    Checks the number of files starting with 'mask' in each 'merge' folder.
    Returns the folder names where the count of 'mask' files is not 1.

    Args:
    target_dir: The directory to search for 'merge' folders.

    Returns:
    A list of folder names where the count of 'mask' files is not 1.
    """
    folders_with_issues = []

    for patient_case in os.listdir(target_dir):
        merge_folder = os.path.join(target_dir, patient_case, 'merge')

        # Check if the merge folder exists
        if os.path.exists(merge_folder) and os.path.isdir(merge_folder):
            mask_files = [f for f in os.listdir(merge_folder) if f.startswith('mask')]

            # Check if the count of mask files is not 1
            if len(mask_files) != 1:
                folders_with_issues.append(patient_case)

    return folders_with_issues

# 使用示例
target_directory = "D:\\projects\\ITHscore\\02_FJMU_cases\\FJMU_GBM_152_beta"
folders_with_extra_mask_files = check_mask_files_count(target_directory)
print("Folders with incorrect number of mask files:", folders_with_extra_mask_files)















# import os
# import shutil
# #
# #
# def extract_and_save_dicoms(source_dir, target_dir):
#     if not os.path.exists(target_dir):
#         os.makedirs(target_dir)
#
#     # 遍历source_dir中的每个子目录（即每个病例编号目录）
#     for patient_case in os.listdir(source_dir):
#         patient_case_path = os.path.join(source_dir, patient_case)
#
#         # 确保它是一个目录
#         if os.path.isdir(patient_case_path):
#             target_patient_dir = os.path.join(target_dir, patient_case)
#             if not os.path.exists(target_patient_dir):
#                 os.makedirs(target_patient_dir)
#
#             # 在该病例编号目录中查找.dcm文件
#             for file in os.listdir(patient_case_path):
#                 if file.endswith('.dcm'):
#                     source_file_path = os.path.join(patient_case_path, file)
#                     target_file_path = os.path.join(target_patient_dir, file)
#                     shutil.copy2(source_file_path, target_file_path)
#
#     print("DICOM files have been copied to the target directory.")
#
#
# # 使用示例
# source_directory = "D:\\projects\\ITHscore\\02_FJMU_cases\\FJMU_GBM_152"
# target_directory = "D:\\projects\\ITHscore\\02_FJMU_cases\\FJMU_GBM_152_beta"
# extract_and_save_dicoms(source_directory, target_directory)
#
#
# def move_dicom_files(target_dir):
#     for patient_case in os.listdir(target_dir):
#         patient_case_path = os.path.join(target_dir, patient_case)
#
#         # 确保它是一个目录
#         if os.path.isdir(patient_case_path):
#             image_folder = os.path.join(patient_case_path, 'image')
#             mask_folder = os.path.join(patient_case_path, 'mask')
#
#             # 创建image和mask文件夹（如果它们不存在）
#             if not os.path.exists(image_folder):
#                 os.makedirs(image_folder)
#             if not os.path.exists(mask_folder):
#                 os.makedirs(mask_folder)
#
#             # 在病例编号目录中查找MR和RS开头的.dcm文件
#             for file in os.listdir(patient_case_path):
#                 if file.startswith('MR'):
#                     source_file_path = os.path.join(patient_case_path, file)
#                     target_file_path = os.path.join(image_folder, file)
#                     shutil.move(source_file_path, target_file_path)
#                 elif file.startswith('RS'):
#                     source_file_path = os.path.join(patient_case_path, file)
#                     target_file_path = os.path.join(mask_folder, file)
#                     shutil.move(source_file_path, target_file_path)
#
#     # 检查image和mask文件夹是否都有文件
#     for patient_case in os.listdir(target_dir):
#         patient_case_path = os.path.join(target_dir, patient_case)
#         image_folder = os.path.join(patient_case_path, 'image')
#         mask_folder = os.path.join(patient_case_path, 'mask')
#
#         if os.path.isdir(patient_case_path):
#             if not os.listdir(image_folder) or not os.listdir(mask_folder):
#                 print(f"Warning: Missing files in image or mask folders for case {patient_case}")
#             else:
#                 print(f"Case {patient_case}: Image and mask folders are populated.")
#
# # 使用示例
# target_directory = "D:\\projects\\ITHscore\\02_FJMU_cases\\FJMU_GBM_152_beta"
# move_dicom_files(target_directory)
#
# # 参数分别为struct文件、病人图像文件夹，输出文件夹
#
# import os
# from dcmrtstruct2nii import dcmrtstruct2nii
#
# def process_all_cases(base_dir):
#     # 遍历base_dir下的所有子文件夹
#     for case in os.listdir(base_dir):
#         case_path = os.path.join(base_dir, case)
#
#         # 确保它是一个目录
#         if os.path.isdir(case_path):
#             image_folder = os.path.join(case_path, 'img')
#             mask_folder = os.path.join(case_path, 'RT')
#             merge_folder = os.path.join(case_path, 'merge')
#
#             # 创建merge文件夹（如果它不存在）
#             if not os.path.exists(merge_folder):
#                 os.makedirs(merge_folder)
#
#             # 查找RTSTRUCT文件
#             rtstruct_file = None
#             for file in os.listdir(mask_folder):
#                 if file.startswith('RS') and file.endswith('.dcm'):
#                     rtstruct_file = os.path.join(mask_folder, file)
#                     break
#
#             # 如果找到了RTSTRUCT文件，就进行处理
#             if rtstruct_file:
#                 dcmrtstruct2nii(rtstruct_file, image_folder, merge_folder)
#                 print(f"Processed case: {case}")
#             else:
#                 print(f"No RTSTRUCT file found for case: {case}")
#
# # 使用示例
# base_directory = "D:\\projects\\ITHscore\\02_FJMU_cases\\FJMU_GBM_152_beta"
# process_all_cases(base_directory)
