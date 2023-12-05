import os
import pandas as pd
import numpy as np
import shutil
from glob import glob
from tqdm import tqdm


def cat_intersection_ids(case_path, label_path):
    """
    Finds the intersection of TCGA case IDs and labels, and saves the intersected IDs to a CSV file.

    Args:
        case_path: + '/GBM_T1C_first_nrrd' Path to the directory with case files named by TCGA ID.
        label_path: Path to the CSV file with columns: ID, OS, and OS.time.

    Returns:
        intersected_ids: A sorted list of intersected IDs.

    Usage:
        intersected_ids = cat_intersection_ids(case_path, label_path)
        print(np.shape(intersected_ids))
    """
    df = pd.read_csv(label_path + '/02_FJMU_label.csv')

    # Ensure the ID column is of type string
    df['ID'] = df['ID'].astype(str)
    # Pad the ID values with leading zeros to make them at least three digits
    df['ID'] = df['ID'].str.zfill(3)
    ids = set(df['ID'].values)

    filenames = set(os.listdir(case_path + '/02_GBM-SEG-FINAL'))
    intersection = ids.intersection(filenames)

    intersected_ids = sorted(intersection)

    intersection_result = os.path.join(label_path, 'intersected_ids.csv')
    pd.Series(intersected_ids).to_csv(intersection_result, index=False, header=['ID'])

    return intersected_ids


def processing_interactions(intersected_ids, label_path, case_path):
    """
    Process and save the intersected rows and folders from label CSV and case folders.

    Args:
        intersection_ids: List of intersected IDs.
        label_path: Path to a CSV file with the original row names. Intersect rows with intersection_ids
                    and save them to a new CSV file '01_TCGA_label_intersection.csv'.
        case_path: Path to the folders containing case folders. Each folder includes an image file (nrrd)
                   and a seg file (nrrd).

    Usage:
        intersected_ids = cat_intersection_ids(case_path, label_path)
        processing_interactions(intersected_ids, label_path, case_path)
    """

    df_labels = pd.read_csv(os.path.join(label_path, '02_FJMU_label.csv'))

    # Ensure the ID column is a string and zero-padded
    df_labels['ID'] = df_labels['ID'].astype(str).str.zfill(3)

    # Also ensure the intersected_ids are strings and zero-padded
    intersected_ids = [str(id).zfill(3) for id in intersected_ids]

    # Filter the rows that match the intersection IDs
    intersected_rows = df_labels[df_labels['ID'].isin(intersected_ids)]

    # Save the intersected rows to a new CSV file in the result directory
    intersected_rows.to_csv(os.path.join(label_path, '02_FJMU_label_intersection.csv'), index=False)

    intersection_folder_path = os.path.join(case_path, 'intersection_folders')
    if not os.path.exists(intersection_folder_path):
        os.makedirs(intersection_folder_path)

    # Process and copy the case folders
    for case_id in intersected_ids:
        # Define the source and destination paths
        source_folder = os.path.join(case_path + '/02_GBM-SEG-FINAL', case_id)
        destination_folder = os.path.join(intersection_folder_path, case_id)

        # Check if the source folder exists
        if os.path.exists(source_folder):
            # Copy the folder to the destination
            shutil.copytree(source_folder, destination_folder, dirs_exist_ok=True)

import os
from glob import glob
from tqdm import tqdm

def achieve_img_and_mask_path(case_path):
    """
    Achieve the mask_path and image_path by 'glob’ function.

    Args:
        case_path: Path including 'intersection_folders'.

    Returns:
        image_path: List of image.nii.gz paths.
        mask_path: List of mask_pre-tumor-T1.nii.gz paths.
        missing_folders: List of folders that are missing image or mask files.

    Usage:
        image_path, mask_path, missing_folders = achieve_img_and_mask_path(case_path)
        print("Image Paths Count:", len(image_paths))
        print("Mask Paths Count:", len(mask_paths))
        print("Missing Folders Count:", len(missing_folders))

    """
    image_path = []
    mask_path = []
    found_folders = []

    # 获取'intersection_folders'中所有数字编码的文件夹路径
    coded_folders_path = os.path.join(case_path, 'intersection_folders', '*')
    coded_folders = sorted(glob(coded_folders_path))

    # 遍历每个数字编码的文件夹
    for coded_folder in tqdm(coded_folders):
        # 在当前数字编码文件夹中搜索对应的MASK文件夹
        mask_folder_path = os.path.join(coded_folder, 't1_fl2d_tra-MASK')
        mask_folder = glob(mask_folder_path)

        # 如果MASK文件夹存在
        if mask_folder:
            # 在MASK文件夹中搜索image和mask文件
            data_paths = glob(os.path.join(mask_folder[0], 'image.nii.gz'))
            mask_paths_glob = glob(os.path.join(mask_folder[0], 'mask_pre-tumor-T1.nii.gz'))

            # 如果找到了图像和掩模文件，将它们的路径添加到列表中
            if data_paths and mask_paths_glob:
                image_path.append(data_paths[0])
                mask_path.append(mask_paths_glob[0])
                found_folders.append(coded_folder)  # 添加找到文件的文件夹路径

    # 所有文件夹的编号
    all_folder_numbers = [os.path.basename(folder) for folder in coded_folders]
    # 找到文件的文件夹的编号
    found_folder_numbers = [os.path.basename(folder) for folder in found_folders]
    # 缺失文件的文件夹的编号
    missing_folder_numbers = list(set(all_folder_numbers) - set(found_folder_numbers))

    return image_path, mask_path, missing_folder_numbers

# # Preprocessing steps
# case_path = 'D:/projects/ITHscore/02_FJMU_cases'
# label_path = 'D:/projects/ITHscore/02_FJMU_labels'
# cat_intersection_ids(case_path, label_path)
# intersected_ids = cat_intersection_ids(case_path, label_path)
# processing_interactions(intersected_ids, label_path, case_path)
#
# # Get image and mask path
# # 使用函数
# image_path, mask_path, missing_folders = achieve_img_and_mask_path(case_path)
# print("Image Paths Count:", len(image_path))
# print("Mask Paths Count:", len(mask_path))
# print("Image Paths:", image_path[:5])
# print("Mask Paths:", mask_path[:5])

