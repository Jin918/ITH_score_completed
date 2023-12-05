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
    df = pd.read_csv(label_path + '/03_ZDFE_label.csv')
    df = pd.read_csv(label_path + '/03_ZDFE_label.csv')

    # Ensure the ID column is of type string
    df['ID'] = df['ID'].astype(str)
    # Pad the ID values with leading zeros to make them at least three digits
    df['ID'] = df['ID'].str.zfill(3)
    ids = set(df['ID'].values)

    filenames = set(os.listdir(case_path + '/03_N4_GBMfix'))
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

    df_labels = pd.read_csv(os.path.join(label_path, '03_ZDFE_label.csv'))

    # Ensure the ID column is a string and zero-padded
    df_labels['ID'] = df_labels['ID'].astype(str).str.zfill(3)

    # Also ensure the intersected_ids are strings and zero-padded
    intersected_ids = [str(id).zfill(3) for id in intersected_ids]

    # Filter the rows that match the intersection IDs
    intersected_rows = df_labels[df_labels['ID'].isin(intersected_ids)]

    # Save the intersected rows to a new CSV file in the result directory
    intersected_rows.to_csv(os.path.join(label_path, '03_ZDFE_label_intersection.csv'), index=False)

    intersection_folder_path = os.path.join(case_path, 'intersection_folders')
    if not os.path.exists(intersection_folder_path):
        os.makedirs(intersection_folder_path)

    # Process and copy the case folders
    for case_id in intersected_ids:
        # Define the source and destination paths
        source_folder = os.path.join(case_path + '/03_N4_GBMfix', case_id)
        destination_folder = os.path.join(intersection_folder_path, case_id)

        # Check if the source folder exists
        if os.path.exists(source_folder):
            # Copy the folder to the destination
            shutil.copytree(source_folder, destination_folder, dirs_exist_ok=True)


def achieve_img_and_mask_path(case_path):
    """
    Achieve the mask_path and image_path by 'glob’ function.

    Args:
        case_path: + intersection_fonders.

    Returns:
        image_path: List of image.nrrd.
        mask_path: List of seg.nrrd.

    Usage:
        image_path, mask_path = achieve_img_and_mask_path(case_path)
        print("Image Paths Count:", len(image_paths))
        print("Mask Paths Count:", len(mask_paths))

    """
    image_path = []
    mask_path = []

    search_path = os.path.join(case_path, 'intersection_folders', '*')

    folders = sorted(glob(search_path))

    for folder in tqdm(folders):
        data_paths = glob(os.path.join(folder, 'T1C.nrrd'))
        mask_paths_glob = glob(os.path.join(folder, 'T1Cseg-label.nrrd'))

        if data_paths and mask_paths_glob:
            image_path.append(data_paths[0])
            mask_path.append(mask_paths_glob[0])

    return image_path, mask_path

# Preprocessing steps
case_path = 'D:/projects/ITHscore/03_ZDFE_cases'
label_path = 'D:/projects/ITHscore/03_ZDFE_labels'
cat_intersection_ids(case_path, label_path)
intersected_ids = cat_intersection_ids(case_path, label_path)
processing_interactions(intersected_ids, label_path, case_path)

# Get image and mask path
image_path, mask_path = achieve_img_and_mask_path(case_path)
print("Image Paths Count:", len(image_path))
print("Mask Paths Count:", len(mask_path))
print("Image Paths:", image_path[:5])
print("Mask Paths:", mask_path[:5])

