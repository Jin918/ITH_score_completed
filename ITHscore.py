import os
import numpy as np
import pydicom as dicom
import functools
import concurrent.futures
from scipy import ndimage
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
import SimpleITK as sitk
from utils import read_dcm_series
from utils import extract_feature_unit
from matplotlib import pyplot as plt
from preprocessing import cat_intersection_ids
from preprocessing import processing_interactions
from preprocessing import achieve_img_and_mask_path
import csv
from tqdm import tqdm
random_seed = 42

def load_image(path):
    """
    Load CT image volume
    Args:
        path: Str. Path to the .nii(.gz) file or dicom series directory
    Returns:
        image: Numpy array. The 3D CT volume.
    """
    if os.path.isdir(path):
        sitk_image = read_dcm_series(path)
    else:
        sitk_image = sitk.ReadImage(path)
    image = sitk.GetArrayFromImage(sitk_image)

    return image

def load_seg(path):
    """
    Load segmentation mask
    Args:
        path: Str. Path to the .nii or .dcm mask.
    Returns:
        seg: Numpy array. The mask of ROI with the same shape of image.
    """
    if path.endswith(".dcm"):
        # RTStruct dcm file sometimes cannot be loaded by SimpleITK
        ds = dicom.read_file(path)
        seg = ds.pixel_array
    else:
        sitk_seg = sitk.ReadImage(path)
        seg = sitk.GetArrayFromImage(sitk_seg)

    return seg

def get_largest_slice(img3d, mask3d):
    """
    Get the slice with largest tumor area
    Args:
        img3d: Numpy array. The whole CT volume (3D)
        mask3d: Numpy array. Same size as img3d, binary mask with tumor area set as 1, background as 0
    Returns:
        img: Numpy array. The 2D image slice with largest tumor area
        mask: Numpy array. The subset of mask in the same position of sub_img
    """
    # 计算每个切片中肿瘤区域的面积。mask3d == 1 产生一个布尔数组，其中肿瘤区域为True。
    # np.sum 对每个切片的肿瘤区域进行计数（沿着第一和第二轴），得到每个切片中肿瘤区域的像素数。
    area = np.sum(mask3d == 1, axis=(1, 2))

    # np.argsort 返回数组排序后的索引。[-1]选择了最大面积的索引。
    area_index = np.argsort(area)[-1]

    # 使用找到的索引从原始3D图像中提取出具有最大肿瘤面积的切片。
    img = img3d[area_index, :, :]

    # 使用相同的索引从肿瘤掩码数组中提取对应的切片。
    mask = mask3d[area_index, :, :]

    # 返回找到的图像切片及其对应的肿瘤掩码切片。
    return img, mask


def locate_tumor(img, mask, padding=2):
    """
    Locate and extract tumor from CT image using mask
    Args:
        img: Numpy array. The whole image
        mask: Numpy array. Same size as img, binary mask with tumor area set as 1, background as 0
        padding: Int. Number of pixels padded on each side after extracting tumor
    Returns:
        sub_img: Numpy array. The tumor area defined by mask
        sub_mask: Numpy array. The subset of mask in the same position of sub_img
    """
    top_margin = min(np.where(mask == 1)[0])
    bottom_margin = max(np.where(mask == 1)[0])
    left_margin = min(np.where(mask == 1)[1])
    right_margin = max(np.where(mask == 1)[1])
    # padding two pixels at each edges for further computation
    sub_img = img[top_margin - padding:bottom_margin + padding + 1, left_margin - padding:right_margin + padding + 1]
    sub_mask = mask[top_margin - padding:bottom_margin + padding + 1,
                    left_margin - padding:right_margin + padding + 1]

    return sub_img, sub_mask


def extract_radiomic_features(sub_img, sub_mask, parallel=False, workers=None):
    '''
    Extract radiomic features for each pixel from its surrounding area
    Args:
        sub_img: Numpy array. The tumor area defined by mask
        sub_mask: Numpy array. Same size as sub_img, binary values, 1:tumor area; 0:background
        parallel: Bool. Whether to process in parallel with multiple workers
        workers: Int. Number of workers used to process. Only works when "parallel" set to be True

    Returns:
        features: Dict. A dictionary contains all the radiomic features with keys used in "pyradiomics"
    '''
    features = {}
    first_features = []
    shape_features = []
    glcm_features = []
    gldm_features = []
    glrlm_features = []
    glszm_features = []
    ngtdm_features = []

    if parallel:
        ps, qs = [], []
        partial_extract_feature_unit = functools.partial(extract_feature_unit, sub_img)
        for p in range(len(sub_img)):
            for q in range(len(sub_img[0])):
                if (sub_mask[p][q] == 1):
                    ps.append(p)
                    qs.append(q)
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
            results = executor.map(partial_extract_feature_unit, ps, qs)
            for result in results:
                first_features.append(result["first"])
                shape_features.append(result["shape"])
                glcm_features.append(result["glcm"])
                gldm_features.append(result["gldm"])
                glrlm_features.append(result["glrlm"])
                glszm_features.append(result["glszm"])
                ngtdm_features.append(result["ngtdm"])
    else:
        for p in range(len(sub_img)):
            for q in range(len(sub_img[0])):
                if (sub_mask[p][q] == 1):
                    features_temp = extract_feature_unit(sub_img, p, q, padding=2)
                    first_features.append(features_temp["first"])
                    shape_features.append(features_temp["shape"])
                    glcm_features.append(features_temp["glcm"])
                    gldm_features.append(features_temp["gldm"])
                    glrlm_features.append(features_temp["glrlm"])
                    glszm_features.append(features_temp["glszm"])
                    ngtdm_features.append(features_temp["ngtdm"])
    features['first'] = MinMaxScaler().fit_transform(first_features)
    features['shape'] = MinMaxScaler().fit_transform(shape_features)
    features['glcm'] = MinMaxScaler().fit_transform(glcm_features)
    features['gldm'] = MinMaxScaler().fit_transform(gldm_features)
    features['glrlm'] = MinMaxScaler().fit_transform(glrlm_features)
    features['glszm'] = MinMaxScaler().fit_transform(glszm_features)
    features['ngtdm'] = MinMaxScaler().fit_transform(ngtdm_features)

    return features

def pixel_clustering(sub_mask, features, cluster=6, random_seed=random_seed):
    """
    Args:
        sub_mask: Numpy array. ROI mask only within the tumor bounding box
        features: Numpy array or dict (output of previous step). Matrix of radiomic features. Rows are pixels and columns are features
        cluster: Int. The cluster number in clustering
    Returns:
        label_map: Numpy array. Labels of pixels within tumor. Same size as tumor_img
    """
    if isinstance(features, dict):
        features = np.hstack((features['first'], features['shape'], features['glcm'], features['gldm'],
                              features['glrlm'], features['glszm'], features['ngtdm']))
    features = MinMaxScaler().fit_transform(features)
    label_map = sub_mask.copy()

    clusters = KMeans(n_clusters=cluster,random_state=random_seed).fit_predict(features)
    cnt = 0
    for i in range(len(sub_mask)):
        for j in range(len(sub_mask[0])):
            if sub_mask[i][j] == 1:
                label_map[i][j] = clusters[cnt] + 1
                cnt += 1
            else:
                label_map[i][j] = 0

    return label_map



def visualize(img, sub_img, mask, sub_mask, features, cluster=6, save_dir="result_fig"):
    """
    Args:
        img: Numpy array. Original whole image, used for display
        sub_img: Numpy array. Tumor image
        mask: Numpy array. Same size as img, 1 for tumor and 0 for background, used for display
        sub_mask: Numpy array. Same size as sub_img, 1 for nodule and 0 for background
        features: Numpy array. Matrix of radiomic features. Rows are pixels and columns are features
        cluster: Int or Str. Integer defines the cluster number in clustering. "all" means iterate clusters from 3 to 9 to generate multiple cluster pattern.
    Returns:
        fig: figure for display
    """
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    if cluster != "all":
        if not isinstance(cluster, int):
            raise Exception("Please input an integer or string 'all'")
        fig = plt.figure()
        label_map = pixel_clustering(sub_mask, features, cluster)
        plt.matshow(label_map, fignum=0)
        plt.xlabel(f"Cluster pattern (K={cluster})", fontsize=15)

        # 保存图像
        plt.savefig(os.path.join(save_dir, f'cluster_{cluster}.png'))
        plt.close()

        return fig


    else:  # generate cluster pattern with multiple resolutions, together with whole lung CT
        max_cluster = 9
        for clu in range(3, max_cluster + 1):
            plt.figure()
            label_map = pixel_clustering(sub_mask, features, clu)
            plt.matshow(label_map)
            plt.xlabel(str(clu) + ' clusters', fontsize=15)

            # 保存图像
            plt.savefig(os.path.join(save_dir, f'cluster_{clu}.png'))
            plt.close()

def calITHscore(label_map, min_area=200, thresh=1):
    """
    Calculate ITHscore from clustering label map
    Args:
        label_map: Numpy array. Clustering label map
        min_area: Int. For tumor area (pixels) smaller than "min_area", we don't consider connected-region smaller than "thresh"
        thresh: Int. The threshold for connected-region's area, only valid when tumor area < min_area
    Returns:
        ith_score: Float. The level of ITH, between 0 and 1
    """
    size = np.sum(label_map > 0)  # Record the number of total pixels
    num_regions_list = []
    max_area_list = []
    for i in np.unique(label_map)[1:]:  # For each gray level except 0 (background)
        flag = 1  # Flag to count this gray level, in case this gray level has only one pixel
        # Find (8-) connected-components. "num_regions" is the number of connected components
        labeled, num_regions = ndimage.label(label_map==i, structure=ndimage.generate_binary_structure(2,2))
        max_area = 0
        for j in np.unique(labeled)[1:]:  # 0 is background (here is all the other regions)
            # Ignore the region with only 1 or "thresh" px
            if size <= min_area:
                if np.sum(labeled == j) <= thresh:
                    num_regions -= 1
                    if num_regions == 0:  # In case there is only one region
                        flag = 0
                else:
                    temp_area = np.sum(labeled == j)
                    if temp_area > max_area:
                        max_area = temp_area
            else:
                if np.sum(labeled == j) <= 1:
                    num_regions -= 1
                    if num_regions == 0:  # In case there is only one region
                        flag = 0
                else:
                    temp_area = np.sum(labeled == j)
                    if temp_area > max_area:
                        max_area = temp_area
        if flag == 1:
            num_regions_list.append(num_regions)
            max_area_list.append(max_area)
    # Calculate the ITH score
    ith_score = 0
    # print(num_regions_list)
    for k in range(len(num_regions_list)):
        ith_score += float(max_area_list[k]) / num_regions_list[k]
    # Normalize each area with total size
    ith_score = ith_score / size
    ith_score = 1 - ith_score

    return ith_score



def batch_process(image_path, mask_path, result_path, case_path, label_path, limit=None):
    """
    Batch processes images and masks, saves visualizations and ITH scores.
    Retrieves intersected IDs from case_path and label_path.

    Args:
        image_path: List of paths to images.
        mask_path: List of paths to masks.
        result_path: Path to save results.
        case_path: Path to case files.
        label_path: Path to label files.
        limit: Optional; limits the number of processed cases.
    """
    # Get intersected IDs
    intersected_ids = cat_intersection_ids(case_path, label_path)

    if limit:
        image_path = image_path[:limit]
        mask_path = mask_path[:limit]
        intersected_ids = intersected_ids[:limit]

    # Ensure result directory exists
    if not os.path.exists(result_path):
        os.makedirs(result_path)

    result_fig_dir = os.path.join(result_path, "result_fig")
    if not os.path.exists(result_fig_dir):
        os.makedirs(result_fig_dir)

    ith_scores = []

    # Only process up to the 'limit' number of cases
    for img_path, msk_path, patient_id in tqdm(zip(image_path, mask_path, intersected_ids), total=min(len(image_path), limit) if limit else len(image_path), desc="Processing"):
        img = load_image(img_path)
        msk = load_seg(msk_path)

        # Get the largest slice
        largest_img, largest_msk = get_largest_slice(img, msk)

        # Locate tumor in the slice
        sub_img, sub_msk = locate_tumor(largest_img, largest_msk)

        # Extract radiomic features
        features = extract_radiomic_features(sub_img, sub_msk)

        # Perform clustering
        label_map = pixel_clustering(sub_msk, features)

        # Calculate ITH score
        ith_score = calITHscore(label_map)
        ith_scores.append((patient_id, ith_score))

        # Create patient-specific subfolder for visualizations
        patient_fig_dir = os.path.join(result_fig_dir, patient_id)
        if not os.path.exists(patient_fig_dir):
            os.makedirs(patient_fig_dir)

        # Visualization
        # visualize(largest_img, sub_img, largest_msk, sub_msk, features, cluster="all", save_dir=patient_fig_dir)
        visualize(largest_img, sub_img, largest_msk, sub_msk, features, cluster=6, save_dir=patient_fig_dir)

    # Save ITH scores to a CSV file
    ith_score_path = os.path.join(result_path, 'ith_scores.csv')
    with open(ith_score_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['ID', 'ITH_Score'])
        for id, score in ith_scores:
            writer.writerow([id, score])



# Preprocessing steps
case_path = 'D:/projects/ITHscore/01_TCGA_cases'
label_path = 'D:/projects/ITHscore/01_TCGA_labels'
intersected_ids = cat_intersection_ids(case_path, label_path)
processing_interactions(intersected_ids, label_path, case_path)

# Get image and mask path
image_path, mask_path = achieve_img_and_mask_path(case_path)


# Run batch processing
result_path = 'D:/projects/ITHscore/01_TCGA_results'
batch_process(image_path, mask_path, result_path, case_path, label_path, limit=None)
# batch_process(image_path, mask_path, result_path, limit=5)  # 只处理前5个图像




# # 加载图像和掩码
# image = load_image(image_path)
# mask = load_seg(mask_path)
#
# # 获取具有最大肿瘤面积的切片
# largest_slice_img, largest_slice_mask = get_largest_slice(image, mask)
#
# # 定位肿瘤并提取
# sub_img, sub_mask = locate_tumor(largest_slice_img, largest_slice_mask)
#
# # 提取放射组学特征
# import csv
# import numpy as np
# from sklearn.preprocessing import MinMaxScaler
# import functools
# import concurrent.futures
#
# features = extract_radiomic_features(sub_img, sub_mask, parallel=False)
#
# # Print extracted features
# for feature_name, feature_values in features.items():
#     print(f"{feature_name}: {feature_values}")
#
# # Save extracted features to a CSV file
# with open('extracted_features.csv', 'w', newline='') as csvfile:
#     feature_writer = csv.writer(csvfile)
#
#     # Writing header
#     headers = ['Feature_Type', 'Feature_Values']
#     feature_writer.writerow(headers)
#
#     # Writing feature values
#     for feature_name, feature_values in features.items():
#         # Feature values are expected to be in list-like structures
#         feature_writer.writerow([feature_name] + list(feature_values))
#
# # 进行像素聚类
# label_map = pixel_clustering(sub_mask, features, cluster=6)
#
# # 可视化
# visualize(image, sub_img, mask, sub_mask, features, cluster="all")
#
# plt.show()
# # 计算ITH分数
# ith_score = calITHscore(label_map)
# print("ITH Score:", ith_score)

