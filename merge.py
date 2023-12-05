import pandas as pd
import os

def merge_csv_file(label_path, result_path, file_name):
    labels_file_path = os.path.join(label_path, '02_FJMU_label_intersection.csv')
    ith_score_file_path = os.path.join(result_path, 'ith_scores.csv')

    labels_df = pd.read_csv(labels_file_path, dtype={'ID': str})
    ith_score_df = pd.read_csv(ith_score_file_path, dtype={'ID': str})

    labels_df['ID'] = labels_df['ID'].str.strip()
    ith_score_df['ID'] = ith_score_df['ID'].str.strip()

    merged_df = labels_df.merge(ith_score_df, on='ID')

    column_order = ['ID', 'ITH_Score', 'OS', 'OS.time']
    merged_df = merged_df[column_order]

    ith_score_path = os.path.join(result_path, file_name)

    merged_df.to_csv(ith_score_path, index=False)

    return merge_csv_file

label_path = 'D:/projects/ITHscore/02_FJMU_labels'
result_path = 'D:/projects/ITHscore/02_FJMU_results'
file_name = 'ith_score_merge.csv'
merge_csv_file(label_path, result_path, file_name)