import os
import json
import pickle
import pandas as pd
import xgboost as xgb
from joblib import dump, load
from skopt import BayesSearchCV
from xgboost import XGBClassifier
from functools import reduce
from sklearn.model_selection import KFold
from skopt.space import Real, Integer, Categorical
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_auc_score, accuracy_score, average_precision_score
from sklearn.metrics import precision_recall_curve, roc_curve
import matplotlib.pyplot as plt
import numpy as np
from functools import partialmethod
from warnings import filterwarnings
from sklearn import metrics
import matplotlib.pyplot as plt


# 设置scorch的输出目录
def scale_multipose_features(df):

    """
    Function: Scale features using prefitted feature scaler

    Input:    Dataframe of features to scale

    Output:   Scaled dataframe of features
    """

    # open the reference features file
    reference_headers = json.load(open(os.path.join('utils','params','features.json')))
    scaler_58 = reference_headers.get('for_scaler_58')
    headers_58 = reference_headers.get('492_models_58')

    # store the ligand and receptor information
    PDB_code, Pose, Label = df['PDBCode'], df['Pose_name'], df['Label']


    # get the missing columns that the scaler expects
    missing_columns = list(set(scaler_58) - set(list(df)))

    # fill in dummy columns for the scaler
    for col in missing_columns:
        df[col] = 0
    df = df[scaler_58]

    # load the scaler and scale the data
    scaler = load(os.path.join('utils','params','58_maxabs_scaler_params.save'))
    scaled = scaler.transform(df)
    df[df.columns] = scaled

    # then only keep the columns we want for the models
    df = df[headers_58]

    # add the receptor and ligand info back in
    df['PDBCode'], df['Pose_name'], df['Label'] =  PDB_code, Pose, Label

    # return the dataframe
    return df

# 设置输出工作目录
file_dir = '0729_test_01_100tree_withoutrmd/'
os.makedirs(file_dir, exist_ok=True)

#1.载入test集
test_df_50 = pd.read_csv('/home/s2331261/Master_Project/8_Model_trian/A1_model_dataset_50/01_ALL/all_3v_test_910.csv')
test_df_100 = pd.read_csv('/home/s2331261/Master_Project/8_Model_trian/A2_model_dataset_100tree/01_ALL/all_3v_test_1166.csv')
test_df_scorch = pd.read_csv('/home/s2331261/Master_Project/7_Feature_Engineer/A1_features_files/05_split_feature/all_dataset/all_TEST.csv')

#2.处理scorch要使用的test集
X_test = test_df_score.drop(['PDBCode', 'Pose_name', 'Label'], axis=1)
y_test = test_df_score.iloc[:,-1]



# 加载模型并在测试集上进行评估
loaded_model = load('/home/s2331261/Master_Project/8_Model_trian/0729_train_01/model.pkl')
# Evaluate on test set
y_pred = loaded_model.predict(X_test)
y_prob = loaded_model.predict_proba(X_test)[:, 1]

print("Test Accuracy: ", accuracy_score(y_test, y_pred))
print("Test AUC: ", roc_auc_score(y_test, y_prob))
print("Test AUCPR: ", average_precision_score(y_test, y_prob))

# 计算模型在测试集上的性能指标
accuracy = accuracy_score(y_test, y_pred)
roc_auc = roc_auc_score(y_test, y_prob)
pr_auc = average_precision_score(y_test, y_prob)

# 获取特征重要性
feature_importances = loaded_model.feature_importances_

# 创建一个DataFrame来存储特征重要性
feature_importances_df = pd.DataFrame({'Feature': X_test.columns, 'Importance': feature_importances})

# 计算各个性能指标，并存储到DataFrame中
metrics_df = pd.DataFrame({
    'Metric': ['Accuracy', 'AUC ROC', 'AUC PR'],
    'Score': [accuracy, roc_auc, pr_auc]
})
# 新建一个DataFrame来存储'PDBCode', 'Pose_name', 'Label'以及预测结果和预测概率
results_df = test_df[['PDBCode', 'Pose_name', 'Label']].copy()
results_df['Prediction'] = y_pred
results_df['Probability'] = y_prob



# 保存两个DataFrame到CSV文件
results_df.to_csv(os.path.join(file_dir, 'predictions.csv'), index=False)
feature_importances_df.to_csv(os.path.join(file_dir, 'feature_importances.csv'), index=False)
metrics_df.to_csv(os.path.join(file_dir, 'model_metrics.csv'), index=False)

# 绘制性能指标的条形图
plt.bar(metrics_df['Metric'], metrics_df['Score'])
plt.ylabel('Score')
plt.title('Performance Metrics')
plt.savefig(os.path.join(file_dir, 'performance_metrics.png'))
plt.show()

# 计算并绘制ROC曲线和PR曲线
fpr, tpr, _ = roc_curve(y_test, y_prob)
plt.figure()
plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")
plt.savefig(os.path.join(file_dir, 'roc_curve.png'))
plt.show()

precision, recall, _ = precision_recall_curve(y_test, y_prob)
plt.figure()
plt.plot(recall, precision, label='PR curve (area = %0.2f)' % pr_auc)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall curve')
plt.legend(loc="lower right")
plt.savefig(os.path.join(file_dir, 'pr_curve.png'))
plt.show()



# # 原来的XGBoost
# dtest = xgb.DMatrix(X_test)
# loaded_model = load('/home/s2331261/Master_Project/8_Model_trian/495_models_58_booster.pkl')
# X_test = X_test[loaded_model.feature_names]

# y_pred = loaded_model.predict(dtest, output_margin=True)
# y_prob = 1.0 / (1.0 + np.exp(-y_pred))  # convert margin to probability


# accuracy = accuracy_score(y_test, (y_prob > 0.5).astype(int))
# roc_auc = roc_auc_score(y_test, y_prob)
# pr_auc = average_precision_score(y_test, y_prob)
# feature_importances = loaded_model.get_score(importance_type='gain')