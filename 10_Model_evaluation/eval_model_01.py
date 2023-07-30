import pandas as pd
from sklearn.model_selection import KFold
from skopt import BayesSearchCV
from xgboost import XGBClassifier
import xgboost as xgb
from skopt.space import Real, Integer, Categorical
from joblib import dump, load
from sklearn.metrics import roc_auc_score, accuracy_score, average_precision_score
from sklearn.metrics import precision_recall_curve, roc_curve
import matplotlib.pyplot as plt
import numpy as np
import os


test_df = pd.read_csv('/home/s2331261/Master_Project/7_Feature_Engineer/A3_model_dataset/01_ALL/all_3v_test_910.csv')
X_test = test_df.drop(['PDBCode', 'Pose_name', 'Label'], axis=1)
y_test = test_df.iloc[:,-1]
file_dir = '0729_test_01_100tree_withoutrmd/'
os.makedirs(file_dir, exist_ok=True)


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