import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, roc_auc_score, precision_recall_curve
import os

# 设定输出的目录
output_dir = "0730_test_03_all_100tree_withoutrmd_mutipose/"
os.makedirs(output_dir, exist_ok=True)

eval_results = pd.read_csv('/home/s2331261/Master_Project/9_Model_evaluation/0729_test_01_100tree_withoutrmd/predictions.csv')

eval_results ['PoseScoreRank'] = eval_results .groupby('PDBCode')['Probability'].rank(method='dense', ascending=False)

# this subsets the dataframe to only the highest scoring pose per pdb complex
best_df  = eval_results .loc[eval_results.PoseScoreRank == 1]

print(best_df)

# 获取标签列和预测得分列
y_true = best_df ['Label']
y_scores = best_df ['Probability']

# 计算ROC曲线和AUC
fpr, tpr, _ = roc_curve(y_true, y_scores)
roc_auc = roc_auc_score(y_true, y_scores)

# 计算PR曲线和AUCPR
precision, recall, _ = precision_recall_curve(y_true, y_scores)
pr_auc = auc(recall, precision)

# 画图
plt.figure(figsize=(15, 5))

plt.subplot(1, 2, 1)
plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")

plt.subplot(1, 2, 2)
plt.plot(recall, precision, color='blue', lw=2, label='PR curve (area = %0.2f)' % pr_auc)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.legend(loc="lower right")

# 保存图表到指定的目录
plt.savefig(os.path.join(output_dir, 'evaluation_curves.png'))
