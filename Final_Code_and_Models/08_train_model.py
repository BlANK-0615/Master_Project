# This script is designed to optimize and train 
# an XGBoost model using Bayesian optimization and then evaluate its performance.

import pandas as pd
import os
import pickle
from sklearn.model_selection import KFold
from skopt import BayesSearchCV
from xgboost import XGBClassifier
from skopt.space import Real, Integer, Categorical
from joblib import dump, load
from sklearn.metrics import roc_auc_score, accuracy_score, average_precision_score
from sklearn.metrics import precision_recall_curve, roc_curve
import matplotlib.pyplot as plt
import timeit
from skopt.callbacks import VerboseCallback
from skopt.utils import point_asdict
from sklearn.metrics import make_scorer, average_precision_score
from sklearn.model_selection import GroupKFold

start=timeit.default_timer()
output_dir = 'path/to/output/dir'
os.makedirs(output_dir, exist_ok=True)

# define the hyperparameter search space
search_space = {
    "colsample_bylevel": Real(0.3, 1),
    "colsample_bytree": Real(0.3, 1),
    "gamma": Real(0.01, 1),
    "eta": Real(0.0001, 1),
    "max_delta_step": Integer(3, 15),
    "min_child_weight": Real(1, 8),
    "reg_alpha": Real(0.1, 100),
    "reg_lambda": Real(0.1, 100),
    "subsample": Real(0.3, 1.0),
    "scale_pos_weight": Integer(5, 15),
    "objective": Categorical(['binary:logistic'])
}


# 1.2 Create AUCPR scorer
aucpr_scorer = make_scorer(average_precision_score, needs_proba=True)

#1.3 print the optimization progress
def print_progress(res):
    """Print the parameters and target values of the current iteration"""
    # Get the parameters of the current iteration
    params = point_asdict(search_space, res.x_iters[-1])
   # Get the target value of the current iteration
    target = res.fun
    # print parameters and target values
    print(f"Iteration {len(res.x_iters)}")
    print(f"Target: {target:.4f}")
    for k, v in params.items():
        print(f"{k}: {v}")
    print()
    # Define the string to write
    output_str = ""
    output_str += f"Iteration {len(res.x_iters)}\n"
    output_str += f"Target: {target:.4f}\n"
    for k, v in params.items():
        output_str += f"{k}: {v}\n"
    output_str += "\n"
    # open the file and write
    with open(os.path.join(output_dir, 'progress.txt'), 'a') as f:
        f.write(output_str)

# 2.Load the training and validation data
train_df = pd.read_csv('path to train data')
val_df = pd.read_csv('path to validation data')

y_train = train_df['Label']
X_train = train_df.drop(['PDBCode', 'Pose_name', 'Label'], axis=1)

y_val = val_df['Label']
X_val = val_df.drop(['PDBCode', 'Pose_name', 'Label'], axis=1)

# 3.Define Bayesian Search CV
opt = BayesSearchCV(
    XGBClassifier(use_label_encoder=False, eval_metric='aucpr', tree_method='gpu_hist',learning_rate=0.1),
    #XGBClassifier(use_label_encoder=False, eval_metric='aucpr',learning_rate=0.1),
    search_space,
    scoring=aucpr_scorer,
    #cv=KFold(n_splits=5, shuffle=True, random_state=42),
    cv=GroupKFold(n_splits=5),
    n_jobs=20,
    n_iter=100,
    return_train_score=True,
    refit=True,
    optimizer_kwargs={'base_estimator': 'GP'},
    random_state=42,
    verbose=1,
)

# 4.Parameter optimization using the training set
print("########################## Search HyperParameters ##########################")
opt.fit(X_train, y_train, callback=[print_progress], groups=train_df['PDBCode']) 


# 5.Train the model on the training and validation sets with the optimal parameters
best_params = opt.best_params_
best_params.update({
    'use_label_encoder': False,
    'tree_method': 'gpu_hist',  # using GPU to train
    'verbosity': 2,
    'n_estimators': 10000,
    'learning_rate': 0.02,
    'random_state': 42
})
with open(os.path.join(output_dir, 'best_params.txt'), 'w') as file:
    file.write(str(best_params))
with open(os.path.join(output_dir, "best_params.pkl"), "wb") as f:
    pickle.dump(best_params, f)


print("########################## Begin to train model ##########################")
xgb = XGBClassifier(**best_params)
progress_dictionary = dict()  # dict to store the evaluation results of each round
xgb.fit(X_train, y_train,
        early_stopping_rounds=100,  # change to your desired early stopping rounds
        eval_set=[(X_train, y_train), (X_val, y_val)],
        eval_metric=['aucpr', 'auc'],
        verbose=True
        )  
progress_dictionary = xgb.evals_result() # store evaluation results

# 6. Save the trained model as a pkl file
dump(xgb, os.path.join(output_dir, 'model.pkl'))

# Save the evaluation results to a CSV file
train_AUCPRs = progress_dictionary['validation_0']['aucpr']
val_AUCPRs = progress_dictionary['validation_1']['aucpr']
train_AUROCs = progress_dictionary['validation_0']['auc']
val_AUROCs = progress_dictionary['validation_1']['auc']
df = pd.DataFrame({
    'train-AUCPR': train_AUCPRs, 
    'val-AUCPR': val_AUCPRs, 
    'train-AUROC': train_AUROCs, 
    'val-AUROC': val_AUROCs
})
# save the evaluation results to a CSV file
df.to_csv(os.path.join(output_dir, 'Model_Training_History.csv'), index=False)
end=timeit.default_timer()
print('Running time: %s Seconds'%(end-start))
