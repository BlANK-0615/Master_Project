# Uses sklearns cross-validated recursive feature elimination function
# with a random forest classifier estimator to determine optimum
# number of features to use based on the highest acheived ROC AUC score. 
# Original script from https://github.com/miles-mcgibbon/XGBScore, modified for this project


import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold
import sys
import pandas as pd
import numpy as np
import timeit
import os

# define the function to parse CLI inputs
def parse_args(args): # parse CLI user inputs
    h5_path = args[args.index('-h5') + 1]
    output_dir =  args[args.index('-out') + 1]
    return h5_path, output_dir

# define main function to perform feature selection
def main(): # run script using CLI
    h5_path, output_dir = parse_args(sys.argv)
    os.makedirs(output_dir, exist_ok=True)

    # load dataset
    df = pd.read_csv(h5_path)
    label_headers = ['PDBCode','Pose_name','Label']
    X = df.drop(label_headers , axis=1)
    headers = X.columns
    y = df['Label'].copy()

    # define model
    model = RandomForestClassifier(n_estimators=50, n_jobs=-1)
    # model = RandomForestClassifier(n_estimators=100, n_jobs=-1)

    # define recursive eliminator
    cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=0)
    selector = RFECV(model, step=1, cv=cv, n_jobs= -1, scoring='roc_auc', verbose=10)

    # run elimination
    selector = selector.fit(X, y)

    # fetch boolean list of features to include
    print('############################# Bool ###############################')
    print('selector.support_list:')
    print(selector.support_)
    choices = list(selector.support_)

    # fetch feature rankings
    print('##########################  Ranking  #############################')
    print('selector.ranking_list:')
    print(selector.ranking_,'\n')
    ranks = list(selector.ranking_)

    # return optimal feature counts
    print('####################  Optimal features  ########################')
    print("Optimal number of features : %d" % selector.n_features_)
    optimal_features = np.array(headers)[selector.support_]
    print("Optimal features:", optimal_features, '\n')

    # save importance, boolean and ranking results to dataframe
    results = pd.DataFrame({'Feature':headers, 'includeBool':choices, 'Rank':ranks})
    results.to_csv(f'{output_dir}RFECVFeatureRanks.csv', index=False)

    # save plot of aurocs to number of features
    print('##########################  Plot  ###############################')
    plt.figure()
    plt.xlabel("Number of features selected")
    plt.ylabel("Cross validation score (nb of correct classifications)")
    plt.plot(np.linspace(1, len(headers), len(selector.grid_scores_)),
             selector.grid_scores_)
    print('the number of grid_scores_:', len(selector.grid_scores_), '\n')
    print('selector.grid_scores_:')
    print(selector.grid_scores_, '\n')
    plt.savefig(f'{output_dir}featureROCAUCPlot.png')

    # save aurocs of feature counts to dataframe
    print('#####################  Saving results  ##########################')
    feature_numbers = np.linspace(1, len(headers), len(selector.grid_scores_))
    auc_scores = selector.grid_scores_
    auc_results = pd.DataFrame({'Num Features':feature_numbers, 'ROC_AUC':auc_scores})
    auc_results.to_csv(f'{output_dir}Ideal{selector.n_features_}ROCAUCResults.csv', index=False)


if __name__ == '__main__':
    start=timeit.default_timer()
    main()
    end=timeit.default_timer()
    print('Running time: %s Seconds'%(end-start))
