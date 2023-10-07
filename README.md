# **Prediction of protein-ligand interactions via modern machine learning methods to accelerate medicine discovery**

## The related scripts, models and environment for dissertation project of B224911

### Absctract of this project
In this study, TocoDecoy (a decoy generation method) is utilized as an alternative to DeepCoy to further reduce dataset bias 
and retrain the XGBoost model in SCORCH in an attempt to improve performance.

The results showed that the decoys created by the new TocoDecoy did not enhance the original SCORCH model due to limited chemical diversity.
Furthermore, TocoDecoy alone is not suitable for model training, but its combination with various decoy sources reveals the potential for model improvement. 

Furthermore, models retrained using the combined TocoDecoy and DeepCoy decoys with tuned hyperparameters show improved performance compared to SCORCH. However, all models exhibit insufficient generalization ability on independent benchmark datasets. The findings highlight the critical importance of the chemical diversity and breadth of the bait dataset for MLSF performance.

At the same time, the hyperparameter optimization step in the training model will also affect the performance, which needs to be paid attention in the follow-up research.

Although this study did not significantly improve the performance of SCORCH's XGBoost model, it provided valuable insights into factors influencing MLSF construction. Continuous optimization of 
data representation and training procedures will further enhance MLSF, advance structure-based virtual screening, and accelerate drug discovery.

### The scripts and modles in repository

- **Models_file** : including the total nine re-trained XGBoost models in this project. The first six models are models trained only with TocoDecoy, and the last three are models trained after combining TocoDecoy and adjusting hyperparameter search. Compared with the first six models, the performance of the last three models has been greatly improved.

- **utils**：including the scripts for feature extraction and the original XGBoost model from SCORCH1.0 related work(https://github.com/SMVDGroup/SCORCH)

- **Master_project.ymal**: incluing the conda environment setting of project.

- **All scripts**: including the key scripts for project. For details, please see the comments in the script. Some other scripts for grid filter can be accessed in TocoDecoy work(https://github.com/5AGE-zhang/TocoDecoy).

- **Notice**: Raw data are not presented as they are unpublished work belonging to the Houston lab.




