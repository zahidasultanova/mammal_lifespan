#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd 
import numpy as np
from catboost import CatBoostRegressor
from sklearn.model_selection import RepeatedKFold, RandomizedSearchCV
from scipy.stats import randint, uniform


os.chdir('/gpfs/home/djs19ctu/ML_project/anage_all_data_nogen')
os.getcwd()


contents = os.listdir(os.getcwd())

print("Contents of the current working directory:")
for item in contents:
    print(item)


train_df_imputed = pd.read_csv('train_mammal_data_imputed_nogen.csv')
test_df_imputed = pd.read_csv('test_mammal_data_imputed_nogen.csv')

try:
    import sklearn
    print(f"scikit-learn version: {sklearn.__version__}")
    print("scikit-learn imported successfully!")
except ImportError as e:
    print("Error importing scikit-learn:", e)

try:
    import catboost
    print(f"CatBoost version: {catboost.__version__}")
    print("CatBoost imported successfully!")
except ImportError as e:
    print("Error importing CatBoost:", e)
    print("To install CatBoost, run: pip install catboost")
    

from catboost import CatBoostRegressor
from sklearn.model_selection import RepeatedKFold, RandomizedSearchCV
from scipy.stats import randint, uniform

# Define the features and target
features = ['order', 'family', 'genus',
            'adult_mass_g', 'adult_brain_mass_g',
            'female_maturity_d', 'gestation_length_d',
            'litter_size_n', 'litters_per_year_n',
            'weaning_age_d',
            'hibernation_torpor',
            'trophic_level', 'activity_cycle',
            'freshwater', 'marine', 'terrestrial_non-volant', 'terrestrial_volant', 'habitat_breadth_n',
           'specimen origin']

target = 'maximum longevity (yrs)'

# Assuming train_df_imputed is defined elsewhere in your code
data_train = train_df_imputed.copy()

X = data_train[features].copy()  # Explicitly make a copy
y = data_train[target]

# Convert categorical features to category type
categorical_features = ['order', 'family', 'genus', 'specimen origin']
for col in categorical_features:
    X.loc[:, col] = X[col].astype('category')  # Use .loc to avoid SettingWithCopyWarning

# Define the model
model = CatBoostRegressor(thread_count=4, verbose=0)  # Set thread_count to 4

# Define the parameter grid for randomized search
param_dist = {
    'iterations': randint(100, 1000),
    'learning_rate': uniform(0.01, 0.3),
    'depth': randint(3, 12),
    'l2_leaf_reg': uniform(0, 20),
    'cat_features': [categorical_features]
}


# Define the repeated k-fold cross-validator
rkf = RepeatedKFold(n_splits=5, n_repeats=10, random_state=42)

# Perform randomized search with cross-validation
random_search = RandomizedSearchCV(model, param_distributions=param_dist, n_iter=100, cv=rkf, scoring='r2', random_state=42, n_jobs=-1)
random_search.fit(X, y)

# Get the best model
catboost_best_model = random_search.best_estimator_
catboost_params = random_search.best_params_

# Save the results to a text file
with open("catboost_param_results_withQuality.txt", "w") as file:
    file.write(f"Best Model: {catboost_best_model}\n")
    file.write(f"Best Parameters: {catboost_params}\n")