from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report
import numpy as np
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt


'''Reading in the csv.files'''

# creating pandas dataframe
dataset = pd.read_csv("DF_With_Distances.csv")
dataset = dataset.dropna(axis=0, how='any')

final_dataset = pd.read_csv("all_features.csv")
final_dataset = final_dataset.dropna(axis=0, how='any')





def random_forest(X,Y, test_s = 0.4):
    '''Input: X = target variable. Y = Prediction variables, percentage_test (standard 0.4),
    Output: Classification report, Feature Importance graph
     '''

    # Split the data into training and testing sets
    X_train, X_test, Y_train, Y_test = train_test_split(
        X, Y, test_size=test_s, random_state=60)

    # Create separate Random Forest classifiers for upstream and downstream
    clf = RandomForestClassifier(n_estimators=300, random_state=60)

    # Train the classifiers
    clf.fit(X_train, Y_train)

    # Make predictions for both upstream and downstream
    Y_pred = clf.predict(X_test)

    # Evaluate the models (you can calculate accuracy, precision, recall, etc.)
    accuracy_upstream = accuracy_score(Y_test, Y_pred)

    # Print the results for both upstream and downstream
    print("Accuracy:", accuracy_upstream)

    # Print the results
    print("Classification Report:")
    print(classification_report(Y_test, Y_pred))
    print("Classification Report:")

    # Get feature importances
    importances = clf.feature_importances_
    feature_names = X_train.columns

    # Sort features by importance
    sorted_indices = importances.argsort()[::-1]

    # Plot feature importances
    plt.figure(figsize=(10, 6))
    plt.title("Feature Importances")
    plt.bar(range(X_train.shape[1]), importances[sorted_indices], align="center")
    plt.xticks(range(X_train.shape[1]), [feature_names[i] for i in sorted_indices], rotation=45)
    plt.xlabel("Feature")
    plt.ylabel("Importance")
    plt.show()

    # Return the accuracy and classification report




''' Test run with up and downstream emthylation values'''
beta_values = dataset["B_Val"]
upstream_methylation = dataset['Upstream_methylation']
downstream_methylation = dataset['Downstream_methylation']
upstream_distance = dataset['Upstream_distance']
downstream_distance = dataset['Downstream_distance']

# Define features and target variables
X = dataset[['Upstream_methylation', 'Downstream_methylation', 'Upstream_distance', 'Downstream_distance']]
y_upstream = dataset['B_Val'].apply(lambda x: 0 if x < 0.5 else 1)  # Y MOET ALLEEN B VALUES ZIJN
random_forest(X, y_upstream)


# final dataset
b_val_final = final_dataset["Beta"]




















