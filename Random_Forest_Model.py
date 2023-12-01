from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report
from sklearn.metrics import precision_recall_curve, auc
import numpy as np
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt


'''Reading in the csv.files'''

# creating pandas dataframe
dataset = pd.read_csv("DF_With_Distances.csv")
dataset = dataset.dropna(axis=0, how='any')

final_dataset = pd.read_csv("it_worked.csv")
final_dataset = final_dataset.dropna(axis=0, how='any')

patient_3 = pd.read_csv("Patient_3.csv")
patient_3 = patient_3.dropna(axis=0, how='any')




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

    # Ensure X_train and X_test have the same column order
    X_train = X_train[X_test.columns]

    # Sort features by importance
    sorted_indices = importances.argsort()[::-1]
    sorted_importances = importances[sorted_indices]

    # Plot feature importances
    plt.figure(figsize=(15, 15))
    plt.title("Feature Importance's RF all features", fontsize=27) ###CHANGE WHEN CHANGING SINGLE PATIENT TO FULL RUN OR REVERSE!

    # Directly specify labels for each bar
    plt.bar(range(X_train.shape[1]), sorted_importances, align="center",
            tick_label=[X_train.columns[i] for i in sorted_indices])
    plt.xticks(rotation=45, fontsize=8, ha="right")  # Adjust rotation and alignment
    plt.xlabel('Feature', fontsize=27)  # Set your desired font size
    plt.ylabel('Importance', fontsize=27)  # Set your desired font size
    # Adjust font size for x and y axis tick labels
    plt.xticks(fontsize=14)  # Set your desired font size
    plt.yticks(fontsize=23)  # Set your desired font size

    plt.tight_layout()
    plt.savefig("Feature importance RF all features.png", dpi = 300)
    plt.show()


    # Precision-Recall Curve
    precision, recall, _ = precision_recall_curve(Y_test, clf.predict_proba(X_test)[:, 1])
    pr_auc = auc(recall, precision)

    plt.figure(figsize=(8, 6))
    plt.step(recall, precision, color='b', alpha=0.2, where='post')
    plt.fill_between(recall, precision, step='post', alpha=0.15, color='b')

    # Adjust font size for title, x-axis label, and y-axis label
    plt.xlabel('Recall', fontsize=14)  # Set your desired font size
    plt.ylabel('Precision', fontsize=14)  # Set your desired font size
    plt.title(f'Precision-Recall Curve RF 4 features (AUC = {pr_auc:.2f})', fontsize=16)  # Set your desired font size  ###CHANGE WHEN CHANGING SINGLE PATIENT TO FULL RUN OR REVERSE!

    # Adjust font size for x and y axis tick labels
    plt.xticks(fontsize=12)  # Set your desired font size
    plt.yticks(fontsize=12)  # Set your desired font size

    plt.show()
    # Return the accuracy and classification report









if __name__ == "__main__":

    Single_patient = False
    Full_run = True
    ''' Test run with up and downstream methylation values'''
    beta_values = dataset["B_Val"]
    upstream_methylation = dataset['Upstream_methylation']
    downstream_methylation = dataset['Downstream_methylation']
    upstream_distance = dataset['Upstream_distance']
    downstream_distance = dataset['Downstream_distance']

    # Define features and target variables
    X = dataset[['Upstream_methylation', 'Downstream_methylation', 'Upstream_distance', 'Downstream_distance']]
    y_upstream = dataset['B_Val'].apply(lambda x: 0 if x < 0.5 else 1)  # Y MOET ALLEEN B VALUES ZIJN

    if Single_patient == True:
        random_forest(X, y_upstream)

    '''Actual full data run'''
    # final dataset
    X_final_columns = final_dataset.columns.tolist()
    X_final_columns = X_final_columns[5:]
    X_final_dataset = final_dataset[X_final_columns]
    Y_final_dataset= final_dataset["Patient_1"]

    if Full_run == True:
        random_forest(X_final_dataset,Y_final_dataset)























