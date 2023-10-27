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


beta_values = dataset["B_Val"]
upstream_methylation = dataset['Upstream_methylation']
downstream_methylation = dataset['Downstream_methylation']
upstream_distance = dataset['Upstream_distance']
downstream_distance = dataset['Downstream_distance']

# Define features and target variables
X = dataset[['Upstream_methylation','Downstream_methylation', 'Upstream_distance', 'Downstream_distance']]
y_upstream=dataset['B_Val'].apply(lambda x: 0 if x < 0.5 else 1) # Y MOET ALLEEN B VALUES ZIJN
# y_downstream = dataset['Downstream_methylation']

# Split the data into training and testing sets
X_train, X_test, y_train_upstream, y_test_upstream = train_test_split(
    X, y_upstream, test_size=0.4, random_state=60)
# _, _, y_train_downstream, y_test_downstream = train_test_split(
#     X, y_downstream, test_size=0.4, random_state=60)



# Create separate Random Forest classifiers for upstream and downstream
clf_upstream = RandomForestClassifier(n_estimators=300, random_state=60)
# clf_downstream = RandomForestClassifier(n_estimators=300, random_state=60)

# Train the classifiers
clf_upstream.fit(X_train, y_train_upstream)
# clf_downstream.fit(X_train, y_train_downstream)

# Make predictions for both upstream and downstream
y_pred_upstream = clf_upstream.predict(X_test)
# y_pred_downstream = clf_downstream.predict(X_test)

# Evaluate the models (you can calculate accuracy, precision, recall, etc.)
accuracy_upstream = accuracy_score(y_test_upstream, y_pred_upstream)
# accuracy_downstream = accuracy_score(y_test_downstream, y_pred_downstream)

# Print the results for both upstream and downstream
print("Accuracy:", accuracy_upstream)
# print("Downstream Accuracy:", accuracy_downstream)



# Print the results
print("Classification Report:")
print(classification_report(y_test_upstream, y_pred_upstream))

print("Downstream Classification Report:")
# print(classification_report(y_test_downstream, y_pred_downstream))




# # Print the results
# print(f'Accuracy: {accuracy}')
# print('Classification Report:\n', classification_report_str)
#
#
# # Get feature importances
# importances = clf.feature_importances_
# feature_names = X_train.columns
#
# # Sort features by importance
# sorted_indices = importances.argsort()[::-1]
#
# # Plot feature importances
# plt.figure(figsize=(10, 6))
# plt.title("Feature Importances")
# plt.bar(range(X_train.shape[1]), importances[sorted_indices], align="center")
# plt.xticks(range(X_train.shape[1]), [feature_names[i] for i in sorted_indices], rotation=90)
# plt.xlabel("Feature")
# plt.ylabel("Importance")
# plt.show()











