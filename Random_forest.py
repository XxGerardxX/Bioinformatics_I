from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report
import numpy as np
import pandas as pd

'''Reading in the csv.files'''

# creating pandas dataframe
dataset = pd.read_csv("DF_With_Distances.csv")
dataset = dataset.dropna(axis=0, how='any')

dataset = dataset[['Upstream_distance', 'Downstream_distance', 'Upstream_methylation',
                   'Downstream_methylation']]

upstream_methylation = dataset['Upstream_distance']
downstream_methylation = dataset['Downstream_methylation']
upstream_distance = dataset['Upstream_distance']
downstream_distance = dataset['Downstream_distance']


# setting up the model
n_samples = len(upstream_distance)
target = np.random.randint(2, size=n_samples)

# Split the data into training and testing sets
X = dataset[['Upstream_distance', 'Downstream_distance', 'Upstream_methylation',
                   'Downstream_methylation']]
y = target
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4, random_state=60)

# Create a Random Forest classifier
clf = RandomForestClassifier(n_estimators=300, random_state=60)

# Train the classifier on the training data
clf.fit(X_train, y_train)

# Make predictions on the test data
y_pred = clf.predict(X_test)

# Evaluate the model's performance
accuracy = accuracy_score(y_test, y_pred)
classification_report_str = classification_report(y_test, y_pred)

# Print the results
print(f'Accuracy: {accuracy}')
print('Classification Report:\n', classification_report_str)
