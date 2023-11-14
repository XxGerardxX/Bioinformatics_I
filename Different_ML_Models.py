from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import MinMaxScaler
from sklearn.svm import SVC, SVR
from sklearn.metrics import accuracy_score, classification_report, mean_squared_error
import numpy as np
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import New_Raw_Data_Sort as NRDS

final_dataset = pd.read_csv("it_worked.csv")
final_dataset = final_dataset.dropna(axis=0, how='any')
non_binary_data_p_1 = pd.read_csv("test_non_binary.csv")

patient_3 = pd.read_csv("Patient_3.csv")
patient_3 = patient_3.dropna(axis=0, how='any')




def baseline(Y_val_non_binary, bp_length_comp = 1500):
    #TODO: import methylation without conversion values
    #TODO: convert to 1 and 0 depending on mean calculation of 1500bp left and right
    # comparing percentage correct to actual values from imported data in function
    Y_values = np.array(Y_val_non_binary.sort_values())
    comparison_df = {"Non_binary": Y_values, "Binary": None, "WinCPG": None}

    Binary_list =[]
    WinCPG_list = []
    middle_of_list = Y_values[Y_values // 2]

    for index,value in enumerate(Y_values):

        # calculating binary values and appending them to a list
        Binary_val = value
        if Binary_val >= 0.5:
            Binary_val = 1
        else:
            Binary_val = 0
        Binary_list.append(Binary_val)

        #calculating upper and lower limit and finding indices

        WinCPG_val = value
        upper_limit = WinCPG_val + 1500
        lower_limit = WinCPG_val - 1500
        if lower_limit < middle_of_list:
                












    return

def random_forest(X,Y, test_s = 0.4):
    '''Input: X = target variable. Y = Prediction variables, percentage_test (standard 0.4),
    Output: Classification report, Feature Importance graph
     '''

    # Split the data into training and testing sets
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_s, random_state=60)

    # Create separate Random Forest classifiers for upstream and downstream
    clf = RandomForestClassifier(n_estimators=300, random_state=60)

    # Train the classifiers
    clf.fit(X_train, Y_train)

    # Make predictions for both upstream and downstream
    Y_pred = clf.predict(X_test)

    # Evaluate the models (you can calculate accuracy, precision, recall, etc.)
    accuracy = accuracy_score(Y_test, Y_pred)

    # Print the results for both upstream and downstream
    print("Accuracy:", accuracy)

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
    plt.title("Feature Importance's")
    plt.bar(range(X_train.shape[1]), importances[sorted_indices], align="center")
    plt.xticks(range(X_train.shape[1]), [feature_names[i] for i in sorted_indices], rotation=45)
    plt.xlabel("Feature")
    plt.ylabel("Importance")
    plt.show()

    # Return the accuracy and classification report
    return accuracy

def random_forest_cross_sample(X_train, Y_train, X_test, Y_test):
    '''Input: X_train, Y_train for training; X_test, Y_test for testing.
       Output: Classification report, Feature Importance graph
        '''

    # Create a Random Forest classifier
    clf = RandomForestClassifier(n_estimators=300, random_state=60)

    # Train the classifier on the training data
    clf.fit(X_train, Y_train)

    # Make predictions on the testing data
    Y_pred = clf.predict(X_test)

    # Evaluate the model (you can calculate accuracy, precision, recall, etc.)
    accuracy = accuracy_score(Y_test, Y_pred)

    # Print the results
    print("Accuracy:", accuracy)
    print("Classification Report:")
    print(classification_report(Y_test, Y_pred))

    # Get feature importances
    importances = clf.feature_importances_
    feature_names = X_train.columns

    # Sort features by importance
    sorted_indices = importances.argsort()[::-1]

    # Plot feature importances
    plt.figure(figsize=(10, 6))
    plt.title("Feature Importance's")
    plt.bar(range(X_train.shape[1]), importances[sorted_indices], align="center")
    plt.xticks(range(X_train.shape[1]), [feature_names[i] for i in sorted_indices], rotation=45)
    plt.xlabel("Feature")
    plt.ylabel("Importance")
    plt.show()

    # Return the accuracy and classification report
    return accuracy

def logistic_regression(X, Y, test_s=0.4):
    # Min-Max Scaling
    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(X)

    # Split the data into training and testing sets
    X_train, X_test, Y_train, Y_test = train_test_split(X_scaled, Y, test_size=test_s, random_state=60)

    # Create and train the logistic regression model
    logistic_regression_model = LogisticRegression(max_iter=1000)
    logistic_regression_model.fit(X_train, Y_train)

    # Make predictions on the testing data
    Y_pred = logistic_regression_model.predict(X_test)

    # Evaluate the model
    accuracy = accuracy_score(Y_test, Y_pred)
    print(f"Accuracy: {accuracy}")
    return accuracy

# def Support_Vector_Machines(X,Y,test_s=0.4):
#
#     # Split the data into training and testing sets
#     X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_s, random_state=60)
#
#     # Create an SVM model for classification
#     cls_model = SVC(kernel='linear', C=1.0)
#
#     # Train the classification model
#     cls_model.fit(X_train, Y_train)
#
#     # Make predictions for classification
#     Y_cls_pred = cls_model.predict(X_test)
#
#     # Evaluate classification model
#     cls_accuracy = accuracy_score(Y_test, Y_cls_pred)
#     print(f"Classification Accuracy: {cls_accuracy}")
#
#     # Create an SVM model for regression
#     reg_model = SVR(kernel='linear', C=1.0)
#
#     # Train the regression model
#     reg_model.fit(X_train, Y_train)
#
#     # Make predictions for regression
#     Y_reg_pred = reg_model.predict(X_test)
#
#     # Evaluate regression model
#     reg_mse = mean_squared_error(Y_test, Y_reg_pred)
#     print(f"Regression Mean Squared Error: {reg_mse}")





if __name__ == "__main__":
    #TODO: PCA and SVM?

    full_Data_run = False

    if full_Data_run == True:

        '''Baseline data creation'''

        # assay = pd.read_csv("assay.csv")
        # LIHC = pd.read_csv("brca.csv")
        # p_1_baseline_model = NRDS.dataframe_filter(assay,chromosome_filter="no", patient_filter="yes", b_values_binary="no", creating_csv="no")
        #
        # # 1) which patient and which chromosome?
        # ## merge LIHC and the filtered assay (csv was created to test)
        # merge_assay_columns = p_1_baseline_model.columns.tolist()
        # merged_LIHC_filtered_assay = NRDS.merge_df([LIHC,p_1_baseline_model],merge_assay_columns,0)
        #
        # # 3) filter on chromosome 1 (the csv was created to test)
        # final_df_chr_1 = NRDS.dataframe_filter(merged_LIHC_filtered_assay, chromosome_filter="yes") #creates unnamed column because of first column 48 features for this drop the unnamed column
        #
        # final_df_chr_1.to_csv("test_non_binary.csv")



        '''Actual full data run'''
        # final dataset
        X_final_columns = final_dataset.columns.tolist()
        X_final_columns = X_final_columns[5:]
        X_final_dataset = final_dataset[X_final_columns]
        Y_final_dataset= final_dataset["Patient_1"]

        '''Full data run cross-section'''
        X_train_set = X_final_dataset
        Y_train_set = Y_final_dataset

        X_test_set_columns = patient_3.columns.tolist()[5:]
        X_test_set = patient_3[X_test_set_columns]

        Y_test_set = patient_3["Patient_3"]



        #running random forest
        # random_forest(X_final_dataset,Y_final_dataset)

        # running random forest cross sample
        # random_forest_cross_sample(X_train_set,Y_train_set,X_test_set,Y_test_set)

        # running logistic regression
        # logistic_regression(X_final_dataset, Y_final_dataset)

        #SVM
        # Support_Vector_Machines(X_final_dataset,Y_final_dataset)


    #Baseline
    data_baseline_comp  =non_binary_data_p_1.iloc[:,4]
    baseline(data_baseline_comp)





