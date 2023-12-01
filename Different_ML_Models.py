from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_curve, auc
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


def baseline(Y_val_non_binary, bp_length_comp=1500):
    # TODO: import methylation without conversion values
    # TODO: convert to 1 and 0 depending on mean calculation of 1500bp left and right
    # TODO: comparing percentage correct to actual values from imported data in function
    # TODO: Implementing binary search
    Y_val_genomic_coord = Y_val_non_binary.iloc[:, 4]
    Y_values = np.array(Y_val_genomic_coord.sort_values())

    comparison_df = {"Non_binary": Y_values, "Binary": None, "WinCPG": None}

    Binary_list = []
    WinCPG_list = []

    for index, value in enumerate(Y_values):

        # calculating upper and lower limit and finding indices (TODO binary search)

        mean_Win_list = []
        WinCPG_val = value
        upper_limit = WinCPG_val + 1500
        lower_limit = WinCPG_val - 1500

        if lower_limit <= WinCPG_val <= upper_limit:
            associated_b_val = Y_val_non_binary.iloc[index, 5]
            mean_Win_list.append(associated_b_val)

        # calculating the mean value over the whole list
        mean_Win_Value = np.mean(mean_Win_list)
        WinCPG_list.append(mean_Win_Value)

        # calculating binary values and appending them to a list
        Binary_val = value
        if Binary_val >= 0.5:
            Binary_val = 1
        else:
            Binary_val = 0
        Binary_list.append(Binary_val)

    # print(WinCPG_list[0:10])
    WinCPG_list = list(map(lambda x: 1 if x >= 0.5 else 0, WinCPG_list))
    # print(WinCPG_list[0:10])

    counter = 0
    for index_2, l in enumerate(Binary_list):
        if l == WinCPG_list[index_2]:
            counter += 1

    percentage_correct = counter / len(Binary_list)
    print(f'The percentage correct by the baseline model is {percentage_correct}')

    # creating df out of dictionary

    return percentage_correct


def random_forest(X, Y, test_s=0.4):
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
    plt.figure(figsize=(15, 8))
    plt.title("Feature Importance's RF all features")  ###CHANGE WHEN CHANGING SINGLE PATIENT TO FULL RUN OR REVERSE!
    plt.bar(range(X_train.shape[1]), importances[sorted_indices], align="center")
    plt.xticks(range(X_train.shape[1]), [feature_names[i] for i in sorted_indices], rotation=45, fontsize=8)

    plt.xlabel("Feature")
    plt.ylabel("Importance")
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
    plt.title(f'Precision-Recall Curve RF all features (AUC = {pr_auc:.2f})',
              fontsize=16)  # Set your desired font size  ###CHANGE WHEN CHANGING SINGLE PATIENT TO FULL RUN OR REVERSE!

    # Adjust font size for x and y axis tick labels
    plt.xticks(fontsize=12)  # Set your desired font size
    plt.yticks(fontsize=12)  # Set your desired font size

    plt.show()
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

    # Ensure X_train and X_test have the same column order
    X_train = X_train[X_test.columns]

    # Sort features by importance
    sorted_indices = importances.argsort()[::-1]
    sorted_importances = importances[sorted_indices]

    # Plot feature importances
    plt.figure(figsize=(15, 15))
    plt.title("Feature Importance's RF Cross-Sampling all features", fontsize = 27)

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
    plt.savefig("Feature importance RF all features cross section.png", dpi = 300)
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
    plt.title(f'Precision-Recall Curve RF Cross-Sampling all features (AUC = {pr_auc:.2f})',
              fontsize=16)  # Set your desired font size  ###CHANGE WHEN CHANGING SINGLE PATIENT TO FULL RUN OR REVERSE!

    # Adjust font size for x and y axis tick labels
    plt.xticks(fontsize=12)  # Set your desired font size
    plt.yticks(fontsize=12)  # Set your desired font size

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

    # Make predictions on the testing data
    Y_pred_proba = logistic_regression_model.predict_proba(X_test)[:, 1]

    # Calculate Precision-Recall curve
    precision, recall, _ = precision_recall_curve(Y_test, Y_pred_proba)
    pr_auc = auc(recall, precision)

    # Get the coefficients (feature importances)
    coefficients = logistic_regression_model.coef_[0]

    # Assuming X is your feature matrix (DataFrame)
    feature_names = X.columns

    # Sort features by importance
    sorted_indices = np.argsort(np.abs(coefficients))[::-1]

    # Plot feature importances
    plt.figure(figsize=(15, 15))
    plt.title("Feature Importance's Logistic Regression", fontsize = 27)
    plt.bar(range(X_train.shape[1]), np.abs(coefficients)[sorted_indices], align="center")
    plt.xticks(range(X_train.shape[1]), [feature_names[i] for i in sorted_indices], rotation=45, ha="right", fontsize=8)

    plt.xlabel('Feature', fontsize=27)  # Set your desired font size
    plt.ylabel('Importance', fontsize=27)  # Set your desired font size
    # Adjust font size for x and y axis tick labels
    plt.xticks(fontsize=14)  # Set your desired font size
    plt.yticks(fontsize=23)  # Set your desired font size

    plt.tight_layout()  # Adjust layout to prevent clipping of labels
    plt.savefig("Feature importance Log regression all features.png", dpi = 300)

    plt.show()


    # Plot Precision-Recall curve
    plt.figure(figsize=(8, 6))
    plt.step(recall, precision, color='b', alpha=0.2, where='post')
    plt.fill_between(recall, precision, step='post', alpha=0.15, color='b')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(f'Precision-Recall Curve Logistic Regression (AUC = {pr_auc:.2f})')

    # Adjust font size for title, x-axis label, and y-axis label
    plt.xlabel('Recall', fontsize=14)  # Set your desired font size
    plt.ylabel('Precision', fontsize=14)  # Set your desired font size
    plt.title(f'Precision-Recall logistic regression all features (AUC = {pr_auc:.2f})',
              fontsize=16)  # Set your desired font size  ###CHANGE WHEN CHANGING SINGLE PATIENT TO FULL RUN OR REVERSE!

    # Adjust font size for x and y axis tick labels
    plt.xticks(fontsize=12)  # Set your desired font size
    plt.yticks(fontsize=12)  # Set your desired font size

    plt.show()


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
    # TODO: PCA and SVM?

    full_Data_run = True
    RF_Cross =False
    RF = False
    Log_Reg = False
    Baseline = True

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
        # print(X_final_columns)
        X_final_dataset = final_dataset[X_final_columns]
        Y_final_dataset = final_dataset["Patient_1"]

        '''Full data run cross-section'''
        X_train_set = X_final_dataset
        Y_train_set = Y_final_dataset

        X_test_set_columns = patient_3.columns.tolist()[5:]
        X_test_set = patient_3[X_test_set_columns]

        Y_test_set = patient_3["Patient_3"]

        # running random forest
        if RF == True:
            random_forest(X_final_dataset, Y_final_dataset)

        # running random forest cross sample
        if RF_Cross == True:
            random_forest_cross_sample(X_train_set, Y_train_set, X_test_set, Y_test_set)

        # running logistic regression
        if Log_Reg:
            logistic_regression(X_final_dataset, Y_final_dataset)

        # SVM
        # Support_Vector_Machines(X_final_dataset,Y_final_dataset)

        # Baseline
        if Baseline == True:
            baseline(non_binary_data_p_1)
