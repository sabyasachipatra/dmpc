# Importing the required packages
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report


# Function to perform training with giniIndex.
def train_using_gini(X_train, y_train):
  
    # Creating the classifier object
    clf_gini = DecisionTreeClassifier(criterion = "gini",
            random_state = 100,max_depth=3, min_samples_leaf=5)
  
    # Performing training
    clf_gini.fit(X_train, y_train)
    return clf_gini
      

# Function to make predictions
def prediction(X_test, clf_object):
  
    # Predicton on test with giniIndex
    y_pred = clf_object.predict(X_test)
    return y_pred

# Function to calculate accuracy
def cal_accuracy(y_test, y_pred):
      
    print("Confusion Matrix: ",
        confusion_matrix(y_test, y_pred))
      
    print ("Accuracy : ",
    accuracy_score(y_test,y_pred)*100)
      
    print("Report : ",
    classification_report(y_test, y_pred))



trainf = pd.read_csv("trainfeatures.csv")
testf = pd.read_csv("testfeatures.csv")

X_train = trainf.values[:,0:21]
y_train = trainf.values[:,22]
X_test = testf.values[:,0:21]
y_test = testf.values[:,22]

clf_gini = train_using_gini(X_train, y_train)
y_pred_gini = prediction(X_test, clf_gini)

with open('out.txt', 'w') as f:
    for line in y_pred_gini:
        f.write(line)
        f.write('\n')






  