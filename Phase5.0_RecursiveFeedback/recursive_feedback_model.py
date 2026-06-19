import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

def simulate_recursive_feedback(data, labels, iterations=10):
    acc_list = []
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    
    for i in range(iterations):
        model.fit(data, labels)
        preds = model.predict(data)
        acc = accuracy_score(labels, preds)
        acc_list.append(acc)
        data = np.hstack([data[:, :3], preds.reshape(-1, 1)])  # Append feedback feature
    
    return acc_list
