import numpy as np
import matplotlib.pyplot as plt
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate
import pandas as pd

# Create a random dataset
# rng = np.random.RandomState(1)
# X = np.sort(200 * rng.rand(100, 1) - 100, axis=0)
# y = np.array([np.pi * np.sin(X).ravel(), np.pi * np.cos(X).ravel()]).T
# y[::5, :] += (0.5 - rng.rand(20, 2))

# Fit regression model
# regr_1 = DecisionTreeRegressor(max_depth=2)
# regr_2 = DecisionTreeRegressor(max_depth=5)
# regr_3 = DecisionTreeRegressor(max_depth=8)
# regr_1.fit(X, y)
# regr_2.fit(X, y)
# regr_3.fit(X, y)

# Predict
# X_test = np.arange(-100.0, 100.0, 0.01)[:, np.newaxis]
# y_1 = regr_1.predict(X_test)
# y_2 = regr_2.predict(X_test)
# y_3 = regr_3.predict(X_test)

train_x_pd = pd.read_table('train_x.txt', delimiter='\t')
train_x = train_x_pd.values

train_y_pd = pd.read_table('train_y.txt', delimiter='\t')
train_y = train_y_pd.values

test_x_pd = pd.read_table('test_x.txt', delimiter='\t')
test_x = test_x_pd.values

test_y_pd = pd.read_table('test_y.txt', delimiter='\t')
test_y = test_y_pd.values

regr_1 = DecisionTreeRegressor(max_depth=20)
print(cross_validate(regr_1, train_x, train_y, cv=10))
print(regr_1.get_params())

regr_1.fit(train_x, train_y)

y_1 = regr_1.predict(test_x)

#print(y_1.shape, test_y.shape)
res=np.corrcoef(y_1.flatten(), test_y.flatten())
print(res)

#loss_values = regr_1.estimator.loss_curve_
#plt.plot(loss_values)
#plt.savefig('demo.pdf')
