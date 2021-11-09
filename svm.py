# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 21:50:34 2021

@author: DELL
"""
import pandas as pd
import numpy as np
pd.set_option('display.max_rows',None)
pd.set_option('display.max_columns',None)
from sklearn import svm
import matplotlib.pyplot as plt
import matplotlib
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn import cross_validation
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
help(auc)
sklearn.metrics._ranking
## 读取
data = pd.read_csv('E:\\sirebrowser\\STAD\\miR\\分析\\26k=3.csv', engine='python')
data = data.iloc[:,:28]
data.iloc[0:2,0:3]
data = data.set_index('Sample')
# data = data.drop('Sample',axis=1)
# 确保在每个group中分成0.7 0.3
g0 = data[data['groups'] == 0]
g0.info()
g1 = data[data['groups'] == 1]
g1.info()
g2 = data[data['groups'] == 2]
g2.info()

x0,y0=np.split(g0,indices_or_sections=(26,),axis=1)
train_data0,test_data0,train_label0,test_label0 =train_test_split(x0,y0, random_state=1, train_size=0.7,test_size=0.3)
x1,y1=np.split(g1,indices_or_sections=(26,),axis=1) 
train_data1,test_data1,train_label1,test_label1 =train_test_split(x1,y1, random_state=1, train_size=0.7,test_size=0.3) 
x2,y2=np.split(g2,indices_or_sections=(26,),axis=1) 
train_data2,test_data2,train_label2,test_label2 =train_test_split(x2,y2, random_state=1, train_size=0.7,test_size=0.3) 

frames = [train_data0, train_data1, train_data2]
train_data = pd.concat(frames)
train_data.info()
frames = [test_data0, test_data1, test_data2]
test_data = pd.concat(frames)
test_data.info()
frames = [train_label0, train_label1, train_label2]
train_label = pd.concat(frames)
train_label.info()
frames = [test_label0, test_label1, test_label2]
test_label = pd.concat(frames)
test_label.info()
frames = [x0, x1, x2]
x = pd.concat(frames)
frames = [y0, y1, y2]
y = pd.concat(frames)

 
#3.训练svm分类器
classifier = svm.SVC(C=1, kernel='rbf', gamma=0.05,decision_function_shape='ovo', probability=True)  # ovr:一对多策略
classifier.fit(train_data, train_label.values.ravel())  # ravel函数在降维时默认是行序优先
#交叉验证
scores = cross_val_score(classifier, x, y, cv=10)
scores.mean()

#4.计算svc分类器的准确率
print("训练集：",classifier.score(train_data,train_label))
print("测试集：",classifier.score(test_data,test_label))

#结果写出
print('train_decision_function:\n',classifier.decision_function(train_data)) # (90,3)
print('predict_result:\n',classifier.predict(train_data))
len(classifier.predict(train_data))

data = pd.read_csv('E:\\sirebrowser\\STAD\\miR\\分析\\26k=3.csv', engine='python')
train_data.info()
train_data.iloc[0:2,0:2]
train_data.insert(26,'Group',classifier.predict(train_data))
train_data.insert(0,'Sample',train_data.index)
train_data.to_csv('E:\\sirebrowser\\STAD\\miR\\分析\\26k=3_train.csv',index=0)

test_data.info()
test_data.insert(26,'Group',classifier.predict(test_data))
test_data.insert(0,'Sample',test_data.index)
test_data.to_csv('E:\\sirebrowser\\STAD\\miR\\分析\\26k=3_test.csv',index=0)

## 画图
# pred = classifier.predict_proba(test_data)[:,1]
# #############画图部分
# fpr, tpr, threshold = roc_curve(test_label.ravel(), pred)
# roc_auc = auc(fpr, tpr)
# plt.figure(figsize=(6,6))
# plt.title('Validation ROC')
# plt.plot(fpr, tpr, 'b', label = 'Val AUC = %0.3f' % roc_auc)
# plt.legend(loc = 'lower right')
# plt.plot([0, 1], [0, 1],'r--')
# plt.xlim([0, 1])
# plt.ylim([0, 1])
# plt.ylabel('True Positive Rate')
# plt.xlabel('False Positive Rate')
# plt.show()
x_nd = x.values
y_nd = y.values
y = label_binarize(y_nd, classes=[0, 1, 2])
n_classes = y.shape[1]
n_samples, n_features = x_nd.shape

train_data_1 = train_data.values
train_label_1 = label_binarize(train_label, classes=[0, 1, 2])
test_label_1 = label_binarize(test_label, classes=[0, 1, 2])

classifier = OneVsRestClassifier(svm.SVC(C=1,kernel='rbf',gamma=0.05,decision_function_shape='ovo',probability=True))
y_score = classifier.fit(train_data,train_label_1).decision_function(test_data)
y_score = classifier.fit(test_data,test_label_1).decision_function(train_data)


fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(train_label_1[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

lw = 2
colors = cycle(['#00A087FF', '#F64B35FF', '#7E6148FF'])
for i, color in zip(range(n_classes), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=lw,
             label='ROC curve of Group {0} (area = {1:0.4f})'
             ''.format(i, roc_auc[i]))
 
plt.plot([0, 1], [0, 1], 'k--', lw=1)
# plt.xlim([0.0, 1.0])
# plt.ylim([0.0, 1.05])
plt.xlabel('Specificity')
plt.ylabel('Sencitivity')
plt.title('k=3 top26 train_data')
plt.legend(loc="lower right")
plt.show()

# 参数选择
clf=RandomForestClassifier(random_state = 42, class_weight="balanced")
g = GridSearchCV(estimator=clf,cv=10, n_jobs=4,scoring = 'accuracy',param_grid={'n_estimators':range(10,300,10)})
g.fit(x,y)
print(g.best_params_)
print(g.best_score_)

# random forest 89 183 113
clf=RandomForestClassifier(n_estimators =130, random_state = 42, class_weight="balanced")
output = cross_validate(clf, x_nd, y_nd, cv=183, scoring = 'accuracy', return_estimator =True,n_jobs=4)

s = feature_importances[0]

for idx,estimator in enumerate(output['estimator']):
    feature_importances[idx] = pd.DataFrame(estimator.feature_importances_,
                                       index = train_data.columns,
                                        columns=['importance'])
    # s[idx] = feature_importances
    # print(feature_importances[idx])
    s = pd.concat([s,feature_importances[idx]],axis=1)
s.head()
s.mean(1).sort_values(ascending=False)

################################################################33
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn.multiclass import OneVsOneClassifier

from scipy import interp
 
# # 加载数据
# iris = datasets.load_iris()
# X = iris.data #array
# y = iris.target
# # 将标签二值化
# y = label_binarize(y, classes=[0, 1, 2])
# # 设置种类
# n_classes = y.shape[1]
 
# # 训练模型并预测
# # random_state = np.random.RandomState(0)
# n_samples, n_features = X.shape
 
# # shuffle and split training and test sets
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.3,random_state=0)
# len(y_test)
# # Learn to predict each class against the other
# classifier = OneVsRestClassifier(svm.SVC(kernel='linear', probability=True))#一对多
# y_score = classifier.fit(X_train, y_train).decision_function(X_test)
# len(y_score)
# # 计算每一类的ROC
# fpr = dict()
# tpr = dict()
# roc_auc = dict()
# for i in range(n_classes):
#     fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
#     roc_auc[i] = auc(fpr[i], tpr[i])
 
# # # Compute micro-average ROC curve and ROC area（方法二）
# # fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
# # roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
 
# # # Compute macro-average ROC curve and ROC area（方法一）
# # # First aggregate all false positive rates
# # all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))
# # # Then interpolate all ROC curves at this points
# # mean_tpr = np.zeros_like(all_fpr)
# # for i in range(n_classes):
# #     mean_tpr += interp(all_fpr, fpr[i], tpr[i])
# # # Finally average it and compute AUC
# # mean_tpr /= n_classes
# # fpr["macro"] = all_fpr
# # tpr["macro"] = mean_tpr
# # roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])
 
# # # Plot all ROC curves
# lw=2
# # plt.figure()
# # plt.plot(fpr["micro"], tpr["micro"],
# #          label='micro-average ROC curve (area = {0:0.2f})'
# #                ''.format(roc_auc["micro"]),
# #          color='deeppink', linestyle=':', linewidth=4)
 
# # plt.plot(fpr["macro"], tpr["macro"],
# #          label='macro-average ROC curve (area = {0:0.2f})'
# #                ''.format(roc_auc["macro"]),
# #          color='navy', linestyle=':', linewidth=4)
 
# colors = cycle(['aqua', 'darkorange', 'cornflowerblue'])
# for i, color in zip(range(n_classes), colors):
#     plt.plot(fpr[i], tpr[i], color=color, lw=lw,
#              label='ROC curve of class {0} (area = {1:0.2f})'
#              ''.format(i, roc_auc[i]))
 
# plt.plot([0, 1], [0, 1], 'k--', lw=lw)
# plt.xlim([0.0, 1.0])
# plt.ylim([0.0, 1.05])
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('Some extension of Receiver operating characteristic to multi-class')
# plt.legend(loc="lower right")
# plt.show()

from sklearn import datasets
from sklearn.model_selection import cross_validate
from sklearn.svm import LinearSVC
from sklearn.ensemble import  RandomForestClassifier
import pandas as pd

# diabetes = datasets.load_diabetes()
# X, y = diabetes.data, diabetes.target

# clf=RandomForestClassifier(n_estimators =10, random_state = 42, class_weight="balanced")
# output = cross_validate(clf, X, y, cv=10, scoring = 'accuracy', return_estimator =True)

# for idx,estimator in enumerate(output['estimator']):
#     print("Features sorted by their score for estimator {}:".format(idx))
#     feature_importances = pd.DataFrame(estimator.feature_importances_,
#                                        index = diabetes.feature_names,
#                                         columns=['importance']).sort_values('importance', ascending=False)
#     print(feature_importances)
# 查看权重
clf=RandomForestClassifier(n_estimators =10, random_state = 42, class_weight="balanced")
output = cross_validate(clf, x_nd, y_nd, cv=10, scoring = 'accuracy', return_estimator =True)

# for idx,estimator in enumerate(output['estimator']):
#     print("Features sorted by their score for estimator {}:".format(idx))
#     feature_importances = pd.DataFrame(estimator.feature_importances_,
#                                        index = train_data.columns,
#                                         columns=['importance']).sort_values('importance', ascending=False)
#     print(feature_importances)

############################################## top15/20...重新km
k = 3 #聚类的类别
iteration = 500 #聚类最大循环次数

data = pd.read_csv('E:\\sirebrowser\\STAD\\miR\\分析\\26k=3.csv', engine='python',index_col = 'Sample') 
data.info()
b = data.iloc[:,0:19] # top数+1
b.info()
model = KMeans(n_clusters = k, n_jobs = 4, max_iter = iteration) 
model.fit(b) 
r = pd.concat([b, pd.Series(model.labels_, index = data.index)], axis = 1)  #详细输出每个样本对应的类别
r.columns = list(b.columns) + [u'kgroups19'] #重命名表头
r.iloc[0:2,0:2]
r.insert(0,'Sample',data.index)
r.to_csv('E:\\sirebrowser\\STAD\\miR\\分析\\26k=3importance_19.csv',index=0) #保存结果
# km图
fontt={'color': 'k',
      'size': 27,
      'family': 'Arial'}
fonty={'color': 'k',
      'size': 25,
      'family': 'Arial'}
font={'color': 'k',
      'size': 10,
      'family': 'Arial'}
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\26k=3importance_19.csv',engine='python')
a.info()
a.insert(20,'time',data['time'].values)
a.insert(21,'event',data['event'].values)
a.info()
b = a.iloc[:,[20,21,22]]
b.head()
kmf = KaplanMeierFitter()
groups = b['kgroups19']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)

T = b['time']
E = b['event']
dem1 = (b['kgroups19'] == 0)
dem2 = (b['kgroups19'] == 1)
dem3 = (b['kgroups19'] == 2)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)

kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(show_censors=True,ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,show_censors=True,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(show_censors=True,ci_show=False,color='#7E6148FF'     )

plt.legend(loc=(0.6,0.65),prop={'family' : 'Arial', 'size'   : 22},handletextpad=0.5,frameon=False,labelspacing=0.1)
plt.tick_params(width=2)
ax.spines['bottom'].set_linewidth('2')
ax.spines['top'].set_linewidth('0')
ax.spines['left'].set_linewidth('2')
ax.spines['right'].set_linewidth('0')
plt.xlim(-230,3900,500)
plt.ylim(-0.08,1.08)

plt.title('k=3 top19', fontdict=fontt)
# plt.text(0, 0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
# plt.text(0, 0, 'Group 0     vs     Group 2'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
# plt.text(0, -0.05, 'Group 1     vs     Group 2'+"      P_value=%.2e"%results_2.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel(' ', fontdict=fonty)
plt.yticks(fontproperties = 'Arial', size = 20)
plt.xticks(fontproperties = 'Arial', size = 20,rotation=20)
# plt.figure(dpi=1000,figsize=(24,20))
plt.savefig(fname='E:\\sirebrowser\\STAD\\miR\\分析\\top26\\19k=3.png',figsize=[20,20],dpi=500)
plt.close('all')
#################################################################################
## svm
data = pd.read_csv('E:\\sirebrowser\\STAD\\miR\\分析\\26k=3importance_26.csv', engine='python')
data.info()
# da = data.iloc[:,:17] 
data = data.set_index('Sample')
# 确保在每个group中分成0.7 0.3
g0 = data[data['kgroups26'] == 0]
g0.info()
g1 = data[data['kgroups26'] == 1]
g1.info()
g2 = data[data['kgroups26'] == 2]
g2.info()

x0,y0=np.split(g0,indices_or_sections=(26,),axis=1) #top
train_data0,test_data0,train_label0,test_label0 =train_test_split(x0,y0, random_state=1, train_size=0.7,test_size=0.3) #sklearn.model_selection.

x1,y1=np.split(g1,indices_or_sections=(26,),axis=1) 
train_data1,test_data1,train_label1,test_label1 =train_test_split(x1,y1, random_state=1, train_size=0.7,test_size=0.3) #sklearn.model_selection.

x2,y2=np.split(g2,indices_or_sections=(26,),axis=1) 
train_data2,test_data2,train_label2,test_label2 =train_test_split(x2,y2, random_state=1, train_size=0.7,test_size=0.3) #sklearn.model_selection.

frames = [train_data0, train_data1, train_data2]
train_data = pd.concat(frames)
train_data.info()
frames = [test_data0, test_data1, test_data2]
test_data = pd.concat(frames)
test_data.info()
frames = [train_label0, train_label1, train_label2]
train_label = pd.concat(frames)
train_label.info()
frames = [test_label0, test_label1, test_label2]
test_label = pd.concat(frames)
test_label.info()
frames = [x0, x1, x2]
x = pd.concat(frames)
frames = [y0, y1, y2]
y = pd.concat(frames) 

 
#3.训练svm分类器
classifier=svm.SVC(C=2,kernel='rbf',gamma=0.05,decision_function_shape='ovr',probability=True) 
classifier.fit(train_data,train_label.values.ravel()) #ravel函数在降维时默认是行序优先
#交叉验证
scores = cross_val_score   (classifier, test_data, test_label, cv=10)
scores.mean()

#4.计算svc分类器的准确率
print("训练集：",classifier.score(train_data,train_label))
print("测试集：",classifier.score(test_data,test_label))

# 写出
train_data.info()
train_data.insert(26,'sGroup',classifier.predict(train_data)) #top
train_data.insert(0,'Sample',train_data.index)
train_data.to_csv('E:\\sirebrowser\\STAD\\miR\\分析\\26k=3_train.csv',index=0)

test_data.info()
test_data.insert(26,'sGroup',classifier.predict(test_data))
test_data.insert(0,'Sample',test_data.index)
test_data.to_csv('E:\\sirebrowser\\STAD\\miR\\分析\\26k=3_test.csv',index=0)
# roc画图
frames = [train_data0, train_data1, train_data2]
train_data = pd.concat(frames)
train_data.info()
frames = [test_data0, test_data1, test_data2]
test_data = pd.concat(frames)
test_data.info()
frames = [train_label0, train_label1, train_label2]
train_label = pd.concat(frames)
train_label.info()
frames = [test_label0, test_label1, test_label2]
test_label = pd.concat(frames)
test_label.info()
frames = [x0, x1, x2]
x = pd.concat(frames)
frames = [y0, y1, y2]
y = pd.concat(frames) 
x_nd = x.values
y_nd = y.values
y = label_binarize(y_nd, classes=[0, 1, 2])
n_classes = y.shape[1]
n_samples, n_features = x_nd.shape


train_label_1 = label_binarize(train_label, classes=[0, 1, 2])
test_label_1 = label_binarize(test_label, classes=[0, 1, 2])

# classifier = OneVsRestClassifier(svm.SVC(C=1,kernel='rbf',gamma=0.05,decision_function_shape='ovo',probability=True))
# classifier = OneVsOneClassifier(svm.SVC(C=1,kernel='rbf',gamma=0.05,decision_function_shape='ovo',probability=True))


#test
y_score = classifier.fit(train_data,train_label).decision_function(test_data)

fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(test_label_1[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])
    
# Compute micro-average ROC curve and ROC area
fpr["micro"], tpr["micro"], _ = roc_curve(test_label_1.ravel(), y_score.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))
# Then interpolate all ROC curves at this points
mean_tpr = np.zeros_like(all_fpr)
for i in range(n_classes):
    mean_tpr += interp(all_fpr, fpr[i], tpr[i])
# Finally average it and compute AUC
mean_tpr /= n_classes
fpr["macro"] = all_fpr
tpr["macro"] = mean_tpr
roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])


lw = 2
colors = cycle(['#00A087FF', '#F64B35FF', '#7E6148FF'])
for i, color in zip(range(n_classes), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=lw,
              label='Group {0} (area = {1:0.4f})'
              ''.format(i, roc_auc[i]))
    
plt.plot(fpr["macro"], tpr["macro"],
         label='area = {0:0.4f}'
               ''.format(roc_auc["macro"]),
         color='navy', linestyle='-', linewidth=4)


plt.plot([0, 1], [0, 1], 'k--', lw=3)
# plt.xlim([0.0, 1.0])
# plt.ylim([0.0, 1.05])
plt.legend(prop={'family' : 'Arial', 'size'   : 22},loc="lower right",handletextpad=0.5,frameon=False,labelspacing=0.1)
plt.tick_params(width=2)
ax=plt.gca()
ax.spines['bottom'].set_linewidth('2')
ax.spines['top'].set_linewidth('0')
ax.spines['left'].set_linewidth('2')
ax.spines['right'].set_linewidth('0')
plt.title('k=3 top15 test_data',fontdict=fontt)
plt.xlabel('False Positive Rate',fontdict=fonty)
plt.ylabel('True Positive Rate',fontdict=fonty)
plt.yticks(fontproperties = 'Arial', size = 25)
plt.xticks(fontproperties = 'Arial', size = 25)
plt.savefig(fname='E:\\sirebrowser\\STAD\\miR\\分析\\top26\\26roc.png',figsize=[8,6],dpi=800, bbox_inches = 'tight')
plt.close('all')

#train    
# y_score = classifier.fit(train_data,train_label).decision_function(train_data)
# fpr = dict()
# tpr = dict()
# roc_auc = dict()
# for i in range(n_classes):
#     fpr[i], tpr[i], _ = roc_curve(train_label_1[:, i], y_score[:, i])
#     roc_auc[i] = auc(fpr[i], tpr[i])
# fpr["micro"], tpr["micro"], _ = roc_curve(train_label_1.ravel(), y_score.ravel())
# roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

# all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))
# # Then interpolate all ROC curves at this points
# mean_tpr = np.zeros_like(all_fpr)
# for i in range(n_classes):
#     mean_tpr += interp(all_fpr, fpr[i], tpr[i])
# # Finally average it and compute AUC
# mean_tpr /= n_classes
# fpr["macro"] = all_fpr
# tpr["macro"] = mean_tpr
# roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])



# lw = 2
# colors = cycle(['#00A087FF', '#F64B35FF', '#7E6148FF'])
# for i, color in zip(range(n_classes), colors):
#     plt.plot(fpr[i], tpr[i], color=color, lw=lw,
#              label='Group {0} (area = {1:0.4f})'
#              ''.format(i, roc_auc[i]))
 
# plt.plot([0, 1], [0, 1], 'k--', lw=1)
# # plt.xlim([0.0, 1.0])
# # plt.ylim([0.0, 1.05])
# plt.legend(prop={'family' : 'Arial', 'size'   : 20},loc="lower right",handletextpad=0.5,frameon=False,labelspacing=0.1)
# plt.tick_params(width=2)
# ax=plt.gca()
# ax.spines['bottom'].set_linewidth('2')
# ax.spines['top'].set_linewidth('0')
# ax.spines['left'].set_linewidth('2')
# ax.spines['right'].set_linewidth('0')
# plt.title('k=3 top26 train_data',fontdict=fontt)
# plt.xlabel('Specificity',fontdict=fonty)
# plt.ylabel('Sencitivity',fontdict=fonty)
# plt.yticks(fontproperties = 'Arial', size = 20)
# plt.xticks(fontproperties = 'Arial', size = 20)
# plt.savefig(fname='E:\\sirebrowser\\STAD\\miR\\分析\\top26\\26roc_1.png',figsize=[10,8],dpi=500)
# plt.close('all')

#combine
data = pd.concat([train_data,test_data],axis=0)
data.info()
label = pd.concat([train_label,test_label],axis=0)
label.info()
label = label_binarize(label, classes=[0, 1, 2])

y_score_c = classifier.fit(train_data,train_label).decision_function(data)

fpr_c = dict()
tpr_c = dict()
roc_auc_c = dict()
for i in range(n_classes):
    fpr_c[i], tpr_c[i], _ = roc_curve(label[:, i], y_score_c[:, i])
    roc_auc_c[i] = auc(fpr_c[i], tpr_c[i])
    
fpr_c["micro"], tpr_c["micro"], _ = roc_curve(label.ravel(), y_score_c.ravel())
roc_auc_c["micro"] = auc(fpr_c["micro"], tpr_c["micro"])

all_fpr_c = np.unique(np.concatenate([fpr_c[i] for i in range(n_classes)]))
# Then interpolate all ROC curves at this points
mean_tpr_c = np.zeros_like(all_fpr_c)
for i in range(n_classes):
    mean_tpr_c += interp(all_fpr_c, fpr_c[i], tpr_c[i])
# Finally average it and compute AUC
mean_tpr_c /= n_classes
fpr_c["macro"] = all_fpr_c
tpr_c["macro"] = mean_tpr_c
roc_auc_c["macro"] = auc(fpr_c["macro"], tpr_c["macro"])

lw = 2
colors = cycle(['#00A087FF', '#F64B35FF', '#7E6148FF'])
for i, color in zip(range(n_classes), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=lw,
              label='Group {0} (area = {1:0.4f})'
              ''.format(i, roc_auc[i]))
plt.figure(figsize=(8,8))
plt.plot(fpr["macro"], tpr["macro"],
         label='               {0:0.4f}'
               ''.format(roc_auc["macro"]),
         color='navy', linestyle='-', linewidth=8)

plt.plot(fpr_c["macro"], tpr_c["macro"],
         label='               {0:0.4f}'
               ''.format(roc_auc_c["macro"]),
         color='#F64B35FF', linestyle='-', linewidth=8)
 
plt.plot([0, 1], [0, 1], 'k--', lw=5)
# plt.xlim([0.0, 1.0])
# plt.ylim([0.0, 1.05])
plt.legend(prop={'family' : 'Arial', 'size'   : 38},handletextpad=0.5,frameon=False,labelspacing=0.1,loc=(0.15,-0.03))
plt.tick_params(width=6)
plt.tick_params(length=6)

ax=plt.gca()
ax.spines['bottom'].set_linewidth('4')
ax.spines['top'].set_linewidth('0')
ax.spines['left'].set_linewidth('4')
ax.spines['right'].set_linewidth('0')
# plt.title('k=3 top15',fontdict=fontt)
# plt.xlabel('False Positive Rate',fontdict=fonty)
# plt.ylabel('True Positive Rate',fontdict=fonty)
plt.yticks(fontproperties = 'Arial', size = 0)
plt.xticks(fontproperties = 'Arial', size = 0)
plt.savefig(fname='E:\\sirebrowser\\STAD\\miR\\分析\\top26\\26roc.png',dpi=1000, bbox_inches = 'tight')
plt.close('all')
######################################################################################
## svm km图
fontt={'color': 'k',
      'size': 25,
      'family': 'Arial'}
fonty={'color': 'k',
      'size': 20,
      'family': 'Arial'}
font={'color': 'k',
      'size': 10,
      'family': 'Arial'}
clin = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\cli_385.csv',engine='python')
clin.info()
c = clin.iloc[:,7:10]
c.head()
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\20result.csv',engine='python')
a.info()
a = a.iloc[:,[1,2,4]]
a = a.iloc[:,2:]
# 匹配
list_custom = c['Sample'].values
a['Sample.1'] = a['Sample.1'].astype('category')
a['Sample.1'].cat.reorder_categories(list_custom, inplace=True)
a.sort_values('Sample.1', inplace=True)
a = a.set_index('Sample.1')
a.head()
# a = a.drop('time',axis=1)
a.insert(0,'time',c['time'].values)
a.insert(1,'event',c['event'].values)
b = a

kmf = KaplanMeierFitter()
groups = b['sGroup']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)

T = b['time']
E = b['event']
dem1 = (b['sGroup'] == 0)
dem2 = (b['sGroup'] == 1)
dem3 = (b['sGroup'] == 2)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)

kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(show_censors=True,ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,show_censors=True,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(show_censors=True,ci_show=False,color='#7E6148FF'     )

plt.legend(loc=(0.6,0.65),prop={'family' : 'Arial', 'size'   : 22},handletextpad=0.5,frameon=False,labelspacing=0.1)
plt.tick_params(width=2)
ax.spines['bottom'].set_linewidth('2')
ax.spines['top'].set_linewidth('0')
ax.spines['left'].set_linewidth('2')
ax.spines['right'].set_linewidth('0')
plt.xlim(-230,3900,500)
plt.ylim(-0.08,1.08)

plt.title('k=3 top20(SVM)', fontdict=fontt)
# plt.text(0, 0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
# plt.text(0, 0, 'Group 0     vs     Group 2'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
# plt.text(0, -0.05, 'Group 1     vs     Group 2'+"      P_value=%.2e"%results_2.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel(' ', fontdict=fonty)
plt.yticks(fontproperties = 'Arial', size = 20)
plt.xticks(fontproperties = 'Arial', size = 20,rotation=20)
# plt.figure(dpi=1000,figsize=(24,20))
plt.savefig(fname='E:\\sirebrowser\\STAD\\miR\\分析\\top26\\20svm.png',figsize=[20,20],dpi=500)
plt.close('all')

## 细分表格
a = pd.read_csv('E:\\sirebrowser\\STAD\\clinical\\km\\20stad_pathologic_stage_1.csv',engine='python')
b = a.replace(np.nan,'null')
a.info()
m0 = b[b['groups.1'] == 'm1']
m0[m0['kgroups'] == 0].info()
m0[m0['kgroups'] == 1].info()
m0[m0['kgroups'] == 2].info()
m0[m0['kgroups'] == 'null'].info()


###############################################################################################
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\26k=3.csv',engine='python')
a.info()
b = a.iloc[:,[27,28,30]]
b.head()
b.info()
kmf = KaplanMeierFitter()
groups = b['groupsc']
ix1 = (groups == 1)
ix2 = (groups == 2)

T = b['time']
E = b['event']
dem1 = (b['groupsc'] == 1)
dem2 = (b['groupsc'] == 2)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)

kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0+1')
ax = kmf.plot(show_censors=True,ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 2')
ax = kmf.plot(ax=ax,show_censors=True,ci_show=False,color='#7E6148FF')

plt.legend(loc=(0.6,0.65),prop={'family' : 'Arial', 'size'   : 22},handletextpad=0.5,frameon=False,labelspacing=0.1)
plt.tick_params(width=2)
ax.spines['bottom'].set_linewidth('2')
ax.spines['top'].set_linewidth('0')
ax.spines['left'].set_linewidth('2')
ax.spines['right'].set_linewidth('0')
plt.xlim(-230,3900,500)
plt.ylim(-0.08,1.08)
plt.axvline(x=1825,c='k',ls='--',lw=0.5)
plt.axvline(x=1095,c='k',ls='--',lw=0.5)
plt.title('k=3 top26(combine)', fontdict=fontt)
# plt.text(0, 0.05, 'Group 0+1     vs     Group 2'+"      P_value=%.2e"%results.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel(' ', fontdict=fonty)
plt.yticks(fontproperties = 'Arial', size = 20)
plt.xticks(fontproperties = 'Arial', size = 20,rotation=20)
# plt.figure(dpi=1000,figsize=(24,20))
plt.savefig(fname='E:\\sirebrowser\\STAD\\miR\\分析\\top26\\26combine.png',figsize=[20,20],dpi=500)
plt.close('all')

        
