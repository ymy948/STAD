# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 20:17:57 2021

@author: DELL
"""
import pandas as pd
import numpy as np
pd.set_option('display.max_rows',None)
pd.set_option('display.max_columns',None)
help(np)
## 保留-01样本 删除-11(normal)样本
rpm = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM.csv',engine='python')
rpm.iloc[0:2,0:2]
#转置
data = rpm.values 
data = list(map(list, zip(*data))) 
data = pd.DataFrame(data) 
data.insert(0,'Gene',rpm.columns)
data.iloc[0:2,0:2]

data['Gene'] = data['Gene'].apply(lambda x:x[:12]).tolist()
rpm_no_dup = data[~data.Gene.duplicated()]
rpm_no_dup.info()
rpm_no_dup.iloc[0:2,0:2]
rpm_no_dup.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_no_dup_sample.csv',index=0,header=0)

## 删除4个无time的样本 通过excel筛选工具

## 删除空值大于样本值20%的gene
rpm_no_dup = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_no_dup_sample.csv',engine='python')
rpm = rpm_no_dup.replace('NA',np.nan)
rpm.info()
rpm.iloc[0:2,0:2]
385*0.2 # 430个样本
nan = rpm.isnull().sum(axis=0)
rpm.drop(columns=nan[nan >= 77].index,inplace=True)
rpm.info()

# 删除空值大于gene20%的样本
463*0.2 # 465个gene
nan = rpm.isnull().sum(axis=1)
nan[nan>92] # 无样本
rpm.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_1.csv',index=0)

## knn填补-1
import sklearn.impute    
from sklearn.impute import KNNImputer
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score as CVS
import matplotlib.pyplot as plt
# t1 = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_1.csv',engine='python')
# k的选取
#直接验证
help(KNNImputer)
sklearn.__version__
t1 = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_1.csv',engine='python')
status = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\knn_test.csv',engine='python')
s = status.iloc[:,463] #提取样本status alive-0,dead-1
t2 = t1.drop('Gene',axis=1  )
for i in range(1,8): 
    imputer = KNNImputer(n_neighbors=i) 
    imputed = imputer.fit_transform(t2) 
    df_imputed = pd.DataFrame(imputed, columns=t2.columns) # knn填补完成
    df_imputed.insert(0,'Gene',t1['Gene']) # 添加样本
    df_imputed.insert(463,'status',s) # 添加status信息
    features = df_imputed.iloc[:,1:462]
    targets = df_imputed.iloc[:,463]
    x_train,x_test,y_train,y_test = train_test_split(features,targets,test_size=0.25)
    knn = KNeighborsClassifier(n_neighbors=i)
    knn.fit(x_train,y_train)
    print(i,knn.score(x_test,y_test))
#k折验证
t1 = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\knn_test.csv',engine='python')
t1.iloc[0:2,0:2]
t1.info()
nan1 = t1.isnull().sum(axis=1)

features = t1.iloc[:,1:462]
targets = t1.iloc[:,463]

nan = features.isnull().sum(axis=0)
features = features.drop(columns=nan[nan != 0].index)
features.info()

x_train,x_test,y_train,y_test = train_test_split(features,targets,test_size=0.25)

for j in range(3,8):
    clf = KNeighborsClassifier(n_neighbors = j)
    clf=clf.fit(x_train,y_train)
    for i in range(2,11):
        cvresult=CVS(clf,features,targets,cv=i)
        print(j,i,cvresult.mean()) # 取得均值
  

    
#处理
t1 = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_1.csv',engine='python')
t1.iloc[0:2,0:2]
t2 = t1.drop('Gene',axis=1)
t2.iloc[0:2,0:2]
nan = t2.isnull().sum(axis=1)
nan[nan==0]
imputer = KNNImputer(n_neighbors=7) 
imputed = imputer.fit_transform(t2) 
df_imputed = pd.DataFrame(imputed, columns=t2.columns) # 此时无行标题
df_imputed.insert(0,'Gene',t1['Gene'])
type(df_imputed)
df_imputed.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_T_1_knn1.csv',index=0)
df_imputed.isnull().sum(axis=0)

## 分位数归一化
rpm_no_nan = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_T_1_knn1.csv',engine='python')
rpm_no_nan.info()
# 转置
data = rpm_no_nan.values # data是数组，直接从文件读出来的数据格式是数组
index1 = list(rpm_no_nan.keys()) # 获取原有csv文件的标题，并形成列表
data = list(map(list, zip(*data))) # map()可以单独列出列表，将数组转换成列表
data = pd.DataFrame(data, index=index1) # 将data的行列转换
data.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_T_1_knn1_T.csv', header=0)
# 处理
rpm_no_nan_t = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_T_1_knn1_T.csv',engine='python')
rpm_no_nan_t.iloc[0:2,0:2]
rpm_no_nan_t_2 = rpm_no_nan_t.drop('Gene',axis=1)
rpm_no_nan_t_2.iloc[0:2,0:2]

rank_mean = rpm_no_nan_t_2.stack().groupby(rpm_no_nan_t_2.rank(method='first').stack().astype(int)).mean()
rpkm_qn = rpm_no_nan_t_2.rank(method='min').stack().astype(int).map(rank_mean).unstack()
rpkm_qn.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_qn.csv',index=0)
#写入列名称
rpkm_qn = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_qn.csv',engine='python')
rpkm_qn.iloc[0:2,0:2]
rpkm_qn.insert(0,'Gene',rpm_no_nan_t['Gene']) # 有小问题 通过excel修改
rpkm_qn.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_qn.csv',index=0)


## knn-2
# 找出knn-1之前的空值并删除
rpm_t = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_1.csv',engine='python') # knn-1之前的数据
rpm_t.iloc[0:3,0:4]
#转置 
data = rpm_t.values # data是数组，直接从文件读出来的数据格式是数组
index1 = list(rpm_t.keys()) # 获取原有csv文件的标题，并形成列表
data = list(map(list, zip(*data))) # map()可以单独列出列表，将数组转换成列表
data = pd.DataFrame(data, index=index1) # 将data的行列转换
data.iloc[0:4,0:2]
data.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_1_T.csv', header=0)
#处理
rpm_t = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_1_T.csv',engine='python')
rpm_t.info()
rpm_qn = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_qn.csv',engine='python')
rpm_qn.info()
for i in range(0,462):
    for j in range(0,386):
        if pd.isnull(rpm_t.iloc[i,j]):
            rpm_qn.iloc[i,j] = rpm_t.iloc[i,j]
        else:
            rpm_qn.iloc[i,j] = rpm_qn.iloc[i,j]
nan = rpm_qn.isnull().sum(axis=0)
rpm_qn.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_qn_NA.csv',index=0)
# 填充
rpm_qn_na = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_qn_NA.csv',engine='python')
rpm_qn_na.iloc[0:2,0:2]
#转置
data = rpm_qn_na.values # data是数组，直接从文件读出来的数据格式是数组
index1 = list(rpm_qn_na.keys()) # 获取原有csv文件的标题，并形成列表
data = list(map(list, zip(*data))) # map()可以单独列出列表，将数组转换成列表
data = pd.DataFrame(data, index=index1) # 将data的行列转换
data.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_qn_NA_T.csv', header=0)
#处理
t1 = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_qn_NA_T.csv', engine='python')
t1.iloc[0:3,0:3]
t2 = t1.drop('Gene',axis=1)
imputer = KNNImputer(n_neighbors=7) 
imputed = imputer.fit_transform(t2) 
df_imputed = pd.DataFrame(imputed, columns=t2.columns) # 此时无行标题
df_imputed.insert(0,'Gene',t1['Gene'])
type(df_imputed)
df_imputed.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_knn2.csv',index=0)

## log2
rpm_knn2 = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_knn2.csv',engine='python')
rpm_knn2_log2 = rpm_knn2.drop('Gene',axis=1)
rpm_knn2_log2 = np.log2(rpm_knn2_log2)
rpm_knn2_log2.insert(0,'Gene',rpm_knn2['Gene']) 
rpm_knn2_log2.iloc[0:2,0:2]
rpm_knn2_log2.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_log2.csv',index=0)

## z-score
from sklearn import preprocessing
rpm_knn2_log2 = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_log2.csv',engine='python')
rpm_knn2_log2.iloc[0:2,0:2]
rpm_knn2_log2 = rpm_knn2_log2.drop('Gene',axis=1)
values = rpm_knn2_log2.values #dataframe转换为array
values = values.astype('float32') #定义数据类型
data = preprocessing.scale(values) 
df=pd.DataFrame(data) #将array还原为dataframe
df.columns=rpm_knn2_log2.columns #命名标题行
df.iloc[0:2,0:2]
df.insert(0,'Gene',rpm_knn2['Gene']) 
df.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_z-score.csv',index=0) #另存为excel，删除索引

## knn2按列计算均值和方差
a = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_knn2.csv',engine='python')
a.info()
mean = a.mean().sort_values(ascending=False)
std = a.std().sort_values(ascending=False)
mean.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\mean.csv',index=1)
std.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\std.csv',index=1)

# 匹配
std = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\std.csv',engine='python')
mean = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\mean.csv',engine='python')
list_custom = std['name'].values
mean['name'] = mean['name'].astype('category')
mean['name'].cat.reorder_categories(list_custom, inplace=True)
mean.sort_values('name', inplace=True)
mean = mean.set_index('name')
mean.head()
mean.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\mean_sort.csv',index=1)
# a = a.drop('time',axis=1)
a.insert(0,'time',c['time'].values)
a.insert(1,'event',c['event'].values)
b = a













