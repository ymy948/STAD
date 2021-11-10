# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 15:17:43 2021

@author: DELL
"""
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from lifelines.datasets import load_rossi
# rossi_dataset = load_rossi()
# rossi_dataset.head()
# a = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\cox.csv',engine='python')
# a.iloc[0:2,0:2]
# b = a.drop('Gene',axis=1)
# b.iloc[0:2,0:2]

# cph = CoxPHFitter(penalizer=0.005,alpha=0.05)
# cph.fit(b, duration_col='time', event_col='status',step_size=2)  
# cph.plot()
# c = cph.summary
# c.head()
# c[c['p'] <= 0.05]
##################################在all_clin中提取mir的样本########################################
c = pd.read_csv('E:\sirebrowser\STAD\clinical\\all-clinical.csv',engine='python')
z = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\STAD.miRseq_mature_RPM_z-score.csv',engine='python')
c.iloc[0:2,0:2]
z.iloc[0:2,0:2]
d.iloc[0:2,0:2]
z.info()
c.info()
d = c.set_index(['ID'],drop=False)
for i in range(0,389): # 根据报错修改数值
    if list(d['ID'])[i] in list(z['Gene']):
        continue
    else:
        d.drop(index=[list(d['ID'])[i]],inplace=True)
      
d.to_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\mir_clin.csv',index=0)

## 准备coxdata直接在excel中合并


targets = pd.read_csv('E:\\sirebrowser\\STAD\\miR\\分析\\New-CESC-COX-diff_list.csv',engine='python')
p_value = pd.read_csv('E:\\sirebrowser\\STAD\\miR\\分析\\All_P_value_CESC-COXPH.csv',engine='python')
hr = pd.read_csv('E:\\sirebrowser\\STAD\\miR\\分析\\All_HR_CESC-COXPH.csv',engine='python')

## 在p_value中筛选targets
p_value.iloc[0:2,0:2]
p_value.info()
targets.iloc[1,0:2]
targets.info()
len(p_value.columns)
len(targets.columns)
for i in range(0,73):
    if p_value.columns[i] in targets.columns:
        continue
    else:
        p_value = p_value.drop(p_value.columns[i],axis=1)

p_value.to_csv('E:\\sirebrowser\\STAD\\miR\\分析\\All_P_value_CESC-COXPH_1.csv',index=0)

## 在hr中筛选targets
hr.iloc[0:2,0:2]
hr.info()
for i in range(0,73): # 根据报错一直修改range的第二个值
    if hr.columns[i] in targets.columns:
        continue
    else:
        hr = hr.drop(hr.columns[i],axis=1)

hr.to_csv('E:\\sirebrowser\\STAD\\miR\\分析\\All_HR_CESC-COXPH_1.csv',index=0)

## 生存分析画图
# 数据准备 在z-score中筛选生存相关gene
z_score = pd.read_csv('E:\sirebrowser\STAD\miR\miRseq_Mature_Preprocess\cox_dot.csv',engine='python')
z_score.info()
z_score.iloc[0:2,0:2]
z_score_1 = z_score.iloc[:,1:467]
z_score_1.info()
z_score_1.iloc[0:2,0:2]
for i in range(0,73): # 根据报错一直修改range的第二个值
    if z_score_1.columns[i] in targets.columns:
        continue
    else:
        z_score_1 = z_score_1.drop(z_score_1.columns[i],axis=1)

z_score_1.insert(0,'Gene',z_score['Gene'])         
z_score_1.to_csv('E:\sirebrowser\STAD\miR\分析\\rpm筛选.csv',index=0)


## km
import lifelines
from lifelines import KaplanMeierFitter
from lifelines.utils import median_survival_times
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
help(KaplanMeierFitter)
lifelines.__version__

from lifelines.statistics import logrank_test
font={'color': 'k',
      'size': 15,
      'family': 'Times New Roman'}

k = 3 #聚类数
iteration = 500 

data = pd.read_csv('E:\\sirebrowser\\STAD\\miR\\分析\\26k=3.csv', engine='python',index_col = 'Sample') 

data.info()
data = data.iloc[:,0:25]
model = KMeans(n_clusters = k, n_jobs = 4, max_iter = iteration)
model.fit(data) 

#输出
r = pd.concat([data, pd.Series(model.labels_, index = data.index)], axis = 1)  

r.columns = list(data.columns) + [u'kgroups25'] #重命名表头
r.iloc[0:2,0:2]
r.insert(0,'Sample',data.index)
r.to_csv('E:\\sirebrowser\\STAD\\miR\\分析\\26k=3importance_25.csv',index=0) 
########################################k=2################################
# 3k=2
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\3k=2.csv',engine='python')
a.info()
b = a.iloc[:,[1,2,3]]             
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(show_censors=True,ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(show_censors=True,ax=ax,ci_show=False,color='#F64B35FF')
plt.legend(prop={'family' : 'Arial', 'size'   : 22},handletextpad=0.5,frameon=False,labelspacing=0.1)
plt.tick_params(width=2)
ax.spines['bottom'].set_linewidth('2')
ax.spines['top'].set_linewidth('0')
ax.spines['left'].set_linewidth('2')
ax.spines['right'].set_linewidth('0')
plt.xlim(-230,3900,500)
plt.ylim(-0.08,1.08)
plt.text(0, 0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
plt.title('k=2 top74', fontdict=fontt)
plt.xlabel(' ')
plt.ylabel('', fontdict=fonty)
plt.yticks(fontproperties = 'Arial', size = 20)
plt.xticks(fontproperties = 'Arial', size = 20,rotation=20)
# plt.figure(dpi=1000,figsize=(10,8))
plt.savefig(fname='E:\\sirebrowser\\STAD\\miR\\分析\\km\\74k=2.png',figsize=[20,20],dpi=500)
plt.close('all')

# 14k=2
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\14k=2.csv',engine='python')
b = a.iloc[:,[14,15,16]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)

kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
plt.title('k=2 top14')
plt.text(0, -0.02 , 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 26k=2
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\26k=2.csv',engine='python')
a.info()
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)

T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)

kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(show_censors=True,ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,show_censors=True,ci_show=False,color='#F64B35FF')

plt.title('k=2 top26')
plt.text(0, -0.02, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 39k=2
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\39k=2.csv',engine='python')
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)


kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')

plt.title('k=2 top39')
plt.text(0, 0.02 , 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 55k=2
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\55k=2.csv',engine='python')
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)

T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)


kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')

plt.title('k=2 top55')
plt.text(0, 0.02, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 64k=2
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\64k=2.csv',engine='python')
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)

T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)


kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')

plt.title('k=2 top64')
plt.text(0, 0.04, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 74k=2
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\74k=2.csv',engine='python')
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)

T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)

kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')

plt.title('k=2 top74')
plt.text(0, 0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

########################################k=3################################
# 3k=3
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\3k=3.csv',engine='python')
a.info()
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)

T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
dem3 = (b['groups'] == 2)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)
results_2.p_value

kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(show_censors=True,ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(show_censors=True,ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(show_censors=True,ci_show=False,color='#7E6148FF'     )
plt.legend(loc=(0.5,0.63),prop={'family' : 'Arial', 'size'   : 22},handletextpad=0.5,frameon=False,labelspacing=0.1)
plt.tick_params(width=2)
ax.spines['bottom'].set_linewidth('2')
ax.spines['top'].set_linewidth('0')     
ax.spines['left'].set_linewidth('2')
ax.spines['right'].set_linewidth('0')
plt.xlim(-230,3900,500)
plt.ylim(-0.08,1.08)
plt.title('k=3 top74', fontdict=fontt)

plt.xlabel(' ')
plt.ylabel(' ', fontdict=fonty)
plt.yticks(fontproperties = 'Arial', size = 20)
plt.xticks(fontproperties = 'Arial', size = 20,rotation=20)
# plt.figure(dpi=1000,figsize=(10,8))
plt.savefig(fname='E:\\sirebrowser\\STAD\\miR\\分析\\km\\74k=3.png',figsize=[20,20],dpi=500)
plt.close('all')
# plt.figure(dpi=1000,figsize=(24,20))

# 14k=3
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\14k=3.csv',engine='python')
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)

T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
dem3 = (b['groups'] == 2)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)


kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(ci_show=False,color='#7E6148FF'     )
plt.title('k=3 top14')
plt.text(0, 0.07, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, 0.01, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results_2.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 21k=3
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\26k=3importance_result_21.csv',engine='python')
b = a.iloc[:,[22,23,24]]
b.head()
kmf = KaplanMeierFitter()
groups = b['kgroups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)

T = b['time']
E = b['event']
dem1 = (b['kgroups'] == 0)
dem2 = (b['kgroups'] == 1)
dem3 = (b['kgroups'] == 2)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)


kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(show_censors=True,ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,show_censors=True,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(show_censors=True,ci_show=False,color='#7E6148FF'     )
plt.title('k=3 top21')
plt.text(0, 0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, 0, 'Group 0     vs     Group 2'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.05, 'Group 1     vs     Group 2'+"      P_value=%.2e"%results_2.p_value, fontdict=font)
plt.axvline(x=1825,c='k',ls='--',lw=0.5)
plt.axhline(y=0.7,c='k',ls='--',lw=0.5)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 26k=3
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\26k=3.csv',engine='python')
a.info()
b = a.iloc[:,[27,28,29]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)

T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
dem3 = (b['groups'] == 2)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)


kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(show_censors=True,ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,show_censors=True,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(show_censors=True,ci_show=False,color='#7E6148FF'     )
plt.title('k=3 top20')
plt.text(0, 0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, 0, 'Group 0     vs     Group 2'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.05, 'Group 1     vs     Group 2'+"      P_value=%.2e"%results_2.p_value, fontdict=font)
plt.axvline(x=1825,c='k',ls='--',lw=0.5)
plt.axhline(y=0.7,c='k',ls='--',lw=0.5)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 39k=3
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\39k=3.csv',engine='python')
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)

T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
dem3 = (b['groups'] == 2)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)


kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(ci_show=False,color='#7E6148FF'     )
plt.title('k=3 top39')
plt.text(0, 0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, 0, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results_2.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 55k=3
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\55k=3.csv',engine='python')
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)

T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
dem3 = (b['groups'] == 2)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)


kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(ci_show=False,color='#7E6148FF'     )
plt.title('k=3 top55')
plt.text(0, 0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, 0, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results_2.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 64k=3
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\64k=3.csv',engine='python')
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)

T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
dem3 = (b['groups'] == 2)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)


kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(ci_show=False,color='#7E6148FF'     )
plt.title('k=3 top64')
plt.text(0, 0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, 0, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results_2.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 74k=3
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\74k=3.csv',engine='python')
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)

T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
dem3 = (b['groups'] == 2)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)


kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(ci_show=False,color='#7E6148FF'     )
plt.title('k=3 top74')
plt.text(0, 0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, 0, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results_2.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

########################################k=4################################
# 3k=4
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\3k=4.csv',engine='python')
a.info()
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)
ix4 = (groups == 3)


T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
dem3 = (b['groups'] == 2)
dem4 = (b['groups'] == 3)

kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(show_censors=True,ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(show_censors=True,ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(show_censors=True,ci_show=False,color='#7E6148FF'     )
kmf.fit(b['time'][ix4], b['event'][ix4], label='Group 3') 
ax = kmf.plot(show_censors=True,ci_show=False,color='#3C5488FF'     )
plt.legend(loc=(0.5,0.55),prop={'family' : 'Arial', 'size'   : 22},handletextpad=0.5,frameon=False,labelspacing=0.1)
plt.title('k=4 top74',fontdict=fontt)
plt.tick_params(width=2)
ax.spines['bottom'].set_linewidth('2')
ax.spines['top'].set_linewidth('0')     
ax.spines['left'].set_linewidth('2')
ax.spines['right'].set_linewidth('0')
plt.xlim(-230,3900,500)
plt.ylim(-0.08,1.08)
plt.xlabel(' ')
plt.ylabel(' ', fontdict=fonty)
plt.yticks(fontproperties = 'Arial', size = 20)
plt.xticks(fontproperties = 'Arial', size = 20,rotation=20)
# plt.figure(dpi=1000,figsize=(10,8))
plt.savefig(fname='E:\\sirebrowser\\STAD\\miR\\分析\\km\\74k=4.png',figsize=[20,20],dpi=500)
plt.close('all')


# 14k=4
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\4k=4.csv',engine='python')
a.info()
# b = a.iloc[:,[40,41,42]]
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)
ix4 = (groups == 3)


T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
dem3 = (b['groups'] == 2)
dem4 = (b['groups'] == 3)

results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem1],T[dem4],E[dem1],E[dem4],alpha=.99)
results_3 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)
results_4 = logrank_test(T[dem2],T[dem4],E[dem2],E[dem4],alpha=.99)
results_5 = logrank_test(T[dem3],T[dem4],E[dem3],E[dem4],alpha=.99)



kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(ci_show=False,color='#7E6148FF'     )
kmf.fit(b['time'][ix4], b['event'][ix4], label='Group 3') 
ax = kmf.plot(ci_show=False,color='#3C5488FF'     )
plt.title('k=4 top14')
plt.text(0, -0.07, 'Group 0     vs     Group 1'+"      P_value=%.6f"%results.p_value, fontdict=font)
plt.text(0, -0.13, 'Group 0     vs     Group 2'+"      P_value=%.6f"%results_1.p_value, fontdict=font)
plt.text(0, -0.19, 'Group 0     vs     Group 3'+"      P_value=%.6f"%results_2.p_value, fontdict=font)
plt.text(0, -0.25, 'Group 1     vs     Group 2'+"      P_value=%.6f"%results_3.p_value, fontdict=font)
plt.text(0, -0.31, 'Group 1     vs     Group 3'+"      P_value=%.6f"%results_4.p_value, fontdict=font)
plt.text(0, -0.37, 'Group 2     vs     Group 3'+"      P_value=%.6f"%results_5.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 26k=4
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\26k=4.csv',engine='python')
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)
ix4 = (groups == 3)


T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
dem3 = (b['groups'] == 2)
dem4 = (b['groups'] == 3)

results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem1],T[dem4],E[dem1],E[dem4],alpha=.99)
results_3 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)
results_4 = logrank_test(T[dem2],T[dem4],E[dem2],E[dem4],alpha=.99)
results_5 = logrank_test(T[dem3],T[dem4],E[dem3],E[dem4],alpha=.99)



kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(ci_show=False,color='#7E6148FF'     )
kmf.fit(b['time'][ix4], b['event'][ix4], label='Group 3') 
ax = kmf.plot(ci_show=False,color='#3C5488FF'     )
plt.title('k=4 top26')
plt.text(0, 0.01, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, -0.05, 'Group 0     vs     Group 2'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.11, 'Group 0     vs     Group 3'+"      P_value=%.2e"%results_2.p_value, fontdict=font)
plt.text(0, -0.17, 'Group 1     vs     Group 2'+"      P_value=%.2e"%results_3.p_value, fontdict=font)
plt.text(0, -0.23, 'Group 1     vs     Group 3'+"      P_value=%.2e"%results_4.p_value, fontdict=font)
plt.text(0, -0.29, 'Group 2     vs     Group 3'+"      P_value=%.2e"%results_5.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 39k=4
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\39k=4.csv',engine='python')
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)
ix4 = (groups == 3)


T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
dem3 = (b['groups'] == 2)
dem4 = (b['groups'] == 3)

results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem1],T[dem4],E[dem1],E[dem4],alpha=.99)
results_3 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)
results_4 = logrank_test(T[dem2],T[dem4],E[dem2],E[dem4],alpha=.99)
results_5 = logrank_test(T[dem3],T[dem4],E[dem3],E[dem4],alpha=.99)



kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(ci_show=False,color='#7E6148FF'     )
kmf.fit(b['time'][ix4], b['event'][ix4], label='Group 3') 
ax = kmf.plot(ci_show=False,color='#3C5488FF'     )
plt.title('k=4 top39')
plt.text(0, -0.03, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, -0.09, 'Group 0     vs     Group 2'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.15, 'Group 0     vs     Group 3'+"      P_value=%.2e"%results_2.p_value, fontdict=font)
plt.text(0, -0.21, 'Group 1     vs     Group 2'+"      P_value=%.2e"%results_3.p_value, fontdict=font)
plt.text(0, -0.27, 'Group 1     vs     Group 3'+"      P_value=%.2e"%results_4.p_value, fontdict=font)
plt.text(0, -0.33, 'Group 2     vs     Group 3'+"      P_value=%.2e"%results_5.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

#combine
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\39k=4.csv',engine='python')
b = a.iloc[:,[1,2,5]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups2']
ix1 = (groups == '0+1')
ix2 = (groups == '2')
ix3 = (groups == '3')

T = b['time']
E = b['event']
dem1 = (b['groups2'] == '0+1')
dem2 = (b['groups2'] == '2')
dem3 = (b['groups2'] == '3')

results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_3 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)

kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0+1')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 2')
ax = kmf.plot(ax=ax,ci_show=False,color='#7E6148FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 3') 
ax = kmf.plot(ci_show=False,color='#F64B35FF'     )
plt.title('k=4 top39')
plt.text(0, -0.04, 'Group 0+1         vs     Group 2'+"          P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, -0.1, 'Group 0+1         vs     Group 3'+"          P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.16, 'Group 2             vs     Group 3'+"          P_value=%.2e"%results_3.p_value, fontdict=font)
plt.axvline(x=365,c='k',ls='--',lw=0.5)
plt.axvline(x=730,c='k',ls='--',lw=0.5)
plt.axvline(x=1095,c='k',ls='--',lw=0.5)
plt.axvline(x=1825,c='k',ls='--',lw=0.5)
plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 55k=4
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\55k=4.csv',engine='python')
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)
ix4 = (groups == 3)


T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
dem3 = (b['groups'] == 2)
dem4 = (b['groups'] == 3)

results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem1],T[dem4],E[dem1],E[dem4],alpha=.99)
results_3 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)
results_4 = logrank_test(T[dem2],T[dem4],E[dem2],E[dem4],alpha=.99)
results_5 = logrank_test(T[dem3],T[dem4],E[dem3],E[dem4],alpha=.99)



kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(ci_show=False,color='#7E6148FF'     )
kmf.fit(b['time'][ix4], b['event'][ix4], label='Group 3') 
ax = kmf.plot(ci_show=False,color='#3C5488FF'     )
plt.title('k=4 top55')
plt.text(0, -0.03, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, -0.09, 'Group 0     vs     Group 2'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.15, 'Group 0     vs     Group 3'+"      P_value=%.2e"%results_2.p_value, fontdict=font)
plt.text(0, -0.21, 'Group 1     vs     Group 2'+"      P_value=%.2e"%results_3.p_value, fontdict=font)
plt.text(0, -0.27, 'Group 1     vs     Group 3'+"      P_value=%.2e"%results_4.p_value, fontdict=font)
plt.text(0, -0.33, 'Group 2     vs     Group 3'+"      P_value=%.2e"%results_5.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

#combine
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\55k=4.csv',engine='python')
b = a.iloc[:,[1,2,5]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups2']
ix1 = (groups == '0+2')
ix2 = (groups == '1')
ix3 = (groups == '3')

T = b['time']
E = b['event']
dem1 = (b['groups2'] == '0+2')
dem2 = (b['groups2'] == '1')
dem3 = (b['groups2'] == '3')

results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_3 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)

kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0+2')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 3') 
ax = kmf.plot(ci_show=False,color='#7E6148FF'     )
plt.title('k=4 top55')
plt.text(0, -0.04, 'Group 0+2        vs     Group 1'+"          P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, -0.1, 'Group 0+2        vs     Group 3'+"          P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.16, 'Group 1            vs     Group 3'+"          P_value=%.2e"%results_3.p_value, fontdict=font)
plt.axvline(x=365,c='k',ls='--',lw=0.5)
plt.axvline(x=730,c='k',ls='--',lw=0.5)
plt.axvline(x=1095,c='k',ls='--',lw=0.5)
plt.axvline(x=1825,c='k',ls='--',lw=0.5)
plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

# 64k=4
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\64k=4.csv',engine='python')
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)
ix4 = (groups == 3)


T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
dem3 = (b['groups'] == 2)
dem4 = (b['groups'] == 3)

results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem1],T[dem4],E[dem1],E[dem4],alpha=.99)
results_3 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)
results_4 = logrank_test(T[dem2],T[dem4],E[dem2],E[dem4],alpha=.99)
results_5 = logrank_test(T[dem3],T[dem4],E[dem3],E[dem4],alpha=.99)



kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(ci_show=False,color='#7E6148FF'     )
kmf.fit(b['time'][ix4], b['event'][ix4], label='Group 3') 
ax = kmf.plot(ci_show=False,color='#3C5488FF'     )
plt.title('k=4 top64')
plt.text(0, 0.07, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, 0.01, 'Group 0     vs     Group 2'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.05, 'Group 0     vs     Group 3'+"      P_value=%.2e"%results_2.p_value, fontdict=font)
plt.text(0, -0.11, 'Group 1     vs     Group 2'+"      P_value=%.2e"%results_3.p_value, fontdict=font)
plt.text(0, -0.17, 'Group 1     vs     Group 3'+"      P_value=%.2e"%results_4.p_value, fontdict=font)
plt.text(0, -0.23, 'Group 2     vs     Group 3'+"      P_value=%.2e"%results_5.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

#combine
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\64k=4.csv',engine='python')
b = a.iloc[:,[1,2,5]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups2']
ix1 = (groups == '0+2')
ix2 = (groups == '1')
ix3 = (groups == '3')

T = b['time']
E = b['event']
dem1 = (b['groups2'] == '0+2')
dem2 = (b['groups2'] == '1')
dem3 = (b['groups2'] == '3')

results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_3 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)

kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0+2')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 3') 
ax = kmf.plot(ci_show=False,color='#3C5488FF'     )
plt.title('k=4 top64')
plt.text(0, 0.05, 'Group 0+2        vs     Group 1'+"              P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, -0.01, 'Group 0+2        vs     Group 3'+"          P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.07, 'Group 1        vs     Group 3'+"          P_value=%.2e"%results_3.p_value, fontdict=font)
plt.axvline(x=365,c='k',ls='--',lw=0.5)
plt.axvline(x=730,c='k',ls='--',lw=0.5)
plt.axvline(x=1095,c='k',ls='--',lw=0.5)
plt.axvline(x=1825,c='k',ls='--',lw=0.5)
plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))


# 74k=4
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\74k=4.csv',engine='python')
b = a.iloc[:,[1,2,3]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)
ix4 = (groups == 3)


T = b['time']
E = b['event']
dem1 = (b['groups'] == 0)
dem2 = (b['groups'] == 1)
dem3 = (b['groups'] == 2)
dem4 = (b['groups'] == 3)

results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem1],T[dem4],E[dem1],E[dem4],alpha=.99)
results_3 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)
results_4 = logrank_test(T[dem2],T[dem4],E[dem2],E[dem4],alpha=.99)
results_5 = logrank_test(T[dem3],T[dem4],E[dem3],E[dem4],alpha=.99)



kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(ci_show=False,color='#7E6148FF'     )
kmf.fit(b['time'][ix4], b['event'][ix4], label='Group 3') 
ax = kmf.plot(ci_show=False,color='#3C5488FF'     )
plt.title('k=4 top74')
plt.text(0, 0.01, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, -0.05, 'Group 0     vs     Group 2'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.11, 'Group 0     vs     Group 3'+"      P_value=%.2e"%results_2.p_value, fontdict=font)
plt.text(0, -0.17, 'Group 1     vs     Group 2'+"      P_value=%.2e"%results_3.p_value, fontdict=font)
plt.text(0, -0.23, 'Group 1     vs     Group 3'+"      P_value=%.2e"%results_4.p_value, fontdict=font)
plt.text(0, -0.29, 'Group 2     vs     Group 3'+"      P_value=%.2e"%results_5.p_value, fontdict=font)

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))

#combine
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\74k=4.csv',engine='python')
b = a.iloc[:,[1,2,4]]
b.head()
kmf = KaplanMeierFitter()
groups = b['groups1']
ix1 = (groups == '0')
ix2 = (groups == '1')
ix3 = (groups == '2+3')

T = b['time']
E = b['event']
dem1 = (b['groups1'] == '0')
dem2 = (b['groups1'] == '1')
dem3 = (b['groups1'] == '2+3')

results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_3 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)

kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2+3') 
ax = kmf.plot(ci_show=False,color='#3C5488FF'     )
plt.title('k=4 top74')
plt.text(0, 0.01, 'Group 0        vs     Group 1'+"              P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, -0.05, 'Group 0        vs     Group 2+3'+"          P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.11, 'Group 1        vs     Group 2+3'+"          P_value=%.2e"%results_3.p_value, fontdict=font)
plt.axvline(x=365,c='k',ls='--',lw=0.5)
plt.axvline(x=730,c='k',ls='--',lw=0.5)
plt.axvline(x=1095,c='k',ls='--',lw=0.5)
plt.axvline(x=1825,c='k',ls='--',lw=0.5)
plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)')
plt.figure(dpi=1000,figsize=(24,20))
##################################################
# 26k=3
a = pd.read_csv('E:\sirebrowser\STAD\miR\分析\\26k=3importance_15.csv',engine='python')
a.info()
b = a.iloc[:,[16,17,18]]
b.head()
kmf = KaplanMeierFitter()
groups = b['kgroups15']
ix1 = (groups == 0)
ix2 = (groups == 1)
ix3 = (groups == 2)

T = b['time']
E = b['event']
dem1 = (b['kgroups15'] == 0)
dem2 = (b['kgroups15'] == 1)
dem3 = (b['kgroups15'] == 2)
results = logrank_test(T[dem1],T[dem2],E[dem1],E[dem2],alpha=.99)
results_1 = logrank_test(T[dem1],T[dem3],E[dem1],E[dem3],alpha=.99)
results_2 = logrank_test(T[dem2],T[dem3],E[dem2],E[dem3],alpha=.99)


kmf.fit(b['time'][ix1], b['event'][ix1], label='Group 0')
ax = kmf.plot(show_censors=True,ci_show=False,color='#00A087FF')
kmf.fit(b['time'][ix2], b['event'][ix2], label='Group 1')
ax = kmf.plot(ax=ax,show_censors=True,ci_show=False,color='#F64B35FF')
kmf.fit(b['time'][ix3], b['event'][ix3], label='Group 2') 
ax = kmf.plot(show_censors=True,ci_show=False,color='#7E6148FF'     )
plt.title('k=3 top15(in top26)', fontdict=fontt)
plt.text(0, 0.05, 'Group 0     vs     Group 1'+"      P_value=%.2e"%results.p_value, fontdict=font)
plt.text(0, 0, 'Group 0     vs     Group 2'+"      P_value=%.2e"%results_1.p_value, fontdict=font)
plt.text(0, -0.05, 'Group 1     vs     Group 2'+"      P_value=%.2e"%results_2.p_value, fontdict=font)
# plt.axvline(x=1825,c='k',ls='--',lw=0.5)
# plt.axhline(y=0.7,c='k',ls='--',lw=0.5)

plt.legend(prop={'family' : 'Arial', 'size'   : 13},handletextpad=0.5,frameon=False,labelspacing=0.2)
plt.tick_params(width=2)
ax.spines['bottom'].set_linewidth('2')
ax.spines['top'].set_linewidth('0')
ax.spines['left'].set_linewidth('2')
ax.spines['right'].set_linewidth('0')

plt.xlabel(' ')
plt.ylabel('Cumulative survival (percentage)', fontdict=fonty)
plt.figure(dpi=1000,figsize=(24,20))
















