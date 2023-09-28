from pylab import *
from numpy import *
import pandas as ps
import seaborn as sns

fileName='ShanksRescueData/LifeHistory_greater_1km_withLifetime.csv'
data=ps.read_csv(fileName)

data['distMean']=0.5*(data['distMin']+data['distMax'])
data['durationMean']=0.5*(data['durationMin']+data['durationMax'])

indx=logical_and(isfinite(data.durationMean),isfinite(data.distMean))
dataFilt=data[indx]

#if true, remove long distance outlier
if False:
    indx=data.durationMean<150.0
    dataFilt=dataFilt[indx]

durationVec=linspace(0.0,amax(dataFilt.distMean),200)


#make a column normalized by reduction of Cinvasion relative to
#Cinvasion for lifespan of 1 and a given Rprime. adjust table below is from calculate_Cinvasion_module.py

if True:
    #for Rprime=2.0
    adjTable=array([
        [ 1 , 1.0 , 1.162962962962963 ],
        [ 2 , 0.7197452229299363 , 0.837037037037037 ],
        [ 3 , 0.5668789808917197 , 0.6592592592592593 ],
        [ 4 , 0.464968152866242 , 0.5407407407407407 ],
        [ 5 , 0.39490445859872614 , 0.4592592592592593 ],
        [ 6 , 0.34394904458598724 , 0.4 ],
        [ 7 , 0.29936305732484075 , 0.34814814814814815 ],
        [ 8 , 0.27388535031847133 , 0.31851851851851853 ],
        [ 9 , 0.2484076433121019 , 0.2888888888888889 ],
        [ 10 , 0.2229299363057325 , 0.2592592592592593 ],
        [ 11 , 0.21019108280254778 , 0.24444444444444446 ],
        [ 12 , 0.19108280254777069 , 0.2222222222222222 ],
        [ 13 , 0.17834394904458598 , 0.2074074074074074 ],
        [ 14 , 0.17197452229299362 , 0.2 ],
        [ 15 , 0.1592356687898089 , 0.1851851851851852 ],
        [ 16 , 0.15286624203821655 , 0.17777777777777778 ],
        [ 17 , 0.1464968152866242 , 0.1703703703703704 ],
    ])
else:
    #for Rprime=1.3
    adjTable=array([
        [ 1 , 1.0 , 0.7037037037037037 ],
        [ 2 , 0.6947368421052632 , 0.48888888888888893 ],
        [ 3 , 0.5157894736842106 , 0.36296296296296293 ],
        [ 4 , 0.43157894736842106 , 0.3037037037037037 ],
        [ 5 , 0.35789473684210527 , 0.2518518518518518 ],
        [ 6 , 0.31578947368421056 , 0.2222222222222222 ],
        [ 7 , 0.28421052631578947 , 0.2 ],
        [ 8 , 0.25263157894736843 , 0.17777777777777778 ],
        [ 9 , 0.23157894736842105 , 0.16296296296296295 ],
        [ 10 , 0.21052631578947367 , 0.14814814814814814 ],
        [ 11 , 0.18947368421052632 , 0.13333333333333333 ],
        [ 12 , 0.17894736842105263 , 0.1259259259259259 ],
        [ 13 , 0.16842105263157894 , 0.11851851851851852 ],
        [ 14 , 0.16842105263157894 , 0.11851851851851852 ],
        [ 15 , 0.15789473684210528 , 0.1111111111111111 ],
        [ 16 , 0.1473684210526316 , 0.1037037037037037 ],
        [ 17 , 0.1473684210526316 , 0.1037037037037037 ],
    ])

lifetime=array(data['lifetime_years'][:])
distMean=array(data['distMean'][:])
distMeanAdj=array(data['distMean'][:]+nan)

for n in arange(len(lifetime)):
    thisLife=min(amax(adjTable[:,0]),lifetime[n])
    if isfinite(thisLife):
        #adjFac=adjTable[int(thisLife)-1,1] #adjust by Ci/Ci[0]
        adjFac=adjTable[int(thisLife)-1,2] #adjust by Ci/sigma
        distMeanAdj[n]=distMean[n]/adjFac
data['distMeanAdj']=distMeanAdj

#now fitler the data
if True:
    indx=data.durationMean<150.0
    dataFilt=data[indx]


#set up figures
#==================================================================================
figure(1,figsize=(8,8))
clf()
style.use('ggplot')
ax2=subplot(1,1,1)


sca(ax2)
sns.scatterplot(ax=ax2,data=dataFilt, x="durationMean", y="distMean", hue="howMeasure",
                style='howMeasure')
axisNow=axis()
plot(durationVec,durationVec*8.64e4*0.30/1e3,'r--',alpha=0.95)
plot(durationVec,durationVec*8.64e4*0.10/1e3,'r--',alpha=0.95)
plot(durationVec,durationVec*8.64e4*0.031/1e3,'k--',alpha=0.5)
axisNow=array(axisNow)
axisNow[0]=-5.0
axis(axisNow)
fontsize='xx-large'
xlabel('PLD, days',fontsize=fontsize)
ylabel("Shank's estimate of dispersal distance, km",fontsize=fontsize)
title('PLD < 150 days; R=%4.4f'%corrcoef(dataFilt.durationMean,dataFilt.distMean)[1,0],fontsize=fontsize)

axis(ymax=1000, ymin=-20)

draw()
show()
savefig('Figure01_ShanksData_raw.svg')

#==================================================================================
figure(2,figsize=(8,8))
clf()
style.use('ggplot')
ax2=subplot(1,1,1)

sca(ax2)
sns.scatterplot(ax=ax2,data=dataFilt, x="durationMean", y="distMeanAdj", hue="howMeasure",
                style='howMeasure')
axisNow=axis()
plot(durationVec,durationVec*8.64e4*0.30/1e3,'r--',alpha=0.95)
plot(durationVec,durationVec*8.64e4*0.10/1e3,'r--',alpha=0.95)
plot(durationVec,durationVec*8.64e4*0.031/1e3,'k--',alpha=0.5)
axisNow=array(axisNow)
axisNow[0]=-5.0
axis(axisNow)
xlabel('PLD, days',fontsize=fontsize)
ylabel("adjusted Shank's estimate of dispersal distance, km",fontsize=fontsize)
title('PLD < 150 days; R=%4.4f'%corrcoef(dataFilt.durationMean,dataFilt.distMeanAdj)[1,0],fontsize=fontsize)

axis(ymax=1000,ymin=-20)

draw()
show()
savefig('Figure01_ShanksData_adjusted.svg')

