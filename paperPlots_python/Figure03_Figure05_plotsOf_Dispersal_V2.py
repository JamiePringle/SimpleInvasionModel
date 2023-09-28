from numpy import *
from pylab import *
import multiprocessing as mp

#model code, to run in parallel
def runModel(ageMax,ageRepro,R0,R1,sigma,Ladv,Ngen,sigmaInit):

    #how big is the domain, and how many generations to run
    N=int(1e4)+1  #should be odd

    #when we calculate the invasion speed, we calculate the speed from
    #deltaMore generations after population exceeds threshold to
    #deltaMore+deltaCalc. So we must run for at least
    #deltaMore+deltaCalc+2 generations after the population excees the
    #threshold. Note that the +2 is because the initial condition
    #might exceed the threshold.
    deltaMore=30*3; deltaCalc=15

    xvec=arange(N)-N//2+0.0
    #sigmaInit=1.0 #shape of initial condition
    
    if ageMax>1:
        #Initialize two populations
        P0=1.0*ones((N,ageMax))
        P1=zeros((N,ageMax))
        P1[:,0]=exp(-xvec**2.0/(2*sigmaInit**2))
        #P0=0.0*P0; P0[5000,:]=1.0/ageMax #HACK DEBUGGING

        #breakpoint()
        
        #normalize so total population of both species is one everywhere
        if True: #if true, normalize so sum is 1
            normSum=sum(P0+P1,axis=1)
            for nGen in range(ageMax):
                indx=normSum>0
                P0[indx,nGen]=P0[indx,nGen]/normSum[indx]
                P1[indx,nGen]=P1[indx,nGen]/normSum[indx]

        #breakpoint()
    else:
        #iteroparity
        #Initialize two populations
        P0=1.0*ones((N,))
        P1=zeros((N,))
        P1=exp(-xvec**2.0/(2*sigmaInit**2))

        #normalize so total population is one everywhere
        if True: #if true, normalize so sum is 1
            normSum=P0+P1
            P0=P0/normSum
            P1=P1/normSum
                
    #make kernel
    #sigma=9.0
    xK=arange(N+0.0)-N/2.0+0.5 #this is symmetric around center
    K=(1.0/(sigma*sqrt(2*pi)))*exp(-0.5*((xK+Ladv)**2.0)/sigma**2.0)

    #what is relative growth rate of R1?
    Rcalc=R1/R0

    threshLocPlus=zeros((Ngen,))+nan
    threshLocMinus=zeros((Ngen,))+nan
    popInt=zeros((Ngen,))+nan

    #if for, run Ngen generations. if until, run until we have enough
    #of the invader above the threshold for invasion. Note that if growth is too small, this could be forever...
    for n in range(Ngen):

        #find where population a point where P1> pThresh; this is invasion front
        pThresh=0.05
        try:
            #track youngest generation
            if ageMax>1:
                threshLocPlus[n]=amax(xK[P1[:,0]>pThresh]) #this is slow and sloppy
                threshLocMinus[n]=amin(xK[P1[:,0]>pThresh]) #this is slow and sloppy
            else:
                threshLocPlus[n]=amax(xK[P1>pThresh]) #this is slow and sloppy
                threshLocMinus[n]=amin(xK[P1>pThresh]) #this is slow and sloppy
        except ValueError:
            threshLocPlus[n]=nan
            threshLocMinus[n]=nan

        #what is total population of P1? assume dx=1.0
        popInt[n]=sum(P1.flatten())

        if ageMax>1: #population with age structure

            #all generations older than ageRepro contribute equally R0
            #or R1 settlers to next generation
            P0next=R0*convolve(sum(P0[:,ageRepro:],axis=1),K,mode='same')
            P1next=R1*convolve(sum(P1[:,ageRepro:],axis=1),K,mode='same')

            if False: #print actual acheived growth
                P0growth=sum(P0next)/sum(sum(P0[:,ageRepro:],axis=1))
                P1growth=sum(P1next)/sum(sum(P1[:,ageRepro:],axis=1))
                print('    for %d, acheived growth for P0 and P1 are',P0growth,P1growth)
                #breakpoint()

            
            #now age advance the populations by one generation
            P0[:,1:]=P0[:,:-1]
            P1[:,1:]=P1[:,:-1]

            #now put new recruits into generation 0
            #now limit population if total population>1.0

            #calculate populations of adults who are still alive (generations 1 to nAge-1)
            totalAlive=P0[:,1:].sum(axis=1)+P1[:,1:].sum(axis=1)

            #calculate number of new kids the habat can support, assuming carrying capacity is 1
            totalEmpty=1.0-totalAlive

            #figure out where the settlers for the next generation are more
            #numerous than the space left in the carrying capacity...
            totalKids=P0next+P1next
            indxAllCanSettle=totalKids<=totalEmpty
            indxPopHabLimited=totalKids>totalEmpty

            #breakpoint()
            
            #settlement where no habitat limitation
            P0[indxAllCanSettle,0]=P0next[indxAllCanSettle]
            P1[indxAllCanSettle,0]=P1next[indxAllCanSettle]

            #setlement where habitat limited
            P0[indxPopHabLimited,0]=P0next[indxPopHabLimited]/totalKids[indxPopHabLimited]*totalEmpty[indxPopHabLimited]
            P1[indxPopHabLimited,0]=P1next[indxPopHabLimited]/totalKids[indxPopHabLimited]*totalEmpty[indxPopHabLimited]

        else: #iteroparity
    
            P0next=R0*convolve(P0,K,mode='same')
            P1next=R1*convolve(P1,K,mode='same')

            #now limit population if total population>1.0
            normSum=P0next+P1next
            indx=normSum<1.0
            normSum[indx]=1.0 # no population limit if normSum<1.0, the carrying capacity
            P0=P0next/normSum
            P1=P1next/normSum

    return P0,P1


#plot as f(R1) for various ageRepro NOTE WELL, CRITICAL GROWTH RATE
#FOR NET POP GROWTH IS (1/(ageMax-ageRepro+1)) because
#(ageMax-ageRepro+1) is the number of generations an individual will grow
if __name__ == "__main__":
    figure(3,figsize=(10,8))
    style.use('ggplot')
    clf()

    margins=0.075
    subplots_adjust(wspace=0.2,hspace=0.2,left=margins,right=1.0-margins,
                    bottom=margins,top=1.0-margins)

    #make subplots
    ax1=subplot(2,2,1)
    ax2=subplot(2,2,2)
    ax3=subplot(2,2,3)
    ax4=subplot(2,2,4)

    figure(4,figsize=(10,8))
    style.use('ggplot')
    clf()
    margins=0.075
    subplots_adjust(wspace=0.2,hspace=0.2,left=margins,right=1.0-margins,
                    bottom=margins,top=1.0-margins)
    #make subplots
    fig2ax1=subplot(2,2,1)

    

    #to plot, we need to specify ageMax,ageRepro,R0,R1,sigma,Ladv,Ngen
    
    #first, plot La=0, sigma=10.0 for 3 ages, ageMax=1
    ageMax=1
    ageRepro=0
    R0=3.0
    R1=R0*3.0
    sigma=15.0
    Ladv=0.0
    Ngen=1
    sigmaInit=1.0

    #=========================================================================
    figure(3)
    lineWidth=0.0
    for Rprime in array([3,4,5])*1.2:
        R1=R0*Rprime
        Ngen=2
        P0,P1=runModel(ageMax,ageRepro,R0,R1,sigma,Ladv,Ngen,sigmaInit)
        Nlen=P0.shape[0];xvec=arange(Nlen)-Nlen//2
        sca(ax1)
        plot(xvec,P0,'r')
        lineWidth=lineWidth+0.75
        plot(xvec,P1,'b',linewidth=lineWidth,label=r"R$'$=%4.2f"%(Rprime,))
        axis(xmin=-110,xmax=110)

        ax=gca()
        ax.set_xticklabels([])
        #ax.set_yticklabels([])
        
        #now find the place where the total population of P1 exceeds
        #Pthresh, and mark those locations
        Pthresh=0.2
        Ptotal=P1
        whereCross=where(((Ptotal[1:]-Pthresh)*(Ptotal[:-1]-Pthresh))<0.0)[0]#<0 where cross
        if len(whereCross)==2: #threshold was crossed
            for w in whereCross:
                plot([xvec[w],xvec[w]],[0,Pthresh],'g-',alpha=0.6)
                plot([xvec[w],xvec[w]],[0,Pthresh],'g-',alpha=0.6)
                plot(xvec[w],Pthresh,'go',markersize=7,alpha=0.6)


    legend(fontsize=11)
    plot(xvec,xvec*0.0+Pthresh,'k-',alpha=0.5)
    title(r"3 Growths R$'$, constant L$_{diff}$",fontsize='large')
    ylabel('Population density')

         
    #=========================================================================
    figure(3)
    lineWidth=0.0
    for sigma in [5.0,10.0,25.0]:
        Rprime=4*1.2
        R1=R0*Rprime
        Ngen=2
        P0,P1=runModel(ageMax,ageRepro,R0,R1,sigma,Ladv,Ngen,sigmaInit)
        Nlen=P0.shape[0];xvec=arange(Nlen)-Nlen//2
        sca(ax2)
        #plot(xvec,P0,'r')
        lineWidth=lineWidth+0.75
        plot(xvec,P1,'b',linewidth=lineWidth,label=r"L$_{diff}$=%4.2f"%(sigma,))
        axis(xmin=-110,xmax=110)

        ax=gca()
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        
        #now find the place where the total population of P1 exceeds
        #Pthresh, and mark those locations
        Pthresh=0.2
        Ptotal=P1
        whereCross=where(((Ptotal[1:]-Pthresh)*(Ptotal[:-1]-Pthresh))<0.0)[0]#<0 where cross
        if len(whereCross)==2: #threshold was crossed
            for w in whereCross:
                plot([xvec[w],xvec[w]],[0,Pthresh],'g-',alpha=0.6)
                plot([xvec[w],xvec[w]],[0,Pthresh],'g-',alpha=0.6)
                plot(xvec[w],Pthresh,'go',markersize=7,alpha=0.6)


    legend(fontsize=11)
    plot(xvec,xvec*0.0+Pthresh,'k-',alpha=0.5)
    title(r"3 Dispersals L$_{diff}$, constant R$'$",fontsize='large')
    #ylabel('Population density')
         

    #=========================================================================
    figure(4)
    lineWidth=0.0
    for jnk in [5.0]:
        ageMax=5
        sigma=5.0*2
        Rprime=4*1.2
        R1=R0*Rprime
        Ngen=4
        P0,P1=runModel(ageMax,ageRepro,R0,R1,sigma,Ladv,Ngen,sigmaInit)
        Nlen=P0.shape[0];xvec=arange(Nlen)-Nlen//2
        sca(fig2ax1)
        #plot(xvec,P0,'r')
        for age in range(ageMax):
            lineWidth=lineWidth+0.5
            plot(xvec,P1[:,age],'b',linewidth=lineWidth,label=r"age=%2d"%(ageMax-age-1,))
            plot(xvec,P0[:,age],'r',linewidth=lineWidth)
            axis(xmin=-110,xmax=110)

        ax=gca()
        #ax.set_xticklabels([])
        #ax.set_yticklabels([])
        
        #now find the place where the total population of P1 exceeds
        #Pthresh, and mark those locations
        Pthresh=0.2
        Ptotal=P1
        whereCross=where(((Ptotal[1:]-Pthresh)*(Ptotal[:-1]-Pthresh))<0.0)[0]#<0 where cross
        if len(whereCross)==2: #threshold was crossed
            for w in whereCross:
                plot([xvec[w],xvec[w]],[0,Pthresh],'g-',alpha=0.6)
                plot([xvec[w],xvec[w]],[0,Pthresh],'g-',alpha=0.6)
                plot(xvec[w],Pthresh,'go',markersize=7,alpha=0.6)


    legend(fontsize=11)
    #plot(xvec,xvec*0.0+Pthresh,'k-',alpha=0.5)
    title(r"Age structure of species with 5 generation lifespan",fontsize='large')
    ylabel('Population density')
    xlabel('Alongshore distance')
         
    #=========================================================================
    #=========================================================================
    figure(3)
    lineWidth=0.0
    for ageMax in [1,3,5]:
        #ageMax=5
        sigma=5.0
        Ladv=5.0
        Rprime=3
        R1=R0*Rprime
        Ngen=4+4
        P0,P1=runModel(ageMax,ageRepro,R0,R1,sigma,Ladv,Ngen,sigmaInit)
        Nlen=P0.shape[0];xvec=arange(Nlen)-Nlen//2
        sca(ax4)
        #plot(xvec,P0,'r')
        lineWidth=lineWidth+0.5

        if ageMax==1:
            Ptotal=P1
        else:
            Ptotal=P1.sum(axis=1)

        plot(xvec,Ptotal,'b',linewidth=lineWidth,label=r"lifespan=%2d"%(ageMax,))
        #plot(xvec,P0[:,age],'r',linewidth=lineWidth)
        axis(xmin=-110,xmax=110)

        ax=gca()
        #ax.set_xticklabels([])
        #ax.set_yticklabels([])

        #now find the place where the total population of P1 exceeds
        #Pthresh, and mark those locations
        Pthresh=0.2
        whereCross=where(((Ptotal[1:]-Pthresh)*(Ptotal[:-1]-Pthresh))<0.0)[0]#<0 where cross
        if len(whereCross)==2: #threshold was crossed
            for w in whereCross:
                plot([xvec[w],xvec[w]],[0,Pthresh],'g-',alpha=0.6)
                plot([xvec[w],xvec[w]],[0,Pthresh],'g-',alpha=0.6)
                plot(xvec[w],Pthresh,'go',markersize=7,alpha=0.6)


    legend(fontsize=11)
    plot(xvec,xvec*0.0+Pthresh,'k-',alpha=0.5)
    title(r"Different lifespans, with L$_{adv}$",fontsize='large')
    ylabel('Population density')
    xlabel('Alongshore distance')

    draw()
    show()

    #=========================================================================================
    figure(3)

    lineWidth=0.0
    for ageMax in [1,3,5]:
        #ageMax=5
        #sigma=5.0
        Ladv=0.0
        Rprime=3
        R1=R0*Rprime
        Ngen=4
        P0,P1=runModel(ageMax,ageRepro,R0,R1,sigma,Ladv,Ngen,sigmaInit)
        Nlen=P0.shape[0];xvec=arange(Nlen)-Nlen//2
        sca(ax3)
        #plot(xvec,P0,'r')
        lineWidth=lineWidth+0.5

        if ageMax==1:
            Ptotal=P1
        else:
            Ptotal=P1.sum(axis=1)
        
        plot(xvec,Ptotal,'b',linewidth=lineWidth,label=r"lifespan=%2d"%(ageMax,))
        #plot(xvec,P0[:,age],'r',linewidth=lineWidth)
        axis(xmin=-110,xmax=110)

        ax=gca()
        #ax.set_xticklabels([])
        #ax.set_yticklabels([])
        
        #now find the place where the total population of P1 exceeds
        #Pthresh, and mark those locations
        Pthresh=0.2
        whereCross=where(((Ptotal[1:]-Pthresh)*(Ptotal[:-1]-Pthresh))<0.0)[0]#<0 where cross
        if len(whereCross)==2: #threshold was crossed
            for w in whereCross:
                plot([xvec[w],xvec[w]],[0,Pthresh],'g-',alpha=0.6)
                plot([xvec[w],xvec[w]],[0,Pthresh],'g-',alpha=0.6)
                plot(xvec[w],Pthresh,'go',markersize=7,alpha=0.6)


    legend(fontsize=11)
    plot(xvec,xvec*0.0+Pthresh,'k-',alpha=0.5)
    title(r"Different lifespans",fontsize='large')
    ylabel('Population density')
    xlabel('Alongshore distance')
         
    

    figure(3)
    draw()
    show()
    savefig('Figure03_figInvasionParams.svg')

    figure(4)
    draw()
    show()
    savefig('Figure05_figMultiGens.svg')

