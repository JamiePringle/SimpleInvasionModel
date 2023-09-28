from numpy import *
from pylab import *
import multiprocessing as mp

#model code, to run in parallel
def runModel(ageMax,ageRepro,R0,R1,sigma,Ladv,Ngen):

    #how big is the domain, and how many generations to run
    N=int(1e4)+1  #should be odd

    #when we calculate the invasion speed, we calculate the speed from
    #deltaMore generations after population exceeds threshold to
    #deltaMore+deltaCalc. So we must run for at least
    #deltaMore+deltaCalc+2 generations after the population excees the
    #threshold. Note that the +2 is because the initial condition
    #might exceed the threshold.
    deltaMore=30*3; deltaCalc=15

    
    if ageMax>1:
        #Initialize two populations
        P0=1.0*ones((N,ageMax))
        P1=zeros((N,ageMax))
        P1[N//2,0]=1.0/ageMax
        #P0=0.0*P0; P0[5000,:]=1.0/ageMax #HACK DEBUGGING

        #breakpoint()
        
        #normalize so total population is one everywhere
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
        P1[N//2]=1.0/2

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
    n=-1; stopLoop=False
    #for n in range(Ngen):
    while not stopLoop:
        n+=1

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

        if sum(isfinite(threshLocPlus))>(deltaMore+deltaCalc+2):
            #we have enough data to quit
            stopLoop=True
            
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

        if False: #debugging plots
            figure(1)
            clf()
            plot(P0,'r')
            plot(sum(P0,axis=1),'r--')
            plot(P1,'b')
            plot(sum(P1,axis=1),'b--')
            title('%d'%(n,))
            draw();show();pause(0.2)
            print('done with',n,sum(P0),sum(P1))
        else:
            #print('done with',n,sum(P0),sum(P1))
            jnk=0

    #draw last time
    if False:
        figure(1)
        clf()
        style.use('ggplot')
        ax1=subplot(3,1,1)
        plot(xK,P0,'r-')
        plot(xK,P1,'-')#many colors?
        plot(threshLocPlus[n],0.0,'g*',markersize=15)
        plot(threshLocMinus[n],0.0,'b*',markersize=15)
        axis(xmin=xK[0],xmax=xK[-1])

        subplot(3,1,2)
        plot(threshLocPlus,label='data')
        plot(threshLocMinus,label='data')
        #plot((sigma*sqrt(2*log(Rcalc)))*arange(Ngen),label='predict') #fisher invasion speed*time for growth Rcalc
        legend()

        subplot(3,1,3)
        plot(popInt)

        draw()
        show()
        pause(0.01)

    #now find the invasion speed from the 10th location after locations become valid.
    threshLocPlus[0]=nan #avoid initial condition
    threshLocMinus[0]=nan #avoid initial condition

    if False:
        minValid=arange(len(threshLocPlus))[(nanmin(threshLocPlus)==threshLocPlus)]
        speedPlus=(threshLocPlus[minValid+deltaMore+deltaCalc]-threshLocPlus[minValid+deltaMore])/deltaCalc
        minValid=arange(len(threshLocMinus))[(nanmax(threshLocMinus)==threshLocMinus)]
        speedMinus=(threshLocMinus[minValid+deltaMore+deltaCalc]-threshLocMinus[minValid+deltaMore])/deltaCalc
    else:
        threshLocValid=threshLocPlus[isfinite(threshLocPlus)]
        speedPlus=(threshLocValid[deltaMore+deltaCalc]-threshLocValid[deltaMore])/deltaCalc

        threshLocValid=threshLocMinus[isfinite(threshLocMinus)]
        speedMinus=(threshLocValid[deltaMore+deltaCalc]-threshLocValid[deltaMore])/deltaCalc

        
    print('For R0,R1 of',R0,R1,'speeds are',speedPlus,speedMinus,flush=True)
    #pause(1.0)
    #breakpoint()
    
    return speedPlus,speedMinus


#plot as f(R1) for various ageRepro NOTE WELL, CRITICAL GROWTH RATE
#FOR NET POP GROWTH IS (1/(ageMax-ageRepro+1)) because
#(ageMax-ageRepro+1) is the number of generations an individual will grow
if __name__ == "__main__":
    figure(3,figsize=(11,8))
    style.use('ggplot')
    clf()

    nTrys=21 #length of parameter vector
    
    #==================================================================================
    subplot(2,1,1)
    colList=['r','g','b','c']
    for ageMax in [1,3,5]: #[0,2,4,5]:
        ageRepro=0
        RVec=logspace(log10(1.02),log10(5.0),nTrys)
        #RVec=array([1.2,1.3,1.4])
        speedVecPlus=zeros((len(RVec),))+nan
        speedVecMinus=zeros((len(RVec),))+nan
        sigma=9.0
        Ladv=sigma*0.0
        R0=2.5
        #R1=2.5*R0
        Ngen=30000 #maximum number of generations to run

        if False:
            print('only do this if P0 is is 0 everywhere')
            print('and so critical growth rate for P1 is not')
            print('R1/R0>1 but is 1/(ageMax-ageRepro+1)')
            RVec=RVec*(1.0/(ageMax-ageRepro+1.0))

        if False: #serial run
            for n in range(len(speedVecMinus)):
                R1=RVec[n]*R0
                speedVecPlus[n],speedVecMinus[n]=runModel(ageMax,ageRepro,R0,R1,sigma,Ladv,Ngen)
        else:
            #parallel run
            nProcs=8
            R1Vec=RVec*R0
            argvec=[(ageMax,ageRepro,R0,R1p,sigma,Ladv,Ngen) for R1p in R1Vec]
            with mp.Pool() as pool:
                output=pool.starmap(runModel,argvec)
            for n in range(len(output)):
                speedVecPlus[n],speedVecMinus[n]=output[n]
            print('DONE with parallel runs',flush=True)


        figure(3)
        col=colList.pop()
        plot(RVec,speedVecPlus/sigma,col+'-',label='lifespan is %d generations'%(ageMax,))
        plot(RVec,speedVecMinus/sigma,col+'-')
        plot(RVec,speedVecMinus/sigma*0,'k-')

        if False:
            #plot estimated speed.
            Linv=sigma*sqrt(2*log(RVec)/ageMax)
            LinvPlus=-(Ladv+Linv) #minus sign because of how I defined Ladv...
            plot(RVec,LinvPlus/sigma,col+'--')

            LinvMinus=-(Ladv-Linv) #minus sign because of how I defined Ladv...
            plot(RVec,LinvMinus/sigma,col+'--')

        
        title(r'$C_{invasion}$/$L_{diff}$')
        #xlabel('relative reproductive success, (R invader)/(R native) ')
        ylabel(r'downstream <== $C_{invasion}$/$L_{diff}$ ==> upstream')
        legend()
        draw(); show(); pause(0.2) #for profiling, comment this out, so it does not stop here

    #==================================================================================
    subplot(2,1,2)
    colList=['r','g','b','c']
    for ageMax in [1,3,5]: #[0,2,4,5]:
        ageRepro=0
        RVec=logspace(log10(1.02),log10(5.0),nTrys)
        #RVec=array([1.2,1.3,1.4])
        speedVecPlus=zeros((len(RVec),))+nan
        speedVecMinus=zeros((len(RVec),))+nan
        sigma=9.0
        Ladv=sigma*1.0
        R0=2.5
        #R1=2.5*R0
        #Ngen=3000 #maximum number of generations to run

        if False:
            print('only do this if P0 is is 0 everywhere')
            print('and so critical growth rate for P1 is not')
            print('R1/R0>1 but is 1/(ageMax-ageRepro+1)')
            RVec=RVec*(1.0/(ageMax-ageRepro+1.0))

        if False: #serial run
            for n in range(len(speedVecMinus)):
                R1=RVec[n]*R0
                speedVecPlus[n],speedVecMinus[n]=runModel(ageMax,ageRepro,R0,R1,sigma,Ladv,Ngen)
        else:
            #parallel run
            nProcs=8
            R1Vec=RVec*R0
            argvec=[(ageMax,ageRepro,R0,R1p,sigma,Ladv,Ngen) for R1p in R1Vec]
            with mp.Pool() as pool:
                output=pool.starmap(runModel,argvec)
            for n in range(len(output)):
                speedVecPlus[n],speedVecMinus[n]=output[n]
            print('DONE with parallel runs')


        figure(3)
        col=colList.pop()
        plot(RVec,speedVecPlus/sigma,col+'-',label='lifespan is %d generations'%(ageMax,))
        plot(RVec,speedVecMinus/sigma,col+'-')
        plot(RVec,speedVecMinus/sigma*0,'k-')

        if False:
            #plot estimated speed.
            Linv=sigma*sqrt(2*log(RVec)/ageMax)
            LinvPlus=-(Ladv+Linv) #minus sign because of how I defined Ladv...
            plot(RVec,LinvPlus/sigma,col+'--')

            LinvMinus=-(Ladv-Linv) #minus sign because of how I defined Ladv...
            plot(RVec,LinvMinus/sigma,col+'--')

        
        #title('Vary relative reproductive success')
        xlabel('relative reproductive success, (R invader)/(R native) ')
        #ylabel(r'downstream <== $C_{invasion}$/$L_{diff}$ ==> upstream')
        draw(); show(); pause(0.01) #for profiling, comment this out, so it does not stop here


        
    draw()
    show()
    pause(0.2)

savefig('Figure04_figInvasionSpeeds.svg')


