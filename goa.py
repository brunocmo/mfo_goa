from math import exp, cos, pi

import matplotlib as plt
import numpy as np

import matplotlib as mpl

def s_function( r ):
    f = 0.5
    l = 1.5
    return f*exp( -r/l ) - exp( -r )

def tracklsq(GDMH, output, input):
    W2 = GDMH[0]
    W = GDMHH[1]
    B = GDMH[3]
    (m,n) = np.shape()
    
    MaximoInt = max(input)
    MaxOut = max(output)
    sumerror = 0

    for i in range(0, round(n*0.7)):
        inputnormal = input[i]
        NNGDMH = W2*inputnormal**2 + W*inputnormal + B
        outDes = NNGDMH
        sumerror = sumerror + (outDes - output[i])**2
    
    return sumerror/round(n*0.7)


def initialization(N, dim, up, down):

    rng = np.random.default_rng(42)
    X = 0

    if up.shape[1-1] == 1:
        X = rng.random(N, dim)*(up-down)+down
    
    if up.shape[1-1] > 1:
        for i in range(1,dim):
            high=up[i]
            low=down[i]
            X[:,i]= rng.random[1,N]*(high-low)+low

    return X

def GOA(N, Max_iter, lb,ub, dim,sp):

    flag=0
    if ub.shape[1-1]==1:
        ub=np.ones((dim,1))*ub
        lb=np.ones((dim,1))*lb

    # this algorithm should be run with a even number of variables. This line is to handle odd number of variables

    if (dim%2)!=0:
        dim = dim+1
        ub = np.concatenate((ub,100))
        lb = np.concatenate((lb,100))
        flag=1

    # Initialize the population of grasshoppers
    GrassHopperPositions=initialization(N,dim,ub,lb)
    GrassHopperFitness = zeros(1,N)

    fitness_history=np.zeros((N,Max_iter))
    position_history=np.zeros((N,Max_iter,dim))
    Convergence_curve=np.zeros((1,Max_iter))
    Trajectories=np.zeros((N,Max_iter))

    cMax=1
    cMin=0.00004
    # Calculate the fitness of initial grasshoppers

    for i in range (1, GrassHopperPositions.shape[1-1] ):
        if flag == 1:
            GrassHopperFitness[1,i]=tracklsq(GrassHopperPositions[i,1:-1],sp)
        else:
            GrassHopperFitness[1,i]=tracklsq(GrassHopperPositions[i,:],sp)
        fitness_history[i,1]=GrassHopperFitness[1,i]
        position_history[i,1,:]=GrassHopperPositions[i,:]
        Trajectories[:,1]=GrassHopperPositions[:,1]

    [sorted_fitness,sorted_indexes]=np.sort(GrassHopperFitness)

    # Find the best grasshopper (target) in the first population 
    for newindex in range(1,N):
        Sorted_grasshopper[newindex,:]=GrassHopperPositions[sorted_indexes[newindex],:]

    TargetPosition=Sorted_grasshopper[1,:]
    TargetFitness=sorted_fitness[1]

    # Main loop
    l=2; # Start from the second iteration since the first iteration was dedicated to calculating the fitness of antlions
    while (l<Max_iter+1):
        
        c=cMax-l*((cMax-cMin)/Max_iter); # Eq. (2.8) in the paper
        
        for i in range(1, GrassHopperPositions.shape(1-1)):
            temp= GrassHopperPositions
            for k in range(1,dim,2):
                S_i=np.zeros(2,1)
                for j in range(1,N):
                    if i!=j:
                        # TODO procurar um jeito de criar ou pegar algo pronto distance
                        Dist=distance(temp[k:k+1,j], temp[k:k+1,i]); # Calculate the distance between two grasshoppers
                        
                        r_ij_vec=(temp[k:k+1,j]-temp[k:k+1,i])/(Dist+np.finfo(float).eps) # xj-xi/dij in Eq. (2.7)
                        xj_xi=2+(Dist%2) # |xjd - xid| in Eq. (2.7) 
                        
                        s_ij=((ub[k:k+1] - lb[k:k+1])*c/2)*S_func(xj_xi)*r_ij_vec # The first part inside the big bracket in Eq. (2.7)
                        S_i=S_i+s_ij
                S_i_total[k:k+1,:] = S_i
            
            X_new = c * S_i_total.conj().transpose()+ (TargetPosition); # Eq. (2.7) in the paper      
            GrassHopperPositions_temp[i,:]=X_new.conj().transpose(); 

        # GrassHopperPositions
        GrassHopperPositions=GrassHopperPositions_temp
        
        for i in range(1, GrassHopperPositions.shape(1-1)):
            # Relocate grasshoppers that go outside the search space 
            Tp=GrassHopperPositions[i,:]>ub.conj().transpose()
            Tm=GrassHopperPositions[i,:]<lb.conj().transpose()
            GrassHopperPositions[i,:]=(GrassHopperPositions[i,:]*(~(Tp+Tm)))+ub.conj().transpose()*Tp+lb.conj().transpose()*Tm
            
            # Calculating the objective values for all grasshoppers
            if flag == 1:
                GrassHopperFitness[1,i]=tracklsq(GrassHopperPositions[i,1:-1],sp)
            else:
                GrassHopperFitness[1,i]=tracklsq(GrassHopperPositions[i,:],sp)
            fitness_history[i,l]=GrassHopperFitness[1,i]
            position_history[i,l,:]=GrassHopperPositions[i,:]
            
            Trajectories[:,l]=GrassHopperPositions[:,1]
            
            # Update the target
            if GrassHopperFitness[1,i]<TargetFitness:
                TargetPosition=GrassHopperPositions[i,:];
                TargetFitness=GrassHopperFitness[1,i];
            
        Convergence_curve[l]=TargetFitness

        print(f'At iteration: {l}, target\'s objective = {TargetFitness}')
        
        l = l + 1
        print(f'Iteration ==> {l}')
        globalmin = Convergence_curve[Max_iter]

    if (flag==1):
        TargetPosition = TargetPosition[1:dim-1]
