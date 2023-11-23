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
    W = GDMH[1]
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


def initialization(SearchAgents_no, dim, ub, lb):

    rng = np.random.default_rng(42)
    Boundary_no = ub.shape[2-1]
    X = 0

    if Boundary_no == 1:
        X = rng.random(SearchAgents_no, 1)*(ub-lb)+lb
        return X 

    if Boundary_no > 1:
        for i in dim:
            ub_i = ub[i]
            lb_i = lb[i]
            X[:,i] = rng.random(SearchAgents_no, 1)*(ub_i-lb_i)+lb_i
        return X

def MFO(N,Max_iteration,lb,ub,dim,sp):

    Moth_pos = initialization(N,dim,ub,lb)
    # stoppMOF = 0
    # lastBestScore = 0
    Convergence_curve = np.zeros((1,Max_iteration))

    rng = np.random.default_rng(42)

    Iteration = 1
    Moth_fitness= np.zeros((1, Moth_pos.shape[1-1]))

    while(Iteration<Max_iteration+1):
        Flame_no=round(N-Iteration*((N-1)/Max_iteration))

        for i in Moth_pos.shape[1,1]:
            Flag4ub = Moth_pos[i,:] > ub
            Flag4lb = Moth_pos[i,:] < ub
            Moth_pos[i,:] = ( Moth_pos[i,:] * ( not( Flag4ub+Flag4lb ) ) ) + ub * Flag4ub + lb * Flag4lb

            # O que é tracklsq?
            Moth_fitness[1,i] = tracklsq(Moth_pos[i,:],sp)
        
        if Iteration==1:
            [fitness_sorted, I] = Moth_fitness.sort()
            sorted_population = Moth_pos[I,:]

            best_flames = sorted_population
            best_flames_fitness = fitness_sorted
        else:
            # Sort the moths
            double_population=np.concatenate((previous_population,best_flames))
            double_fitness= np.concatenate(previous_fitness,best_flame_fitness, axis=None)
            
            [double_fitness_sorted, I]=double_fitness.sort()
            double_sorted_population=double_population[I,:]
            
            fitness_sorted=double_fitness_sorted[1:N]
            sorted_population=double_sorted_population[1:N,:]
            
        # Update the flames
        best_flames=sorted_population
        best_flame_fitness=fitness_sorted

        Best_flame_score = fitness_sorted[0]
        Best_flamse_pos = sorted_population[1,:]

        previous_population=Moth_pos
        previous_fitness=Moth_fitness

        a = 1 + Iteration*(-1/Max_iteration)

        for i in Moth_pos.shape[Moth_pos,1]:
            for j in Moth_pos.shape[Moth_pos,2]:
                if i <= Flame_no:
                    distance_to_flame=abs(sorted_population[i,j]-Moth_pos[i,j])
                    b = 1
                    t = (a-1)*rng.random()+1

                    Moth_pos[i,j] = distance_to_flame*exp(b*t)*cos(t*2*pi) + sorted_population[Flame_no,j]
        
        Convergence_curve[Iteration]=Best_flame_score

        if Iteration % 50 == 0:
            print(f'At iteration: {Iteration}, the best fitness is {Best_flame_score}')
        
        Iteration= Iteration+1
        print(f'Iteration ===> {Iteration}')