#
# STAT525-SISPER, 2016
#
# Problem:
#   Searching for the lowest-energy folded state of protein in a 2D hydrophobic-hydrophilic (HP) lattice model
#
# Algorithm:
#   Sequential Importance Sampling with Pilot-Exploration Resampling (SISPER)
#
# Reference:
#   J. L. Zhang and J. S. Liu, A new sequential importance sampling method and its application to the two-dimensional hydrophobic–hydrophilic model,
# The Journal of Chemical Physics, 117, 3492-3498 (2002).
#
# Developers:
# 	Zhikun Cai, Chenchao Shou, Guanfeng Gao
#


# import modules
from __future__ import division

import numpy as np
import random
import sys


# initialize global parameters
N     = 5000        # num of conformations
delta = 5           # num of steps lookahead in the regular SIS steps
lamb  = 2           # frequency of resampling
m     = 20          # num of independent paths in the resampling steps
Delta = 20          # num of steps explored in the resampling steps
rho   = 1           # num of steps lookahead in the resampling steps
alpha = 0.5         # power of the Boltzmann weight in calculating the resampling probability
tau   = 0.5         # temperature

# define global constants
res_H = 1           # flag indicating the residue type, H or P
res_P = 0

esp_HH = -1         # pairwised energy
esp_HP = 0
esp_PP = 0

left  = 0           # flag indicating the torsion angle
right = 1
ahead = 2


# define functions
def read_input_sequence(file_name):
    """
    read in the input HP sequence from file
    """
    input_HP_sequence = ''
    with open(file_name, 'r') as file_id:
        for line in file_id:
            input_HP_sequence += line.replace('\n', '')

    file_id.close()
    return input_HP_sequence


def compute_conformation_coordinates(x):
    """
    Compute the conformation cooridates from the vector of torsion angles
    """
    # return a list of conformation coordinates (x_i, y_i)
    # feel free to change the function naming
    x_i=np.zeros(x.size+2) #initialize all the points in x and y coordinates
    y_i=np.zeros(x.size+2) #initialize all the points in x and y coordinates
    x_i[0]=0#get the first point
    y_i[0]=0
    x_i[1]=0#get the second point in the positive y direction
    y_i[1]=1
    direction=0 #represent current direction, 0 is +y, 1 is -x, 2 is -y, 3 is +x
    def next_position1(x_c,y_c,x,flag):#compute the next point in direction 0 and 2
            if x==0:
                x_c=x_c+pow(-1,flag)
            elif x==1:
                x_c=x_c+pow(-1,flag+1)
            elif x==2:
                y_c=y_c+pow(-1,flag+1)
            return x_c,y_c 
    def next_position2(x_c,y_c,x,flag):#compute the next point in direction 1 and 3
            if x==0:
                y_c=y_c+pow(-1,flag)
            elif x==1:
                y_c=y_c+pow(-1,flag+1)
            elif x==2:
                x_c=x_c+pow(-1,flag)
            return x_c,y_c 
    for ii in range(x.size):#loop through every x  
        if direction==0:
            x_i[ii+2],y_i[ii+2]=next_position1(x_i[ii+1],y_i[ii+1],x[ii],1)
        elif direction==1:
            x_i[ii+2],y_i[ii+2]=next_position2(x_i[ii+1],y_i[ii+1],x[ii],1)
        elif direction==2:
            x_i[ii+2],y_i[ii+2]=next_position1(x_i[ii+1],y_i[ii+1],x[ii],0)
        elif direction==3:
            x_i[ii+2],y_i[ii+2]=next_position2(x_i[ii+1],y_i[ii+1],x[ii],0)
        #change the direction denpends on x[ii]
        if x[ii]==0:
            direction=direction+1
            direction=np.mod(direction,4)
        elif x[ii]==1:
            direction=direction-1
            direction=np.mod(direction,4)
    return x_i,y_i
                
def compute_couples(input_HP_sequence):#compute all the possible couples in the series
    num=0#initialize the size of couples
    couples=np.zeros([len(input_HP_sequence)**2,2])    
    for ii in range(len(input_HP_sequence)):
        if input_HP_sequence[ii]=='H':   
            for jj in range(ii+3,len(input_HP_sequence),2):#loop through the rest possible nodes
                if input_HP_sequence[jj]=='H':
                    couples[num,0]=ii
                    couples[num,1]=jj
                    num=num+1
    couples=couples[0:num,:]#get the effective couples
    return couples
            
def compute_energy(x, couples):
    """
    Compute the energy of a specific conformation
    """
    # add code here, feel free to change the argument list
    # Given a input HP sequence, we already which points are H's.
    U=0    
    x_i,y_i=compute_conformation_coordinates(x)#compute all the positions needed
    for ii in range(int(couples.size/2)):
        if couples[ii,0]>x.size+1 or couples[ii,1]>x.size+1:#exclude the couples exceed the current set
            continue
        else:
            point1=np.array([x_i[couples[ii,0]],y_i[couples[ii,0]]])
            point2=np.array([x_i[couples[ii,1]],y_i[couples[ii,1]]])
            if np.dot(point1-point2,point1-point2)==1:  #neighbour point 
                U=U+1
    U=np.exp(-U/tau)
    return U


def multi_step_look_ahead(x, steps):
    """
    Collect future information using the multi-step-look-ahead method to bias the movement of the next step
    """
    # return the unnormalized probabilities based on the immediate next step
    return [p_left, p_right, p_ahead]


def resample_conformations(S, w, a, N_star):
    """
    perform standard resampling from the conformation set S and the probability vector a
    """
    # return the resampled conformations S_star and the weights w_star
    return (S_star, w_star)


def resample_with_pilot_exploration(S, w, N_star):
    """
    Resample the current set of conformations using the pilot-exploration resampling scheme
    """
    # put one resampling step in a function, which may call multiple small
    # functions to complete the tasks.
    # for example, it may need a support function to conduct the residual
    # resampling described in the appendix of the paper

    a = []      # resampling probability vector
    # generate the next Delta residues m independent times for each conformation in S
    for j, x in enumerate(S):
        bj = 0
        for l in xrange(m):
            xl = list(x)
            pi_l = 0
            for i in xrange(Delta):
                # compute the probabilities along different directions
                [p_left, p_right, p_ahead] = multi_step_look_ahead(xl, rho)
                p_sum    = p_left + p_right + p_ahead

                # move a new step and update the weight of the conformation
                rand_num = random.random()
                if (rand_num < p_left/p_sum):
                    xl.append(left)
                    pi_l = p_left
                elif (rand_num < (p_left + p_right)/p_sum):
                    xl.append(right)
                    pi_l = p_right
                else:
                    xl.append(ahead)
                    pi_l = p_ahead

            # sum up the Boltzmann weight of each path l
            bj += pi_l          # pi_l value exiting the loop

        # compute the unnormalized resampling probability aj
        bj /= len(S)
        aj = bj^alpha
        a.append(aj)

    # perform a standard resampling step with the probability vector
    S_star, w_star = resample_conformations(S, w, a, N_star)

    # return the resampled conformations S_star and the weights w_star
    return (S_star, w_star)


def main():
    """
    Main function to implement the alogrithm SISPER
    """
    if len(sys.argv) == 2:
        sequence_file_name = str(sys.argv[1])
    else:
        print "ERROR! Missing a sequence file!"

    # read in HP sequence
    input_HP_sequence = read_input_sequence(sequence_file_name)

    # sequence length
    d = len(input_HP_sequence) - 2

    # fix the first step along the horizontal direction
    # this step, x0, doesn't matter at all

    # initialization
    S = []          # conformation set
    w = []          # conformation weights
    U = []          # conformation energy
    for i in xrange(N):
        S.append([ahead])
        w.append(1)
        U.append(0)

    # sequentially generate conformations
    for t in xrange(1, d+1, 1):
        # perform a regular SIS step with multi-step-look-ahead
        for i, x in enumerate(S):
            # compute the probabilities along different directions
            [p_left, p_right, p_ahead] = multi_step_look_ahead(x, delta)
            p_sum    = p_left + p_right + p_ahead
            p_left  /= p_sum
            p_right /= p_sum
            p_ahead /= p_sum

            # move a new step and update the weight of the conformation
            rand_num = random.random()
            if (rand_num < p_left):
                x.append(left)
                U_new = compute_energy(x, input_HP_sequence)
                w[i] *= np.exp(- (U_new - U[i]) / tau) / p_left
            elif (rand_num < (p_left + p_right)):
                x.append(right)
                U_new = compute_energy(x, input_HP_sequence)
                w[i] *= np.exp(- (U_new - U[i]) / tau) / p_right
            else:
                x.append(ahead)
                U_new = compute_energy(x, input_HP_sequence)
                w[i] *= np.exp(- (U_new - U[i]) / tau) / p_ahead

            # save the energy of the current configuration
            U[i] = U_new

        # perform resampling
        if (t % (lamb + 1) == 0):
            S, w = resample_with_pilot_exploration(S, w, N)

            # update the conformation energies
            for i, x in enumerate(S):
                U[i] = compute_energy(x, input_HP_sequence)


if __name__ == '__main__':
    main()
