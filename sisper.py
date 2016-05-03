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
#   J. L. Zhang and J. S. Liu, A new sequential improtance sampling method and its application to the two-dimensional hydrophobic-hydrophilic model,
#   Journal of Chemical Physics, 117, 3492 (2002)
#
# Developers:
# 	Zhikun Cai, Chenchao Shou, Guanfeng Gao
#


# import modules
from __future__ import division
from scipy import stats

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
    return input_HP_sequence


#def compute_conformation_coordinates(x_angle):
#    """
#    Compute the conformation cooridates from the vector of torsion angles
#    """
#    # return a list of conformation coordinates (x_i, y_i)
#    # feel free to change the function naming
#    x = np.asarray(x_angle) #change from list to array by guanfeng on May 2nd
#    x_i=np.zeros(x.size+2) #initialize all the points in x and y coordinates
#    y_i=np.zeros(x.size+2) #initialize all the points in x and y coordinates
#    x_i[0]=0#get the first point
#    y_i[0]=0
#    x_i[1]=0#get the second point in the positive y direction
#    y_i[1]=1
#    direction=0 #represent current direction, 0 is +y, 1 is -x, 2 is -y, 3 is +x
#    print x_angle
#    def next_position1(x_c,y_c,x_tmp,flag):#compute the next point in direction 0 and 2
#        print "loop:", x_tmp, "--"
#        if x_tmp==0:
#            x_c=x_c+pow(-1,flag)
#        elif x_tmp==1:
#            x_c=x_c+pow(-1,flag+1)
#        elif x_tmp==2:
#            y_c=y_c+pow(-1,flag+1)
#        return (x_c,y_c)
#    def next_position2(x_c,y_c,x_tmp,flag):#compute the next point in direction 1 and 3
#        if x_tmp==0:
#            y_c=y_c+pow(-1,flag)
#        elif x_tmp==1:
#            y_c=y_c+pow(-1,flag+1)
#        elif x_tmp==2:
#            x_c=x_c+pow(-1,flag)
#        return (x_c,y_c)
#    for ii in range(x.size):#loop through every x
#        if direction==0:
#            x_i[ii+2],y_i[ii+2]=next_position1(x_i[ii+1],y_i[ii+1],x[ii],1)
#        elif direction==1:
#            x_i[ii+2],y_i[ii+2]=next_position2(x_i[ii+1],y_i[ii+1],x[ii],1)
#        elif direction==2:
#            x_i[ii+2],y_i[ii+2]=next_position1(x_i[ii+1],y_i[ii+1],x[ii],0)
#        elif direction==3:
#            x_i[ii+2],y_i[ii+2]=next_position2(x_i[ii+1],y_i[ii+1],x[ii],0)
#        #change the direction denpends on x[ii]
#        if x[ii]==0:
#            direction=direction+1
#            direction=np.mod(direction,4)
#        elif x[ii]==1:
#            direction=direction-1
#            direction=np.mod(direction,4)
#    coordinate=np.transpose(np.vstack((x_i,y_i))) #change from array to list by guanfeng on March 2nd
#    coordinate=coordinate.tolist()
#    return coordinate


def compute_couples(input_HP_sequence):
    """
    compute all the possible couples in the series
    """
    num=0   #initialize the size of couples
    couples=np.zeros([len(input_HP_sequence)**2,2])
    for ii in range(len(input_HP_sequence)):
        if input_HP_sequence[ii]=='H':
            for jj in range(ii+3,len(input_HP_sequence),2): #loop through the rest possible nodes
                if input_HP_sequence[jj]=='H':
                    couples[num,0]=ii
                    couples[num,1]=jj
                    num=num+1
    couples=couples[0:num,:]    #get the effective couples
    return couples


def compute_energy(x_angle, couples):
    """
    Compute the energy of a specific conformation
    """
    U=0
    coordinate=compute_conformation_coordinates(x_angle)#compute all the positions needed
    x = np.asarray(x_angle)#change from list to array by guanfeng on May 2nd
    #change from list to array
    for ii in range(int(couples.size/2)):
        if couples[ii,0]>x.size+1 or couples[ii,1]>x.size+1:#exclude the couples exceed the current set
            continue
        else:
            point1=np.array(coordinate[int(couples[ii,0])])
            point2=np.array(coordinate[int(couples[ii,1])])
            if np.dot(point1-point2,point1-point2)==1:  #neighbour point
                U=U+1
    return U


def one_step(x):
    """
    return all possible configurations after one step given the current position
    and direction for each configuration

    Input: current position x = [[x1,y1],[x2,y2]...[xn,yn]]
    Output: [a list of configurations, a list of directions]
    """
    assert (len(x)>=2)
    nbrs,dirs = neighbor(x[-1],x[-2]) # find 3 neighbors of the end point and their directions
    x3 = x[:-3] # excluding the last 3 points
    configs = []
    dir_configs = []
    for i,nb in enumerate(nbrs):
        d = dirs[i]
        valid_nb = True
        for pt in x3:
            if pt == nb:
                valid_nb = False
        if valid_nb:
            configs.append(x+[nb])
            dir_configs.append(d)
    return [configs,dir_configs]


def neighbor(pt,cpt):
    """
    return the neighbors of pt other than cpt as well as direction

    Input: pt = [x1,y1], cpt = [x2,y2]
    Output: [a list of 3 points, a list of 3 directions]

    """
    nbrs = []
    dirs = []
    if pt[0] == cpt[0]:
        if pt[1] == cpt[1]+1:
            nbrs = [[pt[0]+1,pt[1]],[pt[0]-1,pt[1]],[pt[0],pt[1]+1]]
            dirs = [right,left,ahead]
        elif pt[1] == cpt[1]-1:
            nbrs = [[pt[0]+1,pt[1]],[pt[0]-1,pt[1]],[pt[0],pt[1]-1]]
            dirs = [left,right,ahead]
    if pt[1] == cpt[1]:
        if pt[0] == cpt[0]+1:
            nbrs = [[pt[0]+1,pt[1]],[pt[0],pt[1]+1],[pt[0],pt[1]-1]]
            dirs = [ahead,left,right]
        elif pt[0] == cpt[0]-1:
            nbrs = [[pt[0]-1,pt[1]],[pt[0],pt[1]+1],[pt[0],pt[1]-1]]
            dirs = [ahead,right,left]
    assert (len(nbrs) == 3 and len(dirs) == 3), 'error! number of neighbors is invalid!' # sanity check
    return [nbrs,dirs]


def multi_step_look_ahead(x, input_HP_sequence, steps_tmp):
    """
    Collect future information using the multi-step-look-ahead method to bias the movement of the next step

    Input:  current position x = [[x1,y1],[x2,y2]...[xn,yn]],
            input_HP_sequence: the protein HP sequence input by user
            step_tmp: num of steps look-ahead in the algorithm
    Output: unormalized probabilities towards three different directions, i.e. left, right, and ahead
    
    """
    steps = min(steps_tmp,len(input_HP_sequence)-len(x))
    input_seq = input_HP_sequence[:(len(x)+steps)]

    couples = compute_couples(input_seq)

    # look ahead to get all possible configurations
    # do the first step
    [c1,d1] = one_step(x)
    # then do the remaining steps
    configs = c1
    dirs = d1 # directions for each configuration
    for s in range(steps-1):
        new_config = [] # a list of new configurations
        new_dirs = [] # a list of directions for new configurations
        for i,cf in enumerate(configs):
            [cfg,d] = one_step(cf)
            new_config += cfg
            new_dirs += [dirs[i]]*len(cfg)
        configs = new_config
        dirs = new_dirs

    # compute unnormalized probability for each configuration and find marginals
    # with respect to the direction
    [p_left, p_right, p_ahead] = [0.0,0.0,0.0]
    for i,cf in enumerate(configs):
        d = dirs[i]
        prob = np.exp(-compute_energy(cf,couples)/tau)
        if d == left:
            p_left += prob
        elif d == right:
            p_right += prob
        elif d == ahead:
            p_ahead += prob
    # return the unnormalized probabilities based on the immediate next step
    return [p_left, p_right, p_ahead]


def compute_next_point(pt1, pt2, move):
    """
    compute the next position of a newly added point
    """
    disp_vec = [pt1[0] - pt2[0], pt1[1] - pt2[1]]
    next_pt = [pt1[0], pt1[1]]
    
    if disp_vec == [0, 1]:
        if move == "left":
            next_pt[0] -= 1
        elif move == "right":
            next_pt[0] += 1
        elif move == "ahead":
            next_pt[1] += 1
        else:
            print "Unrecognized move direction!"
    elif disp_vec == [0, -1]:
        if move == "left":
            next_pt[0] += 1
        elif move == "right":
            next_pt[0] -= 1
        elif move == "ahead":
            next_pt[1] -= 1
        else:
            print "Unrecognized move direction!"
    elif disp_vec == [1, 0]:
        if move == "left":
            next_pt[1] += 1
        elif move == "right":
            next_pt[1] -= 1
        elif move == "ahead":
            next_pt[0] += 1
        else:
            print "Unrecognized move direction!"
    elif disp_vec == [-1, 0]:
        if move == "left":
            next_pt[1] -= 1
        elif move == "right":
            next_pt[1] += 1
        elif move == "ahead":
            next_pt[0] -= 1
        else:
            print "Unrecognized move direction!"
    else:
        print "Wrong displacement bewteen two continuous points!"

    return next_pt


def resample_conformations(S, w, a, N_star):
    """
    perform standard resampling from the conformation set S and the probability vector a
    """
    assert(min(a) >= 0)

    # normalize the probability vector a
    a_normal = np.array(a)
    a_normal = a_normal / np.sum(a_normal)

    # resample conformations S_star with probabilities proportional to a_normal
    S_indice = np.arrange(len(S))
    resample_obj = stats.rv_discrete(name = "resample_obj", values = (S_indice, a_normal))
    S_star_indice = resample_obj.rvs(size = N_star)

    S_star = []
    w_star = []
    for i in S_star_indice:
        S_star.append(S[i])
        w_star.append(w[i])     # TODO: this weight may be wrong

    # return the resampled conformations S_star and the weights w_star
    return (S_star, w_star)


def resample_with_pilot_exploration(S, w, input_HP_sequence, N_star):
    """
    Resample the current set of conformations using the pilot-exploration resampling scheme
    """
    a = []      # resampling probability vector
    # generate the next Delta residues m independent times for each conformation in S
    for j, x in enumerate(S):
        bj = 0
        for l in xrange(m):
            xl = list(x)
            pi_l = 0
            for i in xrange(Delta):
                # compute the probabilities along different directions
                [p_left, p_right, p_ahead] = multi_step_look_ahead(xl, input_HP_sequence, rho)
                p_sum = p_left + p_right + p_ahead

                if p_sum == 0:
                    pi_l = 0
                    break

                # move a new step and update the weight of the conformation
                rand_num = random.random()
                if (rand_num < p_left/p_sum):
                    next_pt = compute_next_point(xl[t], xl[t-1], "left")
                    xl.append(next_pt)
                    pi_l = p_left
                elif (rand_num < (p_left + p_right)/p_sum):
                    next_pt = compute_next_point(xl[t], xl[t-1], "right")
                    xl.append(next_pt)
                    pi_l = p_right
                else:
                    next_pt = compute_next_point(xl[t], xl[t-1], "ahead")
                    xl.append(next_pt)
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
    if len(sys.argv) > 1:
        sequence_file_name = str(sys.argv[1])
    else:
        print "ERROR! Missing a sequence file!"

    # read in HP sequence
    input_HP_sequence = read_input_sequence(sequence_file_name)
    seq_len = len(input_HP_sequence)

    # connectivity couples
    couples = compute_couples(input_HP_sequence)

    # fix the first step along the vertical direction
    # this step, x0, doesn't matter at all

    # initialization
    S = N * [[[0, 0], [0, 1]]]      # conformation set
    w = N * [1]                     # conformation weights
    U = N * [0]                     # conformation energy

    # sequentially generate conformations
    for t in xrange(1, seq_len-1):
        
        # initialize a list to save the indice of incorrectly terminated conformations
        failed_indice = []

        # perform a regular SIS step with multi-step-look-ahead
        for i, x in enumerate(S):
            # compute the probabilities along different directions
            [p_left, p_right, p_ahead] = multi_step_look_ahead(x, input_HP_sequence, delta)
            p_sum = p_left + p_right + p_ahead

            # record the indices of incorrectly terminated conformations
            if p_sum == 0:
                failed_indice.append(i)
                continue

            p_left  /= p_sum
            p_right /= p_sum
            p_ahead /= p_sum

            # move a new step and update the weight of the conformation
            rand_num = random.random()
            if rand_num < p_left:
                next_pt = compute_next_point(x[t], x[t-1], "left")
                x.append(next_pt)
                U_new = compute_energy(x, couples)
                w[i] *= np.exp(- (U_new - U[i]) / tau) / p_left
            elif (rand_num < (p_left + p_right)):
                next_pt = compute_next_point(x[t], x[t-1], "right")
                x.append(next_pt)
                U_new = compute_energy(x, couples)
                w[i] *= np.exp(- (U_new - U[i]) / tau) / p_right
            else:
                next_pt = compute_next_point(x[t], x[t-1], "ahead")
                x.append(next_pt)
                U_new = compute_energy(x, couples)
                w[i] *= np.exp(- (U_new - U[i]) / tau) / p_ahead

            # save the energy of the current configuration
            U[i] = U_new

        # remove incorrectly terminated conformations
        for i, index in enumerate(failed_indice):
            del S[index - i]    # shift the index due to deletion; work for the sorted list
            del w[index - i]
            del U[index - i]

        # perform resampling
        if t % (lamb + 1) == 0:
            S, w = resample_with_pilot_exploration(S, w, input_HP_sequence, len(S))

            # update the conformation energies
            for i, x in enumerate(S):
                U[i] = compute_energy(x, couples)


if __name__ == '__main__':
    main()
