# -*- coding: utf-8 -*-
"""
Stat 525 final project

@author: Chenchao
"""
def compute_energy(x, input_HP_sequence):
    
    """
    Compute the energy of a specific conformation
    """
    # add code here, feel free to change the argument list
    # Given a input HP sequence, we already which points are H's.
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
    
    direction = 0: turn left; 1: turn right; 2: go ahead
    
    Input: pt = [x1,y1], cpt = [x2,y2]
    Output: [a list of 3 points, a list of 3 directions]
    
    """
    nbrs = []
    dirs = []
    if pt[0] == cpt[0]:
        if pt[1] == cpt[1]+1:
            nbrs = [[pt[0]+1,pt[1]],[pt[0]-1,pt[1]],[pt[0],pt[1]+1]]
            dirs = [1,0,2]
        elif pt[1] == cpt[1]-1:
            nbrs = [[pt[0]+1,pt[1]],[pt[0]-1,pt[1]],[pt[0],pt[1]-1]]
            dirs = [0,1,2]
    if pt[1] == cpt[1]:
        if pt[0] == cpt[0]+1:
            nbrs = [[pt[0]+1,pt[1]],[pt[0],pt[1]+1],[pt[0],pt[1]-1]]
            dirs = [2,0,1]
        elif pt[0] == cpt[0]-1:
            nbrs = [[pt[0]-1,pt[1]],[pt[0],pt[1]+1],[pt[0],pt[1]-1]]
            dirs = [2,1,0]
    assert (len(nbrs) == 3 and len(dirs) == 3) # sanity check
    return [nbrs,dirs]
    
def multi_step_look_ahead(x, input_seq, steps):
    """
    Collect future information using the multi-step-look-ahead method to bias the movement of the next step
    
    Input: list of positions: x = [[x1,y1],[x2,y2]...[xn,yn]]
           input sequence: input_seq = [z1,z2...zm] (m=n+steps)
    """
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
    return [configs,dirs]
#    # compute unnormalized probability for each configuration and find marginals
#    # with respect to the direction
#    [p_left, p_right, p_ahead] = [0.0,0.0,0.0]
#    for i,cf in enumerate(configs):
#        d = dirs[i]
#        prob = exp(-compute_energy(cf,input_seq)/tau)
#        if d == 0:
#            p_left += prob
#        elif d == 1:
#            p_right += prob
#        elif d == 2:
#            p_ahead += prob
#    # return the unnormalized probabilities based on the immediate next step
#    return [p_left, p_right, p_ahead]

# test functions
pt1 = [0,0]
pt2 = [1,0]
x = [pt1,pt2]
print '---- test neighbor function ---'
print x
print neighbor(pt2,pt1)
print '---- test one-step function ---'
pt1 = [0,0]
pt2 = [1,0]
pt3 = [1,1]
pt4 = [0,1]
x = [pt1,pt2,pt3,pt4]
print x
print one_step(x)
print '---- test multi_step_look_ahead function ---'
pt1 = [0,0]
pt2 = [1,0]
pt3 = [1,1]
x = [pt1,pt2,pt3]
steps = 2
input_seq = [0]*(len(x)+steps)
print x
print 'step = %d' % steps
config,dirs = multi_step_look_ahead(x,input_seq,steps)
print config
print dirs