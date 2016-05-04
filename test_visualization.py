# -*- coding: utf-8 -*-
"""
Visualize protein configuration (stat 525 final project)

@author: Chenchao
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as tk
import numpy as np


def plot_config(x_coord,input_HP_sequence,title,figname,plot_coord):
    """
    Input: x_coord: coordinates of all the protein residues = [[x1,y1],[x2,y2]...[xn,yn]]
           input_HP_sequence: HP sequence of protein (e.g. 'HPPHHP')
           title: figure title
           figname: figure name that will be saved as
           plot_coord: whether to plot coordinates = True/False. Plotting coordinates is usually for debugging purpose
    """
    
    plt.figure(1)
    plt.clf()
    # figure parameters
    dot_size = 30
    line_wid = 1    
    # first plot all individual residues    
    for i,coord in enumerate(x_coord):
        HP = input_HP_sequence[i]
        if HP == 'H':
            fill_col = 'k' # filled color of the dot is black for 'H' residue
        elif HP == 'P':
            fill_col = 'none' # no fill for 'P' residue      
        plt.scatter(coord[0],coord[1],s = dot_size, facecolor = fill_col, edgecolor = 'k')
    # then plot connected lines between adjacent residues    
    x_coord_np = np.array(x_coord) # convert to numpy array
    plt.plot(x_coord_np[:,0],x_coord_np[:,1],'k-',linewidth=line_wid)
    # add labels and title
    plt.title(title)
    if plot_coord:
        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid(True)
        # only show tickers at integers
        plt.gca().axes.get_yaxis().set_major_locator(tk.MaxNLocator(integer=True))
        plt.gca().axes.get_xaxis().set_major_locator(tk.MaxNLocator(integer=True))
    else:
        # make both axes invisible
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)
        # remove the frame
        plt.gca().set_frame_on(False)
    # set aspect ratios to be equal
    plt.gca().axes.set_aspect('equal')
    # save figure
    plt.savefig(figname)

# main code
figname = 'residue_visual_1.pdf'
title = 'protein residue plot'
input_HP_sequence = 'HPHPHP'
x_coord = [[0,0],[0,1],[1,1],[2,1],[2,0],[1,0]]
plot_coord = False
# sanity check
assert (len(input_HP_sequence) == len(x_coord))
plot_config(x_coord,input_HP_sequence,title,figname,plot_coord)
