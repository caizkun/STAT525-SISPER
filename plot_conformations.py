# STAT5252-SISPER
#
# This script is used to plot protein conformations.
#
# Usage:
#   python plot_conformations.py sequence_file_name conformation_file_name fig_file_keywords num_of_figs


import matplotlib.pyplot as plt
import matplotlib.ticker as tk
import numpy as np
import sys


def read_input_sequence(sequence_file_name):
    """
        read in the input HP sequence from file
        """
    input_HP_sequence = ''
    with open(sequence_file_name, 'r') as file_id:
        for line in file_id:
            input_HP_sequence += line.replace('\n', '')
    return input_HP_sequence


def read_conformations(conformation_file_name):
    """
    read from file the conformations sorted by energy from low to high
    """
    conformations = np.load(conformation_file_name)
    U = conformations['energies'].tolist()
    S = conformations['coordinates'].tolist()
    w = conformations['weights'].tolist()
    return (U, S, w)


def plot_config(x_coord,input_HP_sequence,title,figname,is_grid_plotted):
    """
    Input:  x_coord: coordinates of all the protein residues = [[x1,y1],[x2,y2]...[xn,yn]]
            input_HP_sequence: HP sequence of protein (e.g. 'HPPHHP')
            title: figure title
            figname: figure name that will be saved as
            is_grid_plotted: whether to plot grids = True/False.
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
    if is_grid_plotted:
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


# --- main code ---
def main():
    if len(sys.argv) == 5:
        sequence_file_name = sys.argv[1]
        conformation_file_name = sys.argv[2]
        fig_file_keywords = sys.argv[3]
        num_of_figs = int(sys.argv[4])
    else:
        "Error in input! Please provide: sequence_file_name, conformation_file_bame, num_of plots"

    # read in data
    input_HP_sequence = read_input_sequence(sequence_file_name)
    U, S, w = read_conformations(conformation_file_name)

    # plot conformations
    seq_len = len(input_HP_sequence)
    assert(num_of_figs <= len(S))

    for i in xrange(num_of_figs):
        title = "conformation: length = %d, energy = %g" % (seq_len, U[i])
        fig_name = fig_file_keywords + "_Umin_%d.pdf" % (i+1)
        plot_config(S[i], input_HP_sequence, title, fig_name, True)


if __name__ == '__main__':
    main()







