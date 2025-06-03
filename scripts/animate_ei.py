
import matplotlib.pyplot as plt
from matplotlib import animation, rc
from matplotlib.animation import FuncAnimation, PillowWriter
from IPython.core.display import HTML
from IPython.display import display
import pandas as pd
import argparse
import numpy as np
import os

def plot_ei(ei, coords, playback_interval=50, threshold=1.0):
    SIZE_SCALAR = 10
    fig, ax = plt.subplots(1, 1, figsize=(10,5), dpi=60)
    # n = ax.scatter(coords[:,0],
    #                coords[:,1],c='k',s=0)
    # n = ax.scatter(coords[:-1,0],
    #                coords[:-1,1],c='k',s=0)
    n = ax.scatter(coords[1:,0],
                   coords[1:,1],c='k',s=0)
    ax.axis('off')
    ei = ei / np.std(ei)
    if threshold > 0.0:
        thresh_idx = np.abs(ei) < threshold
        # thresh_idx = np.abs(ei) < threshold
        ei[thresh_idx] = 0
    # ei -= np.min(ei)
    #ei = np.abs(ei) 
    def animate(i):
        n.set_sizes(ei[:,i]* SIZE_SCALAR)
        # n.set_sizes(ei[1:,i]* SIZE_SCALAR)
        # n.set_sizes(ei[:-1,i]* SIZE_SCALAR)
    anim = animation.FuncAnimation(fig, animate, frames=ei.shape[1], interval=playback_interval) #frames max = 100
    return anim

def visualize_rawdata(ei, coords, filepath, threshold=3.0):
    anim = plot_ei(ei, coords, threshold=threshold)
    HTML(anim.to_jshtml())
    anim.save(filepath, writer='pillow', fps=15)

if __name__ == '__main__':
    '''
    Animate all of the EI files for an experiment. Could be modified to animate a subset.
    '''

    parser = argparse.ArgumentParser(description='Loop through the root directory and average all EI files.')

    parser.add_argument('sortDir', type=str, help='Path to sorted spikes and .ei files')
    parser.add_argument('experimentName', type=str, help='Experiment name')

    args = parser.parse_args()

    # Get the path to the file.
    # d = os.path.join(args.sortDir, args.experimentName + '_EI.p')
    d = os.path.join(args.sortDir, 'EI.p')

    # outdir = os.path.join(args.sortDir, args.experimentName+'_EI')
    outdir = args.sortDir

    # Check if the output directory exists, and make it if not.
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    object = pd.read_pickle(d) 
    E = object['E']

    key=6
    visualize_rawdata(E[key], object['electrode_map'], os.path.join(outdir,'ei_'+str(key)+'.gif'))

    # Iterate through the keys and generate the GIF.
    # for key in E:
    #     visualize_rawdata(E[key], object['electrode_map'], os.path.join(outdir,key+'.gif'))


# anim = visualize_rawdata(ei,coords)
# HTML(anim.to_jshtml())
# anim.save('./ei.gif',writer='pillow', fps=15)
# object = pd.read_pickle(r'filepath') 
# E = object['E']
# import animate_ei as a
# anim = a.plot_ei(E['e348'], object['electrode_map'])