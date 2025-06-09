import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from natsort import natsorted
import re
from tqdm import tqdm
import image

plt.ion()
plt.close('all')

working_path = '/home/eric-neiva/TURLIERLAB Dropbox/Eric Neiva/Codes/SurfaceBulkViscousFlows/examples/OocyteAndSpindle'
contours_path = os.path.join(working_path,'results')
files = os.listdir(contours_path)
files = natsorted(files)

common_size = 300 # or set to desired number of points

for i,file in enumerate(files):
    
    # Read the CSV file
    filepath = os.path.join(contours_path,file)
    df = pd.read_csv(filepath, header=None)

    # Get the first and second columns
    first_column  = df.iloc[1:, 0]
    second_column = df.iloc[1:, 1]
    third_column  = df.iloc[1:, 2]
    fourth_column = df.iloc[1:, 3]

    # Convert the pandas Series to numpy arrays
    first_column_np  = np.array(first_column.values)
    second_column_np = np.array(second_column.values)
    third_column_np  = np.array(third_column.values)
    fourth_column_np = np.array(fourth_column.values)

    # Convert the numpy array of strings to double floats
    first_column_np  = first_column_np.astype(np.float64)
    second_column_np = second_column_np.astype(np.float64)
    third_column_np  = third_column_np.astype(np.float64)
    fourth_column_np = fourth_column_np.astype(np.float64)

    x = first_column_np
    y = second_column_np
    v = third_column_np
    theta = fourth_column_np

    # Compute radius and deviation
    r      = np.sqrt(x**2 + y**2)
    r_mean = np.mean(r)
    r_dev  = 100 * (r - r_mean) / r_mean

    # Interpolate to common length
    orig_theta   = theta
    common_theta = np.linspace(0, 2*np.pi, common_size, endpoint=False)
    r_interp     = np.interp(common_theta, orig_theta, r)
    rdev_interp  = np.interp(common_theta, orig_theta, r_dev)
    nvel_interp  = np.interp(common_theta, orig_theta, v)

    # Compute radii
    if i == 0:
        radii = r_interp
        rdevs = rdev_interp
        nvels = nvel_interp
    else:
        radii = np.vstack((radii, r_interp))
        rdevs = np.vstack((rdevs, rdev_interp))
        nvels = np.vstack((nvels, nvel_interp))

# RADII

plt.figure('Kymograph')
fig, axs = plt.subplots(1,1)
fig.suptitle('Kymograph of Simulated Contour Coordinates')

times = np.arange(1100,-1,-1)
angles = np.linspace(0, 360, common_size, endpoint=False)

axs.pcolormesh(angles, times, rdevs, cmap='viridis')
axs.set_ylabel('Time',labelpad=-20)
axs.set_xlabel('2D polar angle (non x-axis-aligned)',labelpad=30)
xticks = np.arange(0,301,75)
xlabels = ['0', '90', '180', '270', '360']
axs.xaxis.set_ticks(xticks,xlabels)
yticks = np.arange(0,1101,1100)
ylabels = ['P1', 'BD+8H']
axs.yaxis.set_ticks(yticks,ylabels)
axs.collections[0].set_clim(-4.75, 4.75)
cbar = fig.colorbar(axs.collections[0], ax=axs)
cbar.ax.set_ylabel('Deviation from mean radius (%)',rotation=-90,labelpad=15)

im = plt.imread(os.path.join(working_path,'Polar_graph_paper.png'), format='png')

newax = fig.add_axes([0.63, 0.01, 0.18, 0.18], anchor='E', zorder=1)
newax.imshow(im)
newax.axis('off')

plt.subplots_adjust(bottom=0.25,right=1.0,top=0.9,left=0.1)

# plt.show()
plt.savefig('kymograph_radii.png',dpi=200)

# VELOCITIES

plt.figure('Kymograph')
fig, axs = plt.subplots(1,1)
fig.suptitle('Kymograph of Simulated Contour Coordinates')

times = np.arange(1100,-1,-1)
angles = np.linspace(0, 360, common_size, endpoint=False)

axs.pcolormesh(angles, times, nvels, cmap='viridis')
axs.set_ylabel('Time',labelpad=-20)
axs.set_xlabel('2D polar angle (non x-axis-aligned)',labelpad=30)
xticks = np.arange(0,301,75)
xlabels = ['0', '90', '180', '270', '360']
axs.xaxis.set_ticks(xticks,xlabels)
yticks = np.arange(0,1101,1100)
ylabels = ['P1', 'BD+8H']
axs.yaxis.set_ticks(yticks,ylabels)
cbar = fig.colorbar(axs.collections[0], ax=axs)
cbar.ax.set_ylabel('Deviation from mean radius (%)',rotation=-90,labelpad=15)

im = plt.imread(os.path.join(working_path,'Polar_graph_paper.png'), format='png')

newax = fig.add_axes([0.63, 0.01, 0.18, 0.18], anchor='E', zorder=1)
newax.imshow(im)
newax.axis('off')

plt.subplots_adjust(bottom=0.25,right=1.0,top=0.9,left=0.1)

# plt.show()
plt.savefig('kymograph_velocities.png',dpi=200)

# References
# https://community.altair.com/community/en/magnitude-and-phase-angle-of-a-discrete-fourier-transform-dft-function?id=kb_article_view&sysparm_article=KB0121083&sys_kb_id=3d0dfd2b97839d50e3b0361e6253af39&spa=1
# https://towardsdatascience.com/the-fourier-transform-1-ca31adbfb9ef
# https://towardsdatascience.com/the-fourier-transform-2-understanding-phase-angle-a85ad40a194e
# https://www.chegg.com/homework-help/questions-and-answers/36-wave-nodes-fourier-series--based-results-problem-32-find-nodes-following-fourier-series-q22269237