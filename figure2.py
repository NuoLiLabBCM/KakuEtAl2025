import numpy as np
import matplotlib as mpl
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42

import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 12
plt.rcParams['svg.fonttype'] = 'none'

import pickle

#%% Figure 2c
fname = r'.\data\figure_2\data_for_fig_2c.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()

region_rgba, y_ref, anno_depth_bins, spk_rates,spike_bins, depth_bins, depth_edges, ton_onset, ton_offset, session_traces_b_y, session_traces_s_y, good_breathing = data
# canvas setup
fig = plt.figure(figsize=(16, 10))
grid = plt.GridSpec(12, 10)

ax_main = plt.subplot(grid[3:, 0:9])
ax_tong = plt.subplot(grid[0:1, 0:9])
ax_lick = plt.subplot(grid[1:2, 0:9])
ax_brth = plt.subplot(grid[2:3, 0:9])
ax_anno = plt.subplot(grid[3:, 9:])

# -- plot main --
im = ax_main.imshow(spk_rates.T, aspect='auto', cmap='gray_r', vmax=200,
                    extent=[spike_bins[0], spike_bins[-1], depth_bins[-1], depth_bins[0]])
# cosmetic
ax_main.invert_yaxis()
ax_main.set_xlabel('Seconds')
ax_main.set_ylabel('Distance from tip sites (um)')
ax_main.set_ylim(depth_edges[0], depth_edges[-1])
ax_main.set_xlim(0, 20)
ax_main.spines['right'].set_visible(False)
ax_main.spines['top'].set_visible(False)

t_x=np.arange(len(session_traces_b_y))*0.0034

for i, t in enumerate(ton_onset):
    ax_tong.plot(t_x[ton_onset[i]+1:ton_offset[i]+1],session_traces_b_y[ton_onset[i]+1:ton_offset[i]+1],'r')
    
ax_tong.spines['right'].set_visible(False)
ax_tong.spines['top'].set_visible(False)
ax_tong.spines['bottom'].set_visible(False)
ax_tong.set_xticks([])
ax_tong.set_xlim(0, 20)
ax_tong.set_ylabel("mm", fontsize=10)

ax_lick.plot(np.arange(len(session_traces_s_y))*0.0034,session_traces_s_y,'b')
ax_lick.spines['right'].set_visible(False)
ax_lick.spines['top'].set_visible(False)
ax_lick.spines['bottom'].set_visible(False)
ax_lick.set_xticks([])
ax_lick.set_xlim(0, 20)
ax_lick.set_ylabel("mm", fontsize=10)


ax_brth.plot(np.arange(len(good_breathing))*0.0034,-good_breathing,'g')
ax_brth.spines['right'].set_visible(False)
ax_brth.spines['top'].set_visible(False)
ax_brth.spines['bottom'].set_visible(False)
ax_brth.spines['left'].set_visible(False)
ax_brth.set_xticks([])
ax_brth.set_yticks([])
ax_brth.set_xlim(0, 20)

# -- plot colored region annotation
ax_anno.imshow(region_rgba, aspect='auto',
               extent=[0, 10, (anno_depth_bins[-1] + y_ref) / 1000, (anno_depth_bins[0] + y_ref) / 1000])

ax_anno.invert_yaxis()

ax_anno.spines['right'].set_visible(False)
ax_anno.spines['top'].set_visible(False)
ax_anno.spines['bottom'].set_visible(False)
ax_anno.spines['left'].set_visible(False)

ax_anno.set_xticks([])
ax_anno.set_yticks([])
ax_anno.yaxis.tick_right()
ax_anno.yaxis.set_label_position('right')

#to save
ax_main.tick_params(labelbottom=False) 
ax_main.tick_params(labelleft=False) 
ax_main.set_xlabel("")
ax_main.set_ylabel("")


#%% Figure 2f
def plot_2f(data, axs=None):
    tofitx, tofity,brt_l,brt_h, AIR, SPIKES = data
    max_fit_y=np.round(np.amax(brt_h),1)

    fig = None
    if axs is None:
        fig, axs = plt.subplots(subplot_kw={'projection': 'polar'})
    axs.plot(np.append(tofitx,tofitx[0]), np.append(tofity,tofity[0]), color='green')
    axs.fill_between(np.append(tofitx, tofitx[0]), np.append(brt_l, brt_l[0]), np.append(brt_h, brt_h[0]), color='green', alpha=0.3)
    axs.set_rmax(max_fit_y)
    axs.set_rticks([])
    axs.grid(False)
    axs.set_xticks([])
    fig =plt.figure(figsize=(10,6))
    for i in range(len(AIR)):
        
        air=AIR[i]
        air=np.max(air)-air
        kernel = np.ones(85) / 85
        air = np.convolve(air, kernel, 'same')
        spikes=SPIKES[i]
        ax =plt.axes((0.1,0.15*(i+1),0.9,0.1))
        ax.plot(np.arange(0,len(air)/25000,1/25000), air,  color='green')
        ax.plot(spikes, np.full_like(spikes,np.max(air)+100), color='k', marker='$I$',linestyle='None', markersize=12)
        ax.set_xlim([0, 2])
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_yticks([])
        ax.set_xticks([])
    
fname = r'.\data\figure_2\data_for_fig_2f.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()
plot_2f(data)

#%% Figure 2g
def plot_2g(data, axs=None):        
    tofitx, tofity, jaw_l, jaw_h, JAW, SPIKES = data
    max_fit_y = np.round(np.amax(jaw_h), 1)

    fig = None
    if axs is None:
        fig, axs = plt.subplots(subplot_kw={'projection': 'polar'})
    axs.plot(np.append(tofitx, tofitx[0]), np.append(tofity, tofity[0]), color='red')
    axs.fill_between(np.append(tofitx, tofitx[0]), np.append(jaw_l, jaw_l[0]), np.append(jaw_h, jaw_h[0]), color='red', alpha=0.3)
    axs.set_rmax(max_fit_y)
    axs.set_rticks([]) 
    axs.grid(False)
    axs.set_xticks([])
 
    
    fig =plt.figure(figsize=(10,6))
    for c, i in enumerate(range(len(JAW))):

        jaw_dv=JAW[i]
        jaw_dv=jaw_dv-np.min(jaw_dv)
        spikes=SPIKES[i]
        
        ax =plt.axes((0.1,0.15*(i+1),0.9,0.1))
        ax.plot(np.arange(0,1471*0.0034,0.0034), jaw_dv,  color='red')
        ax.plot(spikes, np.full_like(spikes,np.max(jaw_dv)+0.1), color='k', marker='$I$',linestyle='None', markersize=12)
        ax.set_xlim([0, 2])
        ax.set_ylim([-0.1, 3.5])
        if c>0:
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.set_xticks([])
        else:
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.set_ylabel("mm")
            ax.set_xticks([])
    

fname = r'.\data\figure_2\data_for_fig_2g.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()
plot_2g(data)

#%% Figure 2h
def plot_psth(psth, bins, ax=None):
    """Plot average event aligned PSTH
    
    Args:
        psth (ndarray):        PSTH (events x time bins)
        bins (list-like):       time bins
        ax (figure ax):         axis to plot on
    
    Returns axes handle
    """ 
    if ax==None:
        fig, ax = plt.subplots(figsize=(11,3))
    ax.plot(bins, psth.mean(axis=0), color='k', linewidth=2)
    ax.set_xlabel("Time from swallow (s)")
    ax.set_ylabel("Spikes/s")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(-0.5, 0.5)
    return fig, ax

def plot_raster(event_aligned_spikes, ax=None):
    """ """
    if ax==None:
        fig, ax = plt.subplots(figsize=(11,3))
    for i, e in enumerate(event_aligned_spikes):
        ax.plot(e, [i+1]*len(e), 'k.', markersize=9)
        ax.set_xlim(-0.5, 0.5)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_ylabel("swallow #")
        ax.set_xlabel("Time from swallow (s)")
    
    return fig, ax


with open(r'.\data\figure_2\data_for_fig_2h.pkl', 'rb') as f: 
    data = pickle.load(f)
f.close()

bins = data['psth_bins']
rasters = data['rasters'] 
psths = data['psths']

neuron_index = 0 
fig, ax = plot_psth(psths[neuron_index], bins)
fig, ax = plot_raster(rasters[neuron_index])

neuron_index = 1
fig, ax = plot_psth(psths[neuron_index], bins)
fig, ax = plot_raster(rasters[neuron_index])