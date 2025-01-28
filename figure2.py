import numpy as np
import matplotlib as mpl
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42

import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 48
plt.rcParams['svg.fonttype'] = 'none'
import datajoint as dj
from scipy import stats

from pipeline import lab, experiment, ephys, tracking, oralfacial_analysis, histology, ccf
v_oralfacial_analysis = dj.create_virtual_module('oralfacial_analysis', 'map_v2_oralfacial_analysis')
v_tracking = dj.create_virtual_module('tracking', 'map_v2_tracking')

from collections import Counter
from PIL import ImageColor

import pickle

#%% Figure 2c
probe_ins={'subject_id': 3452, 'session': 2, 'insertion_number': 1} 
traces_len=1471
bin_width=0.0034

bad_trial_side,bad_trial_bot,miss_trial_side,miss_trial_bot=(oralfacial_analysis.BadVideo & probe_ins).fetch('bad_trial_side','bad_trial_bot','miss_trial_side','miss_trial_bot')
if (bad_trial_side[0] is None):
    bad_trial_side[0]=np.array([0])
if (miss_trial_side[0] is None):
    miss_trial_side[0]=np.array([0])
if (bad_trial_bot[0] is None):
    bad_trial_bot[0]=np.array([0])
if (miss_trial_bot[0] is None):
    miss_trial_bot[0]=np.array([0])    
bad_trials=np.concatenate((bad_trial_side[0],bad_trial_bot[0],miss_trial_side[0],miss_trial_bot[0]))

start_t=199
end_t=203

# tongue
tongue_thr = 0.05
session_traces_b_y=(v_tracking.TongueTracking3DBot & probe_ins & [{'trial': tr} for tr in np.arange(start_t,end_t)]).fetch('tongue_y')
session_traces_s_l = (tracking.Tracking.TongueTracking & {'tracking_device': 'Camera 3'} & probe_ins & [{'trial': tr} for tr in np.arange(start_t,end_t)]).fetch('tongue_likelihood', order_by='trial')
session_traces_b_l = (tracking.Tracking.TongueTracking & {'tracking_device': 'Camera 4'} & probe_ins & [{'trial': tr} for tr in np.arange(start_t,end_t)]).fetch('tongue_likelihood', order_by='trial')
session_traces_b_y = np.vstack(session_traces_b_y)
session_traces_s_l = np.vstack(session_traces_s_l)
session_traces_b_l = np.vstack(session_traces_b_l)
session_traces_t_l = session_traces_b_l
session_traces_t_l[np.where((session_traces_s_l > tongue_thr) & (session_traces_b_l > tongue_thr))] = 1
session_traces_t_l[np.where((session_traces_s_l <= tongue_thr) | (session_traces_b_l <= tongue_thr))] = 0
session_traces_t_l = np.hstack(session_traces_t_l)
session_traces_b_y = np.hstack(session_traces_b_y)
session_traces_b_y=session_traces_b_y*session_traces_t_l
session_traces_b_y = session_traces_b_y - np.min(session_traces_b_y)

session_traces_s_y=(v_tracking.JawTracking3DSid & probe_ins & [{'trial': tr} for tr in np.arange(start_t,end_t)]).fetch('jaw_y')
session_traces_s_y = np.vstack(session_traces_s_y)
session_traces_s_y = np.hstack(session_traces_s_y)
session_traces_s_y = session_traces_s_y - np.min(session_traces_s_y)

# get breathing
breathing, breathing_ts = (experiment.Breathing & probe_ins & [{'trial': tr} for tr in np.arange(start_t,end_t)]).fetch('breathing', 'breathing_timestamps', order_by='trial')
good_breathing = breathing
for i, d in enumerate(breathing):
    good_breathing[i] = d[breathing_ts[i] < traces_len*3.4/1000]
good_breathing = stats.zscore(np.vstack(good_breathing),axis=None)

good_breathing = np.hstack(good_breathing)
# -- moving-average
window_size = int(round(bin_width/(breathing_ts[0][1]-breathing_ts[0][0]),0))  # sample
kernel = np.ones(window_size) / window_size
good_breathing = np.convolve(good_breathing, kernel, 'same')
# -- down-sample
good_breathing = good_breathing[window_size::window_size]
good_breathing = np.reshape(good_breathing,(-1,1))

probe_insertion=ephys.ProbeInsertion & probe_ins

probe_insertion = probe_insertion.proj()

assert histology.InterpolatedShankTrack & probe_insertion

units = (ephys.Unit * lab.ElectrodeConfig.Electrode
         & probe_insertion & 'unit_quality != "all"')
units = (units.proj('spike_times', 'spike_depths', 'unit_posy')
         * ephys.ProbeInsertion.proj()
         * lab.ProbeType.Electrode.proj('shank'))

# ---- ccf region ----
annotated_electrodes = (lab.ElectrodeConfig.Electrode * lab.ProbeType.Electrode
                        * ephys.ProbeInsertion
                        * histology.ElectrodeCCFPosition.ElectrodePosition
                        * ccf.CCFAnnotation * ccf.CCFBrainRegion.proj(..., annotation='region_name')
                        & probe_insertion)
pos_y, ccf_y, color_code = annotated_electrodes.fetch(
    'y_coord', 'ccf_y', 'color_code', order_by='y_coord DESC')

# CCF position of most ventral recording site
last_electrode_site = np.array((histology.InterpolatedShankTrack.DeepestElectrodePoint
                                & probe_insertion).fetch1(
    'ccf_x', 'ccf_y', 'ccf_z'))
# CCF position of the brain surface where this shank crosses
brain_surface_site = np.array((histology.InterpolatedShankTrack.BrainSurfacePoint
                               & probe_insertion).fetch1(
    'ccf_x', 'ccf_y', 'ccf_z'))

# CCF position of most ventral recording site, with respect to the brain surface
y_ref = -np.linalg.norm(last_electrode_site - brain_surface_site)

# ---- spikes ----brain_surface_site
spike_times, spike_depths = units.fetch('spike_times', 'spike_depths', order_by='unit')

trial_start_tim=(experiment.SessionTrial & probe_ins & [{'trial': tr} for tr in np.arange(start_t,end_t+1)]).fetch('start_time')

spike_times = np.hstack(spike_times)
spike_depths = np.hstack(spike_depths)

spike_time_c=[]
spike_depth_c=[]
for i,trial_s in enumerate(trial_start_tim):
    spike_time_c.append(spike_times[(spike_times>trial_s) & (spike_times<(trial_s+5))]-float(trial_s)+5*i-0.6)
    spike_depth_c.append(spike_depths[(spike_times>trial_s) & (spike_times<(trial_s+5))])
spike_times=np.hstack(np.array(spike_time_c,dtype=object))
spike_depths=np.hstack(np.array(spike_depth_c,dtype=object))
# plot
# time-depth 2D histogram
time_bin_count = 10000
depth_bin_count = 200

spike_bins = np.linspace(0, spike_times.max(), time_bin_count)
depth_bins = np.linspace(0, np.nanmax(spike_depths), depth_bin_count)

spk_count, spk_edges, depth_edges = np.histogram2d(spike_times, spike_depths, bins=[spike_bins, depth_bins])
spk_rates = spk_count / np.mean(np.diff(spike_bins))
spk_edges = spk_edges[:-1]
depth_edges = depth_edges[:-1]

# region colorcode, by depths
binned_hexcodes = []

y_spacing = np.abs(np.nanmedian(np.where(np.diff(pos_y)==0, np.nan, np.diff(pos_y))))
anno_depth_bins = np.arange(0, depth_bins[-1], y_spacing)
for s, e in zip(anno_depth_bins[:-1], anno_depth_bins[1:]):
    hexcodes = color_code[np.logical_and(pos_y > s, pos_y <= e)]
    if len(hexcodes):
        binned_hexcodes.append(Counter(hexcodes).most_common()[0][0])
    else:
        binned_hexcodes.append('FFFFFF')

region_rgba = np.array([list(ImageColor.getcolor("#" + chex, "RGBA")) for chex in binned_hexcodes])
region_rgba = np.repeat(region_rgba[:, np.newaxis, :], 10, axis=1)

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

threshold = 0.5 # tongue detection
ton_onset = (session_traces_t_l < threshold) & (np.roll(session_traces_t_l,-1) >= threshold)
ton_offset = (session_traces_t_l > threshold) & (np.roll(session_traces_t_l,-1) <= threshold)
ton_onset=np.concatenate((np.array([0]),np.argwhere(ton_onset)[:-1,0]))
ton_offset=np.argwhere(ton_offset)[:,0]
t_x=np.arange(len(session_traces_b_y))*0.0034

for i, t in enumerate(ton_onset):
    ax_tong.plot(t_x[ton_onset[i]+1:ton_offset[i]+1],session_traces_b_y[ton_onset[i]+1:ton_offset[i]+1],'r')
    
ax_tong.spines['right'].set_visible(False)
ax_tong.spines['top'].set_visible(False)
ax_tong.spines['bottom'].set_visible(False)
ax_tong.set_xticks([])
ax_tong.set_xlim(0, 20)
ax_tong.set_ylabel("mm")

ax_lick.plot(np.arange(len(session_traces_s_y))*0.0034,session_traces_s_y,'b')
ax_lick.spines['right'].set_visible(False)
ax_lick.spines['top'].set_visible(False)
ax_lick.spines['bottom'].set_visible(False)
ax_lick.set_xticks([])
ax_lick.set_xlim(0, 20)
ax_lick.set_ylabel("mm")


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

#%% Figure 2g
def plot_2g(unit_key, axs=None):        
    tofitx, tofity, jaw_l, jaw_h = (oralfacial_analysis.JawTuning() & unit_key).fetch1('jaw_x', 'jaw_y', 'jaw_l', 'jaw_h')
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
    for c, i in enumerate(range(len(trials))):
        trial = trials[i]

        jaw_dv=(v_tracking.JawTracking3DSid & {'trial':trial} & unit_key).fetch('jaw_y')[0][0,:]
        jaw_dv=jaw_dv-np.min(jaw_dv)
        spikes=(ephys.Unit.TrialSpikes & {'trial':trial} & unit_key).fetch('spike_times')[0]
        
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
    

unit_key={'subject_id': 1111,
 'session': 5,
 'insertion_number': 1,
 'clustering_method': 'kilosort2',
 'unit': 100}
trials=[5, 6, 7]

plot_2g(unit_key)

#%% Figure 2f
def plot_2f(unit_key, axs=None):
    tofitx, tofity,brt_l,brt_h = (oralfacial_analysis.BreathingTuning() & unit_key).fetch1('breathing_x','breathing_y','breathing_l','breathing_h')
    max_fit_y=np.round(np.amax(brt_h),1)

    fig = None
    if axs is None:
        fig, axs = plt.subplots(subplot_kw={'projection': 'polar'})
    tofitx=tofitx-np.pi
    axs.plot(np.append(tofitx,tofitx[0]), np.append(tofity,tofity[0]), color='green')
    axs.fill_between(np.append(tofitx, tofitx[0]), np.append(brt_l, brt_l[0]), np.append(brt_h, brt_h[0]), color='green', alpha=0.3)
    axs.set_rmax(max_fit_y)
    axs.set_rticks([])
    axs.grid(False)
    axs.set_xticks([])
    
    fig =plt.figure(figsize=(10,6))
    for i in range(len(trials)):
        trial = trials[i]
        
        # units_mtl=(oralfacial_analysis.JawTuning * ephys.Unit * histology.ElectrodeCCFPosition.ElectrodePosition).fetch('KEY')
        air=(experiment.Breathing & {'trial':trial} & unit_key).fetch('breathing')[0]
        air=np.max(air)-air
        kernel = np.ones(85) / 85
        air = np.convolve(air, kernel, 'same')
        spikes=(ephys.Unit.TrialSpikes & {'trial':trial} & unit_key).fetch('spike_times')[0]
        
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
    

unit_key={'subject_id': 1111,
 'session': 4,
 'insertion_number': 2,
 'clustering_method': 'kilosort2',
 'unit': 215}
trials = [50, 51, 52]

plot_2f(unit_key)

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