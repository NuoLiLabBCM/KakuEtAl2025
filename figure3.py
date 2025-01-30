import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import pickle

#%% Figure 3c
with open(r'.\data\figure_3\data_for_fig_3c.pkl', 'rb') as f: 
    data = pickle.load(f)
f.close()

def plot_3c(jaw_trace, tongue_onsets, tongue_offsets, cue_time, t):
    """ """
    fig, ax = plt.subplots(figsize=(10,2.5))
    ax.plot(t, jaw_trace, 'b', alpha=1)
    ax.axvline(cue_time, linestyle='--', color='k')
    ax.axvline(cue_time+2, linestyle='--', color='k')
    #tongue traces
    for onset, offset in zip(tongue_onsets, tongue_offsets):
        ax.plot(t[onset:offset], jaw_trace[onset:offset], 'r', alpha=1) #to plot tongue on top of jaw trace

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(0, 5)
    ax.set_xlabel("Time (s)")  
    ax.set_ylabel("mm")  


#example trace 1
jaw_trace = data[0][0]
tongue_onsets = data[0][1]
tongue_offsets = data[0][2]
cue_time = data[0][3]
t = data[0][4]
plot_3c(jaw_trace, tongue_onsets, tongue_offsets, cue_time, t)

#example trace 2
jaw_trace = data[1][0]
tongue_onsets = data[1][1]
tongue_offsets = data[1][2]
cue_time = data[1][3]
t = data[1][4]
plot_3c(jaw_trace, tongue_onsets, tongue_offsets, cue_time, t)

#lick psth
binning = data[2][2]
psth_mean = data[2][0]
psth_sem = data[2][1]

fig, ax = plt.subplots(figsize=(7.5,5.5))
ax.plot(binning, psth_mean, 'b', linewidth=4)
ax.fill_between(binning, psth_mean + psth_sem, y2=psth_mean - psth_sem, color='b', alpha=0.1, edgecolor='none')
ax.set_xlim(-0.5, 3)
ax.set_ylim((-0.4310985510795149, 9.053069572669811))
ax.axvline(0, linestyle='--', color='k', linewidth=2)
ax.axvline(2, linestyle='--', color='k', linewidth=2)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Time from stimulation onset (s)')
ax.set_ylabel('Licks/s')


#%% Figure 3e 
fname = r'.\data\figure_3\SC licking psth combined.mat'
data = loadmat(fname, simplify_cells=True)
psth_time_bins = data['psth_time_bins'] - 0.5719 #laser on time #1.8654 laser off time
psth_sem = data['psth_ydata']
psth_mean = data['psth_mean_data']


fig, ax = plt.subplots(figsize=(7.5,5.5))
ax.plot(psth_time_bins, psth_mean, 'b', linewidth=4)
ax.fill_between(psth_time_bins, psth_mean + psth_sem, y2=psth_mean - psth_sem, color='b', alpha=0.1, edgecolor='none')
ax.set_ylim((-0.4310985510795149, 9.053069572669811))
ax.set_xlim(-0.5, 3)
ax.axvline(0, linestyle='--', color='k', linewidth=2)
ax.axvline(1.2935, linestyle='--', color='k', linewidth=2)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Time from stimulation onset (s)')
ax.set_ylabel('Licks/s')


#%% Figure 3g
with open(r'.\data\figure_3\data_for_fig_3g.pkl', 'rb') as f: 
    data = pickle.load(f)
f.close()


def plot_3g(jaw_trace, tongue_onsets, tongue_offsets):
    """ """
    t = np.linspace(-0.5, 2.5, jaw_trace.size)
    fig, ax = plt.subplots(figsize=(16,6))
    ax.plot(t, jaw_trace, 'b')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(-0.5, 2.5)
    ax.set_xlabel("Time from stimulation onset (s)")
    ax.set_ylabel("mm")
    for onset, offset in zip(tongue_onsets, tongue_offsets):
        ax.plot(t[onset:offset], jaw_trace[onset:offset], 'r')

#example trace 1
trace = 0
jaw_trace = data[trace][0]
tongue_onsets = data[trace][1]
tongue_offsets = data[trace][2]
plot_3g(jaw_trace, tongue_onsets, tongue_offsets)

#example trace 2
trace = 1
jaw_trace = data[trace][0]
tongue_onsets = data[trace][1]
tongue_offsets = data[trace][2]
plot_3g(jaw_trace, tongue_onsets, tongue_offsets)

#example trace 3
trace = 2
jaw_trace = data[trace][0]
tongue_onsets = data[trace][1]
tongue_offsets = data[trace][2]
plot_3g(jaw_trace, tongue_onsets, tongue_offsets)

#example trace 4
trace = 3
jaw_trace = data[trace][0]
tongue_onsets = data[trace][1]
tongue_offsets = data[trace][2]
plot_3g(jaw_trace, tongue_onsets, tongue_offsets)