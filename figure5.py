import numpy as np
import matplotlib.pyplot as plt
import pickle

def plot_traces(data):
    """ """  
    jaw_trace, tongue_onsets, tongue_offsets, breath_trace, cue_time, t = data  
    fig, ax = plt.subplots(figsize=(10,5))
    #jaw traces
    stim_start = np.where(t>0)[0][0] -1
    ax.plot(t[stim_start:], jaw_trace[stim_start:], 'b', alpha=1)
    ax.plot(t[:stim_start], jaw_trace[:stim_start], 'b', alpha=1)
    ax.axvline(cue_time, linestyle='--', color='k')
    ax.axvline(cue_time+2, linestyle='--', color='k')
    
    #tongue traces
    for onset, offset in zip(tongue_onsets, tongue_offsets):
        ax.plot(t[onset:offset], jaw_trace[onset:offset], 'r', alpha=1) #to plot tongue on top of jaw trace
    
    #breath traces
    ax.plot(t, -1*breath_trace+4, 'g')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.set_xlim(0, 5)
    ax.set_xlabel("Time (s)")  
    ax.set_ylabel("mm")  

#%% Figure 5b
fname = r'.\data\figure_5\data_for_fig_5b.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close() 

plot_traces(data[0])
plot_traces(data[1])

#%% Figure 5c
fname = r'.\data\figure_5\first_phase_final.pkl'
with open(fname, 'rb') as f:  
    data = pickle.load(f)
f.close()

first_phase = np.hstack(data['first_phase_ctrl'])
phase_b = data['phase_b_ctrl'] 

t = np.arange(-100,50)*0.0034
plot_times = [-124, -73, -50, -20]
for plot_time in plot_times:
    n_bins = 16
    tofity, tofitx = np.histogram(first_phase[plot_time,:]-np.pi, range=(-np.pi,np.pi), bins=n_bins,density=True)
    baseline, tofitx = np.histogram(phase_b-np.pi, range=(-np.pi,np.pi), bins=n_bins,density=True)  
    tofitx = tofitx[:-1] + (tofitx[1] - tofitx[0])/2
    tofity = tofity / baseline
    tofity=tofity/np.sum(tofity)
    
    # phase as a function of time
    tofity_b, _ = np.histogram(first_phase[-124,:]-np.pi, range=(-np.pi,np.pi), bins=n_bins,density=True)
    tofity_b = tofity_b / baseline
    tofity_b=tofity_b/np.sum(tofity_b)
    
    max_fit_y = 0.14
    max_fit_y_b = 0.14
 
    fig, axs = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(3,3)) 
    axs.plot(np.append(tofitx, tofitx[0]), np.append(tofity, tofity[0]), color='k', linewidth=5)
    axs.fill_between(np.append(tofitx, tofitx[0]), np.append(np.ones(n_bins)/n_bins, np.ones(1)/n_bins)-0.007, np.append(np.ones(n_bins)/n_bins, np.ones(1)/n_bins)+0.007, color='k', alpha=0.3)
    axs.set_rmax(max_fit_y)
    axs.set_rticks([0, max_fit_y])  # Less radial ticks
    axs.grid(False)
    axs.set_title("{:.3f}s" .format(t[plot_time]))
    axs.set_xticks([])
    axs.set_yticks([])
    
    
baseline, tofitx = np.histogram(phase_b-np.pi, range=(-np.pi,np.pi), bins=n_bins,density=True)
fig, ax = plt.subplots(1, 1, figsize=(8, 4))
phase_time=np.zeros((n_bins,np.shape(first_phase)[0]))
for i in np.arange(np.shape(first_phase)[0]):
    phase_time[:,i]=np.histogram(first_phase[i,:]-np.pi, range=(-np.pi,np.pi), bins=n_bins,density=True)[0]/baseline
    phase_time[:,i]=phase_time[:,i]/np.sum(phase_time[:,i])
pos=ax.imshow(phase_time,extent=[t[0],t[-1],-np.pi,np.pi], aspect=0.03)
ax.set_xlabel('Time before lick onset (s)')
ax.set_ylabel('Phase')
ax.set_xlim(-0.2, t[-1])
fig.colorbar(pos,ax=ax)
#%% Figure 5e
fname = r'.\data\figure_5\data_for_fig_5e.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close() 

plot_traces(data[0])
plot_traces(data[1])

#%% Figure 5f
fname = r'.\data\figure_5\first_phase_final.pkl'
with open(fname, 'rb') as f:  
    data = pickle.load(f)
f.close()

first_phase = np.hstack(data['first_phase_stim'])
phase_b = data['phase_b_stim'] 

t = np.arange(-100,50)*0.0034
plot_times = [-124, -73, -50, -20]
for plot_time in plot_times:
    n_bins = 16
    tofity, tofitx = np.histogram(first_phase[plot_time,:]-np.pi, range=(-np.pi,np.pi), bins=n_bins,density=True)
    baseline, tofitx = np.histogram(phase_b-np.pi, range=(-np.pi,np.pi), bins=n_bins,density=True)  
    tofitx = tofitx[:-1] + (tofitx[1] - tofitx[0])/2
    tofity = tofity / baseline
    tofity=tofity/np.sum(tofity)
    
    # phase as a function of time
    tofity_b, _ = np.histogram(first_phase[-124,:]-np.pi, range=(-np.pi,np.pi), bins=n_bins,density=True)
    tofity_b = tofity_b / baseline
    tofity_b=tofity_b/np.sum(tofity_b)
    
    max_fit_y = 0.14
    max_fit_y_b = 0.14
 
    fig, axs = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(3,3)) 
    axs.plot(np.append(tofitx, tofitx[0]), np.append(tofity, tofity[0]), color='k', linewidth=5)
    axs.fill_between(np.append(tofitx, tofitx[0]), np.append(np.ones(n_bins)/n_bins, np.ones(1)/n_bins)-0.007, np.append(np.ones(n_bins)/n_bins, np.ones(1)/n_bins)+0.007, color='k', alpha=0.3)
    axs.set_rmax(max_fit_y)
    axs.set_rticks([0, max_fit_y])  # Less radial ticks
    axs.grid(False)
    axs.set_title("{:.3f}s" .format(t[plot_time]))
    axs.set_xticks([])
    axs.set_yticks([])


baseline, tofitx = np.histogram(phase_b-np.pi, range=(-np.pi,np.pi), bins=n_bins,density=True)
fig, ax = plt.subplots(1, 1, figsize=(8, 4))
phase_time=np.zeros((n_bins,np.shape(first_phase)[0]))
for i in np.arange(np.shape(first_phase)[0]):
    phase_time[:,i]=np.histogram(first_phase[i,:]-np.pi, range=(-np.pi,np.pi), bins=n_bins,density=True)[0]/baseline
    phase_time[:,i]=phase_time[:,i]/np.sum(phase_time[:,i])
pos=ax.imshow(phase_time,extent=[t[0],t[-1],-np.pi,np.pi], aspect=0.03)
ax.set_xlabel('Time before lick onset (s)')
ax.set_ylabel('Phase')
ax.set_xlim(-0.2, t[-1])
fig.colorbar(pos,ax=ax)