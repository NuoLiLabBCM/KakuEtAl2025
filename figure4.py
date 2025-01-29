import pickle
import numpy as np
import matplotlib.pyplot as plt
#%% Figure 4b
def plot_4b(data, tr_type):
    """ """
    
    fig, ax = plt.subplots(figsize=(6.5,8.5))
    LICKS, N_LICKS, IBI = [], [], []
    inspir_onset,lick_onset_time,lick_offset_time,ton_onset,all_trials,trial_types,stim_times = data
    
    RESPONSE_PERIOD = 3
    tr_len = 1471/(1000/3.4)
    stim_trials = all_trials[trial_types == 4]
    non_stim_trials = all_trials[trial_types == 0]
    to_plot = stim_trials if tr_type==4 else non_stim_trials
    io, lon, loff, to = [], [], [], []
    for ix in range(1, all_trials.max()+1):
        if ix in to_plot:
            start = (ix-1) * tr_len
            stop = start + tr_len
            io.append(inspir_onset[(inspir_onset>=start)&(inspir_onset<=stop)])
            idx = np.where((lick_onset_time>=(start)+stim_times[ix-1])&(lick_onset_time<=(start)+stim_times[ix-1]+RESPONSE_PERIOD)) #restrict by response period

            lon.append(lick_onset_time[idx])
            loff.append(lick_offset_time[idx])
            to.append(ton_onset[(ton_onset>=start)&(ton_onset<=stop)])
     
    _x = np.array(list(zip(np.hstack(lon), np.hstack(loff))))
    if _x.size:

        lick_onset_time = _x[:,0]
        lick_offset_time = _x[:,1]
        ton_onset = np.hstack(to)        

        inspir_onset_l=[] # restrict by licking bouts
        for i,val in enumerate(lick_onset_time):
            inspir_onset_l.append(inspir_onset[(inspir_onset>(lick_onset_time[i])) & (inspir_onset<(lick_offset_time[i]))])
    
        inspir_onset_l=np.array(inspir_onset_l)
        inspir_onset_l=np.hstack(inspir_onset_l)
        licks = [] # lick times
        n_licks = [] # number of licks in btw breaths
        
        ibi = []
        ibi2 = []
        lick2_time=[]
        lick_bef_time=[]
        
        max_ibi=np.max(np.diff(inspir_onset)) 
        max_ibi = 0.7 if max_ibi>0.7 else max_ibi 
        
        for i,_ in enumerate(inspir_onset_l[2:-2],2):
            
            lick=ton_onset[(ton_onset > inspir_onset_l[i-2]) & (ton_onset<inspir_onset_l[i+2])]
            
            lick_in=lick[(lick > inspir_onset_l[i]) & (lick<inspir_onset_l[i+1])]
            if lick_in.size != 0:  # only include trials with a lick in between the breaths
                lick_bef=lick[(lick < inspir_onset_l[i])] # find the timing of lick before inspiration onset
                if lick_bef.size != 0:
                    lick_bef=lick_bef[-1]
                licks.append(lick - inspir_onset_l[i])
                lick_bef_time.append(lick_bef - inspir_onset_l[i])
                n_licks.append(lick_in.size)
                lick2_time.append(lick_in[-1] - inspir_onset_l[i])
                ibi.append(inspir_onset_l[i+1] - inspir_onset_l[i])
                ibi2.append(inspir_onset_l[i+2] - inspir_onset_l[i])
        
        ibi=np.array(ibi)
        ibi2=np.array(ibi2)
        n_licks=np.array(n_licks)
        lick2_time=np.array(lick2_time)
        lick_bef_time=np.array(lick_bef_time)
        licks[:] = [ele for i, ele in enumerate(licks) if ibi[i]<max_ibi]
        n_licks=n_licks[ibi<max_ibi]
        ibi2=ibi2[ibi<max_ibi]
        lick_bef_time=lick_bef_time[ibi<max_ibi]
        lick2_time=lick2_time[ibi<max_ibi]
        ibi=ibi[ibi<max_ibi]
        
        idx_all=np.arange(0,len(ibi))
        lick1_rem=np.where((n_licks==1) & (lick_bef_time>-0.05) & (lick_bef_time<0))
        lick2_rem=np.where((n_licks==2) & (lick2_time>(ibi-0.05)) & (lick2_time<ibi))
        idx_keep=np.setdiff1d(idx_all,np.concatenate((lick1_rem[0],lick2_rem[0])))
        
        licks = [licks[i] for i in idx_keep]
        n_licks=n_licks[idx_keep]
        ibi2=ibi2[idx_keep]
        lick_bef_time=lick_bef_time[idx_keep]
        lick2_time=lick2_time[idx_keep]
        ibi=ibi[idx_keep]
        
        sorted_indexes=np.argsort(ibi)
        sorted_indexes=sorted_indexes[::-1]
        
        LICKS.append(licks)
        N_LICKS.append(n_licks)
        IBI.append(ibi)
        
    licks = [a for b in LICKS for a in b]
    n_licks = np.hstack(N_LICKS)
    ibi = np.hstack(IBI)
    
    sorted_indexes=np.argsort(ibi)
    sorted_indexes=sorted_indexes[::-1]
    
    d_bound=(np.mean(ibi[n_licks==2]) + np.mean(ibi[n_licks==1]))/2
    
    plot_db=True
    for i,_ in enumerate(licks):
        if n_licks[sorted_indexes[i]]==1:
            ax.plot(licks[sorted_indexes[i]],i*np.ones(len(licks[sorted_indexes[i]])),'ko',markersize=3)
        elif n_licks[sorted_indexes[i]]==2:
            ax.plot(licks[sorted_indexes[i]],i*np.ones(len(licks[sorted_indexes[i]])),'ko',markersize=3)
        elif n_licks[sorted_indexes[i]]==3:
            ax.plot(licks[sorted_indexes[i]],i*np.ones(len(licks[sorted_indexes[i]])),'ko',markersize=3)
        elif n_licks[sorted_indexes[i]]==4:
            ax.plot(licks[sorted_indexes[i]],i*np.ones(len(licks[sorted_indexes[i]])),'ko',markersize=3)
        else:
            ax.plot(licks[sorted_indexes[i]],i*np.ones(len(licks[sorted_indexes[i]])),'ko',markersize=3)
        ax.plot(0,i,'.r',markersize=4)
        ax.plot(ibi[sorted_indexes[i]],i,'.r',markersize=4)
        if (ibi[sorted_indexes[i]]<d_bound) & plot_db:
            plot_db=False
            
    ax.set_xlim([-0.4,0.7]) 
    ax.set_ylim([-0.5,i+0.5])
    ax.set_xlabel('Time from measured inspiration onset (s)')
    ax.set_ylabel('Breath number')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    
fname = r'.\data\figure_4\data_for_fig_4b1.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close() 
plot_4b(data, 0)

fname = r'.\data\figure_4\data_for_fig_4b2.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close() 
plot_4b(data, 4)

fname = r'.\data\figure_4\data_for_fig_4b3.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()

fig, ax = plt.subplots(1, 1,figsize=(8, 8))
# #control
diff = data[0]
x_min=-0.08
x_max=0.08
_ = ax.hist(diff,range=(x_min, x_max), bins=9, color='k', histtype='step', linewidth=3)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Peak difference (s)')
ax.set_ylabel("# of animals")
ax.set_xlim((x_min,x_max))
ax.axvline(0, color='gray', linestyle='--')

# #stim
diff = data[1]
x_min=-0.08
x_max=0.08
_ = ax.hist(diff,range=(x_min, x_max), bins=9, color='darkturquoise', histtype='step', linewidth=3)

#%% Figure 4c

def plot_4c(data, tr_type):
    """ """
    RESPONSE_PERIOD = 3
    BREATHS, ILI, FIRST_BR = [], [], []
    tr_len = 1471/(1000/3.4)
    inspir_onset,lick_onset_time,lick_offset_time,ton_onset,all_trials,trial_types,stim_times = data
    stim_trials = all_trials[trial_types == 4]
    non_stim_trials = all_trials[trial_types == 0]
    to_plot = stim_trials if tr_type==4 else non_stim_trials
    io, lon, loff, to = [], [], [], []
    for ix in range(1, all_trials.max()+1):
        if ix in to_plot: #select trial type here
            io.append(inspir_onset[(inspir_onset>=(ix-1)*tr_len)&(inspir_onset<=ix*tr_len)])
            idx = np.where((lick_onset_time>=((ix-1)*tr_len)+stim_times[ix-1])&(lick_onset_time<=((ix-1)*tr_len)+stim_times[ix-1]+RESPONSE_PERIOD)) #restrict by response period
            lon.append(lick_onset_time[idx])
            loff.append(lick_offset_time[idx])
            to.append(ton_onset[(ton_onset>=(ix-1)*tr_len)&(ton_onset<=ix*tr_len)])
     
    _x = np.array(list(zip(np.hstack(lon), np.hstack(loff))))
    lick_onset_time = _x[:,0]
    lick_offset_time = _x[:,1]
    ton_onset = np.hstack(to)  
         
    ton_onset_l=[] # restrict by licking bouts
    for i,val in enumerate(lick_onset_time):
        ton_onset_l.append(ton_onset[(ton_onset>(lick_onset_time[i]+0.2)) & (ton_onset<(lick_offset_time[i]-0.2))])
    ton_onset_l=np.array(ton_onset_l)
    ton_onset_l=np.hstack(ton_onset_l)
    
    breaths = []
    

    ili = []
    ili2 = []
    first_br=[]
    
    max_ili=np.max(np.diff(inspir_onset))
    max_ili = max_ili if max_ili <= 1/3 else 1/3 #restrict to 3hz
    
    for i,_ in enumerate(ton_onset_l[2:-2],2):
        
        breath=inspir_onset[(inspir_onset > ton_onset_l[i-2]) & (inspir_onset<ton_onset_l[i+2])]
        
        breath_in=breath[(breath > ton_onset_l[i]) & (breath<ton_onset_l[i+1])]
        if breath_in.size != 0:  # only include trials with a breath in between the licks
            breaths.append(breath - ton_onset_l[i])
            ili.append(ton_onset_l[i+1] - ton_onset_l[i])
            ili2.append(ton_onset_l[i+2] - ton_onset_l[i])
            first_br.append(breath_in[0] - ton_onset_l[i])
    
    #ili = ili[:-1]
    ili=np.array(ili)
    ili2=np.array(ili2)
    first_br=np.array(first_br)
    breaths[:] = [ele for i, ele in enumerate(breaths) if ili[i]<max_ili]
    
    first_br=first_br[ili<max_ili]
    ili2=ili2[ili<max_ili]
    ili=ili[ili<max_ili]
    BREATHS.append(breaths)
    ILI.append(ili)
    FIRST_BR.append(first_br)
    
    ili = np.hstack(ILI)
    breaths = np.hstack(BREATHS)
    first_br = np.hstack(FIRST_BR)
    sorted_indexes=np.argsort(ili)
    sorted_indexes=sorted_indexes[::-1]
    
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 8.5))
    
    for i,_ in enumerate(breaths):
        ax.plot(breaths[sorted_indexes[i]],i*np.ones(len(breaths[sorted_indexes[i]])),'ko',markersize=3)
        ax.plot(0,i,'.r',markersize=4)
        ax.plot(ili[sorted_indexes[i]],i,'.r',markersize=4)
    
    ax.set_xlim([-0.1,0.4])
    ax.set_ylim([-0.5,i+0.5])
    ax.set_xlabel('Time from measured tongue onset (s)')
    ax.set_ylabel('Lick number')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)    

fname = r'.\data\figure_4\data_for_fig_4c1.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close() 
plot_4c(data, 0)

fname = r'.\data\figure_4\data_for_fig_4c2.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close() 
plot_4c(data, 4)

fname = r'.\data\figure_4\data_for_fig_4c3.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()

fig, ax = plt.subplots(1, 1,figsize=(8, 8))
# #control
diff = data[0]
x_min=-0.08
x_max=0.08
ax.hist(diff,range=(x_min, x_max), bins=9, color='k', histtype='step', linewidth=3)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Peak difference (s)')
ax.set_ylabel("# of animals")
ax.set_xlim((x_min,x_max))
_ = ax.axvline(0, color='gray', linestyle='--')

# #stim
diff = data[1]
x_min=-0.08
x_max=0.08
_ = ax.hist(diff,range=(x_min, x_max), bins=9, color='darkturquoise', histtype='step', linewidth=3)

#%% Figure 4e
def plot_4e1(data):
    inspir_onset,lick_onset_time,lick_offset_time,ton_onset = data
    inspir_onset_l=[] # restrict by licking bouts
    for i,val in enumerate(lick_onset_time):
        inspir_onset_l.append(inspir_onset[(inspir_onset>(lick_onset_time[i]+0.2)) & (inspir_onset<(lick_offset_time[i]-0.2))])
    inspir_onset_l=np.array(inspir_onset_l)
    inspir_onset_l=np.hstack(inspir_onset_l)
    
    licks = [] # lick times
    n_licks = [] # number of licks in btw breaths
    
    ibi = []
    ibi2 = []
    lick2_time=[]
    lick_bef_time=[]
    
    max_ibi=np.max(np.diff(inspir_onset))
    
    for i,_ in enumerate(inspir_onset_l[2:-2],2):
        
        lick=ton_onset[(ton_onset > inspir_onset_l[i-2]) & (ton_onset<inspir_onset_l[i+2])]
        
        lick_in=lick[(lick > inspir_onset_l[i]) & (lick<inspir_onset_l[i+1])]
        if lick_in.size != 0:  # only include trials with a lick in between the breaths
            lick_bef=lick[(lick < inspir_onset_l[i])] # find the timing of lick before inspiration onset
            if lick_bef.size != 0:
                lick_bef=lick_bef[-1]
            licks.append(lick - inspir_onset_l[i])
            lick_bef_time.append(lick_bef - inspir_onset_l[i])
            n_licks.append(lick_in.size)
            lick2_time.append(lick_in[-1] - inspir_onset_l[i])
            ibi.append(inspir_onset_l[i+1] - inspir_onset_l[i])
            ibi2.append(inspir_onset_l[i+2] - inspir_onset_l[i])
    
    ibi=np.array(ibi)
    ibi2=np.array(ibi2)
    n_licks=np.array(n_licks)
    lick2_time=np.array(lick2_time)
    lick_bef_time=np.array(lick_bef_time)
    licks[:] = [ele for i, ele in enumerate(licks) if ibi[i]<max_ibi]
    n_licks=n_licks[ibi<max_ibi]
    ibi2=ibi2[ibi<max_ibi]
    lick_bef_time=lick_bef_time[ibi<max_ibi]
    lick2_time=lick2_time[ibi<max_ibi]
    ibi=ibi[ibi<max_ibi]
    
    idx_all=np.arange(0,len(ibi))
    lick1_rem=np.where((n_licks==1) & (lick_bef_time>-0.05) & (lick_bef_time<0))
    lick2_rem=np.where((n_licks==2) & (lick2_time>(ibi-0.05)) & (lick2_time<ibi))
    idx_keep=np.setdiff1d(idx_all,np.concatenate((lick1_rem[0],lick2_rem[0])))

    licks = [licks[i] for i in idx_keep]
    n_licks=n_licks[idx_keep]
    ibi2=ibi2[idx_keep]
    lick_bef_time=lick_bef_time[idx_keep]
    lick2_time=lick2_time[idx_keep]
    ibi=ibi[idx_keep]
    
    sorted_indexes=np.argsort(ibi)
    sorted_indexes=sorted_indexes[::-1]
    
    d_bound=(np.mean(ibi[n_licks==2]) + np.mean(ibi[n_licks==1]))/2
    
    fig, ax = plt.subplots(1, 1, figsize=(6.5,8.5))
    plot_db=False
    for i,_ in enumerate(licks):
        if n_licks[sorted_indexes[i]]==1:
            ax.plot(licks[sorted_indexes[i]],i*np.ones(len(licks[sorted_indexes[i]])),'ko',markersize=3)
        elif n_licks[sorted_indexes[i]]==2:
            ax.plot(licks[sorted_indexes[i]],i*np.ones(len(licks[sorted_indexes[i]])),'ko',markersize=3)
        elif n_licks[sorted_indexes[i]]==3:
            ax.plot(licks[sorted_indexes[i]],i*np.ones(len(licks[sorted_indexes[i]])),'ko',markersize=3)
        elif n_licks[sorted_indexes[i]]==4:
            ax.plot(licks[sorted_indexes[i]],i*np.ones(len(licks[sorted_indexes[i]])),'ko',markersize=3)
        else:
            ax.plot(licks[sorted_indexes[i]],i*np.ones(len(licks[sorted_indexes[i]])),'ko',markersize=3)
        ax.plot(0,i,'.r',markersize=4)
        ax.plot(ibi[sorted_indexes[i]],i,'.r',markersize=4)
        if (ibi[sorted_indexes[i]]<d_bound) & plot_db:
            ax.plot([-0.5,1],[i,i],'k')
            plot_db=False
            
    ax.set_xlim([-0.5,1])
    ax.set_ylim([-0.5,i+0.5])
    ax.set_xlabel('Time from measured inspiration onset (s)')
    ax.set_ylabel('Breath number')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    return fig, ax

def plot_4e2(data):
    """ """
    inspir_onset,lick_onset_time,lick_offset_time,ton_onset,good_spikes,inspir_onset_l = data
    licks = [] # lick times
    n_licks = [] # number of licks in btw breaths
    
    spikes = [] # where the spikes occur
    lick2_time=[]
    ibi = []
    ibi2 = []
    lick_bef_time=[]
    
    max_ibi=np.max(np.diff(inspir_onset))
    
    for i,_ in enumerate(inspir_onset_l[2:-2],2):
        
        lick=ton_onset[(ton_onset > inspir_onset_l[i-2]) & (ton_onset<inspir_onset_l[i+2])]
        
        lick_in=lick[(lick > inspir_onset_l[i]) & (lick<inspir_onset_l[i+1])]
        
        if lick_in.size != 0:  # only include trials with a lick in between the breaths
            lick_bef=lick[(lick < inspir_onset_l[i])] # find the timing of lick before inspiration onset
            lick_bef=lick_bef[-1]
            licks.append(lick - inspir_onset_l[i])
            lick_bef_time.append(lick_bef - inspir_onset_l[i])
            n_licks.append(lick_in.size)
            spike_breath=good_spikes-inspir_onset_l[i]
            spike_breath=spike_breath[spike_breath>-0.5]
            spike_breath=spike_breath[spike_breath<1]
            spikes.append(spike_breath)
            lick2_time.append(lick_in[-1] - inspir_onset_l[i])
            ibi.append(inspir_onset_l[i+1] - inspir_onset_l[i])
            ibi2.append(inspir_onset_l[i+2] - inspir_onset_l[i])
    
    ibi=np.array(ibi)
    ibi2=np.array(ibi2)
    n_licks=np.array(n_licks)
    lick2_time=np.array(lick2_time)
    lick_bef_time=np.array(lick_bef_time)
    licks[:] = [ele for i, ele in enumerate(licks) if ibi[i]<max_ibi]
    spikes[:] = [ele for i, ele in enumerate(spikes) if ibi[i]<max_ibi]
    n_licks=n_licks[ibi<max_ibi]
    ibi2=ibi2[ibi<max_ibi]
    lick_bef_time=lick_bef_time[ibi<max_ibi]
    lick2_time=lick2_time[ibi<max_ibi]
    ibi=ibi[ibi<max_ibi]
    
    idx_all=np.arange(0,len(ibi))
    lick1_rem=np.where((n_licks==1) & (lick_bef_time>-0.05) & (lick_bef_time<0))
    lick2_rem=np.where((n_licks==2) & (lick2_time>(ibi-0.05)) & (lick2_time<ibi))
    idx_keep=np.setdiff1d(idx_all,np.concatenate((lick1_rem[0],lick2_rem[0])))
    
    licks = [licks[i] for i in idx_keep]
    spikes = [spikes[i] for i in idx_keep]
    n_licks=n_licks[idx_keep]
    ibi2=ibi2[idx_keep]
    lick_bef_time=lick_bef_time[idx_keep]
    lick2_time=lick2_time[idx_keep]
    ibi=ibi[idx_keep]
    
    sorted_indexes=np.argsort(ibi)
    sorted_indexes=sorted_indexes[::-1]
    
    fig, ax = plt.subplots(1, 1, figsize=(6.5,8.5))
    
    for i,_ in enumerate(licks):    
        if n_licks[sorted_indexes[i]]==1:
            ax.plot(spikes[sorted_indexes[i]],i*np.ones(len(spikes[sorted_indexes[i]])),'ko',markersize=3)
        elif n_licks[sorted_indexes[i]]==2:
            ax.plot(spikes[sorted_indexes[i]],i*np.ones(len(spikes[sorted_indexes[i]])),'ko',markersize=3)
        elif n_licks[sorted_indexes[i]]==3:
            ax.plot(spikes[sorted_indexes[i]],i*np.ones(len(spikes[sorted_indexes[i]])),'ko',markersize=3)
        elif n_licks[sorted_indexes[i]]==4:
            ax.plot(spikes[sorted_indexes[i]],i*np.ones(len(spikes[sorted_indexes[i]])),'ko',markersize=3)
        ax.plot(0,i,'.r',markersize=4)
        ax.plot(ibi[sorted_indexes[i]],i,'.r',markersize=4)
    ax.set_xlim([-0.5,1])
    ax.set_ylim([-0.5,i+0.5])
    ax.set_xlabel('Time from measured inspiration onset (s)')
    ax.set_ylabel('Breath number')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def plot_4e3(data):
    """ """
    peaks_lick_1,peaks_lick_2 = data
    lick_time_diff=[]
    x_min=-0.08
    x_max=0.08
    for i,_ in enumerate(peaks_lick_1):
        if (peaks_lick_1[i].size>0) & (peaks_lick_2[i].size>0):
            lick_time_diff.append(peaks_lick_2[i][0]-peaks_lick_1[i][0])
    fig, ax = plt.subplots(1, 1,figsize=(8,8))
    ax.hist(lick_time_diff,range=(x_min, x_max), bins=9, color='k', histtype='step', linewidth=3)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Peak difference (s)')
    ax.set_ylabel("# of sessions")
    ax.set_xlim((x_min,x_max))
    ax.axvline(0, color='gray', linestyle='--')

fname = r'.\data\figure_4\data_for_fig_4e1.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()
fig,_=plot_4e1(data)

fname = r'.\data\figure_4\data_for_fig_4e2.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()
plot_4e2(data)

fname = r'.\data\figure_4\data_for_fig_4e3.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()
plot_4e3(data)

#%% Figure 4f
def plot_4f1(data):
    inspir_onset,lick_onset_time,lick_offset_time,ton_onset = data
    ton_onset_l=[] # restrict by licking bouts
    for i,val in enumerate(lick_onset_time):
        ton_onset_l.append(ton_onset[(ton_onset>(lick_onset_time[i]+0.2)) & (ton_onset<(lick_offset_time[i]-0.2))])
    ton_onset_l=np.array(ton_onset_l)
    ton_onset_l=np.hstack(ton_onset_l)
    
    breaths = []

    ili = []
    ili2 = []
    first_br=[]
    
    max_ili=np.max(np.diff(inspir_onset))
    
    for i,_ in enumerate(ton_onset_l[2:-2],2):
        
        breath=inspir_onset[(inspir_onset > ton_onset_l[i-2]) & (inspir_onset<ton_onset_l[i+2])]
        
        breath_in=breath[(breath > ton_onset_l[i]) & (breath<ton_onset_l[i+1])]
        if breath_in.size != 0:  # only include trials with a breath in between the licks
            breaths.append(breath - ton_onset_l[i])

            ili.append(ton_onset_l[i+1] - ton_onset_l[i])
            ili2.append(ton_onset_l[i+2] - ton_onset_l[i])
            first_br.append(breath_in[0] - ton_onset_l[i])
    
    ili=np.array(ili)
    ili2=np.array(ili2)
    first_br=np.array(first_br)
    breaths[:] = [ele for i, ele in enumerate(breaths) if ili[i]<max_ili]
    first_br=first_br[ili<max_ili]
    ili2=ili2[ili<max_ili]
    ili=ili[ili<max_ili]
    
    sorted_indexes=np.argsort(ili)
    sorted_indexes=sorted_indexes[::-1]
    
    fig, ax = plt.subplots(1, 1, figsize=(6.5,8.5))
    
    for i,_ in enumerate(breaths):
        ax.plot(breaths[sorted_indexes[i]],i*np.ones(len(breaths[sorted_indexes[i]])),'ko',markersize=3)
        ax.plot(0,i,'.r',markersize=4)
        ax.plot(ili[sorted_indexes[i]],i,'.r',markersize=4)

    ax.set_xlim([-0.1,0.4])
    ax.set_ylim([-0.5,i+0.5])
    ax.set_xlabel('Time from measured tongue onset (s)')
    ax.set_ylabel('Lick number')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    return fig, ax

def plot_4f2(data):
    """ """
    inspir_onset,lick_onset_time,lick_offset_time,ton_onset,good_spikes,ton_onset_l = data
    breaths = []
    
    ili = []
    ili2 = []
    first_br=[]
    spikes = []
    
    max_ili=np.max(np.diff(inspir_onset))
    
    for i,_ in enumerate(ton_onset_l[2:-2],2):
        
        breath=inspir_onset[(inspir_onset > ton_onset_l[i-2]) & (inspir_onset<ton_onset_l[i+2])]
        
        breath_in=breath[(breath > ton_onset_l[i]) & (breath<ton_onset_l[i+1])]
        if breath_in.size != 0:  # only include trials with a breath in between the licks
            breaths.append(breath - ton_onset_l[i])
            ili.append(ton_onset_l[i+1] - ton_onset_l[i])
            ili2.append(ton_onset_l[i+2] - ton_onset_l[i])
            first_br.append(breath_in[0] - ton_onset_l[i])
            spike_breath=good_spikes-ton_onset_l[i]
            spike_breath=spike_breath[spike_breath>-0.5]
            spike_breath=spike_breath[spike_breath<1]
            spikes.append(spike_breath)
    
    ili=np.array(ili)
    ili2=np.array(ili2)
    first_br=np.array(first_br)
    breaths[:] = [ele for i, ele in enumerate(breaths) if ili[i]<max_ili]
    spikes[:] = [ele for i, ele in enumerate(spikes) if ili[i]<max_ili]
    first_br=first_br[ili<max_ili]
    ili2=ili2[ili<max_ili]
    ili=ili[ili<max_ili]
    
    sorted_indexes=np.argsort(ili)
    sorted_indexes=sorted_indexes[::-1]
    
    fig, ax = plt.subplots(1, 1, figsize=(6.5,8.5))
    
    for i,_ in enumerate(spikes):
        ax.plot(spikes[sorted_indexes[i]],i*np.ones(len(spikes[sorted_indexes[i]])),'ko',markersize=3)
        ax.plot(0,i,'.r',markersize=4)
        ax.plot(ili[sorted_indexes[i]],i,'.r',markersize=4)
    
    ax.set_xlim([-0.1,0.4])
    ax.set_ylim([-0.5,i+0.5])
    ax.set_xlabel('Time from measured tongue onset (s)')
    ax.set_ylabel('Lick number')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

def plot_4f3(data):
    """ """
    peaks_breath_1, peaks_breath_2 = data
    breath_time_diff=[]
    x_min=-0.08
    x_max=0.08
    for i,_ in enumerate(peaks_breath_1):
        if (peaks_breath_1[i].size>0) & (peaks_breath_2[i].size>0):
            breath_time_diff.append(peaks_breath_2[i][0]-peaks_breath_1[i][0])
    fig, ax = plt.subplots(1, 1,figsize=(8,8))
    ax.hist(breath_time_diff,range=(x_min, x_max), bins=10, color='k', histtype='step', linewidth=3)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Peak difference (s)')
    ax.set_ylabel("# of sessions")
    ax.set_xlim((x_min,x_max))
    ax.axvline(0, color='gray', linestyle='--')


fname = r'.\data\figure_4\data_for_fig_4f1.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()
fig, _ = plot_4f1(data)
    
fname = r'.\data\figure_4\data_for_fig_4f2.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()
plot_4f2(data)

fname = r'.\data\figure_4\data_for_fig_4f3.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()
plot_4f3(data)

#%% Figure 4g
def plot_4g1(data):
    """ """
    xcorr,n_trial,xcorr_shuf_m, xcorr_shuf, xcorr_x, breath_unit_key,peak_loc, peak_diff = data
    i = 2
    # Ensure there are enough trials to plot
    if n_trial[i] > 15:
        fig, ax = plt.subplots(figsize=(8, 8))  # Single plot
    
        bin_width = xcorr_x[1] - xcorr_x[0]
        perct5 = np.percentile(xcorr_shuf[i], 5, axis=0)
        perct95 = np.percentile(xcorr_shuf[i], 95, axis=0)
    
        # Plot the 3rd cross-correlation and its shuffled statistics
        ax.bar(xcorr_x, xcorr[i], width=bin_width, color='k')
        ax.plot(xcorr_x, xcorr_shuf_m[i], 'grey', linewidth=4)
        ax.plot(xcorr_x, perct5, color='grey', linestyle='--', linewidth=4)
        ax.plot(xcorr_x, perct95, color='grey', linestyle='--', linewidth=4)
        ax.set_xlim([-0.2, 0.2])
        ax.set_xticks([-0.2, 0, 0.2])
        # Customize the axes and labels
        ax.set_xlabel("Lag (s)")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    
        plt.show()

def plot_corr(data):
    """ """
    lick_units = []
    #negative correlations
    _lick_units = np.array(data['neg_sig'])
    sorting_data = np.squeeze(np.array(data['neg_lag']))
    sorted_ix = np.argsort(sorting_data)
    sorted_ix = sorted_ix[::-1]
    lick_units.append(np.squeeze(_lick_units[sorted_ix]))
    
    #positive correlations
    _lick_units = np.array(data['positive_sig'])
    sorting_data = np.squeeze(np.array(data['positive_lag']))
    sorted_ix = np.argsort(sorting_data)
    sorted_ix = sorted_ix[::-1]
    lick_units.append(np.squeeze(_lick_units[sorted_ix])) 
    line2 = _lick_units.shape[0]
    
    #non significant correlations
    _lick_units = np.array(data['non_sig'])
    lick_units.append(_lick_units) 
    line1 = _lick_units.shape[0]
    line2 = line1 + line2
    lick_units = np.vstack(lick_units)
    
    #plot
    xspan_l=-0.4
    xspan_h=0.4
    minval=-0.05
    maxval=0.05
    fig1, ax = plt.subplots(1, 1, figsize=(16, 10))
    pos=ax.imshow(lick_units,extent=[xspan_l,xspan_h,0,np.size(lick_units,axis=0)], aspect=0.009,vmin=minval,vmax=maxval)
    ax.set_ylabel('Non licking and breathing unit pairs')
    ax.set_xlabel('Lag (s)')
    fig1.colorbar(pos)

fname = r'.\data\figure_4\data_for_fig_4g1.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()
plot_4g1(data)


fname = r'.\data\figure_4\data_for_fig_4g2.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()
plot_corr(data)

#%% Figure 4h
def plot_4h1(data):
    """ """
    xcorr,n_trial,xcorr_shuf_m, xcorr_shuf, xcorr_x, breath_unit_key,peak_loc, peak_diff = data
    i = 3
    # Ensure there are enough trials to plot
    if n_trial[i] > 15:
        fig, ax = plt.subplots(figsize=(8, 8))  # Single plot
    
        bin_width = xcorr_x[1] - xcorr_x[0]
        perct5 = np.percentile(xcorr_shuf[i], 5, axis=0)
        perct95 = np.percentile(xcorr_shuf[i], 95, axis=0)
    
        # Plot the 3rd cross-correlation and its shuffled statistics
        ax.bar(xcorr_x, xcorr[i], width=bin_width, color='k')
        ax.plot(xcorr_x, xcorr_shuf_m[i], 'grey', linewidth=4)
        ax.plot(xcorr_x, perct5, color='grey', linestyle='--', linewidth=4)
        ax.plot(xcorr_x, perct95, color='grey', linestyle='--', linewidth=4)
        ax.set_xlim([-0.2, 0.2])
        ax.set_xticks([-0.2, 0, 0.2])
        # Customize the axes and labels
        ax.set_xlabel("Lag (s)")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    
        plt.show()
    
fname = r'.\data\figure_4\data_for_fig_4h1.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()
plot_4h1(data)

fname = r'.\data\figure_4\data_for_fig_4h2.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()
plot_corr(data)