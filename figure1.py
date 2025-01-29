import numpy as np
import scipy
import matplotlib as mpl
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42

import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 48
plt.rcParams['svg.fonttype'] = 'none'
import datajoint as dj
from scipy import stats
import pickle

from pipeline import oralfacial_analysis, experiment, tracking
v_oralfacial_analysis = dj.create_virtual_module('oralfacial_analysis', 'map_v2_oralfacial_analysis')
v_tracking = dj.create_virtual_module('tracking', 'map_v2_tracking')

#%% figure specific functions
def compute_psth(aligned_spikes, bin_size=0.05, window=[-1, 1]):
    """Compute PSTH 
    
    Args:
        aligned_spikes (iterable):    Event aligned spike times. dims: trial x spike times
        bin_size (float):             PSTH binning size in seconds
        window (list):                PSTH start and stop in seconds
    
    Returns time bins, PSTH (trial x bins)
    """
    bins = np.arange(window[0], window[1] + 0.02, 0.02)
    all_histograms = []
    for trial in aligned_spikes:
        hist, _ = np.histogram(trial, bins)
        all_histograms.append(hist)
    all_histograms = np.vstack(all_histograms)/0.02


    return bins[:-1], all_histograms
#%% Figure 1c
def plot_1c(data):
    """ """
    ton_onset = data[0]
    ton_offset = data[1]
    session_traces_b_y = data[2]
    session_traces_s_y = data[3]
    good_breathing = data[4]
    
    # canvas setup
    fig = plt.figure(figsize=(16,5))
    grid = plt.GridSpec(3, 10)
    
    ax_tong = plt.subplot(grid[0:1, 0:9])
    ax_lick = plt.subplot(grid[1:2, 0:9])
    ax_brth = plt.subplot(grid[2:3, 0:9])
      
    t_x=np.arange(len(session_traces_b_y))*0.0034
    
    for i, t in enumerate(ton_onset):
        ax_tong.plot(t_x[ton_onset[i]+1:ton_offset[i]+1],session_traces_b_y[ton_onset[i]+1:ton_offset[i]+1],'r')
        
    ax_tong.spines['right'].set_visible(False)
    ax_tong.spines['top'].set_visible(False)
    ax_tong.spines['bottom'].set_visible(False)
    ax_tong.set_xticks([])
    ax_tong.set_xlim(0, 10)
    
    ax_lick.plot(np.arange(len(session_traces_s_y))*0.0034,session_traces_s_y,'b')
    ax_lick.spines['right'].set_visible(False)
    ax_lick.spines['top'].set_visible(False)
    ax_lick.spines['bottom'].set_visible(False)
    ax_lick.set_xticks([])
    ax_lick.set_xlim(0, 10)
    
    ax_brth.plot(np.arange(len(good_breathing))*0.0034,-good_breathing,'g')
    ax_brth.spines['right'].set_visible(False)
    ax_brth.spines['top'].set_visible(False)
    ax_brth.set_xlim(0, 10)
    ax_brth.set_xlabel("Time (s)")
    return fig
    

fname = r'.\data\figure_1\data_for_fig_1c.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)

#%% Figure 1d
def plot_1d(breathing_data, licking_data):
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.hist(breathing_data, 40,range=[0,12], histtype='step', color='g', linewidth=2, density=True, label='breaths')
    ax.hist(licking_data, 40,range=[0,12], histtype='step', color='r', linewidth=2, density=True, label='licks')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Probability density')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend()
        
    return fig, ax


fname = r'.\data\figure_1\data_for_fig_1d.pkl'
with open(fname, 'rb') as f: 
    data = pickle.load(f)
f.close()
breathing_data = data['breathing_data']
licking_data = data['licking_data']
fig, ax = plot_1d(breathing_data, licking_data)

#%% Figure 1f
def plot_1f1(session_key):
    inspir_onset,lick_onset_time,lick_offset_time,ton_onset=(oralfacial_analysis.MovementTiming & session_key).fetch1('inspiration_onset','lick_onset','lick_offset','tongue_onset')
    
    inspir_onset_l=[] 
    for i,val in enumerate(lick_onset_time):
        inspir_onset_l.append(inspir_onset[(inspir_onset>(lick_onset_time[i]+0.2)) & (inspir_onset<(lick_offset_time[i]-0.2))])
    inspir_onset_l=np.array(inspir_onset_l)
    inspir_onset_l=np.hstack(inspir_onset_l)
    
    licks = [] 
    n_licks = [] 
    
    ibi = []
    ibi2 = []
    lick2_time=[]
    lick_bef_time=[]
    
    max_ibi=np.max(np.diff(inspir_onset))
    
    for i,_ in enumerate(inspir_onset_l[2:-2],2):
        
        lick=ton_onset[(ton_onset > inspir_onset_l[i-2]) & (ton_onset<inspir_onset_l[i+2])]
        
        lick_in=lick[(lick > inspir_onset_l[i]) & (lick<inspir_onset_l[i+1])]
        if lick_in.size != 0: 
            lick_bef=lick[(lick < inspir_onset_l[i])] 
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
    
    return fig



def plot_1f2(session_key):
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    psth_bins,psth_lick_1,psth_lick_2,sem_lick_1,sem_lick_2,peaks_lick_1,peaks_lick_2=(oralfacial_analysis.LickResetBN & session_key).fetch('psth_bins','psth_lick_1','psth_lick_2','sem_lick_1','sem_lick_2','peaks_lick_1','peaks_lick_2')
    ax.plot(psth_bins[0],psth_lick_1[0],'b')
    ax.plot(psth_bins[0],psth_lick_2[0],'k')
    if peaks_lick_1[0].size>0:
        ax.plot(psth_bins[0][psth_bins[0]==peaks_lick_1[0][0]],psth_lick_1[0][psth_bins[0]==peaks_lick_1[0][0]],'b.', markersize=32)
    if peaks_lick_2[0].size>0:
        ax.plot(psth_bins[0][psth_bins[0]==peaks_lick_2[0][0]],psth_lick_2[0][psth_bins[0]==peaks_lick_2[0][0]],'k.', markersize=32)
    ax.fill_between(psth_bins[0], psth_lick_1[0]-sem_lick_1[0],psth_lick_1[0]+sem_lick_1[0], color='b', alpha=0.3)
    ax.fill_between(psth_bins[0], psth_lick_2[0]-sem_lick_2[0],psth_lick_2[0]+sem_lick_2[0], color='k', alpha=0.3)
    ax.set_xlabel('Time from measured inspiration onset (s)')
    ax.set_ylabel('Licks/s')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim([-0.2,0.5])
    y_max=np.max(np.hstack(np.concatenate((psth_lick_1,psth_lick_2))))
    y_max_s=np.max(np.hstack(np.concatenate((sem_lick_1,sem_lick_2))))
    ax.set_ylim([0,y_max+y_max_s])
    return fig



session_key = {'subject_id': 1114, 'session': 3}
fig=plot_1f1(session_key)
fig=plot_1f2(session_key)

#%% Fgiure 1g
def plot_1g(units_mtl_lick):
    """ """
    peaks_lick_1,peaks_lick_2,keys= units_mtl_lick.fetch('peaks_lick_1','peaks_lick_2','KEY')
    lick_time_diff=[]
    key_list=[]
    x_min=-0.08
    x_max=0.08
    for i,_ in enumerate(peaks_lick_1):
        if (peaks_lick_1[i].size>0) & (peaks_lick_2[i].size>0):
            lick_time_diff.append(peaks_lick_2[i][0]-peaks_lick_1[i][0])
            key_list.append(keys[i])
    fig, ax = plt.subplots(1, 1,figsize=(8, 8))
    ax.hist(lick_time_diff,range=(x_min, x_max), bins=9, color='k', histtype='step', linewidth=2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Peak difference (s)')
    ax.set_xlim((x_min,x_max))
    ax.axvline(0, color='gray', linestyle='--')
    
    return fig, ax, lick_time_diff

units_mtl_lick=oralfacial_analysis.LickResetBN
fig, ax, data = plot_1g(units_mtl_lick)
tstat, p = stats.ttest_1samp(data, 0)
print("Lick latency\nT-statistic: {:.3f}\np-value: {:.3f}" .format(tstat, p))

#%% Figure 1h
def plot_1h1(session_key):
    inspir_onset,lick_onset_time,lick_offset_time,ton_onset=(oralfacial_analysis.MovementTiming & session_key).fetch1('inspiration_onset','lick_onset','lick_offset','tongue_onset')
    
    ton_onset_l=[] 
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
        if breath_in.size != 0: 
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

    ax.set_xlim([-0.1,0.5])
    ax.set_ylim([-0.5,i+0.5])
    ax.set_xlabel('Time from measured tongue onset (s)')
    ax.set_ylabel('Lick number')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    first_br=first_br[sorted_indexes]
    mid_pt=round(len(first_br)/2)
    _, pvalue=stats.mannwhitneyu(first_br[:mid_pt],first_br[mid_pt:],alternative='greater')
    
    return fig

def plot_1h2(session_key):
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    psth_bins,psth_breath_1,psth_breath_2,sem_breath_1,sem_breath_2,peaks_breath_1,peaks_breath_2=(oralfacial_analysis.BreathResetB & session_key).fetch('psth_bins','psth_breath_1','psth_breath_2','sem_breath_1','sem_breath_2','peaks_breath_1','peaks_breath_2')
    ax.plot(psth_bins[0],psth_breath_1[0],'b')
    ax.plot(psth_bins[0],psth_breath_2[0],'k')
    if peaks_breath_1[0].size>0:
        ax.plot(psth_bins[0][psth_bins[0]==peaks_breath_1[0][0]],psth_breath_1[0][psth_bins[0]==peaks_breath_1[0][0]],'b.', markersize=32)
    if peaks_breath_2[0].size>0:
        ax.plot(psth_bins[0][psth_bins[0]==peaks_breath_2[0][0]],psth_breath_2[0][psth_bins[0]==peaks_breath_2[0][0]],'k.', markersize=32)
    ax.fill_between(psth_bins[0], psth_breath_1[0]-sem_breath_1[0],psth_breath_1[0]+sem_breath_1[0], color='b', alpha=0.3)
    ax.fill_between(psth_bins[0], psth_breath_2[0]-sem_breath_2[0],psth_breath_2[0]+sem_breath_2[0], color='k', alpha=0.3)
    ax.set_xlabel('Time from measured tongue onset (s)')
    ax.set_ylabel('Breaths/s')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim([0,0.5])
    y_max=np.max(np.hstack(np.concatenate((psth_breath_1,psth_breath_2))))
    y_max_s=np.max(np.hstack(np.concatenate((sem_breath_1,sem_breath_2))))
    ax.set_ylim([0,y_max+y_max_s])
    return fig


session_key = {'subject_id': 1114, 'session': 3}
fig=plot_1h1(session_key)
fig=plot_1h2(session_key)

#%% Figure 1i
def plot_1i(units_mtl_breath):
    """ """
    peaks_breath_1,peaks_breath_2,keys= units_mtl_breath.fetch('peaks_breath_1','peaks_breath_2','KEY')
    breath_time_diff=[]
    key_list=[]
    x_min=-0.08
    x_max=0.08
    for i,_ in enumerate(peaks_breath_1):
        if (peaks_breath_1[i].size>0) & (peaks_breath_2[i].size>0):
            breath_time_diff.append(peaks_breath_2[i][0]-peaks_breath_1[i][0])
            key_list.append(keys[i])
    fig, ax = plt.subplots(1, 1,figsize=(8, 8))
    ax.hist(breath_time_diff,range=(x_min, x_max), bins=10, color='k', histtype='step', linewidth=2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Peak difference (s)')
    ax.set_xlim((x_min,x_max))
    ax.axvline(0, color='gray', linestyle='--')
    
    return fig, ax, breath_time_diff

units_mtl_breath=oralfacial_analysis.BreathResetB
fig, ax, data = plot_1i(units_mtl_breath)
tstat, p = stats.ttest_1samp(data, 0)
print("\n\nLick Breath:\nT-statistic: {:.3f}\np-value: {:.3f}" .format(tstat, p))

#%% Figure 1k
with open(r'.\data\figure_1\data_for_fig_1k.pkl', 'rb') as f: 
    licks = pickle.load(f)
f.close()

def plot_1k1(licks, num_after=2, num_before=4, max_ili=1/3, min_ili=1/15):
    """ """
    ili = []
    filtered_licks = []
    for lick in licks:
        if (lick[lick>=0].size >= num_after) & (lick[lick<0].size >= num_before): 
            before_4, before_3, before_2, before_1 = lick[lick<0][-num_before:]
            after_1, after_2 = lick[lick>=0][:num_after]
            _ili = after_1 - before_1
            if (_ili < max_ili) & (_ili > min_ili): 
                ili.append(_ili) 
                filtered_licks.append(lick)
    ili = np.array(ili)
    sorted_ix = np.argsort(ili) 
    sorted_ix = sorted_ix[::-1] 
    licks_sorted = [filtered_licks[i] for i in sorted_ix]
    
    fig, ax = plt.subplots(figsize=(6.5,8.5))
    for i, lick in enumerate(licks_sorted):
        ax.plot(lick, [i]*lick.size, 'ko', markersize=3)
    ax.set_ylim(-1, len(licks_sorted)+1)
    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel("Time from swallow (s)")
    ax.set_ylabel("Swallow number")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axvline(0, color='r')
    
    return fig, ax

fig, ax = plot_1k1(licks)


def plot_1k2(licks, num_after=2, num_before=4, max_ili=1/3, min_ili=1/15):
    """ """
    ili = []
    filtered_licks = []
    lick_psth = []
    for lick in licks:
        if (lick[lick>=0].size >= num_after) & (lick[lick<0].size >= num_before): 
            before_4, before_3, before_2, before_1 = lick[lick<0][-num_before:]
            after_1, after_2 = lick[lick>=0][:num_after]
            _ili = after_1 - before_1 
            if (_ili < max_ili) & (_ili > min_ili): 
                ili.append(_ili) 
                filtered_licks.append(lick)
                lick_psth.append(np.concatenate([lick[lick<0][-num_before:], lick[lick>=0][:num_after]]))
    ili = np.array(ili)
    sorted_ix = np.argsort(ili) 
    sorted_ix = sorted_ix[::-1] 
    lick_psth_sorted = [lick_psth[i] for i in sorted_ix]
    bins, lick_psth_top = compute_psth(lick_psth_sorted[-len(lick_psth_sorted)//3:], bin_size=0.07, window=[-1, 1])
    bins, lick_psth_bottom = compute_psth(lick_psth_sorted[:(len(lick_psth_sorted)//3)+1], bin_size=0.07, window=[-1, 1])
    
    fig, ax = plt.subplots(figsize=(8,8))
    mean = lick_psth_top.mean(axis=0)
    ax.plot(bins, mean, 'b')
    sem = scipy.stats.sem(lick_psth_top.mean(axis=0))
    ax.fill_between(bins, lick_psth_top.mean(axis=0)-sem, lick_psth_top.mean(axis=0)+sem, color='b', alpha=0.3)

    mean = lick_psth_bottom.mean(axis=0)
    ax.plot(bins, mean, 'k')
    sem = scipy.stats.sem(lick_psth_bottom.mean(axis=0))
    ax.fill_between(bins, lick_psth_bottom.mean(axis=0)-sem, lick_psth_bottom.mean(axis=0)+sem, color='k', alpha=0.3)


    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel("Time measured from swallow (s)")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axvline(0, color='r')
    ax.set_ylabel("Licks/s")
    return fig, ax

fig, ax = plot_1k2(licks, num_after=2, num_before=4, max_ili=1/3, min_ili=1/15)

#%% Figure 1l
def plot_1l(licks, num_after=2, num_before=4, max_ili=1/3, min_ili=1/15, split=3):
    """ """
    ili = []
    filtered_licks = []
    for lick in licks:
        if (lick[lick>=0].size >= num_after) & (lick[lick<0].size >= num_before): 
            before_4, before_3, before_2, before_1 = lick[lick<0][-num_before:]
            after_1, after_2 = lick[lick>=0][:num_after]
            _ili = after_1 - before_1 
            if (_ili < max_ili) & (_ili > min_ili): 
                ili.append(_ili) 
                filtered_licks.append(lick)
    ili = np.array(ili)
    sorted_ix = np.argsort(ili) 
    sorted_ix = sorted_ix[::-1]
    licks_sorted = [filtered_licks[i] for i in sorted_ix]
    licks = licks_sorted.copy()
    
    ili_0, ili_1, ili_neg1,  = [], [], []
    
    
    for lick in licks:
        if (lick[lick>=0].size >= 2) & (lick[lick<0].size >= 3):
            before_3, before_2, before_1 = lick[lick<0][-3:]
            after_1, after_2 = lick[lick>=0][:2]
            ili_0.append(after_1 - before_1) 
            ili_1.append(after_2 - after_1)
            ili_neg1.append(before_1 - before_2)
    ili_0 = np.array(ili_0)
    ili_1 = np.array(ili_1)
    ili_neg1 = np.array(ili_neg1)
    boundary = ili_0.size // split 
    
    start = -0.3
    stop = 0.3
    step = 0.0075 
    bins = np.arange(start, stop + step + step, step) - step/2
    
    del_0neg1 = ili_0[:boundary] - ili_neg1[:boundary]
    del_1neg1 = ili_1[:boundary] - ili_neg1[:boundary]
    
    
    fig2, ax = plt.subplots()
    ax.hist(del_0neg1, bins=bins, range=(bins[0], bins[-1]), density=True, color='black', linewidth=4, histtype='step', label='\u0394(0,-1)')
    ax.hist(del_1neg1, bins=bins, range=(bins[0], bins[-1]), density=True, color='crimson', linewidth=4, histtype='step', label='\u0394(1,-1)')
    ax.set_xlim(-0.1, 0.1)
    ax.axvline(0, color='gray', linestyle='--')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Probability density')
    ax.legend()

    return fig2, del_0neg1, del_1neg1
    
with open(r'.\data\figure_1\data_for_fig_1l.pkl', 'rb') as f: 
    licks = pickle.load(f)
f.close()
fig2, del_0neg1, del_1neg1 = plot_1l(licks, num_after=2, num_before=4, max_ili=1/3, min_ili=1/15, split=3)

#%% Figure 1m
def plot_1m1(licks, num_after=2, num_before=4, max_ili=1/3, min_ili=1/15):
    """ """
    ili = []
    filtered_licks = []
    for lick in licks:
        if (lick[lick>=0].size >= num_after) & (lick[lick<0].size >= num_before): 
            before_4, before_3, before_2, before_1 = lick[lick<0][-num_before:]
            after_1, after_2 = lick[lick>=0][:num_after]
            _ili = after_1 - before_1
            if (_ili < max_ili) & (_ili > min_ili): 
                ili.append(_ili) 
                filtered_licks.append(lick)
    ili = np.array(ili)
    sorted_ix = np.argsort(ili) 
    sorted_ix = sorted_ix[::-1] 
    licks_sorted = [filtered_licks[i] for i in sorted_ix]
    
    fig, ax = plt.subplots(figsize=(6.5,8.5))
    for i, lick in enumerate(licks_sorted):
        ax.plot(lick, [i]*lick.size, 'ko', markersize=3)
    ax.set_ylim(-1, len(licks_sorted)+1)
    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel("Time from swallow (s)")
    ax.set_ylabel("Swallow number")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axvline(0, color='r')
    
    return fig, ax

def plot_1m2(licks, num_after=2, num_before=4, max_ili=1/3, min_ili=1/15):
    """ """
    ili = []
    filtered_licks = []
    lick_psth = []
    for lick in licks:
        if (lick[lick>=0].size >= num_after) & (lick[lick<0].size >= num_before): 
            before_4, before_3, before_2, before_1 = lick[lick<0][-num_before:]
            after_1, after_2 = lick[lick>=0][:num_after]
            _ili = after_1 - before_1 
            if (_ili < max_ili) & (_ili > min_ili): 
                ili.append(_ili) 
                filtered_licks.append(lick)
                lick_psth.append(np.concatenate([lick[lick<0][-num_before:], lick[lick>=0][:num_after]]))
    ili = np.array(ili)
    sorted_ix = np.argsort(ili) 
    sorted_ix = sorted_ix[::-1] 
    lick_psth_sorted = [lick_psth[i] for i in sorted_ix]
    bins, lick_psth_top = compute_psth(lick_psth_sorted[-len(lick_psth_sorted)//3:], bin_size=0.07, window=[-1, 1])
    bins, lick_psth_bottom = compute_psth(lick_psth_sorted[:(len(lick_psth_sorted)//3)+1], bin_size=0.07, window=[-1, 1])
    
    fig, ax = plt.subplots(figsize=(8,8))
    mean = lick_psth_top.mean(axis=0)
    ax.plot(bins, mean, 'b')
    sem = scipy.stats.sem(lick_psth_top.mean(axis=0))
    ax.fill_between(bins, lick_psth_top.mean(axis=0)-sem, lick_psth_top.mean(axis=0)+sem, color='b', alpha=0.3)

    mean = lick_psth_bottom.mean(axis=0)
    ax.plot(bins, mean, 'k')
    sem = scipy.stats.sem(lick_psth_bottom.mean(axis=0))
    ax.fill_between(bins, lick_psth_bottom.mean(axis=0)-sem, lick_psth_bottom.mean(axis=0)+sem, color='k', alpha=0.3)


    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel("Time measured from swallow (s)")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axvline(0, color='r')
    ax.set_ylabel("Breaths/s")
    return fig, ax


with open(r'.\data\figure_1\data_for_fig_1m.pkl', 'rb') as f:  
    breaths = pickle.load(f)
f.close()
fig, ax = plot_1m1(breaths)
fig, ax = plot_1m2(breaths, num_after=2, num_before=4, min_ili=1/15, max_ili=1/2)


#%% Figure 1n
def plot_1n(licks, num_after=2, num_before=4, max_ili=1/2, min_ili=1/15, split=3):
    """ """
    ili = []
    filtered_licks = []
    for lick in licks:
        if (lick[lick>=0].size >= num_after) & (lick[lick<0].size >= num_before):
            before_4, before_3, before_2, before_1 = lick[lick<0][-num_before:]
            after_1, after_2 = lick[lick>=0][:num_after]
            _ili = after_1 - before_1 
            if (_ili < max_ili) & (_ili > min_ili): 
                ili.append(_ili) 
                filtered_licks.append(lick)
    ili = np.array(ili)
    sorted_ix = np.argsort(ili) 
    sorted_ix = sorted_ix[::-1] 
    licks_sorted = [filtered_licks[i] for i in sorted_ix]
    licks = licks_sorted.copy()
    
    ili_0, ili_1, ili_neg1,  = [], [], []
    
    
    for lick in licks:
        if (lick[lick>=0].size >= 2) & (lick[lick<0].size >= 3):
            before_3, before_2, before_1 = lick[lick<0][-3:]
            after_1, after_2 = lick[lick>=0][:2]
            ili_0.append(after_1 - before_1) 
            ili_1.append(after_2 - after_1)
            ili_neg1.append(before_1 - before_2)
    ili_0 = np.array(ili_0)
    ili_1 = np.array(ili_1)
    ili_neg1 = np.array(ili_neg1)

    boundary = ili_0.size // split
    start = -0.3
    stop = 0.3
    step = 0.0075
    bins = np.arange(start, stop + step + step, step) - step/2
    

    del_0neg1 = ili_0[:boundary] - ili_neg1[:boundary]
    del_1neg1 = ili_1[:boundary] - ili_neg1[:boundary]
    
    
    fig2, ax = plt.subplots()
    ax.hist(del_0neg1, bins=bins, range=(bins[0], bins[-1]), density=True, color='black', linewidth=4, histtype='step', label='\u0394(0,-1)')
    ax.hist(del_1neg1, bins=bins, range=(bins[0], bins[-1]), density=True, color='crimson', linewidth=4, histtype='step', label='\u0394(1,-1)')
    ax.set_xlim(-0.2, 0.2)
    ax.axvline(0, color='gray', linestyle='--')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Probability density')
    ax.legend()

    return fig2, del_0neg1, del_1neg1

with open(r'.\data\figure_1\data_for_fig_1n.pkl', 'rb') as f:  
    breaths = pickle.load(f)
f.close()
fig2, del_0neg1, del_1neg1 = plot_1n(breaths, num_after=2, num_before=4, min_ili=1/15, max_ili=1/2, split=3)
