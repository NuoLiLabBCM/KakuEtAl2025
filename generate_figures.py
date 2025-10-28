from scipy.io import loadmat
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

def volcano_plot_licking(mat, max_ibi=None):
    inspir_onset = mat['processed']['inspir_onset']
    lick_onset_time = mat['processed']['lick_bout_onset']
    lick_offset_time = mat['processed']['lick_bout_offset']
    ton_onset = mat['processed']['tongue_onset']
    
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
    
    max_ibi = np.max(np.diff(inspir_onset)) if max_ibi is None else max_ibi
    
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
    
    return fig, ax


def volcano_plot_breathing(mat):
    inspir_onset = mat['processed']['inspir_onset']
    lick_onset_time = mat['processed']['lick_bout_onset']
    lick_offset_time = mat['processed']['lick_bout_offset']
    ton_onset = mat['processed']['tongue_onset']
    
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

    ax.set_xlim([-0.1,0.4])
    ax.set_ylim([-0.5,i+0.5])
    ax.set_xlabel('Time from measured tongue onset (s)')
    ax.set_ylabel('Lick number')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    return fig, ax

def volcano_plot_licking_unit(mat, unit_key, max_ibi=None):

    trial_len = mat['processed']['trial_len']
    indexes = [i for i, d in enumerate(mat['processed']['movement_timing_unit_keys']) if d == unit_key]
    idx = indexes[0]
    
    good_spikes = mat['processed']['spike_times'][idx].copy()
    for i, d in enumerate(good_spikes):
        good_spikes[i] = d[d<trial_len] + trial_len*i
    good_spikes = np.hstack(good_spikes)
    
    inspir_onset = mat['processed']['inspir_onset']
    lick_onset_time = mat['processed']['lick_bout_onset']
    lick_offset_time = mat['processed']['lick_bout_offset']
    ton_onset = mat['processed']['tongue_onset']    
    
    inspir_onset_l=[]
    for i,val in enumerate(lick_onset_time):
        inspir_onset_l.append(inspir_onset[(inspir_onset>(lick_onset_time[i]+0.2)) & (inspir_onset<(lick_offset_time[i]-0.2))])
    inspir_onset_l=np.array(inspir_onset_l)
    inspir_onset_l=np.hstack(inspir_onset_l)
    
    licks = []
    n_licks = [] 
    
    spikes = [] 
    lick2_time=[]
    ibi = []
    ibi2 = []
    lick_bef_time=[]
    
    max_ibi = np.max(np.diff(inspir_onset)) if max_ibi is None else max_ibi
    
    for i,_ in enumerate(inspir_onset_l[2:-2],2):
        
        lick=ton_onset[(ton_onset > inspir_onset_l[i-2]) & (ton_onset<inspir_onset_l[i+2])]
        
        lick_in=lick[(lick > inspir_onset_l[i]) & (lick<inspir_onset_l[i+1])]
        
        if lick_in.size != 0: 
            lick_bef=lick[(lick < inspir_onset_l[i])] 
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
            ax.plot(spikes[sorted_indexes[i]],i*np.ones(len(spikes[sorted_indexes[i]])),'ko',markersize=2)
        elif n_licks[sorted_indexes[i]]==2:
            ax.plot(spikes[sorted_indexes[i]],i*np.ones(len(spikes[sorted_indexes[i]])),'ko',markersize=2)
        elif n_licks[sorted_indexes[i]]==3:
            ax.plot(spikes[sorted_indexes[i]],i*np.ones(len(spikes[sorted_indexes[i]])),'ko',markersize=2)
        elif n_licks[sorted_indexes[i]]==4:
            ax.plot(spikes[sorted_indexes[i]],i*np.ones(len(spikes[sorted_indexes[i]])),'ko',markersize=2)
        ax.plot(0,i,'.r',markersize=4)
        ax.plot(ibi[sorted_indexes[i]],i,'.r',markersize=4)
    ax.set_xlim([-0.5,1])
    ax.set_ylim([-0.5,i+0.5])
    ax.set_xlabel('Time from measured inspiration onset (s)')
    ax.set_ylabel('Breath number')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    return fig, ax 



def volcano_plot_breath_unit(mat, unit_key):

    trial_len = mat['processed']['trial_len']
    
    indexes = [i for i, d in enumerate(mat['processed']['movement_timing_unit_keys']) if d == unit_key]
    idx = indexes[0]
    good_spikes = mat['processed']['spike_times'][idx].copy()
    for i, d in enumerate(good_spikes):
        good_spikes[i] = d[d<trial_len] + trial_len*i
    good_spikes = np.hstack(good_spikes)
    
    inspir_onset = mat['processed']['inspir_onset']
    lick_onset_time = mat['processed']['lick_bout_onset']
    lick_offset_time = mat['processed']['lick_bout_offset']
    ton_onset = mat['processed']['tongue_onset'] 
    
    ton_onset_l=[] # restrict by licking bouts
    for i,val in enumerate(lick_onset_time):
        ton_onset_l.append(ton_onset[(ton_onset>(lick_onset_time[i]+0.2)) & (ton_onset<(lick_offset_time[i]-0.2))])
    ton_onset_l=np.array(ton_onset_l)
    ton_onset_l=np.hstack(ton_onset_l)
    
    breaths = []
    ili = []
    ili2 = []
    first_br=[]
    spikes = []
    
    max_ili=np.max(np.diff(inspir_onset))
    
    for i,_ in enumerate(ton_onset_l[2:-2],2):
        
        breath=inspir_onset[(inspir_onset > ton_onset_l[i-2]) & (inspir_onset<ton_onset_l[i+2])]
        
        breath_in=breath[(breath > ton_onset_l[i]) & (breath<ton_onset_l[i+1])]
        if breath_in.size != 0: 
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
        ax.plot(spikes[sorted_indexes[i]],i*np.ones(len(spikes[sorted_indexes[i]])),'ko',markersize=2)
        ax.plot(0,i,'.r',markersize=3)
        ax.plot(ili[sorted_indexes[i]],i,'.r',markersize=3)
    
    ax.set_xlim([-0.1,0.4])
    ax.set_ylim([-0.5,i+0.5])
    ax.set_xlabel('Time from measured tongue onset (s)')
    ax.set_ylabel('Lick number')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    return fig, ax

def plot_licking_swallow_summary(licks, num_after=2, num_before=4, max_ili=1/3, min_ili=1/15, split=3):
    """ """
    ili = []
    filtered_licks = []
    for lick in licks:
        lick = np.atleast_1d(lick)
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
    
    
    #bottom 1/3rd
    del_0neg1 = ili_0[:boundary] - ili_neg1[:boundary]
    del_1neg1 = ili_1[:boundary] - ili_neg1[:boundary]
    
    
    fig2, ax = plt.subplots()
    ax.hist(del_0neg1, bins=bins, range=(bins[0], bins[-1]), density=False, color='black', linewidth=4, histtype='step', label='\u0394(0,-1)')
    ax.hist(del_1neg1, bins=bins, range=(bins[0], bins[-1]), density=False, color='crimson', linewidth=4, histtype='step', label='\u0394(1,-1)')
    ax.set_xlim(-0.1, 0.1)
    ax.axvline(0, color='gray', linestyle='--')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('# of licks')
    ax.legend()

def plot_breath_swallow_summary(licks, num_after=2, num_before=4, max_ili=1/2, min_ili=1/15, split=3):
    """ """
    ili = []
    filtered_licks = []
    for lick in licks:
        lick = np.atleast_1d(lick)
        if (lick[lick>=0].size >= num_after) & (lick[lick<0].size >= num_before): 
            before_4, before_3, before_2, before_1 = lick[lick<0][-num_before:]
            after_1, after_2 = lick[lick>=0][:num_after]
            _ili = after_1 - before_1 
            if (_ili < max_ili) & (_ili > min_ili): 
                ili.append(_ili) 
                filtered_licks.append(lick)
    ili = np.array(ili)
    sorted_ix = np.argsort(ili) #sort
    sorted_ix = sorted_ix[::-1] #sort descending
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
    
    #bottom 1/3rd
    del_0neg1 = ili_0[:boundary] - ili_neg1[:boundary]
    del_1neg1 = ili_1[:boundary] - ili_neg1[:boundary]
    
    
    fig2, ax = plt.subplots()
    ax.hist(del_0neg1, bins=bins, range=(bins[0], bins[-1]), density=False, color='black', linewidth=4, histtype='step', label='\u0394(0,-1)')
    ax.hist(del_1neg1, bins=bins, range=(bins[0], bins[-1]), density=False, color='crimson', linewidth=4, histtype='step', label='\u0394(1,-1)')
    ax.set_xlim(-0.2, 0.2)
    ax.axvline(0, color='gray', linestyle='--')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('# of breaths')
    ax.legend()

def plot_summary(export_path):
    """ """
    session_keys = []
    peak_diff_br_behv, peak_diff_br_unit = [], []
    peak_diff_lick_behv, peak_diff_lick_unit = [], []
    
    swallow_aligned_licks, swallow_aligned_breaths = [], []
    
    mat_files = glob.glob(os.path.join(export_path, "*processed.mat"))
    for i, mat_file in enumerate(mat_files):
        print("Working on {} of {}" .format(i+1, len(mat_files)))
        mat = loadmat(mat_file, simplify_cells=True)
        breath_tunings = mat['processed']['breath_tuning']
        indexes = np.where(breath_tunings>0.9)[0]
        if (indexes.size>0) & (mat['processed']['peak_diff_breath']['unit']['peaks_breath_1'].size>0) & (mat['processed']['peak_diff_breath']['unit']['peaks_breath_2'].size>0):
            for _i in indexes:
                peaks_breath_1 = np.atleast_1d(np.array(mat['processed']['peak_diff_breath']['unit']['peaks_breath_1'][_i]))
                peaks_breath_2 = np.atleast_1d(np.array(mat['processed']['peak_diff_breath']['unit']['peaks_breath_2'][_i]))
                if (peaks_breath_1.size>0) & (peaks_breath_2.size>0):
                    peak_diff_br_unit.append(peaks_breath_2[0]-peaks_breath_1[0])
        
        jaw_tunings = mat['processed']['jaw_tuning']
        indexes = np.where(jaw_tunings>0.9)[0]
        if (indexes.size>0) & (mat['processed']['peak_diff_lick']['unit']['peaks_lick_1'].size>0) & (mat['processed']['peak_diff_lick']['unit']['peaks_lick_2'].size>0):
            for _i in indexes:
                peaks_lick_1 = np.atleast_1d(np.array(mat['processed']['peak_diff_lick']['unit']['peaks_lick_1'][_i]))
                peaks_lick_2 = np.atleast_1d(np.array(mat['processed']['peak_diff_lick']['unit']['peaks_lick_2'][_i]))
                if (peaks_lick_1.size>0) & (peaks_lick_2.size>0):
                    peak_diff_lick_unit.append(peaks_lick_2[0]-peaks_lick_1[0])
        
        session_key = mat['processed']['session_key']
        if session_key not in session_keys:
            session_keys.append(session_key)
            peaks_breath_1 = np.atleast_1d(np.array(mat['processed']['peak_diff_breath']['behavior']['peaks_breath_1']))
            peaks_breath_2 = np.atleast_1d(np.array(mat['processed']['peak_diff_breath']['behavior']['peaks_breath_2']))
            if (peaks_breath_1.size>0) & (peaks_breath_2.size>0):
                peak_diff_br_behv.append(peaks_breath_2[0]-peaks_breath_1[0])
            peaks_lick_1 = np.atleast_1d(np.array(mat['processed']['peak_diff_lick']['behavior']['peaks_lick_1']))
            peaks_lick_2 = np.atleast_1d(np.array(mat['processed']['peak_diff_lick']['behavior']['peaks_lick_2']))
            if (peaks_lick_1.size>0) & (peaks_lick_2.size>0):
                peak_diff_lick_behv.append(peaks_lick_2[0]-peaks_lick_1[0])
                
            if 'swallow_aligned_licks' in mat['processed'].keys():
                swallow_aligned_licks.append(mat['processed']['swallow_aligned_licks'])
                swallow_aligned_breaths.append(mat['processed']['swallow_aligned_breaths'])
    licks = [r for row in swallow_aligned_licks for r in row]
    plot_licking_swallow_summary(licks)
    breaths = [r for row in swallow_aligned_breaths for r in row]
    plot_breath_swallow_summary(breaths)
    peak_diff_br_behv = np.hstack(peak_diff_br_behv)
    peak_diff_br_unit = np.hstack(peak_diff_br_unit)
    peak_diff_lick_behv = np.hstack(peak_diff_lick_behv)
    peak_diff_lick_unit = np.hstack(peak_diff_lick_unit)
    x_min=-0.08
    x_max=0.08
    fig, ax = plt.subplots(1, 1,figsize=(8, 8))
    ax.hist(peak_diff_lick_behv,range=(x_min, x_max), bins=9, color='k', histtype='step', linewidth=2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Peak difference (s)')
    ax.set_xlim((x_min,x_max))
    ax.axvline(0, color='gray', linestyle='--')
    ax.set_title("Lick peak difference (Behavior)")
    
    fig, ax = plt.subplots(1, 1,figsize=(8, 8))
    ax.hist(peak_diff_br_behv,range=(x_min, x_max), bins=10, color='k', histtype='step', linewidth=2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Peak difference (s)')
    ax.set_xlim((x_min,x_max))
    ax.axvline(0, color='gray', linestyle='--')
    ax.set_title("Breath peak difference (Behavior)")
    
    fig, ax = plt.subplots(1, 1,figsize=(8,8))
    ax.hist(peak_diff_lick_unit,range=(x_min, x_max), bins=9, color='k', histtype='step', linewidth=2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Peak difference (s)')
    ax.set_ylabel("# of units")
    ax.set_xlim((x_min,x_max))
    ax.axvline(0, color='gray', linestyle='--')
    ax.set_title("Lick peak difference (unit)")
    
    fig, ax = plt.subplots(1, 1,figsize=(8,8))
    ax.hist(peak_diff_br_unit,range=(x_min, x_max), bins=10, color='k', histtype='step', linewidth=2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Peak difference (s)')
    ax.set_ylabel("# of units")
    ax.set_xlim((x_min,x_max))
    ax.axvline(0, color='gray', linestyle='--')
    ax.set_title("Breath peak difference (unit)")

def plot_example_traces(mat, trials):
    """ """
    traces_len = mat['processed']['traces_len']
    good_trials = mat['processed']['trials']
    trial_idx = [np.where(good_trials == val)[0][0] for val in trials]
    bin_width=0.0034
    # tongue
    tongue_thr = 0.05
    session_traces_b_y = mat['processed']['TongueTracking3DBot_y'][trial_idx]
    session_traces_s_l_master = mat['processed']['tracking']['camera_3_side']['tongue_likelihood'][trial_idx]
    session_traces_s_l = []
    for i in range(session_traces_s_l_master.shape[0]):
        arr =session_traces_s_l_master[i].copy()
        arr[0] = 0
        arr[-1] = 0
        session_traces_s_l.append(arr)
    session_traces_b_l_master = mat['processed']['tracking']['camera_4_bottom']['tongue_likelihood'][trial_idx]
    session_traces_b_l = []
    for i in range(session_traces_b_l_master.shape[0]):
        arr =session_traces_b_l_master[i].copy()
        arr[0] = 0
        arr[-1] = 0
        session_traces_b_l.append(arr)
    
    session_traces_b_y = np.vstack(session_traces_b_y)
    session_traces_s_l = np.vstack(session_traces_s_l)
    session_traces_b_l = np.vstack(session_traces_b_l)
    session_traces_t_l = session_traces_b_l
    session_traces_t_l[np.where((session_traces_s_l > tongue_thr) & (session_traces_b_l > tongue_thr))] = 1
    session_traces_t_l[np.where((session_traces_s_l <= tongue_thr) | (session_traces_b_l <= tongue_thr))] = 0
    session_traces_t_l = np.hstack(session_traces_t_l)
    session_traces_b_y = np.hstack(session_traces_b_y)
    session_traces_b_y=session_traces_b_y*session_traces_t_l
    session_traces_b_y=stats.zscore(session_traces_b_y,axis=None)
    
    session_traces_s_y = mat['processed']['JawTracking3DSid_y'][trial_idx]
    session_traces_s_y=stats.zscore(np.vstack(session_traces_s_y),axis=None)
    session_traces_s_y = np.hstack(session_traces_s_y)
    
    # get breathing
    breathing = mat['processed']['breathing'][trial_idx]
    breathing_ts = mat['processed']['breathing_ts'][trial_idx]
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
    
    # canvas setup
    fig = plt.figure(figsize=(16,5))
    grid = plt.GridSpec(3, 10)
    
    ax_tong = plt.subplot(grid[0:1, 0:9])
    ax_lick = plt.subplot(grid[1:2, 0:9])
    ax_brth = plt.subplot(grid[2:3, 0:9])
      
    threshold = 0.5 # tongue detection
    ton_onset = (session_traces_t_l < threshold) & (np.roll(session_traces_t_l,-1) >= threshold)
    ton_offset = (session_traces_t_l > threshold) & (np.roll(session_traces_t_l,-1) <= threshold)
    ton_onset=np.argwhere(ton_onset)[:,0]
    ton_offset=np.argwhere(ton_offset)[:,0]
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


def plot_volcano_licking_swallow(input_dir='.', num_after=2, num_before=4, max_ili=1/3, min_ili=1/15):
    """ """
    licks = []
    fpath = os.path.join(input_dir, 'map-export_DL006_20210329_182702_s1_p1_nwb_processed.mat')
    mat = loadmat(fpath, simplify_cells=True)
    _licks = mat['processed']['swallow_aligned_licks']
    licks.append(_licks)
    
    fpath = os.path.join(input_dir, 'map-export_DL006_20210330_162605_s2_p1_nwb_processed.mat')
    mat = loadmat(fpath, simplify_cells=True)
    _licks = mat['processed']['swallow_aligned_licks']
    licks.append(_licks)
    
    fpath = os.path.join(input_dir, 'map-export_DL006_20210331_165022_s3_p1_nwb_processed.mat')
    mat = loadmat(fpath, simplify_cells=True)
    _licks = mat['processed']['swallow_aligned_licks']
    licks.append(_licks)
    
    fpath = os.path.join(input_dir, 'map-export_DL006_20210404_190316_s4_p1_nwb_processed.mat')
    mat = loadmat(fpath, simplify_cells=True)
    _licks = mat['processed']['swallow_aligned_licks']
    licks.append(_licks)
    
    licks = [r for row in licks for r in row]
    
    
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

def plot_volcano_breathing_swallow(input_dir='.', num_after=2, num_before=4, min_ili=1/15, max_ili=1/2):
    """ """
    breaths = []
    fpath = os.path.join(input_dir, 'map-export_DL006_20210329_182702_s1_p1_nwb_processed.mat')
    mat = loadmat(fpath, simplify_cells=True)
    _breaths = mat['processed']['swallow_aligned_breaths']
    breaths.append(_breaths)
    
    fpath = os.path.join(input_dir, 'map-export_DL006_20210330_162605_s2_p1_nwb_processed.mat')
    mat = loadmat(fpath, simplify_cells=True)
    _breaths = mat['processed']['swallow_aligned_breaths']
    breaths.append(_breaths)
    
    fpath = os.path.join(input_dir, 'map-export_DL006_20210331_165022_s3_p1_nwb_processed.mat')
    mat = loadmat(fpath, simplify_cells=True)
    _breaths = mat['processed']['swallow_aligned_breaths']
    breaths.append(_breaths)
    
    fpath = os.path.join(input_dir, 'map-export_DL006_20210404_190316_s4_p1_nwb_processed.mat')
    mat = loadmat(fpath, simplify_cells=True)
    _breaths = mat['processed']['swallow_aligned_breaths']
    breaths.append(_breaths)
    
    breaths = [r for row in breaths for r in row]
    
    
    ili = []
    filtered_licks = []
    for breath in breaths:
        if (breath[breath>=0].size >= num_after) & (breath[breath<0].size >= num_before): 
            before_4, before_3, before_2, before_1 = breath[breath<0][-num_before:]
            after_1, after_2 = breath[breath>=0][:num_after]
            _ili = after_1 - before_1 
            if (_ili < max_ili) & (_ili > min_ili): 
                ili.append(_ili) 
                filtered_licks.append(breath)
    ili = np.array(ili)
    sorted_ix = np.argsort(ili) 
    sorted_ix = sorted_ix[::-1] 
    breaths_sorted = [filtered_licks[i] for i in sorted_ix]
    
    fig, ax = plt.subplots(figsize=(6.5,8.5))
    for i, breath in enumerate(breaths_sorted):
        ax.plot(breath, [i]*breath.size, 'ko', markersize=3)
    ax.set_ylim(-1, len(breaths_sorted)+1)
    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel("Time from swallow (s)")
    ax.set_ylabel("Swallow number")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.axvline(0, color='r')
    
    return fig, ax

def plot_jaw_tuning_red(mat, unit_key, trials, axs=None):
    
    indexes = [i for i, d in enumerate(mat['processed']['movement_timing_unit_keys']) if d == unit_key]
    idx = indexes[0]
    good_trials = mat['processed']['trials']
    trial_idx = [np.where(good_trials == val)[0][0] for val in trials]
    
    tofitx = mat['processed']['jaw_x'][idx]
    tofity = mat['processed']['jaw_y'][idx]
    jaw_l = mat['processed']['jaw_l'][idx]
    jaw_h = mat['processed']['jaw_h'][idx]
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
    for i in range(len(trial_idx)):
        trial = trial_idx[i]

        jaw_dv = mat['processed']['JawTracking3DSid_y'][trial]
        jaw_dv=jaw_dv-np.min(jaw_dv)
        spikes = mat['processed']['spike_times'][idx][trial]
        
        ax =plt.axes((0.1,0.15*(i+1),0.9,0.1))
        ax.plot(np.arange(0,1471*0.0034,0.0034), jaw_dv,  color='red')
        ax.plot(spikes, np.full_like(spikes,np.max(jaw_dv)+0.1), color='k', marker='$I$',linestyle='None', markersize=12)
        ax.set_xlim([0, 2])
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        
def plot_breathing_tuning_green(mat, unit_key, trials, axs=None):
    indexes = [i for i, d in enumerate(mat['processed']['movement_timing_unit_keys']) if d == unit_key]
    idx = indexes[0]
    good_trials = mat['processed']['trials']
    trial_idx = [np.where(good_trials == val)[0][0] for val in trials]
    
    tofitx = mat['processed']['breath_x'][idx]
    tofitx=tofitx-np.pi
    tofity = mat['processed']['breath_y'][idx]
    jaw_l = mat['processed']['breath_l'][idx]
    jaw_h = mat['processed']['breath_h'][idx]
    max_fit_y = np.round(np.amax(jaw_h), 1)

    fig = None
    if axs is None:
        fig, axs = plt.subplots(subplot_kw={'projection': 'polar'})
    axs.plot(np.append(tofitx, tofitx[0]), np.append(tofity, tofity[0]), color='green')
    axs.fill_between(np.append(tofitx, tofitx[0]), np.append(jaw_l, jaw_l[0]), np.append(jaw_h, jaw_h[0]), color='green', alpha=0.3)
    axs.set_rmax(max_fit_y)
    axs.set_rticks([])  
    axs.grid(False)
    axs.set_xticks([])
    
    fig =plt.figure(figsize=(10,6))
    for i in range(len(trial_idx)):
        trial = trial_idx[i]

        air = mat['processed']['breathing'][trial]
        air=np.max(air)-air
        kernel = np.ones(85) / 85
        air = np.convolve(air, kernel, 'same')
        spikes = mat['processed']['spike_times'][idx][trial]
        
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

def plot_swallow_psth_raster(mat, unit_key):
    indexes = [i for i, d in enumerate(mat['processed']['movement_timing_unit_keys']) if d == unit_key]
    idx = indexes[0]
    fig, ax = plt.subplots(figsize=(11,3))
    bins = mat['processed']['swallow_psth_bins']
    psth = mat['processed']['swallow_psth'][idx]
    ax.plot(bins, psth.mean(axis=0), color='k', linewidth=2)
    ax.set_xlabel("Time from swallow (s)")
    ax.set_ylabel("Spikes/s")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(-0.5, 0.5)
    
    event_aligned_spikes = mat['processed']['swallow_raster'][idx]
    fig, ax = plt.subplots(figsize=(11,3))
    for i, e in enumerate(event_aligned_spikes):
        e = np.atleast_1d(e)
        ax.plot(e, [i+1]*len(e), 'k.', markersize=9)
        ax.set_xlim(-0.5, 0.5)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_ylabel("swallow #")
        ax.set_xlabel("Time from swallow (s)")

#%% figure 1 and 4 driver

###### USER INPUT
export_path = r'D:/mat/' #directory containing *_processed MAT files
#####

fpath = os.path.join(export_path, 'map-export_DL025_20210622_150251_s3_p2_nwb_processed.mat')
mat = loadmat(fpath, simplify_cells=True)
fig, ax = volcano_plot_licking(mat)
fig, ax = volcano_plot_breathing(mat)
fig, ax = plot_volcano_licking_swallow(export_path)
fig, ax = plot_volcano_breathing_swallow(export_path)
fpath = os.path.join(export_path, 'map-export_DL006_20210330_162605_s2_p1_nwb_processed.mat')
mat = loadmat(fpath, simplify_cells=True)
trials = np.array([227, 288]) 
fig = plot_example_traces(mat, trials)


fpath = os.path.join(export_path, 'map-export_DL021_20210531_181515_s6_p1_nwb_processed.mat')
mat = loadmat(fpath, simplify_cells=True)
unit_key = {'subject_id': 3452, 'session': 6, 'insertion_number': 1, 'unit': 268, 'clustering_method': 'kilosort2'}
fig, ax = volcano_plot_licking(mat, max_ibi=0.6)
fig, ax = volcano_plot_licking_unit(mat, unit_key, max_ibi=0.6)

fpath = os.path.join(export_path, 'map-export_DL025_20210625_132657_s6_p2_nwb_processed')
mat = loadmat(fpath, simplify_cells=True)
unit_key = {'subject_id': 1114, 'session': 6, 'insertion_number': 2, 'unit': 25, 'clustering_method': 'kilosort2'}
fig, ax = volcano_plot_breathing(mat)
fig, ax = volcano_plot_breath_unit(mat, unit_key)

plot_summary(export_path)


#%% figure 2 driver
unit_key={'subject_id': 1111, 'session': 5, 'insertion_number': 1, 'clustering_method': 'kilosort2', 'unit': 100}
fpath = os.path.join(export_path, 'map-export_DL014_20210515_185216_s5_p1_nwb_processed.mat')
mat = loadmat(fpath, simplify_cells=True)
trials= np.array([5, 6, 7])
plot_jaw_tuning_red(mat, unit_key, trials)
        
unit_key={'subject_id': 1111, 'session': 4, 'insertion_number': 2, 'clustering_method': 'kilosort2', 'unit': 215}
fpath = os.path.join(export_path, 'map-export_DL014_20210514_122757_s4_p2_nwb_processed.mat')
mat = loadmat(fpath, simplify_cells=True)
trials= np.array([50, 51, 52])
plot_breathing_tuning_green(mat, unit_key, trials)

fpath = os.path.join(export_path, 'map-export_DL026_20210628_193921_s5_p2_nwb_processed.mat')
mat = loadmat(fpath, simplify_cells=True)
unit_key = {'subject_id': 3774, 'session': 5, 'insertion_number': 2, 'clustering_method': 'kilosort2', 'unit': 474}
plot_swallow_psth_raster(mat, unit_key)

fpath = os.path.join(export_path, 'map-export_DL021_20210527_154224_s2_p2_nwb_processed.mat')
mat = loadmat(fpath, simplify_cells=True)
unit_key = {'subject_id': 3452, 'session': 2, 'insertion_number': 2, 'clustering_method': 'kilosort2', 'unit': 61}
plot_swallow_psth_raster(mat, unit_key)