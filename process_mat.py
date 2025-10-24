# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 09:59:06 2025

@author: NuoLiLabTower2
"""
import scipy
from scipy.io import loadmat, savemat
from scipy import stats, signal
import numpy as np
import glob
import os
import time


# def get_good_trials(mat):
#     """ Original version (tested and works)"""
#     tracking = mat['tracking']
#     side = tracking['camera_3_side']['jaw_x']
#     bottom = tracking['camera_4_bottom']['jaw_x']
#     trials_side = tracking['camera_3_side']['trialNum']
#     trials_bottom = tracking['camera_4_bottom']['trialNum']
#     trials = np.array(list(set(trials_side) & set(trials_bottom)))
    
#     good_trials, bad_trial_bot, bad_trial_side = [], list(set(trials_side) - set(trials_bottom)), list(set(trials_bottom) - set(trials_side))
#     side = np.array([arr for i, arr in enumerate(side) if i+1 not in bad_trial_bot], dtype=object)
#     bottom = np.array([arr for i, arr in enumerate(bottom) if i+1 not in bad_trial_side], dtype=object)
#     for s, b, tr in zip(side, bottom, trials):
#         if s.size != 1471:
#             bad_trial_side.append(tr)
#         if b.size != 1471:
#             bad_trial_bot.append(tr)
#         if (s.size == 1471) & (b.size == 1471):
#             good_trials.append(tr)
#     return np.array(good_trials), np.array(bad_trial_side), np.array(bad_trial_bot)

def get_good_trials(mat):
    """ To work with Jian's NWB to MAT export"""
    
    tracking = mat['tracking']
    side = tracking['camera_3_side']['jaw_likelihood']
    bottom = tracking['camera_4_bottom']['jaw_likelihood']
    trials_side = tracking['camera_3_side']['trialNum']
    trials_bottom = tracking['camera_4_bottom']['trialNum']
    trials = np.array(list(set(trials_side) & set(trials_bottom)))
    
    good_trials, bad_trial_bot, bad_trial_side = [], list(set(trials_side) - set(trials_bottom)), list(set(trials_bottom) - set(trials_side))
    side = np.array([arr for i, arr in enumerate(side) if i+1 not in bad_trial_bot], dtype=object)
    bottom = np.array([arr for i, arr in enumerate(bottom) if i+1 not in bad_trial_side], dtype=object)
    for s, b, tr in zip(side, bottom, trials):
        side_flag, bottom_flag = True, True
        if (s<0).any():
            side_flag = False
            bad_trial_side.append(tr)
        if (b<0).any():
            bottom_flag = False
            bad_trial_bot.append(tr)
        if (side_flag) & (bottom_flag):
            good_trials.append(tr)
    return np.array(good_trials), np.array(bad_trial_side), np.array(bad_trial_bot)


def get_good_units(mat):
    """ """
    presence_ratio = mat['neuron_unit_quality_control']['presence_ratio']
    amplitude_cutoff = mat['neuron_unit_quality_control']['amplitude_cutoff']
    avg_firing_rate = mat['neuron_unit_quality_control']['avg_firing_rate']
    isi_violation = mat['neuron_unit_quality_control']['isi_violation']
    unit_amp = mat['neuron_unit_quality_control']['unit_amp']
    
    indices = (presence_ratio>0.9) & (amplitude_cutoff<0.15) & (avg_firing_rate>0.2) &  (isi_violation<10) & (unit_amp>150)
    return np.where(indices==True)[0]

def compute_insta_phase_amp(data, fs, freq_band=(5, 15)):
    """
    :param data: trial x time
    :param fs: sampling rate
    :param freq_band: frequency band for bandpass
    """
    orig_ndim = data.ndim
    if orig_ndim > 1:
        trial_count, time_count = data.shape
        # flatten
        data = data.reshape(-1)

    # band pass
    b, a = signal.butter(5, freq_band, btype='band', fs=fs)
    data = signal.filtfilt(b, a, data)
    # hilbert
    analytic_signal = signal.hilbert(data)
    insta_amp = np.abs(analytic_signal)
    insta_phase = np.angle(analytic_signal)

    if orig_ndim > 1:
        return insta_amp.reshape((trial_count, time_count)), insta_phase.reshape((trial_count, time_count))
    else:
        return insta_amp, insta_phase

def get_movement_timing(mat):
    """ """
    bin_width = 0.0034

    # from the cameras
    tongue_thr = 0.95
    session_traces_s_l = mat[foi]['tracking']['camera_3_side']['tongue_likelihood']
    session_traces_b_l = mat[foi]['tracking']['camera_4_bottom']['tongue_likelihood']

    session_traces_s_l = np.vstack(session_traces_s_l)
    session_traces_b_l = np.vstack(session_traces_b_l)
    session_traces_t_l = session_traces_b_l
    session_traces_t_l[np.where((session_traces_s_l > tongue_thr) & (session_traces_b_l > tongue_thr))] = 1
    session_traces_t_l[np.where((session_traces_s_l <= tongue_thr) | (session_traces_b_l <= tongue_thr))] = 0
    session_traces_t_l = np.hstack(session_traces_t_l)

    session_traces_s_l_f = np.vstack(session_traces_s_l)
    session_traces_b_l_f = np.vstack(session_traces_b_l)
    session_traces_t_l_f = session_traces_b_l_f
    session_traces_t_l_f[np.where((session_traces_s_l_f > tongue_thr) & (session_traces_b_l_f > tongue_thr))] = 1
    session_traces_t_l_f[np.where((session_traces_s_l_f <= tongue_thr) | (session_traces_b_l_f <= tongue_thr))] = 0

    # from 3D calibration
    session_traces_s_y = mat[foi]['JawTracking3DSid_y']
    session_traces_s_x = mat[foi]['JawTracking3DSid_x']
    session_traces_s_z = mat[foi]['JawTracking3DSid_z']
    session_traces_b_y = mat[foi]['TongueTracking3DBot_y']
    session_traces_b_x = mat[foi]['TongueTracking3DBot_x']
    session_traces_b_z = mat[foi]['TongueTracking3DBot_z']
    session_traces_s_y = stats.zscore(np.vstack(session_traces_s_y),axis=None)
    session_traces_s_x = stats.zscore(np.vstack(session_traces_s_x),axis=None)
    session_traces_s_z = stats.zscore(np.vstack(session_traces_s_z),axis=None)
    session_traces_b_y = np.vstack(session_traces_b_y)
    traces_y_mean=np.mean(session_traces_b_y[np.where(session_traces_t_l_f == 1)])
    traces_y_std=np.std(session_traces_b_y[np.where(session_traces_t_l_f == 1)])
    session_traces_b_y = (session_traces_b_y - traces_y_mean)/traces_y_std
    session_traces_b_x = np.vstack(session_traces_b_x)
    traces_x_mean=np.mean(session_traces_b_x[np.where(session_traces_t_l_f == 1)])
    traces_x_std=np.std(session_traces_b_x[np.where(session_traces_t_l_f == 1)])
    session_traces_b_x = (session_traces_b_x - traces_x_mean)/traces_x_std
    session_traces_b_z = np.vstack(session_traces_b_z)
    traces_z_mean=np.mean(session_traces_b_z[np.where(session_traces_t_l_f == 1)])
    traces_z_std=np.std(session_traces_b_z[np.where(session_traces_t_l_f == 1)])
    session_traces_b_z = (session_traces_b_z - traces_z_mean)/traces_z_std

    traces_len = np.size(session_traces_b_z, axis = 1)

    # format the video data
    session_traces_s_y = np.hstack(session_traces_s_y)
    session_traces_s_x = np.hstack(session_traces_s_x)
    session_traces_s_z = np.hstack(session_traces_s_z)
    session_traces_b_y = np.hstack(session_traces_b_y)
    session_traces_b_x = np.hstack(session_traces_b_x)
    session_traces_b_z = np.hstack(session_traces_b_z)
    # -- moving-average and down-sample
    window_size = int(bin_width/0.0034)  # sample
    kernel = np.ones(window_size) / window_size
    session_traces_s_x = np.convolve(session_traces_s_x, kernel, 'same')
    session_traces_s_x = session_traces_s_x[window_size::window_size]
    session_traces_s_y = np.convolve(session_traces_s_y, kernel, 'same')
    session_traces_s_y = session_traces_s_y[window_size::window_size]
    session_traces_s_z = np.convolve(session_traces_s_z, kernel, 'same')
    session_traces_s_z = session_traces_s_z[window_size::window_size]
    session_traces_b_x = np.convolve(session_traces_b_x, kernel, 'same')
    session_traces_b_x = session_traces_b_x[window_size::window_size]
    session_traces_b_y = np.convolve(session_traces_b_y, kernel, 'same')
    session_traces_b_y = session_traces_b_y[window_size::window_size]
    session_traces_b_z = np.convolve(session_traces_b_z, kernel, 'same')
    session_traces_b_z = session_traces_b_z[window_size::window_size]
    session_traces_t_l = np.convolve(session_traces_t_l, kernel, 'same')
    session_traces_t_l = session_traces_t_l[window_size::window_size]
    session_traces_t_l[np.where(session_traces_t_l < 1)] = 0
    session_traces_s_x = np.reshape(session_traces_s_x,(-1,1))
    session_traces_s_y = np.reshape(session_traces_s_y,(-1,1))
    session_traces_s_z = np.reshape(session_traces_s_z,(-1,1))
    session_traces_b_x = np.reshape(session_traces_b_x * session_traces_t_l, (-1,1))
    session_traces_b_y = np.reshape(session_traces_b_y * session_traces_t_l, (-1,1))
    session_traces_b_z = np.reshape(session_traces_b_z * session_traces_t_l, (-1,1))

    # get breathing
    breathing = mat[foi]['breathing'].copy()
    breathing_ts = mat[foi]['breathing_ts'].copy()
    good_breathing = breathing.copy()
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

    
    # coordination of movements
    amp_b, phase_b = compute_insta_phase_amp(good_breathing, 1/bin_width, freq_band=(1, 15)) # breathing
    phase_b = phase_b + np.pi

    threshold = np.pi
    cond = (phase_b < threshold) & (np.roll(phase_b,-1) >= threshold)
    inspir_onset=np.argwhere(cond)[:,0]*bin_width # get onset of breath
    a_threshold = -0.5 # amplitude threshold
    a_cond = (good_breathing > a_threshold) & (np.roll(good_breathing,-1) <= a_threshold)
    inspir_amp=np.argwhere(a_cond)[:,0]*bin_width # amp threshold
    inspir_onset_a = [] # only take inspir with a amp crossing
    for i, inspir_value in enumerate(inspir_onset[:-1]):
        if any((inspir_amp>inspir_value) & (inspir_amp<inspir_onset[i+1])):
            inspir_onset_a.append(inspir_value)
    inspir_onset=np.array(inspir_onset_a)

    # licking epochs
    threshold = 0.5 # tongue detection
    a_cond = (session_traces_t_l < threshold) & (np.roll(session_traces_t_l,-1) >= threshold)
    ton_onset=np.argwhere(a_cond)[:,0]*bin_width # get onset of breath
    a_cond = (session_traces_t_l > threshold) & (np.roll(session_traces_t_l,-1) <= threshold)
    ton_offset=np.argwhere(a_cond)[:,0]*bin_width # get onset of breath
    ilf=1/np.diff(ton_onset)

    ton_onset=ton_onset[:-1] 
    ton_offset=ton_offset[:-1]
    f_cond=(ilf>3) & (ilf<9) # lick freq > 3 & < 9
    ton_onset_idx=np.argwhere(f_cond)[:,0] # index of tongue appearance
    lick_onset_idx=[]
    next_lick=np.diff(ton_onset_idx)
    for i,tongue in enumerate(ton_onset_idx[:-2]):
        if (next_lick[i]==1): # num licks > 3
            lick_onset_idx.append(tongue) # index of tongue
    lick_onset_idx=np.array(lick_onset_idx)
    lick_onset_d=np.diff(lick_onset_idx)
    lick_cond_on = np.roll(lick_onset_d,1) >= 2
    lick_cond_off = lick_onset_d >= 2
    lick_bout_onset=np.argwhere(lick_cond_on)[:,0]
    lick_bout_offset=np.argwhere(lick_cond_off)[:,0]
    if lick_bout_onset[0]!=0:
        lick_bout_onset=np.concatenate((np.array([0]),lick_bout_onset)) # add first lick
        lick_bout_offset=np.concatenate((lick_bout_offset,np.array([len(lick_onset_idx)-1]))) # add last lick

    lick_onset_time=ton_onset[lick_onset_idx[lick_bout_onset]] # get onset of licks
    lick_offset_time=ton_onset[lick_onset_idx[lick_bout_offset]+2]
    
    return ton_onset, ton_offset, lick_onset_time, lick_offset_time, inspir_onset

def kuiper_two(data1, data2):
    """Compute the Kuiper statistic to compare two samples.

    Parameters
    ----------
    data1 : array-like
        The first set of data values.
    data2 : array-like
        The second set of data values.

    Returns
    -------
    D : float
        The raw test statistic.
    fpp : float
        The probability of obtaining two samples this different from
        the same distribution.

    .. warning::
        The fpp is quite approximate, especially for small samples.

    """
    data1 = np.sort(data1)
    data2 = np.sort(data2)
    n1, = data1.shape
    n2, = data2.shape
    common_type = np.find_common_type([], [data1.dtype, data2.dtype])
    if not (np.issubdtype(common_type, np.number)
            and not np.issubdtype(common_type, np.complexfloating)):
        raise ValueError('kuiper_two only accepts real inputs')
    # nans, if any, are at the end after sorting.
    if np.isnan(data1[-1]) or np.isnan(data2[-1]):
        raise ValueError('kuiper_two only accepts non-nan inputs')
    D, _ = stats.ks_2samp(np.asarray(data1, common_type),
                        np.asarray(data2, common_type))
    Ne = len(data1) * len(data2) / float(len(data1) + len(data2))
    return D, kuiper_false_positive_probability(D, Ne)


def kuiper_false_positive_probability(D, N):
    """Compute the false positive probability for the Kuiper statistic.

    Uses the set of four formulas described in Paltani 2004; they report
    the resulting function never underestimates the false positive
    probability but can be a bit high in the N=40..50 range.
    (They quote a factor 1.5 at the 1e-7 level.)

    Parameters
    ----------
    D : float
        The Kuiper test score.
    N : float
        The effective sample size.

    Returns
    -------
    fpp : float
        The probability of a score this large arising from the null hypothesis.

    Notes
    -----
    Eq 7 of Paltani 2004 appears to incorrectly quote the original formula
    (Stephens 1965). This function implements the original formula, as it
    produces a result closer to Monte Carlo simulations.

    References
    ----------

    .. [1] Paltani, S., "Searching for periods in X-ray observations using
           Kuiper's test. Application to the ROSAT PSPC archive",
           Astronomy and Astrophysics, v.240, p.789-790, 2004.

    .. [2] Stephens, M. A., "The goodness-of-fit statistic VN: distribution
           and significance points", Biometrika, v.52, p.309, 1965.

    """
    try:
        from scipy.special import factorial, comb
    except ImportError:
        # Retained for backwards compatibility with older versions of scipy
        # (factorial appears to have moved here in 0.14)
        from scipy.misc import factorial, comb

    if D < 0. or D > 2.:
        raise ValueError("Must have 0<=D<=2 by definition of the Kuiper test")

    if D < 2. / N:
        return 1. - factorial(N) * (D - 1. / N)**(N - 1)
    elif D < 3. / N:
        k = -(N * D - 1.) / 2.
        r = np.sqrt(k**2 - (N * D - 2.)**2 / 2.)
        a, b = -k + r, -k - r
        return 1 - (factorial(N - 1) * (b**(N - 1) * (1 - a) - a**(N - 1) * (1 - b))
                    / N**(N - 2) / (b - a))
    elif (D > 0.5 and N % 2 == 0) or (D > (N - 1.) / (2. * N) and N % 2 == 1):
        # NOTE: the upper limit of this sum is taken from Stephens 1965
        t = np.arange(np.floor(N * (1 - D)) + 1)
        y = D + t / N
        Tt = y**(t - 3) * (y**3 * N
                           - y**2 * t * (3 - 2 / N)
                           + y * t * (t - 1) * (3 - 2 / N) / N
                           - t * (t - 1) * (t - 2) / N**2)
        term = Tt * comb(N, t) * (1 - D - t / N)**(N - t - 1)
        return term.sum()
    else:
        z = D * np.sqrt(N)
        # When m*z>18.82 (sqrt(-log(finfo(double))/2)), exp(-2m**2z**2)
        # underflows.  Cutting off just before avoids triggering a (pointless)
        # underflow warning if `under="warn"`.
        ms = np.arange(1, 18.82 / z)
        S1 = (2 * (4 * ms**2 * z**2 - 1) * np.exp(-2 * ms**2 * z**2)).sum()
        S2 = (ms**2 * (4 * ms**2 * z**2 - 3) * np.exp(-2 * ms**2 * z**2)).sum()
        return S1 - 8 * D / 3 * S2

def min_dist(x1, x2):
    minD = np.abs(x1 - x2)
   
    minD1 = np.mod(minD, 2*np.pi)
    
    minD1 = np.array([minD1]) if isinstance(minD1, (float, np.float64)) else minD1

    minD1[minD1 > np.pi] = 2*np.pi - minD1[minD1 > np.pi]
    return minD1

def r_squared(datay, residuals):
    ss_res = np.sum(residuals**2) 
    ss_tot = np.sum((datay-np.mean(datay))**2)   
    r_2 = 1 - (ss_res / ss_tot)
    return r_2

def vonMise_f(x, std, mean, amp, baseline):
    return amp * np.exp(-0.5 * (min_dist(x, np.full_like(x, mean))/std)**2) + baseline

def compute_phase_tuning(datax, datay):
    from scipy import optimize
    max_fit_y=np.amax(datay)
    min_fit_y=np.amin(datay)
    max_idx=np.where(datay == max_fit_y)
    if np.size(max_idx)>1:
        max_idx=max_idx[0][0]
    # fit curve
    try:
        params_h, pcov = optimize.curve_fit(vonMise_f, datax, datay, p0 = [1, datax[max_idx], max_fit_y-min_fit_y, min_fit_y], bounds=(0, [np.pi/2, 2*np.pi, max_fit_y+min_fit_y, max_fit_y]))
    except:
        print('fitting error')
        params_h = [1, datax[max_idx], max_fit_y-min_fit_y, min_fit_y]
    residuals = datay - vonMise_f(datax, *params_h) 
    r_2_h = r_squared(datay, residuals)

    try:
        params_0, pcov = optimize.curve_fit(vonMise_f, datax, datay, p0 = [1, 0, max_fit_y-min_fit_y, min_fit_y], bounds=(0, [np.pi/2, 2*np.pi, max_fit_y+min_fit_y, max_fit_y]))
    except:
        print('fitting error')
        params_0 = [1, datax[max_idx], max_fit_y-min_fit_y, min_fit_y]
    residuals = datay - vonMise_f(datax, *params_0) 
    r_2_0 = r_squared(datay, residuals)
    
    if r_2_0>r_2_h:
        params = params_0
    else:
        params = params_h
        
    preferred_phase=params[1]
    
    r_max=vonMise_f(params[1], params[0], params[1], params[2], params[3])
    r_min=vonMise_f(params[1]+np.pi, params[0], params[1], params[2], params[3])
    # get MI
    modulation_index=(r_max-r_min)/r_max
        
    return preferred_phase, modulation_index[0]    

def get_unit_jaw_tuning(mat):
    """ """
    modulation_index, preferred_phase, tofitx, tofity, jaw_l,  jaw_h, kuiper_test, di_perm = [], [], [], [], [], [], [], []
    num_frame = 1470   
    num_units = mat[foi]['spike_times'].shape[0]
    session_traces = mat[foi]['tracking']['camera_3_side']['jaw_y'].copy()
    good_traces = np.vstack(session_traces)
    
    fs = mat[foi]['video_fs']
    
    amp, phase = compute_insta_phase_amp(good_traces, float(fs), freq_band=(3, 15))
    phase = phase + np.pi
    phase_s = np.hstack(phase)
    
    # compute phase and MI
    for u in range(num_units):
        all_spikes = mat[foi]['spike_times'][u].copy()
        good_spikes = np.array(all_spikes*float(fs)) # get good spikes and convert to indices
        good_spikes = [np.array(d).astype(int) for d in good_spikes] # convert to intergers
    
        for i, d in enumerate(good_spikes):
            good_spikes[i] = d[d < num_frame]
    
        all_phase = []
        for trial_idx in range(len(good_spikes)):
            all_phase.append(phase[trial_idx][good_spikes[trial_idx]])
    
        all_phase=np.hstack(all_phase)
        _, _kuiper_test = kuiper_two(phase_s, all_phase)
                    
        n_bins = 20
        _tofity, _tofitx = np.histogram(all_phase, bins=n_bins)
        baseline, _tofitx = np.histogram(phase_s, bins=n_bins)  
        _tofitx = _tofitx[:-1] + (_tofitx[1] - _tofitx[0])/2
        _tofity = _tofity / baseline * float(fs)
                       
        _preferred_phase,_modulation_index = compute_phase_tuning(_tofitx, _tofity)
        
        n_perm = 100
        n_spk = len(all_phase)
        di_distr = np.zeros(n_perm)
        di_distr_b = np.zeros(n_perm)
        bootstrap = np.zeros([n_perm,n_bins])
        for i_perm in range(n_perm):
            tofity_p, _ = np.histogram(np.random.choice(all_phase, n_spk), bins=n_bins) 
            tofity_p = tofity_p / baseline * float(fs)
            bootstrap[i_perm,:] = tofity_p
            _, di_distr[i_perm] = compute_phase_tuning(_tofitx, tofity_p)
            
            tofity_b, _ = np.histogram(np.random.choice(phase_s, n_spk), bins=n_bins) 
            tofity_b = tofity_b / baseline * float(fs)
            _, di_distr_b[i_perm] = compute_phase_tuning(_tofitx, tofity_b)
            
        _, _di_perm = stats.mannwhitneyu(di_distr,di_distr_b,alternative='greater')
        _jaw_l=np.percentile(bootstrap,5,axis=0)
        _jaw_h=np.percentile(bootstrap,95,axis=0)
        
        modulation_index.append(_modulation_index)
        preferred_phase.append(_preferred_phase)
        tofitx.append(_tofitx)
        tofity.append(_tofity)
        jaw_l.append(_jaw_l)
        jaw_h.append(_jaw_h)
        kuiper_test.append(_kuiper_test)
        di_perm.append(_di_perm)
    
    return np.array(modulation_index), np.array(preferred_phase), np.array(tofitx), np.array(tofity), np.array(jaw_l), np.array(jaw_h), np.array(kuiper_test), np.array(di_perm)
        
def get_unit_breath_tuning(mat):
    """ """
    modulation_index, preferred_phase, tofitx, tofity, breathing_l,  breathing_h, mi_perm = [], [], [], [], [], [], []
    num_units = mat[foi]['spike_times'].shape[0]
    session_traces = mat[foi]['breathing'].copy()
    breathing_ts = mat[foi]['breathing_ts'].copy()
    fs=25000
    ds=100
    good_traces = session_traces.copy()
    for i, d in enumerate(session_traces):
        good_traces[i] = d[breathing_ts[i] < 5][::ds]
    traces_length = [len(d) for d in good_traces]
    good_trial_ind = np.where(np.array(traces_length) == 5*fs/ds)[0]
    good_traces = good_traces[good_trial_ind]
    good_traces = np.vstack(good_traces)
    
    amp, phase = compute_insta_phase_amp(good_traces, float(fs/ds), freq_band=(1, 10))
    phase = phase + np.pi
    phase_s=np.hstack(phase)
    
    # compute phase and MI
    for u in range(num_units):
        all_spikes = mat[foi]['spike_times'][u].copy()
        good_spikes = np.array(all_spikes[good_trial_ind]*float(fs/ds)) # get good spikes and convert to indices
        good_spikes = [np.array(d).astype(int) for d in good_spikes] # convert to intergers
    
        for i, d in enumerate(good_spikes):
            good_spikes[i] = d[d < int(5*fs/ds)]
    
        all_phase = []
        for trial_idx in range(len(good_spikes)):
            all_phase.append(phase[trial_idx][good_spikes[trial_idx]])
        all_phase=np.hstack(all_phase)
        n_bins = 20
        _tofity, _tofitx = np.histogram(all_phase, bins=n_bins)
        baseline, _tofitx = np.histogram(phase_s, bins=n_bins)  
        _tofitx = _tofitx[:-1] + (_tofitx[1] - _tofitx[0])/2
        _tofity = _tofity / baseline * float(fs/ds)
        
        _preferred_phase,_modulation_index = compute_phase_tuning(_tofitx, _tofity)

        n_perm = 100
        n_spk = len(all_phase)
        di_distr = np.zeros(n_perm)
        di_distr_b = np.zeros(n_perm)
        bootstrap = np.zeros([n_perm,n_bins])
        for i_perm in range(n_perm):
            tofity_p, _ = np.histogram(np.random.choice(all_phase, n_spk), bins=n_bins) 
            tofity_p = tofity_p / baseline * float(fs/ds)
            bootstrap[i_perm,:] = tofity_p
            _, di_distr[i_perm] = compute_phase_tuning(_tofitx, tofity_p)
            
            tofity_b, _ = np.histogram(np.random.choice(phase_s, n_spk), bins=n_bins) 
            tofity_b = tofity_b / baseline * float(fs/ds)
            _, di_distr_b[i_perm] = compute_phase_tuning(_tofitx, tofity_b)
            
        _, _mi_perm = stats.mannwhitneyu(di_distr,di_distr_b,alternative='greater')
        _breathing_l=np.percentile(bootstrap,5,axis=0)
        _breathing_h=np.percentile(bootstrap,95,axis=0)

        modulation_index.append(_modulation_index)
        preferred_phase.append(_preferred_phase)
        tofitx.append(_tofitx)
        tofity.append(_tofity)
        breathing_l.append(_breathing_l)
        breathing_h.append(_breathing_h)
        mi_perm.append(_mi_perm)         
    
    return np.array(modulation_index), np.array(preferred_phase), np.array(tofitx), np.array(tofity), np.array(breathing_l), np.array(breathing_h), np.array(mi_perm)

def get_peak_diff_breath_behv(mat):
    """ """
    
    min_trial=100
    psth_s=-0.1
    psth_e=0.6
    psth_bin=np.arange(psth_s,psth_e,0.02)
    
    inspir_onset = mat[foi]['inspir_onset']
    lick_onset_time = mat[foi]['lick_bout_onset']
    lick_offset_time = mat[foi]['lick_bout_offset']
    ton_onset = mat[foi]['tongue_onset']
    
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
            
    first_br=first_br[sorted_indexes]
    breaths=[breaths[i] for i in sorted_indexes]
    ili=ili[sorted_indexes]

    d_bound=np.median(ili)

    psth_1_i= np.where(ili<d_bound)[0]
    psth_2_i= np.where(ili>d_bound)[0]
    
    if (len(psth_1_i)<min_trial) | (len(psth_2_i)<min_trial):
        print(f'Less than {min_trial} trials')
        return np.array([]), np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]) 
    
    psth_breath_1=[breaths[i] for i in psth_1_i]
    psth_breath_1=[np.histogram(i_p,psth_bin)[0] for i_p in psth_breath_1]
    psth_breath_1=np.array(psth_breath_1)
    psth_breath_2=[breaths[i] for i in psth_2_i]
    psth_breath_2=[np.histogram(i_p,psth_bin)[0] for i_p in psth_breath_2]
    psth_breath_2=np.array(psth_breath_2)
    half_bin=(psth_bin[1]-psth_bin[0])/2
    psth_bins=psth_bin[1:]-half_bin
    psth_breath_2=psth_breath_2/(half_bin*2)
    psth_breath_1=psth_breath_1/(half_bin*2)
    
    psth_breath_1_sem=stats.sem(psth_breath_1)
    psth_breath_2_sem=stats.sem(psth_breath_2)
    
    psth_breath_1=np.mean(psth_breath_1,axis=0)
    psth_breath_2=np.mean(psth_breath_2,axis=0)
    
    # peaks_lick_1=np.mean(lick2_time[psth_1_i])
    # peaks_lick_2=np.mean(lick2_time[psth_2_i])
    
    min_height=7
    psth_1_thr=np.mean(psth_breath_1)+np.std(psth_breath_1)
    if psth_1_thr<min_height:
        psth_1_thr=min_height
    peaks_1=signal.find_peaks(psth_breath_1, height=psth_1_thr, distance=0.24/(half_bin*2))[0]
    peaks_breath_1=psth_bins[peaks_1]
    peaks_breath_1=peaks_breath_1[(peaks_breath_1>0) & (peaks_breath_1<0.22)]
    psth_2_thr=np.mean(psth_breath_2)+np.std(psth_breath_2)
    if psth_2_thr<min_height:
        psth_2_thr=min_height
    peaks_2=signal.find_peaks(psth_breath_2, height=psth_2_thr, distance=0.24/(half_bin*2))[0]
    peaks_breath_2=psth_bins[peaks_2]
    peaks_breath_2=peaks_breath_2[(peaks_breath_2>0) & (peaks_breath_2<0.24)]
    
    return psth_breath_1, psth_breath_2, psth_breath_1_sem, psth_breath_2_sem, psth_bins, peaks_breath_1, peaks_breath_2

def get_peak_diff_breath_unit(mat):
    """ """
    
    min_trial=100
    traces_len=1471
    psth_s=-0.1
    psth_e=0.6
    psth_bin=np.arange(psth_s,psth_e,0.02)
    
    inspir_onset = mat[foi]['inspir_onset']
    lick_onset_time = mat[foi]['lick_bout_onset']
    lick_offset_time = mat[foi]['lick_bout_offset']
    ton_onset = mat[foi]['tongue_onset']
    
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
    
    #ili = ili[:-1]
    ili=np.array(ili)
    ili2=np.array(ili2)
    first_br=np.array(first_br)
    breaths[:] = [ele for i, ele in enumerate(breaths) if ili[i]<max_ili]
    first_br=first_br[ili<max_ili]
    ili2=ili2[ili<max_ili]
    ili=ili[ili<max_ili]
    
    sorted_indexes=np.argsort(ili)
    sorted_indexes=sorted_indexes[::-1]
            
    first_br=first_br[sorted_indexes]
    breaths=[breaths[i] for i in sorted_indexes]
    ili=ili[sorted_indexes]

    d_bound=np.median(ili)

    psth_1_i= np.where(ili<d_bound)[0]
    psth_2_i= np.where(ili>d_bound)[0]
    
    if (len(psth_1_i)<min_trial) | (len(psth_2_i)<min_trial):
        print(f'Less than {min_trial} trials')
        return np.array([]), np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]) 
    

    
    PSTH_BREATH_1, PSTH_BREATH_2, PSTH_BREATH_1_SEM, PSTH_BREATH_2_SEM, PSTH_BINS, PEAKS_BREATH_1, PEAKS_BREATH_2 = [], [], [], [], [], [], []
    num_units = mat[foi]['spike_times'].shape[0]
    for u in range(num_units):          
        spikes = [] # where the spikes occur
        breaths = []        
        ili = []
        ili2 = []
        first_br=[]
        
        all_spikes = mat[foi]['spike_times'][u].copy()
        good_spikes = all_spikes.copy() # get good spikes
        for i, d in enumerate(good_spikes):
            d = np.atleast_1d(d)
            good_spikes[i] = d[d < traces_len*3.4/1000]+traces_len*3.4/1000*i
        good_spikes = np.hstack(good_spikes)
        
        for i,_ in enumerate(ton_onset_l[2:-2],2):
            
            breath=inspir_onset[(inspir_onset > ton_onset_l[i-2]) & (inspir_onset<ton_onset_l[i+2])]
            
            breath_in=breath[(breath > ton_onset_l[i]) & (breath<ton_onset_l[i+1])]
            if breath_in.size != 0:  # only include trials with a breath in between the licks
                breaths.append(breath - ton_onset_l[i])
                ili.append(ton_onset_l[i+1] - ton_onset_l[i])
                ili2.append(ton_onset_l[i+2] - ton_onset_l[i])
                first_br.append(breath_in[0] - ton_onset_l[i])
                spike_breath=good_spikes-ton_onset_l[i]
                spike_breath=spike_breath[spike_breath>psth_s]
                spike_breath=spike_breath[spike_breath<psth_e]
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
                          
        
        first_br=first_br[sorted_indexes]
        breaths=[breaths[i] for i in sorted_indexes]
        ili=ili[sorted_indexes]
        spikes=[spikes[i] for i in sorted_indexes]
        
        psth_breath_1=[spikes[i] for i in psth_1_i]   
        psth_breath_1=[np.histogram(i_p,psth_bin)[0] for i_p in psth_breath_1]
        psth_breath_1=np.array(psth_breath_1)
        psth_breath_2=[spikes[i] for i in psth_2_i]
        psth_breath_2=[np.histogram(i_p,psth_bin)[0] for i_p in psth_breath_2]
        psth_breath_2=np.array(psth_breath_2)
        half_bin=(psth_bin[1]-psth_bin[0])/2
        psth_bins=psth_bin[1:]-half_bin
        psth_breath_2=psth_breath_2/(half_bin*2)
        psth_breath_1=psth_breath_1/(half_bin*2)
        
        psth_breath_1_sem=stats.sem(psth_breath_1)
        psth_breath_2_sem=stats.sem(psth_breath_2)
        
        psth_breath_1=np.mean(psth_breath_1,axis=0)
        psth_breath_2=np.mean(psth_breath_2,axis=0)
        
        min_height=10
        psth_1_thr=np.mean(psth_breath_1)+np.std(psth_breath_1)
        if psth_1_thr<min_height:
            psth_1_thr=min_height
        peaks_1=signal.find_peaks(psth_breath_1, height=psth_1_thr, distance=0.24/(half_bin*2))[0]
        peaks_breath_1=psth_bins[peaks_1]
        peaks_breath_1=peaks_breath_1[(peaks_breath_1>0) & (peaks_breath_1<0.24)]
        psth_2_thr=np.mean(psth_breath_2)+np.std(psth_breath_2)
        if psth_2_thr<min_height:
            psth_2_thr=min_height
        peaks_2=signal.find_peaks(psth_breath_2, height=psth_2_thr, distance=0.24/(half_bin*2))[0]
        peaks_breath_2=psth_bins[peaks_2]
        peaks_breath_2=peaks_breath_2[(peaks_breath_2>0.02) & (peaks_breath_2<0.26)]
        
                
        PSTH_BREATH_1.append(psth_breath_1)
        PSTH_BREATH_2.append(psth_breath_2)
        PSTH_BREATH_1_SEM.append(psth_breath_1_sem)
        PSTH_BREATH_2_SEM.append(psth_breath_2_sem)
        PSTH_BINS.append(psth_bins)
        PEAKS_BREATH_1.append(peaks_breath_1)
        PEAKS_BREATH_2.append(peaks_breath_2)
    
    return np.array(PSTH_BREATH_1, dtype=object), np.array(PSTH_BREATH_2, dtype=object), np.array(PSTH_BREATH_1_SEM, dtype=object),  np.array(PSTH_BREATH_2_SEM, dtype=object), np.array(PSTH_BINS, dtype=object),  np.array(PEAKS_BREATH_1, dtype=object), np.array(PEAKS_BREATH_2, dtype=object) 

def get_peak_diff_lick_behv(mat):
    """ """
    min_trial=40
    psth_s=-0.4 #psth start
    psth_e=1 #psth end
    psth_bin=np.arange(psth_s,psth_e,0.02) #psth time bins
    
    inspir_onset = mat[foi]['inspir_onset']
    lick_onset_time = mat[foi]['lick_bout_onset']
    lick_offset_time = mat[foi]['lick_bout_offset']
    ton_onset = mat[foi]['tongue_onset']
    
    inspir_onset_l=[] # restrict by licking bouts
    for i,val in enumerate(lick_onset_time): #for every lick bout onset...
        inspir_onset_l.append(inspir_onset[(inspir_onset>(lick_onset_time[i]+0.2)) & (inspir_onset<(lick_offset_time[i]-0.2))]) #find every breath occuring within a lick bout
    
    inspir_onset_l=np.array(inspir_onset_l)
    inspir_onset_l=np.hstack(inspir_onset_l)
    
    licks = [] # lick times
    n_licks = [] # number of licks in btw breaths
    ibi = []
    lick_bef_time=[]
    lick2_time = []
    
    max_ibi=np.max(np.diff(inspir_onset))
    
    for i,_ in enumerate(inspir_onset_l[2:-2],2): #for each one of those breaths
        
        lick=ton_onset[(ton_onset > inspir_onset_l[i-2]) & (ton_onset<inspir_onset_l[i+2])] #licks within 2 past and 2 future breaths
        
        lick_in=lick[(lick > inspir_onset_l[i]) & (lick<inspir_onset_l[i+1])] #lick between current breath and next breath
        
        if lick_in.size != 0:  # only include trials with a lick in between the breaths
            lick_bef=lick[(lick < inspir_onset_l[i])] # find the timing of lick before inspiration onset
            if len(lick_bef)>0:
                lick_bef=lick_bef[-1] #take last lick breath inspir onset
                licks.append(lick - inspir_onset_l[i]) #align licks to inspir onset
                lick2_time.append(lick_in[-1] - inspir_onset_l[i]) #next lick time
                lick_bef_time.append(lick_bef - inspir_onset_l[i]) #preceding lick time
                n_licks.append(lick_in.size)                
                ibi.append(inspir_onset_l[i+1] - inspir_onset_l[i])
    
    ibi=np.array(ibi)
    n_licks=np.array(n_licks)
    lick_bef_time=np.array(lick_bef_time)
    lick2_time=np.array(lick2_time)
    
    licks[:] = [ele for i, ele in enumerate(licks) if ibi[i]<max_ibi]
    n_licks=n_licks[ibi<max_ibi]
    lick_bef_time=lick_bef_time[ibi<max_ibi]
    lick2_time=lick2_time[ibi<max_ibi]
    ibi=ibi[ibi<max_ibi]
    
    idx_all=np.arange(0,len(ibi))
    lick1_rem=np.where((n_licks==1) & (lick_bef_time>-0.05) & (lick_bef_time<0))
    lick2_rem=np.where((n_licks==2) & (lick2_time>(ibi-0.05)) & (lick2_time<ibi))
    idx_keep=np.setdiff1d(idx_all,np.concatenate((lick1_rem[0],lick2_rem[0])))

    n_licks=n_licks[idx_keep]
    lick_bef_time=lick_bef_time[idx_keep]
    ibi=ibi[idx_keep]
    licks = [licks[i] for i in idx_keep]
    
    sorted_indexes=np.argsort(ibi)
    sorted_indexes=sorted_indexes[::-1]
    
    licks = [licks[i] for i in sorted_indexes]
    n_licks=n_licks[sorted_indexes]
    lick_bef_time=lick_bef_time[sorted_indexes]
    ibi=ibi[sorted_indexes]
    
    d_bound=np.median(ibi[n_licks==2])

    psth_1_i= np.where((ibi<d_bound) & (n_licks==2))[0]
    psth_2_i= np.where((ibi>d_bound) & (n_licks==2))[0]

    
    if (len(psth_1_i)<min_trial) | (len(psth_2_i)<min_trial):
        print(f'Less than {min_trial} trials')
        return np.array([]), np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]) 
    
    psth_lick_1=[licks[i] for i in psth_1_i]
    psth_lick_1=[np.histogram(i_p,psth_bin)[0] for i_p in psth_lick_1]
    psth_lick_1=np.array(psth_lick_1)
    psth_lick_2=[licks[i] for i in psth_2_i]
    psth_lick_2=[np.histogram(i_p,psth_bin)[0] for i_p in psth_lick_2]
    psth_lick_2=np.array(psth_lick_2)
    half_bin=(psth_bin[1]-psth_bin[0])/2
    psth_bins=psth_bin[1:]-half_bin
    psth_lick_2=psth_lick_2/(half_bin*2)
    psth_lick_1=psth_lick_1/(half_bin*2)
    
    psth_lick_1_sem=stats.sem(psth_lick_1)
    psth_lick_2_sem=stats.sem(psth_lick_2)
    
    psth_lick_1=np.mean(psth_lick_1,axis=0)
    psth_lick_2=np.mean(psth_lick_2,axis=0)
    
    
    min_height=5
    psth_1_thr=np.mean(psth_lick_1)+np.std(psth_lick_1)
    if psth_1_thr<min_height:
        psth_1_thr=min_height
    peaks_1=signal.find_peaks(psth_lick_1, height=psth_1_thr, distance=0.14/(half_bin*2))[0]
    peaks_lick_1=psth_bins[peaks_1]
    peaks_lick_1=peaks_lick_1[(peaks_lick_1>0.12) & (peaks_lick_1<0.32)]
    psth_2_thr=np.mean(psth_lick_2)+np.std(psth_lick_2)
    if psth_2_thr<min_height:
        psth_2_thr=min_height
    peaks_2=signal.find_peaks(psth_lick_2, height=psth_2_thr, distance=0.14/(half_bin*2))[0]
    peaks_lick_2=psth_bins[peaks_2]
    peaks_lick_2=peaks_lick_2[(peaks_lick_2>0.16) & (peaks_lick_2<0.36)]
    
    return psth_lick_1, psth_lick_2, psth_lick_1_sem, psth_lick_2_sem, psth_bins, peaks_lick_1, peaks_lick_2


def get_peak_diff_lick_unit(mat):
    """ """
    n_trial=40
    traces_len=1471
    min_trial=40
    psth_s=-0.4
    psth_e=1
    psth_bin=np.arange(psth_s,psth_e,0.02)

    inspir_onset = mat[foi]['inspir_onset']
    lick_onset_time = mat[foi]['lick_bout_onset']
    lick_offset_time = mat[foi]['lick_bout_offset']
    ton_onset = mat[foi]['tongue_onset']

    inspir_onset_l=[] # restrict by licking bouts
    for i,val in enumerate(lick_onset_time):
        inspir_onset_l.append(inspir_onset[(inspir_onset>(lick_onset_time[i]+0.2)) & (inspir_onset<(lick_offset_time[i]-0.2))])
    inspir_onset_l=np.array(inspir_onset_l)
    inspir_onset_l=np.hstack(inspir_onset_l)

    licks = [] # lick times
    n_licks = [] # number of licks in btw breaths
    ibi = []
    lick_bef_time=[]
    lick2_time = []

    max_ibi=np.max(np.diff(inspir_onset))

    for i,_ in enumerate(inspir_onset_l[2:-2],2):
        
        lick=ton_onset[(ton_onset > inspir_onset_l[i-2]) & (ton_onset<inspir_onset_l[i+2])]
        
        lick_in=lick[(lick > inspir_onset_l[i]) & (lick<inspir_onset_l[i+1])]
        
        if lick_in.size != 0:  # only include trials with a lick in between the breaths
            lick_bef=lick[(lick < inspir_onset_l[i])] # find the timing of lick before inspiration onset
            if len(lick_bef)>0:
                lick_bef=lick_bef[-1]
                licks.append(lick - inspir_onset_l[i])
                lick2_time.append(lick_in[-1] - inspir_onset_l[i])
                lick_bef_time.append(lick_bef - inspir_onset_l[i])
                n_licks.append(lick_in.size)                
                ibi.append(inspir_onset_l[i+1] - inspir_onset_l[i])

    ibi=np.array(ibi)
    n_licks=np.array(n_licks)
    lick_bef_time=np.array(lick_bef_time)
    lick2_time=np.array(lick2_time)

    licks[:] = [ele for i, ele in enumerate(licks) if ibi[i]<max_ibi]
    n_licks=n_licks[ibi<max_ibi]
    lick_bef_time=lick_bef_time[ibi<max_ibi]
    lick2_time=lick2_time[ibi<max_ibi]
    ibi=ibi[ibi<max_ibi]
    
    idx_all=np.arange(0,len(ibi))
    lick1_rem=np.where((n_licks==1) & (lick_bef_time>-0.05) & (lick_bef_time<0))
    lick2_rem=np.where((n_licks==2) & (lick2_time>(ibi-0.05)) & (lick2_time<ibi))
    idx_keep=np.setdiff1d(idx_all,np.concatenate((lick1_rem[0],lick2_rem[0])))

    n_licks=n_licks[idx_keep]
    lick_bef_time=lick_bef_time[idx_keep]
    ibi=ibi[idx_keep]
            
    sorted_indexes=np.argsort(ibi)
    sorted_indexes=sorted_indexes[::-1]

    licks = [licks[i] for i in sorted_indexes]
    n_licks=n_licks[sorted_indexes]
    lick_bef_time=lick_bef_time[sorted_indexes]
    ibi=ibi[sorted_indexes]

    d_bound=np.median(ibi[n_licks==2])

    psth_1_i= np.where((ibi<d_bound) & (n_licks==2))[0]
    psth_2_i= np.where((ibi>d_bound) & (n_licks==2))[0]
    if len(psth_1_i)>n_trial:
        psth_1_i=psth_1_i[:n_trial]
    if len(psth_2_i)>n_trial:
        psth_2_i=psth_2_i[-n_trial:]

    if (len(psth_1_i)<min_trial) | (len(psth_2_i)<min_trial):
        print(f'Less than {min_trial} trials')
        return np.array([]), np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]) 

    PSTH_LICK_1, PSTH_LICK_2, PSTH_LICK_1_SEM, PSTH_LICK_2_SEM, PSTH_BINS, PEAKS_LICK_1, PEAKS_LICK_2 = [], [], [], [], [], [], []
    num_units = mat[foi]['spike_times'].shape[0]
    for u in range(num_units):          
        spikes = [] # where the spikes occur
        ibi = []
        lick2_time = []
        n_licks = [] # number of licks in btw breaths
        lick_bef_time=[]
        all_spikes = mat[foi]['spike_times'][u].copy()
        good_spikes = all_spikes # get good spikes
        for i, d in enumerate(good_spikes):
            d = np.atleast_1d(d)
            good_spikes[i] = d[d < traces_len*3.4/1000]+traces_len*3.4/1000*i
        good_spikes = np.hstack(good_spikes)
        
        for i,_ in enumerate(inspir_onset_l[2:-2],2):
            
            lick=ton_onset[(ton_onset > inspir_onset_l[i-2]) & (ton_onset<inspir_onset_l[i+2])]
            
            lick_in=lick[(lick > inspir_onset_l[i]) & (lick<inspir_onset_l[i+1])]
            
            if lick_in.size != 0:  # only include trials with a lick in between the breaths
                lick_bef=lick[(lick < inspir_onset_l[i])] # find the timing of lick before inspiration onset
                if len(lick_bef)>0:
                    lick_bef=lick_bef[-1]
                    lick2_time.append(lick_in[-1] - inspir_onset_l[i])
                    lick_bef_time.append(lick_bef - inspir_onset_l[i])
                    n_licks.append(lick_in.size)
                    spike_breath=good_spikes-inspir_onset_l[i]
                    spike_breath=spike_breath[spike_breath>psth_s]
                    spike_breath=spike_breath[spike_breath<psth_e]
                    spikes.append(spike_breath)
                    ibi.append(inspir_onset_l[i+1] - inspir_onset_l[i])

        lick2_time=np.array(lick2_time)    
        ibi=np.array(ibi)
        lick_bef_time=np.array(lick_bef_time)
        n_licks=np.array(n_licks)
                
        spikes[:] = [ele for i, ele in enumerate(spikes) if ibi[i]<max_ibi]
        lick2_time=lick2_time[ibi<max_ibi]
        n_licks=n_licks[ibi<max_ibi]
        lick_bef_time=lick_bef_time[ibi<max_ibi]
        ibi=ibi[ibi<max_ibi]
        
        idx_all=np.arange(0,len(ibi))
        lick1_rem=np.where((n_licks==1) & (lick_bef_time>-0.05) & (lick_bef_time<0))
        lick2_rem=np.where((n_licks==2) & (lick2_time>(ibi-0.05)) & (lick2_time<ibi))
        idx_keep=np.setdiff1d(idx_all,np.concatenate((lick1_rem[0],lick2_rem[0])))
    
        spikes = [spikes[i] for i in idx_keep]
        n_licks=n_licks[idx_keep]
        lick_bef_time=lick_bef_time[idx_keep]
        ibi=ibi[idx_keep]
        
        sorted_indexes=np.argsort(ibi)
        sorted_indexes=sorted_indexes[::-1]
        
        n_licks=n_licks[sorted_indexes]
        lick_bef_time=lick_bef_time[sorted_indexes]
        ibi=ibi[sorted_indexes]
        spikes=[spikes[i] for i in sorted_indexes]
        
        psth_lick_1=[spikes[i] for i in psth_1_i]   
        psth_lick_1=[np.histogram(i_p,psth_bin)[0] for i_p in psth_lick_1]
        psth_lick_1=np.array(psth_lick_1)
        psth_lick_2=[spikes[i] for i in psth_2_i]
        psth_lick_2=[np.histogram(i_p,psth_bin)[0] for i_p in psth_lick_2]
        psth_lick_2=np.array(psth_lick_2)
        half_bin=(psth_bin[1]-psth_bin[0])/2
        psth_bins=psth_bin[1:]-half_bin
        psth_lick_2=psth_lick_2/(half_bin*2)
        psth_lick_1=psth_lick_1/(half_bin*2)
        
        psth_lick_1_sem=stats.sem(psth_lick_1)
        psth_lick_2_sem=stats.sem(psth_lick_2)
        
        psth_lick_1=np.mean(psth_lick_1,axis=0)
        psth_lick_2=np.mean(psth_lick_2,axis=0)
        
        min_height=10
        psth_1_thr=np.mean(psth_lick_1)+np.std(psth_lick_1)
        if psth_1_thr<min_height:
            psth_1_thr=min_height
        peaks_1=signal.find_peaks(psth_lick_1, height=psth_1_thr, distance=0.14/(half_bin*2))[0]
        peaks_lick_1=psth_bins[peaks_1]
        peaks_lick_1=peaks_lick_1[(peaks_lick_1>0.1) & (peaks_lick_1<0.3)]
        psth_2_thr=np.mean(psth_lick_2)+np.std(psth_lick_2)
        if psth_2_thr<min_height:
            psth_2_thr=min_height
        peaks_2=signal.find_peaks(psth_lick_2, height=psth_2_thr, distance=0.14/(half_bin*2))[0]
        peaks_lick_2=psth_bins[peaks_2]
        peaks_lick_2=peaks_lick_2[(peaks_lick_2>0.16) & (peaks_lick_2<0.36)]
        
        
        PSTH_LICK_1.append(psth_lick_1)
        PSTH_LICK_2.append(psth_lick_2)
        PSTH_LICK_1_SEM.append(psth_lick_1_sem)
        PSTH_LICK_2_SEM.append(psth_lick_2_sem)
        PSTH_BINS.append(psth_bins)
        PEAKS_LICK_1.append(peaks_lick_1)
        PEAKS_LICK_2.append(peaks_lick_2)
    
    return np.array(PSTH_LICK_1, dtype=object), np.array(PSTH_LICK_2, dtype=object), np.array(PSTH_LICK_1_SEM, dtype=object),  np.array(PSTH_LICK_2_SEM, dtype=object), np.array(PSTH_BINS, dtype=object),  np.array(PEAKS_LICK_1, dtype=object), np.array(PEAKS_LICK_2, dtype=object) 

# def detect_swallows(session_key, trial_key, breathing, phase, phase_criticals=[1.4, 2.3], min_vel_duration=0.1):
def detect_swallows(mat, phase_criticals=[1.2, 1.7], min_vel_duration=0.05):
    """Detect swallow events
    
    Args:
        session_key (dict or key):    session of interest
        trial_key (iterable):         trials to use
        breathing (ndarray):          breathing traces (trials x 1471)
        phase (ndarray):              breathing phase (trials x 1471)
        phase_criticals (list):       upper and lower bounds of phase limits for detection
        min_vel_duration (float):     minimum time signal velocity should be within bounds (in seconds)
        
    Returns (lists):
        all_peaks:                       list of all peak indices (len = trials)
        phase_restricted_peaks:          peak indices retricted to phase criticals
        phase_lick_restricted_peaks:     peak indices restricted to phase criticals and licks
        phase_lick_vel_restricted_peaks: peak indices restricted to phase criticals, licks and velocity
        peak_offs:                       peak offset indices for velocity restricted peaks
    """
    fs = mat[foi]['video_fs']
    tr_len = mat[foi]['trial_len']
    
    tongue_on = mat[foi]['tongue_onset']  
    breathing, _, phase = get_breathing(mat)
   
    a, b = phase_criticals
    all_peaks, phase_restricted_peaks, phase_lick_restricted_peaks, phase_lick_vel_restricted_peaks, peak_offs = [], [], [], [], []
    for i in range(breathing.shape[0]):
        breath = np.convolve(breathing[i], np.ones(5)/5, mode='same')
        _peaks, _ = scipy.signal.find_peaks(breath, distance=10)
        all_peaks.append(_peaks)
        
        #restrict by phase
        peak_phase = phase[i][_peaks]
        _peaks = _peaks[(peak_phase>a) & (peak_phase<b)] 
        
        phase_restricted_peaks.append(_peaks)

        #restrict by licks
        tr_start_time = tr_len*i
        ix = np.where((tongue_on>=tr_start_time) & (tongue_on<=tr_start_time+tr_len))[0]
        _tongue_on = tongue_on[ix]

        ix = []
        if _tongue_on.size == 0:
            phase_lick_restricted_peaks.append(np.array([]))
        else:
            _tongue_on = _tongue_on - tr_start_time
            for _peak in _peaks/fs:
                _bef = np.where((_tongue_on <= _peak) & (_tongue_on >= _peak-1.5))[0]
                _after = np.where((_tongue_on >= _peak) & (_tongue_on <= _peak+0.25))[0]
                
                if (_bef.size >= 5) & (_after.size >= 1): #at least 5 licks before and 1 lick after the event
                    if (_peak - _tongue_on[_bef[-1]] < 0.5): #the last lick before event should be within 500ms before event.
                        ilf = 1/np.diff(_tongue_on[_bef])[-4:]
                        if ((ilf > 3.5) & (ilf < 13)).all():
                            ix.append(int(_peak*fs))
            _peaks = np.array(ix)
            phase_lick_restricted_peaks.append(_peaks)
        
        #restrict by velocity
        velocity = np.diff(breath)
        results = []
        count = 0
        start_index = None
        for j in range(len(velocity)):
            if (velocity[j] >= -0.07) & (velocity[j] <= 0.07):
                count += 1
                if start_index is None:
                    start_index = j
            else:
                if count >= 10:
                    results.append((start_index, start_index+count, count/fs))
                count = 0
                start_index = None
        
        # Check if the last elements in the list are zeros
        if count >= 2:
            results.append((start_index, start_index+count, count/fs))
        results = np.array(results)
        results = results[results[:,2] >= min_vel_duration]
        results[:,0] = results[:,0]+1
        
        val, val_off = [], []
        for peak_on, peak_off in results[:, :2]:
            ix = np.where((phase_lick_restricted_peaks[i] >= peak_on) & (phase_lick_restricted_peaks[i] <= peak_off))[0]
            if ix.size:
                val.append(int(peak_on))
                val_off.append(int(peak_off))
                
        val = np.hstack(val) if len(val) else np.array([])
        val_off = np.hstack(val_off) if len(val_off) else np.array([])
        phase_lick_vel_restricted_peaks.append(val)
        peak_offs.append(val_off)
    return all_peaks, phase_restricted_peaks, phase_lick_restricted_peaks, phase_lick_vel_restricted_peaks, peak_offs

def get_breathing(mat):
    """Fetch breathing traces (1471 samples from trial start)
    
    Args:
        session_key (dict or key):    session of interest
        trial_key (iterable):         trials to fetch
    
    Returns breathing traces, breathing amplitude, breathing phase (dim: trials x 1471)
    """
    breathing = mat[foi]['breathing'].copy()
    breathing_ts = mat[foi]['breathing_ts'].copy()
    num_trials = breathing.size
    good_breathing = breathing.copy()
    for i, d in enumerate(breathing):
        good_breathing[i] = d[breathing_ts[i] < 1471*3.4/1000]
    good_breathing = stats.zscore(np.vstack(good_breathing),axis=None)
    
    good_breathing = np.hstack(good_breathing)
    # -- moving-average
    window_size = int(round(0.0034/(breathing_ts[0][1]-breathing_ts[0][0]),0))  #sample
    kernel = np.ones(window_size) / window_size
    good_breathing = np.convolve(good_breathing, kernel, 'same')
    # -- down-sample
    good_breathing = good_breathing[window_size::window_size] * -1
    result = np.zeros(good_breathing.shape[0]+1)
    result[:-1] = good_breathing
    result[-1] = good_breathing[-1]
    hilbert = scipy.signal.hilbert(result)
    
    return result.reshape(num_trials,-1), np.abs(hilbert).reshape(num_trials,-1), np.angle(hilbert).reshape(num_trials,-1)

def get_licks_around_event(mat):
    """ """
    events = mat[foi]['swallow_times']
    result, indices = [], []
    fs = 1000/3.4
    tr_len = 1471/fs
    tongue_on = mat[foi]['tongue_onset']
    for i, event in enumerate(events):
        if event.size > 0:
            indices.append(i)
            for el in event:
                start = i*tr_len
                event_time = el/fs
                stop = start + tr_len
                ton_on = tongue_on[(tongue_on>=start) & (tongue_on<=stop)] - event_time - start
                result.append(ton_on)
    return np.array(result, dtype=object)

def get_breaths_around_event(mat):
    """ """
    events = mat[foi]['swallow_times']
    result, indices = [], []
    fs = 1000/3.4
    tr_len = 1471/fs
    inspir_on = mat[foi]['inspir_onset']
    for i, event in enumerate(events):
        if event.size > 0:
            indices.append(i)
            for el in event:
                start = i*tr_len
                event_time = el/fs
                stop = start + tr_len
                i_on = inspir_on[(inspir_on>=start) & (inspir_on<=stop)] - event_time - start
                result.append(i_on)
    return np.array(result, dtype=object)

def compute_psth(aligned_spikes, bin_size=0.05, window=[-1, 1]):
    """Compute PSTH 
    
    Args:
        aligned_spikes (iterable):    Event aligned spike times. dims: trial x spike times
        bin_size (float):             PSTH binning size in seconds
        window (list):                PSTH start and stop in seconds
    
    Returns time bins, PSTH (trial x bins)
    """
    # Define the edges of the bins (1ms bins; then later convolve with 50ms window to smoothen)
    bins = np.arange(window[0], window[1] + 0.001, 0.001)
    # Initialize the histogram
    all_histograms = []

    # Compute the histogram for each trial
    for trial in aligned_spikes:
        # Compute histogram for the current trial
        hist, _ = np.histogram(trial, bins)
        all_histograms.append(hist)
    all_histograms = np.vstack(all_histograms)
    kernel = np.ones(int(bin_size*1000))/bin_size
    for i, row in enumerate(all_histograms):
        row = np.convolve(row, kernel, mode='same')
        all_histograms[i] = row

    return bins[:-1], all_histograms

def get_swallow_aligned_spikes(mat):
    """ """
    rasters, psths = [], []
    num_frame = 1471
    fs = 1000/3.4
    num_units = mat[foi]['spike_times'].shape[0]
    for u in range(num_units):  
        all_spikes = mat[foi]['spike_times'][u].copy()
        # all_spikes = (ephys.Unit.TrialSpikes & [{'trial': tr} for tr in trial_key] & unit_key).fetch('spike_times', order_by='trial')
        good_spikes = np.array(all_spikes*fs)
        good_spikes = [np.array(d).astype(int) for d in good_spikes]
        
        for i, d in enumerate(good_spikes):
            good_spikes[i] = d[d < num_frame]
        
        # get event aligned firings
        event_aligned_spikes = []
        psth = []
        for i, el in enumerate(phase_lick_vel_restricted_peaks):
            if el.size > 0:
                for event in el:
                    spikes = good_spikes[i]
                    spikes = spikes - event
                    event_aligned_spikes.append(spikes/fs)
        
        rasters.append(event_aligned_spikes)
        if len(event_aligned_spikes)>0:
            bins, psth = compute_psth(event_aligned_spikes, bin_size=0.05, window=[-1, 1])
        else:
            bins = np.arange(-1, 1 + 0.001, 0.001)[:-1]
            psth = np.array([])
        psths.append(psth)
    return np.array(psths), np.array(rasters), bins

def get_swallow_tuning(mat, n_bootstrap=1000):
    """    
    Parameters
    ----------
    rasters : list of list of np.array
        rasters[i] = list of trials for neuron i,
        each trial is a 1D array of spike times (relative to event).
    n_bootstrap : int, optional
        Number of bootstrap iterations (default=1000).
    
    Returns
    -------
    modulation_indexNEW : np.ndarray
        Shape (n_neurons, 2). Each row has [MI_pre, MI_post].
    pvalsNEW : np.ndarray
        Shape (n_neurons, 2). Bootstrap p-values.
    """
    rasters = mat[foi]['swallow_raster']
    n_neurons = len(rasters)
    modulation_indexNEW = np.full((n_neurons, 2), np.nan)
    pvalsNEW = np.ones((n_neurons, 2))
    count = [len(raster) for raster in rasters]
    for i, neuron_rasters in enumerate(rasters):
        # if isinstance(neuron_rasters, list) and len(neuron_rasters) >= 20:
        if len(neuron_rasters) >= 20:
            # Compute firing rates for each trial
            r = []
            for trial in neuron_rasters:
                r1 = np.sum((trial > -0.3) & (trial < -0.1)) / 0.2
                r2 = np.sum((trial > -0.05) & (trial < 0.05)) / 0.1
                r3 = np.sum((trial > 0.1) & (trial < 0.3)) / 0.2
                r.append([r1, r2, r3])
            r = np.array(r)

            # Bootstrap mean firing rates
            r_btstrp = np.zeros((n_bootstrap, 3))
            n_trials = r.shape[0]
            for b in range(n_bootstrap):
                idx = np.random.randint(0, n_trials, n_trials)
                r_btstrp[b, :] = r[idx, :].mean(axis=0)

            # p-values (bootstrap test)
            mean_diff_21 = r[:, 1].mean() - r[:, 0].mean()
            if mean_diff_21 > 0:
                pvalsNEW[i, 0] = np.mean((r_btstrp[:, 1] - r_btstrp[:, 0]) < 0)
            elif mean_diff_21 < 0:
                pvalsNEW[i, 0] = np.mean((r_btstrp[:, 1] - r_btstrp[:, 0]) > 0)
            else:
                pvalsNEW[i, 0] = 1.0

            mean_diff_23 = r[:, 1].mean() - r[:, 2].mean()
            if mean_diff_23 > 0:
                pvalsNEW[i, 1] = np.mean((r_btstrp[:, 1] - r_btstrp[:, 2]) < 0)
            elif mean_diff_23 < 0:
                pvalsNEW[i, 1] = np.mean((r_btstrp[:, 1] - r_btstrp[:, 2]) > 0)
            else:
                pvalsNEW[i, 1] = 1.0

            # Modulation indices
            mi1 = (r[:, 1].mean() - r[:, 0].mean()) / max(r[:, 1].mean(), r[:, 0].mean(), 1e-12)
            mi2 = (r[:, 1].mean() - r[:, 2].mean()) / max(r[:, 1].mean(), r[:, 2].mean(), 1e-12)
            modulation_indexNEW[i, :] = [mi1, mi2]

    return modulation_indexNEW, pvalsNEW, np.array(count)

def preprocess(mat):
    """ """
    unit_keys = []
    _insertion_num = mat['processed']['insertion_key']['insertion_number']
    _session = mat['processed']['insertion_key']['session']
    _subject_id = mat['processed']['insertion_key']['subject_id']
    _clust_method = 'kilosort2'
    for i in range(mat['neuron_unit_info'].shape[0]):
        _res = {"subject_id": _subject_id, "session": _session, "insertion_number": _insertion_num, "clustering_method": _clust_method, "unit": mat['neuron_unit_info'][i,0]}
        unit_keys.append(_res)
    return np.array(unit_keys, dtype=object)
    
    

#%% driver
############## USER INPUT
foi = 'processed' #change later after veryfing 
input_dir = r'D:/mat/'
# fname = r'map-export_DL004_20210308_180033_s1_p1_nwb.mat'
##############

loop_start = time.perf_counter()

mat_files = glob.glob(os.path.join(input_dir, "*nwb.mat"))
filenames = [os.path.basename(f) for f in mat_files]
for _c, fname in enumerate(filenames):
    start_time = time.perf_counter()
    print("\n\nWorking on {} ({} of {})" .format(fname, _c+1, len(filenames)))
    
    save_name = fname[:-4] + "_processed.mat"
    fpath = input_dir + fname
    save_fpath = input_dir + save_name
    
    if os.path.isfile(save_fpath):
        print("\nFile exists. Skipping...")
    
    else:
        print("\nLoading data...")
        mat = loadmat(fpath, simplify_cells=True)
        if foi not in mat.keys():
            mat[foi] = {}
        
        # pre-process
        print("\nPreprocessing...")
        unit_keys = preprocess(mat)
        mat['unit_keys'] = unit_keys
        
        # good trials
        print("\nFetching good trials...")
        trials, bad_trial_side, bad_trial_bot = get_good_trials(mat)
        mat[foi]['trials'] = trials
        mat[foi]['bad_trial_side'] = bad_trial_side
        mat[foi]['bad_trial_bot'] = bad_trial_bot
        
        # good units
        print("\nFetching spike times...")
        good_unit_ix = get_good_units(mat)
        mat[foi]['ccf_x'] = mat['histology']['ccf_x'][good_unit_ix]
        mat[foi]['ccf_y'] = mat['histology']['ccf_y'][good_unit_ix]
        mat[foi]['ccf_z'] = mat['histology']['ccf_z'][good_unit_ix]
        spike_times = mat['neuron_single_units'][good_unit_ix]
        spike_times = np.array([s[trials-1] for s in spike_times], dtype=object)
        mat[foi]['spike_times'] = spike_times
    
        mat[foi]['movement_timing_unit_keys'] = [mat['unit_keys'][_i] for _i in good_unit_ix]
        mat[foi]['swallow_unit_keys'] = mat[foi]['movement_timing_unit_keys'].copy()
        
        # breathing
        print("\nFetching breathing...")
        mat[foi]['breathing'] = np.array([mat['breathing'][i-1] for i in trials], dtype=object)
        mat[foi]['breathing_ts'] = np.array([mat['breathing_ts'][i-1] for i in trials], dtype=object)
        
        # tracking
        print("\nFetching tracking...")
        mat[foi]['tracking'] = {}
        mat[foi]['tracking']['camera_3_side'] = {}
        side_trials = mat['tracking']['camera_3_side']['trialNum']
        side_trials_idx = np.array([np.where(side_trials == t)[0][0] for t in trials])
        mat[foi]['tracking']['camera_3_side']['jaw_likelihood'] = mat['tracking']['camera_3_side']['jaw_likelihood'][side_trials_idx]
        mat[foi]['tracking']['camera_3_side']['jaw_x'] = mat['tracking']['camera_3_side']['jaw_x'][side_trials_idx]
        mat[foi]['tracking']['camera_3_side']['jaw_y'] = mat['tracking']['camera_3_side']['jaw_y'][side_trials_idx]
        mat[foi]['tracking']['camera_3_side']['lickport_likelihood'] = mat['tracking']['camera_3_side']['lickport_likelihood'][side_trials_idx]
        mat[foi]['tracking']['camera_3_side']['lickport_x'] = mat['tracking']['camera_3_side']['lickport_x'][side_trials_idx]
        mat[foi]['tracking']['camera_3_side']['lickport_y'] = mat['tracking']['camera_3_side']['lickport_y'][side_trials_idx]
        mat[foi]['tracking']['camera_3_side']['tongue_likelihood'] = mat['tracking']['camera_3_side']['tongue_likelihood'][side_trials_idx]
        mat[foi]['tracking']['camera_3_side']['tongue_x'] = mat['tracking']['camera_3_side']['tongue_x'][side_trials_idx]
        mat[foi]['tracking']['camera_3_side']['tongue_y'] = mat['tracking']['camera_3_side']['tongue_y'][side_trials_idx]
        
        mat[foi]['tracking']['camera_4_bottom'] = {}
        bottom_trials = mat['tracking']['camera_4_bottom']['trialNum']
        bottom_trials_idx = np.array([np.where(bottom_trials == t)[0][0] for t in trials])
        mat[foi]['tracking']['camera_4_bottom']['jaw_likelihood'] = mat['tracking']['camera_4_bottom']['jaw_likelihood'][bottom_trials_idx]
        mat[foi]['tracking']['camera_4_bottom']['jaw_x'] = mat['tracking']['camera_4_bottom']['jaw_x'][bottom_trials_idx]
        mat[foi]['tracking']['camera_4_bottom']['jaw_y'] = mat['tracking']['camera_4_bottom']['jaw_y'][bottom_trials_idx]
        mat[foi]['tracking']['camera_4_bottom']['lickport_likelihood'] = mat['tracking']['camera_4_bottom']['lickport_likelihood'][bottom_trials_idx]
        mat[foi]['tracking']['camera_4_bottom']['lickport_x'] = mat['tracking']['camera_4_bottom']['lickport_x'][bottom_trials_idx]
        mat[foi]['tracking']['camera_4_bottom']['lickport_y'] = mat['tracking']['camera_4_bottom']['lickport_y'][bottom_trials_idx]
        mat[foi]['tracking']['camera_4_bottom']['tongue_likelihood'] = mat['tracking']['camera_4_bottom']['tongue_likelihood'][bottom_trials_idx]
        mat[foi]['tracking']['camera_4_bottom']['tongue_x'] = mat['tracking']['camera_4_bottom']['tongue_x'][bottom_trials_idx]
        mat[foi]['tracking']['camera_4_bottom']['tongue_y'] = mat['tracking']['camera_4_bottom']['tongue_y'][bottom_trials_idx]
        
        # movement timing
        print("\nFetching movement timing...")
        tongue_onset, tongue_offset, lick_bout_onset, lick_bout_offset, inspir_onset = get_movement_timing(mat)
        mat[foi]['tongue_onset'] = tongue_onset
        mat[foi]['tongue_offset'] = tongue_offset
        mat[foi]['lick_bout_onset'] = lick_bout_onset
        mat[foi]['lick_bout_offset'] = lick_bout_offset
        mat[foi]['inspir_onset'] = inspir_onset
        
        # unit jaw and breath tuning
        print("\nFetching jaw tuning...")
        jaw_tuning, jaw_preferred_phase, jaw_x, jaw_y, jaw_l, jaw_h, jaw_kuiper_test, jaw_di_perm = get_unit_jaw_tuning(mat)
        mat[foi]['jaw_tuning'] = jaw_tuning
        mat[foi]['jaw_preferred_phase'] = jaw_preferred_phase
        mat[foi]['jaw_x'] = jaw_x
        mat[foi]['jaw_y'] = jaw_y
        mat[foi]['jaw_l'] = jaw_l
        mat[foi]['jaw_h'] = jaw_h
        mat[foi]['jaw_kuiper_test'] = jaw_kuiper_test
        mat[foi]['jaw_di_perm'] = jaw_di_perm
        
        print("\nFetching breath tuning...")
        breath_tuning, breath_preferred_phase, breath_x, breath_y, breath_l, breath_h, breath_mi_perm = get_unit_breath_tuning(mat)
        mat[foi]['breath_tuning'] = breath_tuning
        mat[foi]['breath_preferred_phase'] = breath_preferred_phase
        mat[foi]['breath_x'] = breath_x
        mat[foi]['breath_y'] = breath_y
        mat[foi]['breath_l'] = breath_l
        mat[foi]['breath_h'] = breath_h
        mat[foi]['breath_mi_perm'] = breath_mi_perm
        
        # peak difference breath and lick (behavior and unit)
        # breath
        print("\nFetching breath peak difference...")
        mat[foi]['peak_diff_breath'] = {}
        mat[foi]['peak_diff_breath']['behavior'] = {}
        mat[foi]['peak_diff_breath']['unit'] = {}
        psth_breath_1, psth_breath_2, psth_breath_1_sem, psth_breath_2_sem, psth_bins, peaks_breath_1, peaks_breath_2 = get_peak_diff_breath_behv(mat)
        mat[foi]['peak_diff_breath']['behavior']['psth_breath_1'] = psth_breath_1
        mat[foi]['peak_diff_breath']['behavior']['psth_breath_2'] = psth_breath_2
        mat[foi]['peak_diff_breath']['behavior']['sem_breath_1'] = psth_breath_1_sem
        mat[foi]['peak_diff_breath']['behavior']['sem_breath_2'] = psth_breath_2_sem
        mat[foi]['peak_diff_breath']['behavior']['psth_bins'] = psth_bins
        mat[foi]['peak_diff_breath']['behavior']['peaks_breath_1'] = peaks_breath_1
        mat[foi]['peak_diff_breath']['behavior']['peaks_breath_2'] = peaks_breath_2
        
        psth_breath_1, psth_breath_2, psth_breath_1_sem, psth_breath_2_sem, psth_bins, peaks_breath_1, peaks_breath_2 = get_peak_diff_breath_unit(mat)
        mat[foi]['peak_diff_breath']['unit']['psth_breath_1'] = psth_breath_1
        mat[foi]['peak_diff_breath']['unit']['psth_breath_2'] = psth_breath_2
        mat[foi]['peak_diff_breath']['unit']['sem_breath_1'] = psth_breath_1_sem
        mat[foi]['peak_diff_breath']['unit']['sem_breath_2'] = psth_breath_2_sem
        mat[foi]['peak_diff_breath']['unit']['psth_bins'] = psth_bins
        mat[foi]['peak_diff_breath']['unit']['peaks_breath_1'] = peaks_breath_1
        mat[foi]['peak_diff_breath']['unit']['peaks_breath_2'] = peaks_breath_2
        
        # lick
        print("\nFetching lick peak difference...")
        mat[foi]['peak_diff_lick'] = {}
        mat[foi]['peak_diff_lick']['behavior'] = {}
        mat[foi]['peak_diff_lick']['unit'] = {}
        psth_lick_1, psth_lick_2, psth_lick_1_sem, psth_lick_2_sem, psth_bins, peaks_lick_1, peaks_lick_2 = get_peak_diff_lick_behv(mat)
        mat[foi]['peak_diff_lick']['behavior']['psth_lick_1'] = psth_lick_1
        mat[foi]['peak_diff_lick']['behavior']['psth_lick_2'] = psth_lick_2
        mat[foi]['peak_diff_lick']['behavior']['sem_lick_1'] = psth_lick_1_sem
        mat[foi]['peak_diff_lick']['behavior']['sem_lick_2'] = psth_lick_2_sem
        mat[foi]['peak_diff_lick']['behavior']['psth_bins'] = psth_bins
        mat[foi]['peak_diff_lick']['behavior']['peaks_lick_1'] = peaks_lick_1
        mat[foi]['peak_diff_lick']['behavior']['peaks_lick_2'] = peaks_lick_2
        
        psth_lick_1, psth_lick_2, psth_lick_1_sem, psth_lick_2_sem, psth_bins, peaks_lick_1, peaks_lick_2 = get_peak_diff_lick_unit(mat)
        mat[foi]['peak_diff_lick']['unit']['psth_lick_1'] = psth_lick_1
        mat[foi]['peak_diff_lick']['unit']['psth_lick_2'] = psth_lick_2
        mat[foi]['peak_diff_lick']['unit']['sem_lick_1'] = psth_lick_1_sem
        mat[foi]['peak_diff_lick']['unit']['sem_lick_2'] = psth_lick_2_sem
        mat[foi]['peak_diff_lick']['unit']['psth_bins'] = psth_bins
        mat[foi]['peak_diff_lick']['unit']['peaks_lick_1'] = peaks_lick_1
        mat[foi]['peak_diff_lick']['unit']['peaks_lick_2'] = peaks_lick_2
        
        # swallow times
        print("\nFetching swallow times...")
        all_peaks, phase_restricted_peaks, phase_lick_restricted_peaks, phase_lick_vel_restricted_peaks, peak_offs = detect_swallows(mat)
        mat[foi]['swallow_times'] = np.array(phase_lick_vel_restricted_peaks, dtype=object)
        mat[foi]['swallow_aligned_licks'] = get_licks_around_event(mat)
        mat[foi]['swallow_aligned_breaths'] = get_breaths_around_event(mat)
        
        # swallow aligned psth and raster
        print("\nFetching swallow psth & raster...")
        psth, raster, psth_bins = get_swallow_aligned_spikes(mat)
        mat[foi]['swallow_psth_bins'] = psth_bins
        mat[foi]['swallow_raster'] = raster
        mat[foi]['swallow_psth'] = psth
        
        # swallow tuning
        print("\nFetching swallow tuning...")
        swallow_tuning, unit_swallow_pval, unit_swallow_count = get_swallow_tuning(mat)
        mat[foi]['swallow_tuning'] = swallow_tuning
        mat[foi]['unit_swallow_pval'] = unit_swallow_pval
        mat[foi]['unit_swallow_count'] = unit_swallow_count
        
        #save
        print("\nsaving...")
        savemat(save_fpath, mat, do_compression=True)
        
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        print("\nElapsed time: {:.2f} minutes" .format(elapsed_time/60))
        
loop_end = time.perf_counter()
loop_elapsed_time = end_time - start_time
print("\n\nTotal Elapsed time: {:.2f} minutes" .format(loop_elapsed_time/60))