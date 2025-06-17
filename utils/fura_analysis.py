# --- Step 2: Analyze a single trace with multiple peaks ---
def analyze_single_trace(time, trace, stim_regions=None, plot=True):
    """ Analyze a single calcium trace to extract peak characteristics and other metrics."
    
    Parameters:
    time (array-like): Time points corresponding to the calcium trace.
    trace (array-like): Calcium trace data.
    stim_regions (list of tuples): Optional list of stimulus regions as (start, end) tuples.
    plot (bool): Whether to plot the trace and detected peaks.      
    Returns:
    pd.DataFrame: DataFrame containing the analysis results for each detected peak.
    
    This function analyzes a single calcium trace to extract peak characteristics including:
    - Peak index and time
    - Peak amplitude
    - Delta Ca²⁺ (change from baseline)
    - Time to peak
    - Area Under the Curve (AUC) during stimulus region
    - Decay time constant (tau) using exponential decay fitting
    - Signal-to-Noise Ratio (SNR)
    - Duration of the peak
    - Maximum upstroke and downstroke slopes
    - Repolarisation slope
    It also plots the trace with detected peaks and stimulus regions if `plot` is set to True.
    The function uses scipy's `find_peaks` to detect peaks, `curve_fit` for exponential decay fitting, and `simpson` for AUC calculation.
    The results are returned as a pandas DataFrame with columns for each metric.
    The function also handles baseline estimation and noise calculation from the initial part of the trace.
    If `stim_regions` is provided, it uses these regions to calculate AUC; otherwise, it defaults to a fixed region around each peak.
    If `plot` is True, it generates a plot showing the calcium trace, detected peaks, and stimulus regions.
    The function is designed to handle multiple peaks in the trace and returns a DataFrame with results for each peak detected.
    The analysis includes:
    - Baseline estimation from the first 10% of the trace
    - Noise estimation as the standard deviation of the baseline
    - Peak detection using a threshold based on baseline and noise
    - Exponential decay fitting to estimate decay time constant (tau)
    - Calculation of AUC using Simpson's rule
    - Calculation of signal-to-noise ratio (SNR)
    - Duration, maximum upstroke, downstroke slopes, and repolarisation slope calculations
    The function is robust to handle cases where no peaks are detected or fitting fails, returning NaN for those metrics.
    The results DataFrame contains the following columns:
    - 'Peak_Index': Index of the detected peak in the trace
    - 'Peak_Time': Time at which the peak occurs
    - 'Peak': Amplitude of the peak
    - 'Delta_Ca': Change in calcium concentration from baseline
    - 'Time_to_Peak': Time elapsed from the start to the peak
    - 'AUC': Area Under the Curve during the stimulus region
    - 'Decay_Tau': Decay time constant estimated from exponential fitting
    - 'SNR': Signal-to-Noise Ratio calculated as Delta Ca / noise standard deviation
    - 'Duration': Duration of the peak (time from start to end of the peak)
    - 'Max_Upstroke': Maximum upstroke slope before the peak
    - 'Max_Downstroke': Maximum downstroke slope after the peak
    - 'Repolarisation_Slope': Mean slope after the peak over a specified number of points
    The function is designed to be flexible and can handle different types of calcium imaging data, making it suitable for various experimental setups.
    It is particularly useful for analyzing Fura-2 AM calcium imaging data, but can be adapted for other types of calcium indicators as well.
    The function uses several libraries including pandas, numpy, matplotlib, scipy, and seaborn for data manipulation, numerical analysis, and plotting.
    The function is intended to be used in a Jupyter notebook or similar environment where data visualization is important for analysis.
    The function can be extended or modified to include additional metrics or analysis steps as needed for specific research questions or experimental designs. 

    """
    results = []
    baseline_end_idx = int(0.1 * len(trace))
    baseline = np.mean(trace[:baseline_end_idx])
    noise_std = np.std(trace[:baseline_end_idx])

    peaks, _ = find_peaks(trace, height=baseline + 3 * noise_std, distance=10)

    def exp_decay(t, A, tau, C):
        return A * np.exp(-t / tau) + C

    for i, peak_idx in enumerate(peaks):
        peak_time = time[peak_idx]
        peak = trace[peak_idx]
        delta_ca = peak - baseline/baseline
        time_to_peak = peak_time - time[0]

        if stim_regions is not None:
            region = next(((start, end) for start, end in stim_regions if start <= peak_time <= end), (time[0], time[-1]))
            stim_start, stim_end = region
        else:
            stim_start = time[max(0, peak_idx - 10)]
            #stim_end = time[min(len(time)-1, peak_idx + 30)]
            stim_end = time[min(len(time)-1, peak_idx)]

        stim_mask = (time >= stim_start) & (time <= stim_end)
        auc = simpson(trace[stim_mask] - baseline, dx=(time[1] - time[0]))

        decay_tau = np.nan
        if peak_idx < len(trace) - 5:
            decay_time = time[peak_idx:] - time[peak_idx]
            decay_trace = trace[peak_idx:]
            if np.any(np.diff(decay_trace) < 0):
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", OptimizeWarning)
                        popt, _ = curve_fit(
                            exp_decay,
                            decay_time,
                            decay_trace,
                            p0=(delta_ca, 1, baseline),
                            bounds=([0, 0.001, 0], [np.inf, 100, np.max(decay_trace)])
                        )
                        decay_tau = popt[1]
                except (RuntimeError, ValueError):
                    pass
    snr = delta_ca / noise_std if noise_std else np.nan

        duration = np.nan
        max_upstroke = np.nan
        max_downstroke = np.nan
        repolarisation_slope = np.nan

        if i == 0:
            window_size = 20
            start_idx = max(0, peak_idx - window_size)
            end_idx = min(len(trace) - 1, peak_idx + window_size)
            sub_time = time[start_idx:end_idx+1]
            sub_trace = trace[start_idx:end_idx+1]

            half_max = baseline + 0.5 * delta_ca
            above_half = np.where(sub_trace >= half_max)[0]
            if len(above_half) >= 2:
                duration = sub_time[above_half[-1]] - sub_time[above_half[0]]

            dt = np.diff(sub_time)
            dy = np.diff(sub_trace)
            slopes = dy / dt
            peak_pos = peak_idx - start_idx
            max_upstroke = np.max(slopes[:peak_pos]) if peak_pos > 0 else np.nan
            max_downstroke = np.min(slopes[peak_pos:]) if peak_pos < len(slopes) else np.nan

            # Estimate repolarisation slope (mean slope after the peak over 10 points)
            if peak_pos + 10 < len(slopes):
                repolarisation_slope = np.mean(slopes[peak_pos:peak_pos + 10])

        result = {
            'Peak_Index': peak_idx,
            'Peak_Time': peak_time,
            'Peak': peak,
            'Delta_Ca': delta_ca,
            'Time_to_Peak': time_to_peak,
            'AUC': auc,
            'Decay_Tau': decay_tau,
            'SNR': snr,
            'Duration': duration,
            'Max_Upstroke': max_upstroke,
            'Max_Downstroke': max_downstroke,
            'Repolarisation_Slope': repolarisation_slope
        }
        results.append(result)
    
            if plot:
            plt.figure(figsize=(10, 4))
            plt.plot(time, trace, label='Calcium Trace')
            plt.axhline(baseline, color='gray', linestyle='--', label='Baseline')
            plt.scatter(time[peak_idx], peak, color='red', label='Peak')
            plt.axvspan(stim_start, stim_end, color='yellow', alpha=0.2, label='Stimulus Region')
            plt.xlabel('Time (s)')
            plt.ylabel('Fura-2 Ratio')
            plt.title(f'Peak {i+1} Analysis')
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.show()

    return pd.DataFrame(results)

# --- Step 3: Analyze all cells ---
def analyze_all_cells(df, stim_regions=None, plot_each=False):
    """ Analyze all cells in the DataFrame to extract peak characteristics and other metrics."
    
    Parameters:
    df (pd.DataFrame): DataFrame containing the calcium imaging data with time and cell traces.
    stim_regions (list of tuples): Optional list of stimulus regions as (start, end) tuples.
    plot_each (bool): Whether to plot each cell's trace and detected peaks. Default is False.'
    Returns:
    pd.DataFrame: DataFrame containing the analysis results for each cell and its detected peaks.   
    This function iterates over each cell in the DataFrame, extracts the time and trace data, and calls the `analyze_single_trace` function to analyze each trace.
    It collects the results for each cell and returns a concatenated DataFrame with all results.
    The function assumes the first column of the DataFrame contains time data and the subsequent columns contain traces for different cells.    
    The `stim_regions` parameter allows specifying regions of interest for calculating the Area Under the Curve (AUC) during stimulus periods.
    If `plot_each` is set to True, it generates plots for each cell's trace with detected peaks and stimulus regions.'
    '    The function is designed to handle multiple cells in a single DataFrame, making it suitable for batch processing of calcium imaging data.
    The results DataFrame contains the following columns for each cell:
        
    - 'Cell': Name of the cell (column name in the original DataFrame)
    - 'Peak_Index': Index of the detected peak in the trace
    - 'Peak_Time': Time at which the peak occurs
    - 'Peak': Amplitude of the peak
    - 'Delta_Ca': Change in calcium concentration from baseline
    - 'Time_to_Peak': Time elapsed from the start to the peak
    - 'AUC': Area Under the Curve during the stimulus region            
    - 'Decay_Tau': Decay time constant estimated from exponential fitting
    - 'SNR': Signal-to-Noise Ratio calculated as Delta Ca / noise standard deviation
    - 'Duration': Duration of the peak (time from start to end of the peak)
    - 'Max_Upstroke': Maximum upstroke slope before the peak
    - 'Max_Downstroke': Maximum downstroke slope after the peak
    - 'Repolarisation_Slope': Mean slope after the peak over a specified number of points
    The function is intended to be used in a Jupyter notebook or similar environment where data visualization is important for analysis.
    The function can be extended or modified to include additional metrics or analysis steps as needed for specific research questions or experimental designs.
    The function uses several libraries including pandas, numpy, matplotlib, scipy, and seaborn for data manipulation, numerical analysis, and plotting.
    The function is designed to be flexible and can handle different types of calcium imaging data, making it suitable for various experimental setups. 
  """  
    time = df.iloc[:, 0].values
    all_results = []

    for cell in df.columns[1:]:
        trace = df[cell].values
        cell_results = analyze_single_trace(time, trace, stim_regions=stim_regions, plot=plot_each)
        cell_results['Cell'] = cell
        all_results.append(cell_results)

    return pd.concat(all_results, ignore_index=True)
