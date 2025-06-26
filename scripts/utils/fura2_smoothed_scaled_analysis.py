import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter
from scipy.optimize import curve_fit, OptimizeWarning
from scipy.integrate import simpson
import warnings
import os

'''
def exp_decay(t, A, tau, C):
    return A * np.exp(-t / tau) + C


def smooth_and_scale_trace(trace, window_length=7, polyorder=2):
    trace_smooth = savgol_filter(trace, window_length=window_length, polyorder=polyorder)
    trace_min = np.min(trace_smooth)
    trace_max = np.max(trace_smooth)
    if trace_max - trace_min == 0:
        trace_scaled = np.zeros_like(trace_smooth)
    else:
        trace_scaled = (trace_smooth - trace_min) / (trace_max - trace_min)
    return trace_smooth, trace_scaled


def plot_trace_scaling(time, raw_trace, scaled_trace, cell_id, output_dir=None):
    fig, ax1 = plt.subplots(figsize=(10, 4))
    ax1.plot(time, raw_trace, color='tab:blue', label='Raw Trace')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Raw Fura-2 Ratio', color='tab:blue')
    ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax2 = ax1.twinx()
    ax2.plot(time, scaled_trace, color='tab:red', linestyle='--', label='Smoothed + Scaled Trace')
    ax2.set_ylabel('Scaled Signal [0–1]', color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red')
    fig.suptitle(f'Cell {cell_id}: Raw vs Smoothed & Scaled')
    fig.tight_layout()
    fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.9))
    plt.grid(True)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(os.path.join(output_dir, f'cell_{cell_id}_scaling_plot.png'))
        plt.close()
    else:
        plt.show()


def analyze_trace(time, trace, cell_id, stim_regions=None, smoothing_window=7, smoothing_poly=2):
    trace_smooth, trace_scaled = smooth_and_scale_trace(trace, smoothing_window, smoothing_poly)
    baseline_end_idx = int(0.1 * len(trace_scaled))
    baseline = np.mean(trace_scaled[:baseline_end_idx])
    noise_std = np.std(trace_scaled[:baseline_end_idx])
    peaks, _ = find_peaks(trace_scaled, prominence = 1)#, height=baseline + 3 * noise_std, distance=10)

    results = []
    for peak_idx in peaks:
        peak_time = time[peak_idx]
        peak = trace_scaled[peak_idx]
        delta_ca = (peak - baseline)/baseline
        time_to_peak = peak_time - time[0]

        # Stim region
        if stim_regions is not None:
            region = next(((start, end) for start, end in stim_regions if start <= peak_time <= end),
                          (time[0], time[-1]))
            stim_start, stim_end = region
        else:
            stim_start = time[max(0, peak_idx - 10)]
            stim_end = time[min(len(time) - 1, peak_idx + 30)]

        stim_mask = (time >= stim_start) & (time <= stim_end)
        auc = simpson(trace_scaled[stim_mask] - baseline, dx=(time[1] - time[0]))

        # Decay tau
        decay_tau = np.nan
        if peak_idx < len(trace_scaled) - 5:
            decay_time = time[peak_idx:] - time[peak_idx]
            decay_trace = trace_scaled[peak_idx:]
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

        # Kinetics
        window_size = 20
        start_idx = max(0, peak_idx - window_size)
        end_idx = min(len(trace_scaled) - 1, peak_idx + window_size)
        sub_time = time[start_idx:end_idx + 1]
        sub_trace = trace_scaled[start_idx:end_idx + 1]
        duration = np.nan
        max_upstroke = np.nan
        max_downstroke = np.nan
        repolarisation_slope = np.nan

        half_max = baseline + 0.5 * delta_ca
        above_half = np.where(sub_trace >= half_max)[0]
        if len(above_half) >= 2:
            duration = sub_time[above_half[-1]] - sub_time[above_half[0]]

        dt = np.diff(sub_time)
        dy = np.diff(sub_trace)
        slopes = dy / dt
        peak_pos = peak_idx - start_idx
        if peak_pos > 0:
            max_upstroke = np.max(slopes[:peak_pos])
        if peak_pos < len(slopes):
            max_downstroke = np.min(slopes[peak_pos:])
        if peak_pos + 10 <= len(slopes):
            repolarisation_slope = np.mean(slopes[peak_pos:peak_pos + 10])

        results.append({
            'Cell_ID': cell_id,
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
        })

    return pd.DataFrame(results), trace, trace_scaled


def batch_process_csv(df, output_dir=None, plot_scaling=True):
    #df = pd.read_csv(file_path)
    time = df.iloc[:, 0].values
    all_results = []

    for i, col in enumerate(df.columns[1:], start=1):
        trace = df[col].values
        cell_id = col
        print(f"Processing {cell_id}...")

        cell_results, raw_trace, scaled_trace = analyze_trace(time, trace, cell_id)

        if plot_scaling:
            plot_trace_scaling(time, raw_trace, scaled_trace, cell_id, output_dir)

        all_results.append(cell_results)

    combined_df = pd.concat(all_results, ignore_index=True)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        out_csv = os.path.join(output_dir, "calcium_trace_peak_metrics.csv")
        combined_df.to_csv(out_csv, index=False)
        print(f"\nResults saved to: {out_csv}")

    return combined_df
'''
# === Utility functions ===

def exp_decay(t, A, tau, C):
    return A * np.exp(-t / tau) + C

def smooth_and_scale_trace(trace, window_length=7, polyorder=2):
    trace_smooth = savgol_filter(trace, window_length=window_length, polyorder=polyorder)
    trace_min = np.min(trace_smooth)
    trace_max = np.max(trace_smooth)
    trace_scaled = (trace_smooth - trace_min) / (trace_max - trace_min) if (trace_max - trace_min) != 0 else np.zeros_like(trace_smooth)
    return trace_smooth, trace_scaled

def plot_trace_scaling(time, raw_trace, scaled_trace, cell_id, peak_indices=None, output_dir=None):
    fig, ax1 = plt.subplots(figsize=(10, 4))
    ax1.plot(time, raw_trace, color='tab:blue', label='Raw Trace')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Raw Fura-2 Ratio', color='tab:blue')
    ax1.tick_params(axis='y', labelcolor='tab:blue')

    ax2 = ax1.twinx()
    ax2.plot(time, scaled_trace, color='tab:red', linestyle='--', label='Smoothed + Scaled Trace')
    ax2.set_ylabel('Scaled Signal [0–1]', color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red')

    if peak_indices is not None and len(peak_indices) > 0:
        ax2.scatter(time[peak_indices], scaled_trace[peak_indices], color='black', s=40, label='Detected Peaks', zorder=5)

    fig.suptitle(f'Cell {cell_id}: Trace with Detected Peaks')
    fig.tight_layout()
    fig.legend(loc='upper right')
    plt.grid(True)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(os.path.join(output_dir, f'{cell_id}_trace_plot.png'))
        plt.close()
    else:
        plt.show()

# === Analysis logic ===

def find_best_peak_window(time, trace_scaled, window_size=30, step_size=5, min_peaks=1):
    """
    Slide a window and find the segment with the most peaks.
    Returns: (start_time, end_time)
    """
    max_peaks = 0
    best_range = (time[0], time[-1])

    for start in np.arange(time[0], time[-1] - window_size, step_size):
        end = start + window_size
        mask = (time >= start) & (time <= end)
        sub_trace = trace_scaled[mask]
        if len(sub_trace) < 10:
            continue
        peaks, _ = find_peaks(sub_trace, prominence=0.1)
        if len(peaks) > max_peaks and len(peaks) >= min_peaks:
            max_peaks = len(peaks)
            best_range = (start, end)

    return best_range

def analyze_trace(time, trace, cell_id, stim_regions=None, region_range=None,
                  smoothing_window=7, smoothing_poly=2):
    if region_range:
        mask = (time >= region_range[0]) & (time <= region_range[1])
        time = time[mask]
        trace = trace[mask]

    trace_smooth, trace_scaled = smooth_and_scale_trace(trace, smoothing_window, smoothing_poly)
    baseline_end_idx = int(0.1 * len(trace_scaled))
    baseline = np.mean(trace_scaled[:baseline_end_idx])
    noise_std = np.std(trace_scaled[:baseline_end_idx])
    peaks, _ = find_peaks(trace_scaled, prominence=0.2)

    results = []
    for peak_idx in peaks:
        peak_time = time[peak_idx]
        peak = trace_scaled[peak_idx]
        delta_ca = (peak - baseline) / baseline if baseline != 0 else np.nan
        time_to_peak = peak_time - time[0]

        stim_start = time[max(0, peak_idx - 10)]
        stim_end = time[min(len(time) - 1, peak_idx + 30)]
        stim_mask = (time >= stim_start) & (time <= stim_end)
        auc = simpson(trace_scaled[stim_mask] - baseline, dx=(time[1] - time[0]))

        decay_tau = np.nan
        if peak_idx < len(trace_scaled) - 5:
            decay_time = time[peak_idx:] - time[peak_idx]
            decay_trace = trace_scaled[peak_idx:]
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
        window_size = 20
        start_idx = max(0, peak_idx - window_size)
        end_idx = min(len(trace_scaled) - 1, peak_idx + window_size)
        sub_time = time[start_idx:end_idx + 1]
        sub_trace = trace_scaled[start_idx:end_idx + 1]

        duration = max_upstroke = max_downstroke = repolarisation_slope = np.nan
        half_max = baseline + 0.5 * delta_ca
        above_half = np.where(sub_trace >= half_max)[0]
        if len(above_half) >= 2:
            duration = sub_time[above_half[-1]] - sub_time[above_half[0]]

        dt = np.diff(sub_time)
        dy = np.diff(sub_trace)
        slopes = dy / dt
        peak_pos = peak_idx - start_idx
        if peak_pos > 0:
            max_upstroke = np.max(slopes[:peak_pos])
        if peak_pos < len(slopes):
            max_downstroke = np.min(slopes[peak_pos:])
        if peak_pos + 10 <= len(slopes):
            repolarisation_slope = np.mean(slopes[peak_pos:peak_pos + 10])

        results.append({
            'Cell_ID': cell_id,
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
        })

    return pd.DataFrame(results), trace, trace_scaled, peaks

# === Main driver ===

def batch_process_csv(df, output_dir=None, plot_scaling=True, auto_window=True,
                      window_size=30, step_size=5):
    time = df.iloc[:, 0].values
    all_results = []

    for i, col in enumerate(df.columns[1:], start=1):
        trace = df[col].values
        cell_id = col
        print(f"\n--- Processing {cell_id} ---")

        _, trace_scaled = smooth_and_scale_trace(trace)
        region_range = None
        if auto_window:
            region_range = find_best_peak_window(time, trace_scaled,
                                                 window_size=window_size,
                                                 step_size=step_size)
            print(f"   ⤷ Using window {region_range[0]:.2f}s to {region_range[1]:.2f}s")

        cell_results, raw_trace, scaled_trace, peak_indices = analyze_trace(
            time, trace, cell_id, region_range=region_range
        )

        if plot_scaling:
            plot_trace_scaling(time, raw_trace, scaled_trace, cell_id,
                               peak_indices=peak_indices, output_dir=output_dir)

        all_results.append(cell_results)

    combined_df = pd.concat(all_results, ignore_index=True)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        out_csv = os.path.join(output_dir, "calcium_trace_peak_metrics.csv")
        combined_df.to_csv(out_csv, index=False)
        print(f"\n✅ Results saved to: {out_csv}")

    return combined_df