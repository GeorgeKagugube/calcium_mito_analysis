
import matplotlib.pyplot as plt
import seaborn as sns

# Set plotting style
sns.set_theme(style='whitegrid', palette='muted')
plt.rcParams["font.family"] = "DejaVu Sans"
# Set default figure size
plt.rcParams['figure.figsize'] = (10, 6)    
# Set default font size
plt.rcParams['font.size'] = 15  
# Set default line width
plt.rcParams['lines.linewidth'] = 2.5

def plot_fura2_rhodamine(time, fura2_trace, rhodamine_trace, position='upper left'):
    """ Plot Fura-2 AM and Rhodamine 123 traces on the same figure with dual y-axes.
    Parameters:
    time (array-like): Time points for the x-axis.
    fura2_trace (array-like): Fura-2 AM trace for the first condition.
    fura2_trace2 (array-like): Fura-2 AM trace for the second condition.
    rhodamine_trace (array-like): Rhodamine 123 trace for the first condition.
    rhodamine_trace2 (array-like): Rhodamine 123 trace for the second condition.
    position (str): Position of the legend in the plot. Default is 'upper left'.
    """
    fig, ax1 = plt.subplots(figsize=(7, 5))

    # Plot Fura-2 AM on the left y-axis
    ax1.plot(time, fura2_trace, color='blue', label='Fura-2 AM (Ca²⁺)', linestyle='-')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Fura-2 AM Ratio (Ca²⁺)', color='black')
    ax1.tick_params(axis='y', labelcolor='black')


    # Create a second y-axis for Rhodamine 123
    ax2 = ax1.twinx()
    ax2.plot(time, rhodamine_trace, color='black', label='Rhodamine 123 (ΔΨm)')
    ax2.set_ylabel('Rhodamine 123 (ΔΨm)', color='black')
    ax2.tick_params(axis='y', labelcolor='black')

    # Title and combined legend
    #plt.title(title)
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc=position)

    plt.grid(False)
    plt.tight_layout()
    #plt.savefig(path+'/calcium_rhodamine123_mut.png')
    plt.show()