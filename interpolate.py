# interpolate the cross section vs mass
# the data is from a text file, the first column is the xs and the second column is the mass point
# plot the xs vs mass and do interpolation to get the xs at a given mass
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.ROOT)
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline

def load_interpolate(filename, m_min=8, m_max=11):
    """
    Load data from a file and perform cubic spline interpolation.
    The first column is the x values (mass) and the second column is the y values (cross section).
    """
    # Load the data from the text file
    data = np.loadtxt(filename)
    xs = data[:, 0]  # Cross section values
    mass = data[:, 1]  # Mass points
    # Create an interpolation function
    interp_func = CubicSpline(mass, xs)
    # Define a range of mass values for interpolation
    mass_range = np.linspace(m_min, m_max, 200)
    # Calculate the interpolated cross section values
    interpolated_xs = interp_func(mass_range)
    result = {'mass': mass, 'xs': xs, 'mass_range': mass_range, 'interpolated_xs': interpolated_xs}
    return result

m_min = 8
m_max = 11

results_dict = {}

# results_dict['CT14nnlo from theorist'] = load_interpolate('data/theory_xs.txt', m_min, m_max)
results_dict['CT10 (NtPoints=1,000) from theorist'] = load_interpolate('data/theory_CT10_steps1000.txt', m_min, m_max)
results_dict['CT10 (NtPoints=10,000) from theorist'] = load_interpolate('data/theory_CT10_steps10000.txt', m_min, m_max)
# results_dict['CT14nnlo from Danyi'] = load_interpolate('data/sphaleron_xs_CT14nnlo_0.txt', m_min, m_max)
# results_dict[r'CT14nnlo ($x_1$-$x_2$, $\Delta N = -1$) from Danyi'] = load_interpolate('data/sphaleron_xs_CT14nnlo_0_pNCS0.txt', m_min, m_max)
# results_dict[r'CT14nnlo ($\tau$-y, $\Delta N = -1$) from Danyi'] = load_interpolate('data/sphaleron_xs_tau_y_CT14nnlo_0_pNCS0.txt', m_min, m_max)
# results_dict[r'CT14nnlo ($\Delta N = 1$) from Danyi'] = load_interpolate('data/sphaleron_xs_CT14nnlo_0_pNCS1.txt', m_min, m_max)
# results_dict['CT14nnlo (Q = 100 GeV) from Danyi'] = load_interpolate('data/sphaleron_xs_CT14nnlo_0_Q100.txt', m_min, m_max)
# results_dict['CT10 from Danyi'] = load_interpolate('data/sphaleron_xs_CT10_0.txt', m_min, m_max)
results_dict['CT10 ($\tau$-y, $\Delta N = -1$) from Danyi'] = load_interpolate('data/sphaleron_xs_tau_y_CT10_0_pNCS0.txt', m_min, m_max)
# results_dict['NNPDF31_nnlo_as_0118 from Danyi'] = load_interpolate('data/sphaleron_xs_NNPDF31_nnlo_as_0118_0.txt', m_min, m_max)
# results_dict['NNPDF31_nnlo_as_0118 ($\tau$-y, $\Delta N = -1$) from Danyi'] = load_interpolate('data/sphaleron_xs_tau_y_NNPDF31_nnlo_as_0118_0_pNCS0.txt', m_min, m_max)

colors = ['blue', 'red', 'green', 'orange', 'purple']

# ------ Plotting ------
plt.figure()
for i, (label, result) in enumerate(results_dict.items()):
    mass = result['mass']
    xs = result['xs']
    mass_range = result['mass_range']
    interpolated_xs = result['interpolated_xs']
    # Plot the original data points
    marker_style = '^' if 'from theorist' in label else 'o'
    plt.plot(mass, xs, marker_style, color=colors[i])
    # Plot the interpolated data
    line_style = '--' if 'from theorist' in label else '-'
    plt.plot(mass_range, interpolated_xs, line_style, label=label, color=colors[i])

# log scale on y-axis
plt.yscale('log')

# Add a text with PDF version
plt.text(0.65, 0.95, r'$\sqrt{s} = 13$ TeV', transform=plt.gca().transAxes,
         fontsize=20, verticalalignment='top')
plt.text(0.65, 0.90, r'$p_{sph}=1$', transform=plt.gca().transAxes,
         fontsize=20, verticalalignment='top')
# plt.text(0.65, 0.84, 'CT14nnlo PDF', transform=plt.gca().transAxes,
#          fontsize=20, verticalalignment='top')

plt.xlabel(r'$E_{sph}$ [TeV]')
plt.ylabel('Cross Section [fb]')
# plt.title('Cross Section vs Mass with Interpolation')
plt.legend(loc='lower left')
plt.grid(True)
plt.savefig('xs_vs_mass_interpolate_cubicSpline_log_compare_method.png')

# print the interpolated xs at 8.5, 8.75, 9, 9.25, 9.5
# masses_to_check = [8.5, 8.75, 9, 9.25, 9.5]
# for m in masses_to_check:
#     print(f"Interpolated cross section at {m} TeV: {interp_func(m)} fb")
# Output:
# Interpolated cross section at 8.5 TeV: 28.3304 fb
# Interpolated cross section at 8.75 TeV: 14.512099114949768 fb
# Interpolated cross section at 9 TeV: 7.11033 fb
# Interpolated cross section at 9.25 TeV: 3.4119388499679584 fb
# Interpolated cross section at 9.5 TeV: 1.55813 fb


# Draw a ratio plot between two curves and original curve on the same plot
# Calculate the ratio of the two curves
def draw_ratio_plot(x, y1, y2, label1, label2):
    """
    Draw a ratio plot between two curves and the original curve on the same plot.
    """
    # Calculate the ratio of the two curves
    ratio = y1 / y2

    # Create a new figure
    fig, ax = plt.subplots()

    # Plot the original curve
    # ax.plot(x, y1, label=label1, color='blue')

    # Plot the second curve
    # ax.plot(x, y2, label=label2, color='red')

    # Plot the ratio
    ax.plot(x, ratio, label='Ratio', color='green')

    # Set labels and title
    ax.set_xlabel('Mass [TeV]')
    ax.set_ylabel('Ratio')
    # ax.set_title('Cross Section Ratio Plot')
    ax.legend()
    ax.grid(True)
    # Set y-axis to log scale
    # ax.set_yscale('log')
    # Add a horizontal line at y=1
    ax.axhline(y=1, color='black', linestyle='--')
    # Show the plot
    plt.savefig('xs_ratio_plot_CT10.png')
# Example usage
draw_ratio_plot(results_dict['CT10 ($\tau$-y, $\Delta N = -1$) from Danyi']['mass_range'],
                results_dict['CT10 ($\tau$-y, $\Delta N = -1$) from Danyi']['interpolated_xs'],
                results_dict['CT10 (NtPoints=10,000) from theorist']['interpolated_xs'],
                'CT10 ($\tau$-y, $\Delta N = -1$) from Danyi',
                'CT10 (NtPoints=10,000) from theorist')


