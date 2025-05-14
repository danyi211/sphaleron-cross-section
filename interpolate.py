# interpolate the cross section vs mass
# the data is from a text file, the first column is the xs and the second column is the mass point
# plot the xs vs mass and do interpolation to get the xs at a given mass
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.ROOT)
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline

# Load the data from the text file
data = np.loadtxt('theory_xs.txt')
xs = data[:, 0]  # Cross section values
mass = data[:, 1]  # Mass points

# Create an interpolation function
# interp_func = interp1d(mass, xs, kind='linear', fill_value='extrapolate')
interp_func = CubicSpline(mass, xs)

# Define a range of mass values for interpolation
m_min = 8
m_max = 11
mass_range = np.linspace(m_min, m_max, 200)
# Calculate the interpolated cross section values
interpolated_xs = interp_func(mass_range)

# ------ import another data ------
# Load member 0
data0 = np.loadtxt('sphaleron_xs_CT14nnlo_0.txt')
xs0 = data0[:, 0]  # Cross section values
mass0 = data0[:, 1]  # Mass points
# Create an interpolation function
interp_func0 = CubicSpline(mass0, xs0)
# Define a range of mass values for interpolation
m_min0 = 8
m_max0 = 11
mass_range0 = np.linspace(m_min0, m_max0, 200)
# Calculate the interpolated cross section values
interpolated_xs0 = interp_func0(mass_range0)

# Load member 1
data1 = np.loadtxt('sphaleron_xs_CT14nnlo_1.txt')
xs1 = data1[:, 0]  # Cross section values
mass1 = data1[:, 1]  # Mass points
# Create an interpolation function
interp_func1 = CubicSpline(mass1, xs1)
# Define a range of mass values for interpolation
m_min1 = 8
m_max1 = 11
mass_range1 = np.linspace(m_min1, m_max1, 200)
# Calculate the interpolated cross section values
interpolated_xs1 = interp_func1(mass_range1)

# Load member 2
data2 = np.loadtxt('sphaleron_xs_CT14nnlo_2.txt')
xs2 = data2[:, 0]  # Cross section values
mass2 = data2[:, 1]  # Mass points
# Create an interpolation function
interp_func2 = CubicSpline(mass2, xs2)
# Define a range of mass values for interpolation
m_min2 = 8
m_max2 = 11
mass_range2 = np.linspace(m_min2, m_max2, 200)
# Calculate the interpolated cross section values
interpolated_xs2 = interp_func2(mass_range2)

# ------ Plotting ------
plt.figure()
plt.plot(mass, xs, 'o', label='CT14nnlo from theorist')
plt.plot(mass_range, interpolated_xs, '-')
plt.plot(mass0, xs0, 'o', label='CT14nnlo member 0 from Danyi')
plt.plot(mass_range0, interpolated_xs0, '-')
plt.plot(mass1, xs1, 'o', label='CT14nnlo member 1 from Danyi')
plt.plot(mass_range1, interpolated_xs1, '-')
plt.plot(mass2, xs2, 'o', label='CT14nnlo member 2 from Danyi')
plt.plot(mass_range2, interpolated_xs2, '-')

# log scale on y-axis
plt.yscale('log')

# Add a text with PDF version
plt.text(0.65, 0.95, r'$\sqrt{s} = 13$ TeV', transform=plt.gca().transAxes,
         fontsize=20, verticalalignment='top')
plt.text(0.65, 0.90, r'$p_{sph}=1$', transform=plt.gca().transAxes,
         fontsize=20, verticalalignment='top')
plt.text(0.65, 0.84, 'CT14nnlo PDF', transform=plt.gca().transAxes,
         fontsize=20, verticalalignment='top')

plt.xlabel(r'$E_{sph}$ [TeV]')
plt.ylabel('Cross Section [fb]')
# plt.title('Cross Section vs Mass with Interpolation')
plt.legend(loc='lower left')
plt.grid(True)
plt.savefig('xs_vs_mass_interpolate_cubicSpline_log.png')

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
