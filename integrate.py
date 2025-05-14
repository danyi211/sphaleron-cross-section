import lhapdf
import numpy as np
from scipy.integrate import dblquad
import sys


# --- Configuration ---
# arguments: E_thr in TeV
if len(sys.argv) != 2:
    print("Usage: python integrate.py <E_thr_TEV>")
    print("E_thr_TEV: Threshold energy in TeV (e.g., 9 for 9 TeV)")
    exit(1)
E_THR_TEV = float(sys.argv[1])  # Threshold energy in TeV

PDF_SET_NAME = "CT14nnlo"  # PDF set
P_SPH = 1.0  # Dimensionless constant for sphaleron probability (user defined)
M_W_GEV = 80.379  # W boson mass in GeV (PDG 2023)
E_THR_GEV = E_THR_TEV * 1e3  # Threshold energy in GeV (e.g., 9 TeV, user defined)
S_PP_GEV = 13000.0  # Proton-proton center-of-mass energy in GeV (e.g., 13 TeV)

# --- Load PDF Set ---
try:
    pdf = lhapdf.mkPDF(PDF_SET_NAME, 0) # Member 0 is usually the central PDF
except lhapdf.LhapdfException as e:
    print(f"Error loading PDF set {PDF_SET_NAME}: {e}")
    print("Please ensure LHAPDF is configured and the PDF set is installed.")
    print("You can install PDF sets using: lhapdf install PDF_SET_NAME")
    exit()

# Get PDF operational range
X_MIN_PDF = pdf.xMin
X_MAX_PDF = pdf.xMax
Q2_MIN_PDF = pdf.q2Min
Q2_MAX_PDF = pdf.q2Max

# --- Constants and Derived Values ---
s_pp_val = S_PP_GEV**2
E_thr_sq_val = E_THR_GEV**2
sigma_hat_const = P_SPH / (M_W_GEV**2)

# Conversion factor: 1 GeV^-2 = 0.389379 mb = 3.89379e8 pb = 3.89379e11 fb
GEV_SQ_TO_PB = 3.89379e8
GEV_SQ_TO_FB = GEV_SQ_TO_PB * 1e3

# --- Integrand Function ---
# This function will be called by dblquad.
# Note: dblquad passes arguments as func(y, x), so here func(x2, x1)
def integrand(x2, x1, parton_id1, parton_id2):
    """
    Calculates the value f1(x1,Q^2) * f2(x2,Q^2) * sigma_hat_factor
    for the given x1, x2 and parton IDs.
    The sigma_hat_factor is p_sph / m_W^2.
    The Heaviside function is implicitly handled by the s_hat_val check.
    """
    # Ensure x1 and x2 are within the valid range of the PDF
    if not (X_MIN_PDF <= x1 <= X_MAX_PDF and X_MIN_PDF <= x2 <= X_MAX_PDF):
        return 0.0

    s_hat_val = x1 * x2 * s_pp_val

    # Apply Heaviside Theta function condition
    if s_hat_val < E_thr_sq_val:
        return 0.0

    # Set the PDF scale Q^2 = s_hat
    current_Q2 = s_hat_val
    if current_Q2 < Q2_MIN_PDF or current_Q2 > Q2_MAX_PDF:
        print(f"Warning: Q2={current_Q2} out of PDF range [{Q2_MIN_PDF}, {Q2_MAX_PDF}] for x1={x1}, x2={x2}")
        return 0.0

    try:
        # LHAPDF returns x*f(x, Q^2), so divide by x to get f(x, Q^2)
        # Ensure x1 and x2 are not zero (they shouldn't be if x_min_pdf > 0)
        f1 = pdf.xfxQ2(parton_id1, x1, current_Q2) / x1
        f2 = pdf.xfxQ2(parton_id2, x2, current_Q2) / x2
    except ZeroDivisionError:
        return 0.0
    except lhapdf.LhapdfException as e: # Catch errors from LHAPDF (e.g., out of range)
        print(f"LHAPDF error for x1={x1}, x2={x2}, Q2={current_Q2}: {e}")
        return 0.0

    return f1 * f2 * sigma_hat_const

# --- Perform Numerical Integration ---
total_hadronic_cross_section_gev_sq = 0.0

# Define parton PDG IDs to sum over.
# PDG IDs: d=1, u=2, s=3, c=4, b=5
parton_ids_to_sum = [
    1, 2, 3, 4, 5,  # d, u, s, c, b
    -1, -2, -3, -4, -5 # dbar, ubar, sbar, cbar, bbar
]

print(f"Calculating hadronic cross-section for proton-proton at sqrt(s) = {S_PP_GEV} GeV")
print(f"Sphaleron energy threshold E_thr = {E_THR_GEV} GeV")
print(f"Using PDF set: {PDF_SET_NAME}")
print(f"sigma_hat constant factor (p_sph/m_W^2) = {sigma_hat_const:.4e} GeV^-2")
print(f"Integrating for x1, x2 from {X_MIN_PDF} to {X_MAX_PDF}")

# Integration limits for x1
x1_lower_limit = X_MIN_PDF
x1_upper_limit = X_MAX_PDF

for id1 in parton_ids_to_sum:
    for id2 in parton_ids_to_sum:
        if id1 * id2 < 0:
            # Skip combinations of quark and antiquark
            continue
        print(f"Integrating for initial state: parton1_id={id1}, parton2_id={id2}")

        x2_lower_limit_func = lambda x1_val: X_MIN_PDF # Min x for PDF
        x2_upper_limit_func = lambda x1_val: X_MAX_PDF # Max x for PDF

        integral_val, integral_err = dblquad(
            integrand,
            x1_lower_limit,
            x1_upper_limit,
            x2_lower_limit_func, # x2 lower limit (can be a function of x1)
            x2_upper_limit_func, # x2 upper limit (can be a function of x1)
            args=(id1, id2),
            epsabs=1.49e-08, # Default absolute error
            epsrel=1.49e-08  # Default relative error
        )
        
        # Only left-handed quarks contribute
        integral_val /= 4.0
        
        # For identical partons, divide by 3, as the colour of the pair must be different
        if id1 == id2:
            integral_val /= 3.0

        print(f"  Contribution from ({id1}, {id2}): {integral_val:.4e} GeV^-2 (error: {integral_err:.2e})")
        total_hadronic_cross_section_gev_sq += integral_val
    

# --- Print Results ---
total_hadronic_cross_section_pb = total_hadronic_cross_section_gev_sq * GEV_SQ_TO_PB
total_hadronic_cross_section_fb = total_hadronic_cross_section_gev_sq * GEV_SQ_TO_FB

print(f"\n--- Total Hadronic Cross-Section ---")
print(f"Sum over specified parton combinations:")
print(f"Cross-section (GeV^-2): {total_hadronic_cross_section_gev_sq:.6e}")
# print(f"Sigma (pb): {total_hadronic_cross_section_pb:.6e}")
print(f"Cross-section (fb): {total_hadronic_cross_section_fb:.6e}")