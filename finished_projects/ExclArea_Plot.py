import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erfc
from matplotlib.backends.backend_pdf import PdfPages
from scipy.special import logsumexp

sigKa = np.array([13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 1, 0.75, 0.5, 0.25, 0])
AreaKa =     np.array([0.86, 1.07, 1.52, 2.50, 3.23, 4.12, 5.13, 6.50, 9.23, 14.34, 24.94, 29.88, 35.77, 44.57, 54.62, 67.28, 87.89, 112.10, 141.43, 175.75, 219.02, 263.88, 308.81])
AreaKa_err = np.array([1.31, 2.18, 2.11, 2.32, 2.29, 1.95, 2.21, 2.30, 2.30, 2.37, 2.30, 2.52, 2.63, 2.74, 2.78, 2.97, 3.21, 3.53, 3.75, 4.15, 4.58, 4.90, 5.42])
chi2Ka = np.array([1.66247, 1.66475, 1.65947, 1.65515, 1.65082, 1.6473, 1.63918, 1.63777, 1.63336, 1.62737, 1.6334, 1.63111, 1.61772, 1.61145, 1.60784, 1.61561, 1.6073, 1.64207, 1.70633, 1.80732, 1.90965, 2.00597, 2.12548])
AICKa = np.array([3.9461e+06, 3.94688e+06, 3.94745e+06, 3.94817e+06, 3.94874e+06, 3.94922e+06, 3.94982e+06, 3.95038e+06, 3.95113e+06, 3.95215e+06, 3.95394e+06, 3.95468e+06, 3.95563e+06, 3.95687e+06, 3.9587e+06, 3.9608e+06, 3.96356e+06, 3.96715e+06, 3.97159e+06, 3.97693e+06, 3.98318e+06, 3.9894e+06, 3.99589e+06])
BICKa = np.array([3.94631e+06, 3.94709e+06, 3.94766e+06, 3.94838e+06, 3.94895e+06, 3.94943e+06, 3.95003e+06, 3.95059e+06, 3.95134e+06, 3.95236e+06, 3.95415e+06, 3.95489e+06, 3.95584e+06, 3.95708e+06, 3.95891e+06, 3.96101e+06, 3.96377e+06, 3.96736e+06, 3.9718e+06, 3.97714e+06, 3.98339e+06, 3.98961e+06, 3.9961e+06])
DNKa= np.array([4.52325, 4.52327, 4.52326, 4.52325, 4.52323, 4.52321, 4.52322, 4.52322, 4.52323, 4.52326, 4.52331, 4.52334, 4.52341, 4.52344, 4.52362, 4.52379, 4.5239, 4.52417, 4.5245, 4.52496, 4.52539, 4.52576, 4.52622])
#const Double_t pStart = 0.45, pEnd = 0.55, step = 0.1;
#const Bool_t FitKaonExclComp = true;
#const Bool_t FitProtonExclComp = false;
#manualNGauss = 4;
#manualMeans = {-6, -5, 0, 1.5};
#manualSigmas = {1, 1, 1, 1};
#manualAmps = {1e4, 1e4, 1e2, 1e1};
#par[offM + 2] = -0.271677;
#par[offS + 2] = 1.08472;
#par[offA1 + 2] = 78.0623;
#par[offA2 + 2] = par[offA1 + 2];
#//par[offA2 + 3] = 0;
#std::cout << "Mean: " << par[offM + 2] << ", Sigma: " << par[offS + 2] << ", Amp2: "<< par[offA2 + 2] << std::endl;

sigPr = np.array([13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 1, 0.75, 0.5, 0.25, 0])
AreaPr = np.array([0, 0, 0, 0.09, 0.08, 0.16, 0.65, 0.41, 0.82, 1.51, 4.92, 6.89, 9.30, 12.81, 16.63, 23.85, 34.82, 49.95, 69.80, 93.46, 123.34, 157.39, 189.79])
AreaPr_err = np.array([0, 0, 0.03, 1.49, 1.55, 1.50, 0.52, 1.50, 1.56, 1.66, 1.64, 1.67, 1.67, 1.97, 1.84, 1.97, 2.10, 2.36, 2.69, 3.08, 3.44, 3.85, 3.99])
chi2Pr = np.array([1.43889, 1.43374, 1.43729, 1.4323, 1.43458, 1.43238, 1.43212, 1.43347, 1.43785, 1.42486, 1.41977, 1.43227, 1.41255, 1.40898, 1.4024, 1.39478, 1.39959, 1.38814, 1.38673, 1.36941, 1.35767, 1.36301, 1.37432])
AICPr = np.array([1.02512e+06, 1.02593e+06, 1.02621e+06, 1.02632e+06, 1.02648e+06, 1.02661e+06, 1.02673e+06, 1.0268e+06, 1.02697e+06, 1.02718e+06, 1.02755e+06, 1.02771e+06, 1.02803e+06, 1.02849e+06, 1.02908e+06, 1.02992e+06, 1.03131e+06, 1.03316e+06, 1.03547e+06, 1.03813e+06, 1.04158e+06, 1.04554e+06, 1.04927e+06])
BICPr = np.array([1.02526e+06, 1.02607e+06, 1.02636e+06, 1.02646e+06, 1.02662e+06, 1.02675e+06, 1.02687e+06, 1.02695e+06, 1.02711e+06, 1.02732e+06, 1.02769e+06, 1.02786e+06, 1.02818e+06, 1.02864e+06, 1.02922e+06, 1.03006e+06, 1.03145e+06, 1.03331e+06, 1.03561e+06, 1.03827e+06, 1.04172e+06, 1.04569e+06, 1.04941e+06])
DNPr = np.array([4.53047, 4.53047, 4.53051, 4.53045, 4.53048, 4.53051, 4.53066, 4.53052, 4.53056, 4.53063, 4.53048, 4.53041, 4.53036, 4.53039, 4.53053, 4.5305, 4.53064, 4.53064, 4.53062, 4.53055, 4.53043, 4.53043, 4.5305])
#const Double_t pStart = 0.9, pEnd = 1.0, step = 0.1;
#const Bool_t FitKaonExclComp = false;
#const Bool_t FitProtonExclComp = true;
#manualNGauss = 3;
#manualMeans = {-5, -0.5, 1};
#manualSigmas = {1, 1, 1};
#manualAmps = {1e3, 1e2, 1e2};
#par[offM + 1] = -1;
#par[offS + 1] = 1.46314;
#par[offA1 + 1] = 22.3524;
#par[offA2 + 1] = par[offA1 + 1];
#//par[offA2 + 2] = 0;
#std::cout << "Mean: " << par[offM + 1] << ", Sigma: " << par[offS + 1] << ", Amp2: "<< par[offA2 + 1] << std::endl;


def tail_area(n, A, sigma0):
    return A * erfc(n / (np.sqrt(2) * sigma0))

poptKa, pcovKa = curve_fit(tail_area, sigKa, AreaKa, sigma=AreaKa_err, absolute_sigma=True, p0=[AreaKa.max(), 1.0])
Aka, sigma0_ka = poptKa

poptPr, pcovPr = curve_fit(tail_area, sigPr, AreaPr, sigma=AreaPr_err, absolute_sigma=True, p0=[AreaPr.max(), 1.0])
Apr, sigma0_pr = poptPr

def dn_weights(D_over_N, N_ref=None, target_ratio=0.1, quantile=0.9):
    r = np.asarray(D_over_N, float)
    dr = r - r.min()
    if N_ref is None:
        pos = dr[dr > 0]
        if pos.size == 0:
            return np.ones_like(r) / r.size, 0.0
        dq = np.quantile(pos, quantile)
        N_ref = 2.0 * np.log(1.0/target_ratio) / dq
    logw = -0.5 * N_ref * dr
    logw -= logsumexp(logw)
    return np.exp(logw), N_ref


with PdfPages('PeakAreas_Kaon_Proton.pdf') as pdf:
    fig1, ax1 = plt.subplots()
    n_fit_ka = np.linspace(sigKa.min(), sigKa.max(), 200)
    ax1.errorbar(sigKa, AreaKa, yerr=AreaKa_err, fmt='o',markersize=3.5, color="blue", label='Kaon-Area ± error', capsize=3)
    ax1.plot(n_fit_ka, tail_area(n_fit_ka, Aka, sigma0_ka), '-', label=f'Fit: A={Aka:.1f}, $\sigma_{0}$={sigma0_ka:.2f}')
    ax1.set_xlabel('Sigma')
    ax1.set_ylabel('Kaon-Peak-area', rotation=90)
    ax1.tick_params(axis='y')

    ax1b = ax1.twinx()
    ax1b.plot(sigKa, chi2Ka, marker='o', markersize=3.5, color="red", linestyle='', label='$\chi^{2}$ / NDF')
    ax1b.set_ylabel('$\chi^{2}$ / NDF', rotation=90)
    ax1b.set_ylim(0,5)
    ax1b.tick_params(axis='y')
    plt.title('Kaon-Peak-Area vs Excl-Sigma')
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1b.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1+ labels2, loc= "best")
    pdf.savefig(fig1)
    plt.close(fig1)

    fig2, ax2 = plt.subplots()
    n_fit_pr = np.linspace(sigPr.min(), sigPr.max(), 200)
    ax2.errorbar(sigPr, AreaPr, yerr=AreaPr_err, fmt='o',markersize=3.5, color="blue", label='Proton-Area ± error', capsize=3)
    ax2.plot(n_fit_pr, tail_area(n_fit_pr, Apr, sigma0_pr), '-', label=f'Fit: A={Apr:.1f}, $\sigma_{0}$={sigma0_pr:.2f}')
    ax2.set_xlabel('Sigma')
    ax2.set_ylabel('Proton-Peak-area', rotation=90)
    ax2.tick_params(axis='y')

    ax2b = ax2.twinx()
    ax2b.plot(sigPr, chi2Pr, marker='o', markersize=3.5, color="red", linestyle='', label='$\chi^{2}$ / NDF')
    ax2b.set_ylabel('$\chi^{2}$ / NDF', rotation=90)
    ax2b.set_ylim(0,5)
    ax2b.tick_params(axis='y')
    plt.title('Proton-Peak-Area vs Excl-Sigma')
    lines3, labels3 = ax2.get_legend_handles_labels()
    lines4, labels4 = ax2b.get_legend_handles_labels()
    ax2.legend(lines3 + lines4, labels3+ labels4, loc= "best")
    pdf.savefig(fig2)
    plt.close(fig2)
        
    fig3, ax3 = plt.subplots()
    w_ka, Nref_ka = dn_weights(DNKa)  
    ax3.plot(sigKa, w_ka, 'o', markersize=3.5)
    ax3.set_xlabel('Sigma-Exclusion')
    ax3.set_ylabel('D/N-weight w($\sigma$)')
    ax3.set_title(f'Kaon: Poisson-Deviance/Number-of-Events-weights vs Excl-Sigma')
    ax3.set_ylim(0, 1)
    ax3.grid(True, alpha=0.2)
    pdf.savefig(fig3, bbox_inches='tight', pad_inches=0.1); plt.close(fig3)

    fig4, ax4 = plt.subplots()
    w_pr, Nref_pr = dn_weights(DNPr)
    ax4.plot(sigPr, w_pr, 'o', markersize=3.5)
    ax4.set_xlabel('Sigma-Exclusion')
    ax4.set_ylabel('D/N-weight w($\sigma$)')
    ax4.set_title(f'Proton: Poisson-Deviance/Number-of-Events-weights vs Excl-Sigma')
    ax4.set_ylim(0, 1)
    ax4.grid(True, alpha=0.2)
    pdf.savefig(fig4, bbox_inches='tight', pad_inches=0.1); plt.close(fig4)
    

