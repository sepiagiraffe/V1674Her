from astropy.modeling import models, fitting

from scipy.optimize import minimize

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import emcee
import corner
import random
import warnings
warnings.filterwarnings("ignore");

#read in data
lc_c = pd.read_csv('./lc_all.csv',delimiter=',');
VLBA = pd.read_csv('./VLBA.csv', delimiter=',');
lc_c =  lc_c.iloc[:, 1:-1];
x = np.linspace(lc_c.iloc[0,0], lc_c.iloc[-1,0], num=100);

#get initial model using astropy
p_init = models.BrokenPowerLaw1D(amplitude=2, x_break=19, alpha_1=0.7, alpha_2=-1); ## init power law
fit_p = fitting.LevMarLSQFitter(); ## which fitter to use
p = fit_p(p_init, lc_c.iloc[:,0], lc_c.iloc[:,2], maxiter=10000); ## fit the data


#make data nicer for emcee
t = lc_c.iloc[:,0].to_numpy(dtype=np.float128);
y = lc_c.iloc[:,2].to_numpy(dtype=np.float128);
yerr = lc_c.iloc[:,3].to_numpy(dtype=np.float128);


# t,y,yerr = [np.array(x[:-2]) for x in [t,y,yerr]];
# t,y,yerr = [np.array(x[1:]) for x in [t,y,yerr]];
x0 = np.linspace(t[0], t[-1], num=100);

#Maximum Likelihood Estimate
# def broken_powerlaw(t, t_b, beta_1, beta_2, F_c, s = 1.0):
#     #Granot & sari 2002
    
#     a = (t/t_b)**(-s*beta_1);
#     b = (t/t_b)**(-s*beta_2);

#     F_nu = F_c*(a+b)**(-1/s);
    
#     return F_nu


# def log_likelihood(theta, t, y, yerr, s = 1):
#     #Granot & sari 2002
#     t_b, beta_1, beta_2, F_c = theta;

#     model = broken_powerlaw(t,t_b,beta_1,beta_2,F_c)

#     return -0.5 * np.sum((y - model) ** 2 / yerr + np.log(yerr))

# nll = lambda *args: -log_likelihood(*args);
# okay now what if we put s back in it
def broken_powerlaw(t, t_b, s, beta_1, beta_2, F_c):
    #Granot & sari 2002
    
    a = (t/t_b)**(-s*beta_1);
    b = (t/t_b)**(-s*beta_2);

    F_nu = F_c*(a+b)**(-1/s);
    
    return F_nu

def log_likelihood(theta, t, y, yerr):
    t_b, s, beta_1, beta_2, F_c = theta;
    #Granot & sari 2002

    model = broken_powerlaw(t, t_b, s, beta_1, beta_2, F_c);
    
    return -0.5 * np.sum((y - model) ** 2 / yerr + np.log(yerr))

nll = lambda *args: -log_likelihood(*args);

initial = np.array([p.x_break.value, 1., p.alpha_1.value, p.alpha_2.value, p.amplitude.value]);
soln = minimize(nll, initial, args=(t, y, yerr),  method='Nelder-Mead', options={'disp':True, 'maxiter':10000});

t_bml, s_ml, beta_1ml, beta_2ml, F_cml = soln.x;


newy = broken_powerlaw(x0, t_bml, s_ml, beta_1ml, beta_2ml, F_cml);
# print(t_bml, s_ml, beta_1ml, beta_2ml, F_cml);


# initial = np.array([p.x_break.value, p.alpha_1.value, p.alpha_2.value, p.amplitude.value]);
# initial = np.array([5, 0.5, -7, 1, 9]);

# soln = minimize(nll, initial, args=(t, y, yerr),  method='Nelder-Mead', options={'disp':True});
#print(soln);
#t_bml, s_ml, beta_1ml, beta_2ml, F_cml = soln.x;
# t_bml, beta_1ml, beta_2ml, F_cml = soln.x;



#plot intial fits
# plt.style.use('light');
fig, ax = plt.subplots(figsize=(17,11), facecolor='white');

#generate some data
t_line = np.linspace(10,100,100);
# newy = broken_powerlaw(t_line, t_bml, beta_1ml, beta_2ml, F_cml);
# VLBA_plt = broken_powerlaw(VLBA['Day'], t_bml, beta_1ml, beta_2ml, F_cml);

newy = broken_powerlaw(t_line, t_bml, s_ml, beta_1ml, beta_2ml, F_cml);
VLBA_plt = broken_powerlaw(VLBA['Day'], s_ml, t_bml, beta_1ml, beta_2ml, F_cml);

#actual plots
ax.plot(VLBA['Day'], VLBA['Flux'], linestyle='dashed', marker='*', markersize=9, color='crimson');
ax.scatter(VLBA['Day'], VLBA_plt, label='Guesstimate, using ML', color='red', s=17, marker='x');
ax.plot(x0, p(x0), marker='.', label='astropy broken powerlaw fit');
ax.plot(t_line, newy, ":k", label="ML broken powerlaw fit");

#format the plots
# plt.rcParams.update({'lines.markersize': 9});
# ax.set_xscale('log');
# ax.set_yscale('log');
# ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'));
# ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'));
# ax.tick_params(which='both', direction='in', length=13, width=2);
# ax.grid(True, alpha=0.25);
# ax.set_ylim(0.1, 10);
# ax.set_xlim(10, 90);
# ax.set_title('V1674Her Cute Radio Light Curve + VLBA', fontsize=27, weight='bold');
# ax.set_ylabel('Flux Density (mJy)', fontsize=17, weight='bold');
# ax.set_xlabel('Day After Detection', fontsize=17, weight='bold');
# ax.legend();
# plt.show();
# plt.savefig(fname='./fits/initial_plots.png', format='png');


#now for the emcee stuff
#first up, the prior and log prob
# def log_prior(theta):
#     t_b, beta_1, beta_2, F_c = theta;
#     # print(theta);
#     if 10 < t_b < 21  and 0 < beta_2 < 10 and -3 < beta_1 < 0 and 2 < F_c < 8:
#         return 0.0
#     return -np.inf


def log_probability(theta, x, y, yerr):
    lp = log_prior(theta);
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr);


def log_prior(theta):
    t_b, s, beta_1, beta_2, F_c = theta;
    # print(theta);
    if 16 < t_b < 19  and 0.1 < s < 4 and -3 < beta_1 < 0 and 3 < beta_2 < 10 and 2 < F_c < 6:
        return 0.0
    return -np.inf

#walkers and sampler
pos = soln.x + 0.1*np.random.randn(128, 5);

nwalkers, ndim = pos.shape;

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(t, y, yerr));

sampler.run_mcmc(pos, int(1.5e6), progress=True);

#the walker plots
plt.style.use('default');
fig, axes = plt.subplots(5, figsize=(20, 7), sharex=True);

samples = sampler.get_chain();

# labels = ["t_b", "beta_1", "beta_2", "F_c"];

labels = ["t_b", "s", "beta_1", "beta_2", "F_c"];


for i in range(ndim):
    ax = axes[i]
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'));
    ax.plot(samples[:, :, i], "k", alpha=0.2)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number");

plt.savefig(fname='./fits/walker_plot.png', format='png');

#now for the corner plots
# tau = sampler.get_autocorr_time();
# discard_num = int(tau[0]*2.5);
# thin_num = int(tau[0]*0.5);

# flat_samples = sampler.get_chain(discard=discard_num, thin=thin_num, flat=True);
flat_samples = sampler.get_chain(discard=664, thin=300, flat=True);
fig = corner.corner(flat_samples, labels=labels);

plt.savefig(fname='./fits/corner_plots.png', format='png');

t_plot = np.linspace(8,500,10000);
fig.clear(True);


plt.errorbar(t,y,yerr,fmt = 'H');

for i in range(1000):
    theta = random.choice(flat_samples)
    t_b, s, beta_1, beta_2, F_c = theta
    model = broken_powerlaw(t_plot, t_b, s, beta_1, beta_2, F_c);
    plt.plot(t_plot,model,alpha = 0.1, color = 'black');


plt.loglog();

plt.savefig('./fits/sample_plots.png', format='png');   

from IPython.display import display, Math

for i in range(ndim):
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
    txt = txt.format(mcmc[1], q[0], q[1], labels[i])
    print(Math(txt))

print('beep boop, it is done');
