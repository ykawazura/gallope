# -*- coding: utf-8 -*-
from load import *
from fft import *
import sys
sys.path.append('../')
from plots import *

print('\nplotting energy\n')
outdir = './fig_energy/'
import os
if not os.path.isdir(outdir):
  os.makedirs(outdir)

#--------------------------------------------------------#
#          Alfvenic (EM) power balance & energy          #
#--------------------------------------------------------#
# plot energy balance (dW/dt = P - D, so dW/dt + D - P should vanish)
ys = [
       upe2dot_sum + bpe2dot_sum ,
       upe2dissip_sum + bpe2dissip_sum,
       -p_phi_sum,
       -p_psi_sum,
       -p_phi_sum - p_psi_sum,
       upe2dot_sum + bpe2dot_sum + upe2dissip_sum + bpe2dissip_sum - p_phi_sum - p_psi_sum,
     ]
xs = [
       tt,
       tt,
       tt,
       tt,
       tt,
       tt,
     ]
ls = [
        '',
        '',
        '',
        '',
        '',
        'k--',
     ]
legends = [
       r'$\rmd W/\rmd t$',
       r'$D_\mr{AW}$',
       r'$-P_{\Phi}$',
       r'$-P_{\Psi}$',
       r'$-P_\mr{total}$',
       r'balance',
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_all.pdf')


# plot energy change (EM + compressive free energy)
ys = [
       upe2_sum,
       bpe2_sum,
       W_free,
     ]
ls = [
        '',
        '',
        '',
     ]
xs = [
       tt,
       tt,
       tt,
     ]
legends = [
       r'$W_{u_\+}$',
       r'$W_{\delta B_\+}$',
       r'$W_\mr{free}$',
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'energy_all.pdf')

# plot helicity change
ys = [
       (zppe2_sum - zmpe2_sum)/(zppe2_sum + zmpe2_sum),
       (upe2_sum - bpe2_sum)/(upe2_sum + bpe2_sum),
     ]
ls = [
       '', ''
     ]
xs = [
       tt, tt
     ]
legends = [
       r'$H := \f{\int\rmd^3\bm{x}[(Z_\+^+)^2 - (Z_\+^-)^2]}{\int\rmd^3\bm{x}[(Z_\+^+)^2 + (Z_\+^-)^2]}$',
       r'$R := \f{\int\rmd^3\bm{x}[(u_\+)^2 - (b_\+^-)^2]}{\int\rmd^3\bm{x}[(u_\+)^2 + (b_\+)^2]}$',
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'helicity_resiene.pdf')


ys = [
       p_phi_sum + p_psi_sum,
       p_xhl_sum,
       p_phi_sum - p_psi_sum,
     ]
ls = [
       '', '', '',
     ]
xs = [
       tt, tt, tt,
     ]
legends = [
       r'$\int \rmd^3\bm{x}(-\zeta_\+^+\nbl^2 f_{\zeta_\+^+} - \zeta_\+^-\nbl^2 f_{\zeta_\+^-})$',
       r'$\int \rmd^3\bm{x}(-\zeta_\+^+\nbl^2 f_{\zeta_\+^+} + \zeta_\+^-\nbl^2 f_{\zeta_\+^-})$',
       r'$\int \rmd^3\bm{x}(-\Phi\nbl^2 f_{\Phi} + \Psi\nbl^2 f_{\Psi})$',
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, ymin=0.0, save=outdir + 'helicity_resiene_injection.pdf')


#--------------------------------------------------------#
#          compressive (g) free energy power balance      #
#--------------------------------------------------------#
# g has its OWN forcing (g_1 injection, force.F90) and dissipation
# (mu*kperp^(2*nexp_perp) + nu*m^nexp_m), independent of the EM fields, so
#   dW_free/dt = P_g - D_g          (Meyrand et al. 2019, free energy eq 6).
# P_g (p_g_sum) and D_g (Dg_sum) are now output by the Fortran diagnostics
# (diagnostics.F90:get_g_power_balance). dW_free/dt is formed here as a
# finite difference of the stored W_free(t) time series; hence the balance
#   dW_free/dt + D_g - P_g
# should vanish (the coarser the diagnostic cadence, the noisier the FD term).
Wfree_dot = np.gradient(W_free, tt)

# plot g power balance (dW_free/dt = P_g - D_g, so dW_free/dt + D_g - P_g -> 0)
ys = [
       Wfree_dot,
       Dg_sum,
       -p_g_sum,
       Wfree_dot + Dg_sum - p_g_sum,
     ]
xs = [ tt, tt, tt, tt ]
ls = [ '', '', '', 'k--' ]
legends = [
       r'$\rmd W_\mr{free}/\rmd t$',
       r'$D_g$',
       r'$-P_g$',
       r'balance',
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'balance_g_all.pdf')

# plot injection vs dissipation (steady state: <P_g> ~ <D_g>)
ys = [ p_g_sum, Dg_sum ]
xs = [ tt, tt ]
ls = [ '', '' ]
legends = [ r'$P_g$', r'$D_g$' ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'Pg_Dg.pdf')

# cumulative balance (resolution-independent): int_0^t (P_g - D_g) dt' must
# overlay W_free(t) - W_free(0) exactly (energy conservation, trapezoid rule).
net_g    = p_g_sum - Dg_sum
cum_netg = np.concatenate(([0.0], np.cumsum(0.5*(net_g[1:] + net_g[:-1])*np.diff(tt))))
ys = [ W_free - W_free[0], cum_netg ]
xs = [ tt, tt ]
ls = [ '', 'k--' ]
legends = [
       r'$W_\mr{free}(t) - W_\mr{free}(0)$',
       r'$\int_0^t (P_g - D_g)\,\rmd t^\prime$',
     ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'Wfree_cumulative.pdf')

# W_free(t)
ys = [ W_free ]
xs = [ tt ]
ls = [ '' ]
legends = [ r'$W_\mr{free}$' ]
plot_1d_many(xs, ys, xlab=tlab, legends=legends, ls=ls, legendloc='upper left', title='', ylab='', term=True, save=outdir + 'Wfree.pdf')


# ascii output
np.savetxt(outdir + 'energies.txt' , np.column_stack((tt, upe2_sum, bpe2_sum,
                                                          upe2dot_sum, bpe2dot_sum,
                                                          upe2dissip_sum, bpe2dissip_sum,
                                                          p_phi_sum, p_psi_sum,
                                                          zppe2_sum, zmpe2_sum,
                                                          W_free, Wfree_dot,
                                                          p_g_sum, Dg_sum,
                                                        )), fmt='%E')
np.savetxt(outdir + 'balance.txt' , np.column_stack((tt, upe2dot_sum + bpe2dot_sum + upe2dissip_sum + bpe2dissip_sum - p_phi_sum - p_psi_sum
                                                      )), fmt='%E')
np.savetxt(outdir + 'energy_dot.txt' , np.column_stack((tt, upe2dot_sum + bpe2dot_sum
                                                      )), fmt='%E')



# calculate balance
balance = upe2dot_sum + bpe2dot_sum + upe2dissip_sum + bpe2dissip_sum - p_phi_sum - p_psi_sum

print ('|balance| > 1e1 at')
print (np.where(abs(balance) > 1e1)[0])
