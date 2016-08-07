import sys
import os
from scipy.optimize import minimize, curve_fit
import numpy as np
import matplotlib.pyplot as plt
import random
import misc
import intensity_array


class KdFitIA(object):
    """
    A class to fit Kd Values using four schemes:

        Simple:         Const Imin/Imax from lowest neg control / perfect fit
        Imax_fixed:     Imin(x) from neg control, Imax from perfect fit
        Imax_adjusted:  Imin(x) from neg control, Imax(x) > thresh adjusted from fit
        Imax_ML:        Imin(x) from neg control, Imax(x) from ML of model

    Input:
        Intensity Array
    """
    def __init__(self, IA, max_clust=2000):
        self.IA = IA
        assert self.IA.course_trait_name == 'concentration_pM', self.IA.course_trait_name
        self.concentrations = self.IA.course_trait_list
        self.nM_concentrations = [conc/1000.0 for conc in self.concentrations]
        self.target = self.IA.target
        self.neg_control_target = self.IA.neg_control_target
        self.max_clust = max_clust

    def find_Imin_and_background_noise(self):
        """
        Find Imin and Imin_stdev, the mode and stdev of the negative control intensities at each
        concentration.
        """
        self.Imin_neg_cont = self.IA.modes_given_seq(self.neg_control_target)
        self.Imin_const = self.Imin_neg_cont[0]
        self.Imin_given_conc = {
            conc: Imin for conc, Imin in zip(self.concentrations, self.Imin_neg_cont)
        }
        self.Imin_stdev = self.IA.stdevs_given_seq(self.neg_control_target)
        self.Imin_stdev_given_conc = {
            conc: Imin_stdev for conc, Imin_stdev in zip(self.concentrations, self.Imin_stdev)
        }

    def find_Imax(self, ML_seqs=None):
        def Iobs(x, Kd, Imax):
            return (Imax - self.Imin_const)/(1.0 + (float(Kd)/x)) + self.Imin_const

        all_concentrations, all_intensities = self.IA.all_trait_and_inten_vals_given_seq(
            self.target,
            max_clust=self.max_clust
        )

        popt, pcov = curve_fit(Iobs, all_concentrations, all_intensities)
        perfect_Kd, self.Imax_const = popt

        perfect_medians = self.IA.medians_given_seq(self.target)
        self.Imax_adjusted = []
        for conc, Imin, med in zip(self.concentrations, self.Imin_neg_cont, perfect_medians):
            fit_Iobs = Iobs(conc, perfect_Kd, self.Imax_const)
            if fit_Iobs < 0.90 * self.Imax_const:
                self.Imax_adjusted.append(self.Imax_const)
            else:
                # Using the median value as real Iobs, solve for Imax_adjusted
                self.Imax_adjusted.append(Imin + (1 + (perfect_Kd/conc)) * (med - Imin))
        self.Imax_given_conc = {
            conc: Imax for conc, Imax in zip(self.concentrations, self.Imax_adjusted)
        }

        if ML_seqs:
            self.fit_Imax_ML(ML_seqs)

    def model_logL(self, seqs, Kds, Imax_list, sigma_consts, sigI, max_clust=None, bootstrap=False):
        """
        Returns model logL probability. See documentation.
        """
        if max_clust is None:
            max_clust = self.max_clust
        assert len(Imax_list) == len(self.concentrations), Imax_list
        Imax_arr = np.array(Imax_list)
        def theta(x, Kd):
            return 1.0/(1.0 + (float(Kd)/x))
        thetas = np.empty((len(Kds), len(Imax_arr)))
        for i, Kd in enumerate(Kds):
            thetas[i] = [theta(conc, Kd) for conc in self.concentrations]
        sigma_clusts = (
            np.tile(self.Imin_stdev, (len(Kds), 1))
            + np.multiply(thetas, np.tile(sigma_consts, (len(Kds), 1)))
        )
        yhat = (
            np.tile(self.Imin_neg_cont, (len(Kds), 1))
            + np.multiply(thetas, np.tile(Imax_arr, (len(Kds), 1)))
        )
        logL = 0
        for sidx, seq in enumerate(seqs):
            loarr = self.IA.intensity_loarr_given_seq[seq]
            for cidx, inten_arr in enumerate(loarr):
                if bootstrap:
                    inten_arr = np.random.choice(
                        inten_arr,
                        size=min(max_clust, len(inten_arr)),
                        replace=True
                    )
                logL += (
                    -len(inten_arr) * np.log(sigma_clusts[sidx, cidx])
                    - 1.0/(2.0 * sigma_clusts[sidx, cidx]**2) 
                    * np.square(inten_arr - yhat[sidx, cidx]).sum()
                )
        logL += (
            - len(Imax_arr) * np.log(sigI)
            - 1.0/(2.0 * sigI**2) * np.square(Imax_arr - self.Imax_const).sum()
        )
        return logL
        
    def fit_Imax_ML(self, ML_seqs):
        idx1 = len(ML_seqs)
        idx2 = idx1 + self.IA.course_len
        idx3 = idx2 + self.IA.course_len
        def neg_log_L(params):
            params = map(abs, params)
            Kds = params[:idx1]
            Imax_list = params[idx1:idx2]
            sig_cs = params[idx2:idx3]
            sigI = (self.Imax_const - self.Imin_const)/100.0 #params[idx3]
            return -self.model_logL(ML_seqs, Kds, Imax_list, sig_cs, sigI)

        x0 = list(
            [self.simple_fit_Kd(seq) for seq in ML_seqs]
            + self.Imax_adjusted
            + [1]*self.IA.course_len
            + [(self.Imax_const - self.Imin_const)/10]
        )

        assert len(x0) == idx3 + 1
        res = minimize(neg_log_L, x0=x0, method='powell', options=dict(maxiter=1000000, disp=True))
        assert res.success, res
        self.Imax_ML = map(abs, res.x[idx1:idx2])
        self.sigma_consts = map(abs, res.x[idx2:idx3])
        self.sigma_I = abs(res.x[idx3])

    def ML_fit_Kd(self, seq, *args, **kw_args):
        def neg_log_L(Kd):
            return -self.model_logL([seq], [Kd], self.Imax_ML, self.sigma_consts, self.sigma_I)
        res = minimize(neg_log_L, x0=20, method='powell', options=dict(maxiter=1000000))
        assert res.success, res
        return float(res.x)

    def simple_fit_Kd(self, seq, max_clust=None, bootstrap=False, *args, **kwargs):
        if max_clust is None:
            max_clust = self.max_clust
        def Iobs(x, Kd):
            return (self.Imax_const - self.Imin_const)/(1 + (float(Kd)/x)) + self.Imin_const
        all_concentrations, all_intensities = self.IA.all_trait_and_inten_vals_given_seq(
            seq,
            max_clust=max_clust,
            bootstrap=bootstrap
        )
        popt, pcov = curve_fit(Iobs, all_concentrations, all_intensities)
        return popt[0]

    def Imax_of_choice_curve_fit_Kd(self, seq, Imax, max_clust=None, bootstrap=False):
        if max_clust is None:
            max_clust = self.max_clust
        def Iobs(x, Kd):
            return 1.0/(1 + (float(Kd)/x))
        all_concentrations, all_intensities = self.IA.all_normalized_trait_and_inten_vals_given_seq(
            seq,
            self.Imin_neg_cont,
            Imax,
            max_clust=max_clust,
            bootstrap=bootstrap
        )
        popt, pcov = curve_fit(Iobs, all_concentrations, all_intensities)
        return popt[0]

    def fit_all_Kds(self, num_bootstraps=50):
        self.fit_types = ['Simple', 'Imax fixed', 'Imax adjusted', 'Imax ML']
        self.fit_funcs = [
            self.simple_fit_Kd,
            self.Imax_of_choice_curve_fit_Kd,
            self.Imax_of_choice_curve_fit_Kd,
            self.ML_fit_Kd
        ]
        self.fit_Imins = [
            self.Imin_const,
            self.Imin_neg_cont,
            self.Imin_neg_cont,
            self.Imin_neg_cont,
        ]
        self.fit_Imaxs = [
            self.Imax_const,
            self.Imax_const,
            self.Imax_adjusted,
            self.Imax_ML
        ]

        self.neg_control_Kds = {}
        print 'Negative controls'
        for fit_type, fit_func, fit_Imax in zip(self.fit_types, self.fit_funcs, self.fit_Imaxs):
            self.neg_control_Kds[fit_type] = fit_func(self.neg_control_target, Imax=fit_Imax)

        def ABA(Kd, neg_cont_Kd):
            return np.log(neg_cont_Kd) - np.log(Kd)

        self.Kds = {name: [] for name in self.fit_types}
        self.Kd_errors = {name: [] for name in self.fit_types}
        self.ABAs = {name: [] for name in self.fit_types}
        self.ABA_errors = {name: [] for name in self.fit_types}
        dot_val = 10
        print '{} Seqs, \'.\'={}\n'.format(self.IA.nseqs, dot_val)
        for fit_type, fit_func, fit_Imax in zip(self.fit_types, self.fit_funcs, self.fit_Imaxs):
            print '\n', fit_type
            for i, seq in enumerate(self.IA.seqs):
                if i % dot_val == 0:
                    sys.stdout.write('.')
                    sys.stdout.flush()
                Kd = fit_func(seq, Imax=fit_Imax)
                bs_Kds = [fit_func(seq, Imax=fit_Imax, bootstrap=True) for _ in range(num_bootstraps)]
                self.Kds[fit_type].append(Kd)
                self.Kd_errors[fit_type].append(np.std(bs_Kds))
                self.ABAs[fit_type].append(ABA(Kd, self.neg_control_Kds[fit_type]))
                bs_ABAs = [ABA(kkd, self.neg_control_Kds[fit_type]) for kkd in bs_Kds]
                self.ABA_errors[fit_type].append(np.std(bs_ABAs))

    def write_results(self, out_dir, bname):
        for fit_type, fit_func, fit_Imin, fit_Imax in zip(self.fit_types,
                                                          self.fit_funcs,
                                                          self.fit_Imins,
                                                          self.fit_Imaxs):
            fit_type_fname = fit_type.replace(' ', '_')
            fit_Imin = misc.list_if_scalar(fit_Imin, len(self.concentrations))
            fit_Imax = misc.list_if_scalar(fit_Imax, len(self.concentrations))

            out_fpath = os.path.join(out_dir, '{}_{}_Kds_and_ABAs.txt'.format(bname, fit_type_fname))
            with open(out_fpath, 'w') as out:
                out.write('# Target: {}\n'.format(self.target))
                out.write('# Neg Control: {}\n'.format(self.neg_control_target))
                out.write('# Concentration\tImin\tImax\n')
                for conc, Imin, Imax in zip(self.concentrations, fit_Imin, fit_Imax):
                    out.write('\t'.join(map(str, map(float, (conc, Imin, Imax)))) + '\n')
                out.write('\t'.join(['# Seq', 'Kd (pM)', 'Kd error', 'ABA (kB T)', 'ABA error']) + '\n')
                out_zipper = zip(
                    self.IA.seqs,
                    self.Kds[fit_type],
                    self.Kd_errors[fit_type],
                    self.ABAs[fit_type],
                    self.ABA_errors[fit_type]
                )
                out.write(
                    '\n'.join(
                        '\t'.join(map(str, [seq, Kd, Kd_err, ABA, ABA_err]))
                        for seq, Kd, Kd_err, ABA, ABA_err in out_zipper
                    )
                )

    def plot_raw_fit(self, ax, seq, Kd, Imin, Imax):
        self.IA.plot_raw_intensities(ax, seq, xvals=self.nM_concentrations)
        Imin = misc.list_if_scalar(Imin, self.IA.course_len)
        Imax = misc.list_if_scalar(Imax, self.IA.course_len)
        nM_Kd = Kd/1000.0

        ax.plot(self.nM_concentrations, Imin, 'ko', alpha=0.8, label='$I_{min}$')
        ax.plot(self.nM_concentrations, Imax, 's', color='darkgoldenrod', alpha=0.8, label='$I_{max}$')

        def Iobs(x, Kd, Imin, Imax):
            return (Imax - Imin) / (1.0 + (float(Kd)/x)) + Imin

        fit_path = [Iobs(conc, nM_Kd, imn, imx)
                    for conc, imn, imx in zip(self.nM_concentrations, Imin, Imax)]
        ax.plot(self.nM_concentrations, fit_path, 'k--')

        x = np.linspace(self.nM_concentrations[0], self.nM_concentrations[-1], 200)
        y = [Iobs(xx, nM_Kd, self.Imin_const, self.Imax_const) for xx in x]
        ax.plot(x, y, 'r')
        ax.set_xscale('log')

        ax.set_xlabel('Concentration (nM)', fontsize=18)
        ax.set_axis_bgcolor('white')
        ax.grid(False)
        ax.set_ylabel('Intensity', fontsize=18)

    def plot_normalized_fit(self, ax, seq, Kd, Imin, Imax):
        self.IA.plot_normalized_intensities(ax, seq, Imin, Imax, xvals=self.nM_concentrations)
        Imin = misc.list_if_scalar(Imin, self.IA.course_len)
        Imax = misc.list_if_scalar(Imax, self.IA.course_len)
        nM_Kd = Kd/1000.0

        def Iobs(x, Kd):
            return 1.0 / (1.0 + (float(Kd)/x))

        x = np.linspace(self.nM_concentrations[0], self.nM_concentrations[-1], 200)
        y = [Iobs(xx, nM_Kd) for xx in x]
        ax.plot(x, y, 'r')
        ax.set_xscale('log')

        ax.set_xlabel('Concentration (nM)', fontsize=18)
        ax.set_axis_bgcolor('white')
        ax.grid(False)
        ax.set_ylabel('Intensity', fontsize=18)

    def example_plots(self, seqs, labels=None):
        if labels is None:
            labels = [None]*len(seqs)
        for seq, label in zip(seqs, labels):
            fig, ax = plt.subplots(figsize=(12, 0.8))
            ax.text(0, 0, label, fontsize=20, ha='center', va='center')
            ax.set_xlim((-1, 1))
            ax.set_ylim((-1, 1))
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)

            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            Kd = self.simple_fit_Kd(seq)
            self.plot_raw_fit(axes[0], seq, Kd, self.Imin_const, self.Imax_const)
            self.plot_normalized_fit(axes[1], seq, Kd, self.Imin_const, self.Imax_const)
            axes[0].set_title('Simple Fit, Kd = %.2f' % (Kd/1000.0))
            axes[1].set_title(label)

            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            Kd = self.Imax_of_choice_curve_fit_Kd(seq, self.Imax_const)
            self.plot_raw_fit(axes[0], seq, Kd, self.Imin_neg_cont, self.Imax_const)
            self.plot_normalized_fit(axes[1], seq, Kd, self.Imin_neg_cont, self.Imax_const)
            axes[0].set_title('$I_{max}$ Fixed, Kd = %.2f' % (Kd/1000.0))
            axes[1].set_title(label)

            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            Kd = self.Imax_of_choice_curve_fit_Kd(seq, self.Imax_adjusted)
            self.plot_raw_fit(axes[0], seq, Kd, self.Imin_neg_cont, self.Imax_adjusted)
            self.plot_normalized_fit(axes[1], seq, Kd, self.Imin_neg_cont, self.Imax_adjusted)
            axes[0].set_title('$I_{max}$ Adjusted, Kd = %.2f' % (Kd/1000.0))
            axes[1].set_title(label)

            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            Kd = self.ML_fit_Kd(seq)
            self.plot_raw_fit(axes[0], seq, Kd, self.Imin_neg_cont, self.Imax_ML)
            self.plot_normalized_fit(axes[1], seq, Kd, self.Imin_neg_cont, self.Imax_ML)
            axes[0].set_title('$I_{max}$ ML, Kd = %.2f' % (Kd/1000.0))
            axes[1].set_title(label)

    def error_analysis_and_figs(self, seq, fit_func, conf_pct=90):
        """
        For the given sequence, performs bootstrap analysis of errors for all numbers of clusters
        from 3 to 100.
        """
        ref_Kd = ff(self.read_names_given_seq[seq])
        ref_dG = np.log(ref_Kd)
        
        read_names = self.read_names_given_seq[seq]
        nclusters = range(3, min(100, len(read_names)))
        Kd_avg_errors, dG_avg_errors = [], []
        Kd_conf_errors, dG_conf_errors = [], []
        for n in nclusters:
            bs_Kds = self.bootstrap_Kds_given_read_names(read_names, 
                                                         sample_size=n,
                                                         num_bootstraps=50,
                                                         fit_func=fit_func)
            bs_Kd_errors = [abs(ref_Kd - Kd)/1000.0 for Kd in bs_Kds]
            bs_dG_errors = [abs(ref_dG - np.log(Kd)) for Kd in bs_Kds]
            Kd_avg_errors.append(np.average(bs_Kd_errors))
            Kd_conf_errors.append(np.percentile(bs_Kd_errors, conf_pct))
            dG_avg_errors.append(np.average(bs_dG_errors))
            dG_conf_errors.append(np.percentile(bs_dG_errors, conf_pct))

        def c_over_sqrt_n(n, c):
            return c / np.sqrt(n)

        def fit_c_over_sqrt_n(ns, data):
            popt, pcov = curve_fit(c_over_sqrt_n, ns, data, maxfev=10000)
            return popt[0]

        for ax, label, avg_errors, conf_errors in zip(axes,
                                                      ('$K_d$', '$\Delta G$'),
                                                      (Kd_avg_errors, dG_avg_errors),
                                                      (Kd_conf_errors, dG_conf_errors)):
            fit_ns = np.linspace(1, max(nclusters), 300)
            c_avg = fit_c_over_sqrt_n(nclusters, avg_errors)
            c_conf = fit_c_over_sqrt_n(nclusters, conf_errors)
            avg_fit_vals = [c_over_sqrt_n(n, c_avg) for n in fit_ns]
            conf_fit_vals = [c_over_sqrt_n(n, c_conf) for n in fit_ns]

            min_reads_avg_fit = c_over_sqrt_n(self.min_reads, c_avg)

            ax.plot(nclusters, avg_errors, '.', label='Average Error')
            ax.plot(fit_ns, avg_fit_vals, label='Average Fit = $%.2f / \sqrt{n}$' % c_avg)
            ax.plot(nclusters, conf_errors, '.', label='90% Confidence Interval')
            ax.plot(fit_ns, conf_fit_vals, '--', label='90%% Conf Interval Fit = $%.2f / \sqrt{n}$' % c_conf)
            ax.plot([0, self.min_reads, self.min_reads], [min_reads_avg_fit, min_reads_avg_fit, 0], ':k')
            ax.set_xlim((0, 100))
            ax.set_xlabel('Number of clusters', fontsize=18)
            ax.set_ylabel('{} Error'.format(label), fontsize=18)
            ax.legend(fontsize=14)
            ax.get_legend().get_frame().set_facecolor('white')
            ax.set_axis_bgcolor('white')
            ax.grid(False)
                                                    
            for item in ax.get_xticklabels() + ax.get_yticklabels():
                item.set_fontsize(16)
