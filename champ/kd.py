import sys
import os
from scipy.optimize import minimize, curve_fit
import numpy as np
import matplotlib.pyplot as plt
import misc
import itertools
from champ import seqtools
from biofits import hyperbola
import lomp

BOOTSTRAP_ROUNDS = 20
MAX_BOOTSTRAP_SAMPLE_SIZE = 2000
MINIMUM_READ_COUNT = 5
TUKEY_CONSTANT = 1.5


class KdFitIA(object):
    """
    A class to fit Kd Values.

    Input:
        Intensity Array:    Data structure with all intensity values
        [max_clust=2000]:   max clusters to use per fit
    """

    def __init__(self, IA, max_clust=2000):
        self.IA = IA
        assert self.IA.course_trait_name == 'concentration_pM', self.IA.course_trait_name
        self.concentrations = self.IA.course_trait_list
        self.nM_concentrations = self.concentrations
        self.nM_concentrations = [conc / 1000.0 for conc in self.concentrations]
        self.target = self.IA.target
        self.neg_control_target = self.IA.neg_control_target
        self.max_clust = max_clust
        self.Imin_names = ['Imin_const']
        self.Imax_names = ['Imax_adjusted']

    def add_Imin_type(self, Imin_name):
        """
        Optionally add more Imin-determining strategies. Default: Imin_const.
        """
        assert Imin_name in ['Imin_const', 'Imin_neg_cont'], Imin_name
        if Imin_name not in self.Imin_names:
            self.Imin_names.append(Imin_name)

    def add_Imax_type(self, Imax_name):
        """
        Optionally add more Imax-determining strategies. Default: Imax_adjusted.
        """
        assert Imax_name in ['Imax_const', 'Imax_adjusted', 'Imax_ML'], Imax_name
        if Imax_name not in self.Imax_names:
            self.Imax_names.append(Imax_name)

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
        """
        Find Imax values according to selected strategies. If ML_seqs are included, finds Imax_ML.
        """

        def Iobs(x, Kd, Imax):
            return (Imax - self.Imin_const) / (1.0 + (float(Kd) / x)) + self.Imin_const

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
                self.Imax_adjusted.append(Imin + (1 + (perfect_Kd / conc)) * (med - Imin))
        self.Imax_given_conc = {
            conc: Imax for conc, Imax in zip(self.concentrations, self.Imax_adjusted)
            }

        if ML_seqs:
            self.fit_Imax_ML(ML_seqs)

    def model_logL(self,
                   seqs,
                   Kds,
                   Imin_list,
                   Imax_list,
                   sigma_consts,
                   sigI,
                   max_clust=None,
                   bootstrap_idxs=None):
        """
        Returns model logL probability. See documentation.
        """
        if bootstrap_idxs is not None:
            assert len(seqs) == 1, 'Bootstrapping can only be performed on one seq at a time.'
        if max_clust is None:
            max_clust = self.max_clust
        assert len(Imax_list) == len(self.concentrations), Imax_list
        Imax_arr = np.array(Imax_list)

        def theta(x, Kd):
            return 1.0 / (1.0 + (float(Kd) / x))

        thetas = np.empty((len(Kds), len(Imax_arr)))
        for i, Kd in enumerate(Kds):
            thetas[i] = [theta(conc, Kd) for conc in self.concentrations]
        sigma_clusts = (
            np.tile(self.Imin_stdev, (len(Kds), 1))
            + np.multiply(thetas, np.tile(sigma_consts, (len(Kds), 1)))
        )
        yhat = (
            np.tile(Imin_list, (len(Kds), 1))
            + np.multiply(thetas, np.tile(Imax_arr, (len(Kds), 1)))
        )
        logL = 0
        for sidx, seq in enumerate(seqs):
            loarr = self.IA.intensity_loarr_given_seq[seq]
            lol = self.IA.intensity_lol_given_seq[seq]
            for cidx, (inten_arr, inten_list) in enumerate(zip(loarr, lol)):
                if bootstrap_idxs is not None:
                    inten_arr = np.array([inten_list[idx] for idx in bootstrap_idxs
                                          if inten_list[idx] is not None])
                else:
                    inten_arr = inten_arr[:max_clust]
                logL += (
                    -len(inten_arr) * np.log(sigma_clusts[sidx, cidx])
                    - 1.0 / (2.0 * sigma_clusts[sidx, cidx] ** 2)
                    * np.square(inten_arr - yhat[sidx, cidx]).sum()
                )
        logL += (
            - len(Imax_arr) * np.log(sigI)
            - 1.0 / (2.0 * sigI ** 2) * np.square(Imax_arr - self.Imax_const).sum()
        )
        return logL

    def fit_Imax_ML_given_Imin(self, ML_seqs, Imin):
        idx1 = len(ML_seqs)
        idx2 = idx1 + self.IA.course_len
        idx3 = idx2 + self.IA.course_len
        Imin_list = misc.list_if_scalar(Imin, len(self.concentrations))

        def neg_log_L(params):
            params = map(abs, params)
            Kds = params[:idx1]
            Imax_list = params[idx1:idx2]
            sig_cs = params[idx2:idx3]
            sigI = params[idx3]
            return -self.model_logL(ML_seqs, Kds, Imin_list, Imax_list, sig_cs, sigI)

        x0 = list(
            [self.curve_fit_Kd(seq, self.Imin_const, self.Imax_const) for seq in ML_seqs]
            + self.Imax_adjusted
            + [1] * self.IA.course_len
            + [(self.Imax_const - self.Imin_const) / 10]
        )

        assert len(x0) == idx3 + 1
        res = minimize(neg_log_L, x0=x0, method='powell', options=dict(maxiter=1000000,
                                                                       maxfev=1000000,
                                                                       disp=True))
        if not res.success:
            print '\nWarning: Failure on {} ({})'.format(seq, seqtools.mm_names(self.target, seq))
        params = map(abs, res.x)
        Kds = params[:idx1]
        Imax_list = params[idx1:idx2]
        sig_cs = params[idx2:idx3]
        sigI = params[idx3]
        return Kds, Imax_list, sig_cs, sigI

    def fit_Imax_ML(self, ML_seqs):
        self.Imax_ML, self.sigma_consts, self.sigma_I = {}, {}, {}
        for Imin_name in self.Imin_names:
            Imin = getattr(self, Imin_name)
            Kds, Imax_list, sig_cs, sigI = self.fit_Imax_ML_given_Imin(ML_seqs, Imin)
            self.Imax_ML[Imin_name] = Imax_list
            self.sigma_consts[Imin_name] = sig_cs
            self.sigma_I[Imin_name] = sigI

    def ML_fit_Kd(self, seq, Imin_name, max_clust=None, bootstrap=False, *args, **kw_args):
        if max_clust is None:
            max_clust = self.max_clust
        if bootstrap:
            len_inten_list = len(self.IA.intensity_lol_given_seq[seq][0])
            bootstrap_idxs = np.random.choice(
                np.arange(len_inten_list),
                size=min(len_inten_list, max_clust),
                replace=True
            )
        else:
            bootstrap_idxs = None

        Imin = getattr(self, Imin_name)
        Imin_list = misc.list_if_scalar(Imin, len(self.concentrations))

        def neg_log_L(Kd):
            return -self.model_logL([seq],
                                    [Kd],
                                    Imin_list,
                                    self.Imax_ML[Imin_name],
                                    self.sigma_consts[Imin_name],
                                    self.sigma_I[Imin_name],
                                    bootstrap_idxs=bootstrap_idxs)

        res = minimize(neg_log_L, x0=20, method='powell', options=dict(maxiter=1000000,
                                                                       maxfev=1000000))
        if not res.success:
            print '\nWarning: Failure on {} ({})'.format(seq, seqtools.mm_names(self.target, seq))
        return float(res.x)

    def curve_fit_Kd(self, seq, Imin, Imax, max_clust=None, bootstrap=False, *args, **kw_args):
        """
        Least square curve fit to values normalized by Imin/Imax.
        """
        if max_clust is None:
            max_clust = self.max_clust

        def Iobs(x, Kd):
            return 1.0 / (1 + (float(Kd) / x))

        all_concentrations, all_intensities = self.IA.all_normalized_trait_and_inten_vals_given_seq(
            seq,
            Imin,
            Imax,
            max_clust=max_clust,
            bootstrap=bootstrap
        )
        popt, pcov = curve_fit(Iobs, all_concentrations, all_intensities, maxfev=100000)
        return popt[0]

    def setup_for_fit(self, force=False):
        if hasattr(self, 'fit_func_given_Imin_max_names') and not force:
            return
        self.fit_func_given_Imin_max_names = {}
        self.Imin_max_pairs_given_names = {}
        for Imin_name, Imax_name in itertools.product(self.Imin_names, self.Imax_names):
            Imin = getattr(self, Imin_name)
            Imax = getattr(self, Imax_name)
            self.fit_func_given_Imin_max_names[(Imin_name, Imax_name)] = self.curve_fit_Kd
            self.Imin_max_pairs_given_names[(Imin_name, Imax_name)] = (Imin, Imax)

        if 'Imax_ML' in self.Imax_names:
            for Imin_name in self.Imin_names:
                Imin = getattr(self, Imin_name)
                Imax = self.Imax_ML[Imin_name]
                Imax_name = 'Imax_ML'
                self.fit_func_given_Imin_max_names[(Imin_name, Imax_name)] = self.ML_fit_Kd
                self.Imin_max_pairs_given_names[(Imin_name, Imax_name)] = (Imin, Imax)

    def fit_all_Kds(self, num_bootstraps=20):
        self.setup_for_fit()

        def ABA(Kd, neg_cont_Kd):
            return np.log(neg_cont_Kd) - np.log(Kd)

        self.Kds = {name_tup: [] for name_tup in self.fit_func_given_Imin_max_names.keys()}
        self.Kd_errors = {name_tup: [] for name_tup in self.fit_func_given_Imin_max_names.keys()}
        self.ABAs = {name_tup: [] for name_tup in self.fit_func_given_Imin_max_names.keys()}
        self.ABA_errors = {name_tup: [] for name_tup in self.fit_func_given_Imin_max_names.keys()}
        dot_val = 100
        print '{} Seqs, \'.\'={}\n'.format(self.IA.nseqs, dot_val)
        for names_tup, (Imin, Imax) in sorted(self.Imin_max_pairs_given_names.items()):
            print '\n', names_tup
            Imin_name = names_tup[0]
            fit_func = self.fit_func_given_Imin_max_names[names_tup]
            neg_control_Kd = fit_func(self.neg_control_target,
                                      Imin=Imin,
                                      Imax=Imax,
                                      Imin_name=Imin_name)
            for i, seq in enumerate(self.IA.seqs):
                if i % dot_val == 0:
                    sys.stdout.write('.')
                    sys.stdout.flush()
                Kd = fit_func(seq, Imin=Imin, Imax=Imax, Imin_name=Imin_name)
                bs_Kds = [fit_func(seq, Imin=Imin, Imax=Imax, Imin_name=Imin_name, bootstrap=True)
                          for _ in range(num_bootstraps)]
                self.Kds[names_tup].append(Kd)
                self.Kd_errors[names_tup].append(np.std(bs_Kds))
                self.ABAs[names_tup].append(ABA(Kd, neg_control_Kd))
                bs_ABAs = [ABA(kkd, neg_control_Kd) for kkd in bs_Kds]
                self.ABA_errors[names_tup].append(np.std(bs_ABAs))

    def write_results(self, out_dir, bname):
        for names_tup, (Imin, Imax) in sorted(self.Imin_max_pairs_given_names.items()):
            Imin_name, Imax_name = names_tup
            Imin = misc.list_if_scalar(Imin, len(self.concentrations))
            Imax = misc.list_if_scalar(Imax, len(self.concentrations))

            out_fname = '{}_{}_{}_Kds_and_ABAs.txt'.format(bname,
                                                           Imin_name.replace(' ', '_'),
                                                           Imax_name.replace(' ', '_'))
            out_fpath = os.path.join(out_dir, out_fname)
            with open(out_fpath, 'w') as out:
                out.write('# Target: {}\n'.format(self.target))
                out.write('# Neg Control: {}\n'.format(self.neg_control_target))
                out.write('# Concentration\tImin\tImax\n')
                for conc, imin, imax in zip(self.concentrations, Imin, Imax):
                    out.write('\t'.join(map(str, map(float, (conc, imin, imax)))) + '\n')
                out.write('\t'.join(['# Seq', 'Kd (pM)', 'Kd error', 'ABA (kB T)', 'ABA error']) + '\n')
                out_zipper = zip(
                    self.IA.seqs,
                    self.Kds[names_tup],
                    self.Kd_errors[names_tup],
                    self.ABAs[names_tup],
                    self.ABA_errors[names_tup]
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
        nM_Kd = Kd / 1000.0

        ax.plot(self.nM_concentrations, Imin, 'ko', alpha=0.8, label='$I_{min}$')
        ax.plot(self.nM_concentrations, Imax, 's', color='darkgoldenrod', alpha=0.8, label='$I_{max}$')

        def Iobs(x, Kd, Imin, Imax):
            return (Imax - Imin) / (1.0 + (float(Kd) / x)) + Imin

        fit_path = [Iobs(conc, nM_Kd, imn, imx)
                    for conc, imn, imx in zip(self.nM_concentrations, Imin, Imax)]
        ax.plot(self.nM_concentrations, fit_path, 'k--')

        x = np.logspace(np.log10(self.nM_concentrations[0]),
                        np.log10(self.nM_concentrations[-1]),
                        200)
        y = [Iobs(xx, nM_Kd, self.Imin_const, self.Imax_const) for xx in x]
        ax.plot(x, y, 'r')
        ax.set_xscale('log')

        ax.set_xlabel('Concentration (nM)', fontsize=18)
        ax.set_axis_bgcolor('white')
        ax.grid(False)
        ax.set_ylabel('Intensity', fontsize=18)

    def plot_normalized_fit(self, ax, seq, Kd, Imin, Imax):
        self.IA.plot_normalized_intensities(ax, seq, Imin, Imax, xvals=self.nM_concentrations)
        nM_Kd = Kd / 1000.0

        def Iobs(x, Kd):
            return 1.0 / (1.0 + (float(Kd) / x))

        x = np.logspace(np.log10(self.nM_concentrations[0]),
                        np.log10(self.nM_concentrations[-1]),
                        200)
        y = [Iobs(xx, nM_Kd) for xx in x]
        ax.plot(x, y, 'r')
        ax.set_xscale('log')

        ax.set_xlabel('Concentration (nM)', fontsize=18)
        ax.set_axis_bgcolor('white')
        ax.grid(False)
        ax.set_ylabel('Intensity', fontsize=18)

    def example_plots(self, seqs, labels=None):
        self.setup_for_fit()
        if labels is None:
            labels = [None] * len(seqs)
        for seq, label in zip(seqs, labels):
            fig, ax = plt.subplots(figsize=(12, 0.8))
            ax.text(0, 0, label, fontsize=20, ha='center', va='center')
            ax.set_xlim((-1, 1))
            ax.set_ylim((-1, 1))
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)

            for names_tup, fit_func in self.fit_func_given_Imin_max_names.items():
                Imin_name = names_tup[0]
                fig, axes = plt.subplots(1, 2, figsize=(12, 5))
                Imin, Imax = self.Imin_max_pairs_given_names[names_tup]
                Kd = fit_func(seq, Imin=Imin, Imax=Imax, Imin_name=Imin_name)
                self.plot_raw_fit(axes[0], seq, Kd, Imin, Imax)
                self.plot_normalized_fit(axes[1], seq, Kd, Imin, Imax)
                axes[0].set_title('%s, Kd = %.2f' % (names_tup, Kd / 1000.0))
                axes[1].set_title(label)

    def all_error_analysis_and_figs(self, *args, **kw_args):
        self.setup_for_fit()
        for names_tup in self.Imin_max_pairs_given_names.keys():
            print names_tup
            sys.stdout.flush()
            self.error_analysis_and_figs(names_tup, *args, **kw_args)

    def error_analysis_and_figs(self,
                                Imin_max_names_tup,
                                seq=None,
                                num_bootstraps=100,
                                conf_pct=90,
                                min_reads=5,
                                out_dir=None,
                                out_bname=None):
        """
        For the given sequence, performs bootstrap analysis of errors for all numbers of clusters
        from 3 to 100.
        """
        if seq is None:
            seq = self.target

        fit_func = self.fit_func_given_Imin_max_names[Imin_max_names_tup]
        Imin, Imax = self.Imin_max_pairs_given_names[Imin_max_names_tup]
        Imin_name, Imax_name = Imin_max_names_tup

        ref_Kd = fit_func(seq, Imin=Imin, Imax=Imax, Imin_name=Imin_name)
        ref_dG = np.log(ref_Kd)

        read_names = self.IA.read_names_given_seq[seq]
        nclusters = range(3, min(100, len(read_names)))
        Kd_avg_errors, dG_avg_errors = [], []
        Kd_conf_errors, dG_conf_errors = [], []
        for n in nclusters:
            sys.stdout.write('.')
            sys.stdout.flush()
            bs_Kds = [fit_func(seq, Imin=Imin, Imax=Imax, Imin_name=Imin_name, max_clust=n, bootstrap=True)
                      for _ in range(num_bootstraps)]
            bs_Kd_errors = [abs(ref_Kd - Kd) / 1000.0 for Kd in bs_Kds]
            bs_dG_errors = [abs(ref_dG - np.log(Kd)) for Kd in bs_Kds]
            Kd_avg_errors.append(np.average(bs_Kd_errors))
            Kd_conf_errors.append(np.percentile(bs_Kd_errors, conf_pct))
            dG_avg_errors.append(np.average(bs_dG_errors))
            dG_conf_errors.append(np.percentile(bs_dG_errors, conf_pct))
        print

        def c_over_sqrt_n(n, c):
            return c / np.sqrt(n)

        def fit_c_over_sqrt_n(ns, data):
            new_ns = [n for n, dd in zip(ns, data) if np.isfinite(n) and np.isfinite(dd) and n > 10]
            new_data = [dd for n, dd in zip(ns, data) if np.isfinite(n) and np.isfinite(dd) and n > 10]
            popt, pcov = curve_fit(c_over_sqrt_n, new_ns, new_data, maxfev=10000)
            return popt[0]

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        for ax, label, units, avg_errors, conf_errors in zip(axes,
                                                             ('$K_d$', 'ABA'),
                                                             ('nM', '$k_B T$'),
                                                             (Kd_avg_errors, dG_avg_errors),
                                                             (Kd_conf_errors, dG_conf_errors)):
            fit_ns = np.linspace(1, max(nclusters), 300)
            c_avg = fit_c_over_sqrt_n(nclusters, avg_errors)
            c_conf = fit_c_over_sqrt_n(nclusters, conf_errors)
            avg_fit_vals = [c_over_sqrt_n(n, c_avg) for n in fit_ns]
            conf_fit_vals = [c_over_sqrt_n(n, c_conf) for n in fit_ns]

            min_reads_avg_fit = c_over_sqrt_n(min_reads, c_avg)

            ax.plot(nclusters, avg_errors, '.', label='Average Error')
            ax.plot(fit_ns, avg_fit_vals, label='Average Fit = $%.2f / \sqrt{n}$' % c_avg)
            ax.plot(nclusters, conf_errors, '.', label='90% Confidence Interval')
            ax.plot(fit_ns, conf_fit_vals, '--', label='90%% Conf Interval Fit = $%.2f / \sqrt{n}$' % c_conf)
            ax.plot([0, min_reads, min_reads], [min_reads_avg_fit, min_reads_avg_fit, 0], ':k')
            ax.set_xlim((0, 100))
            ax.set_xlabel('Number of clusters', fontsize=18)
            ax.set_ylabel('{} Error ({})'.format(label, units), fontsize=18)
            ax.legend(fontsize=14)
            ax.get_legend().get_frame().set_facecolor('white')
            ax.set_axis_bgcolor('white')
            ax.grid(False)

            for item in ax.get_xticklabels() + ax.get_yticklabels():
                item.set_fontsize(16)

        if out_dir:
            out_bname = '{}_{}_{}_error_analysis'.format(out_bname,
                                                         Imin_name.replace(' ', '_'),
                                                         Imax_name.replace(' ', '_'))
            fig.savefig(os.path.join(out_dir, out_bname + '.png'), dpi=300)
            fig.savefig(os.path.join(out_dir, out_bname + '.eps'))


class IAKdData(object):
    def __init__(self, Kd_fpath):
        self.concentrations, self.Imin, self.Imax, = [], [], []
        self.Kd, self.Kd_error, self.ABA, self.ABA_error = {}, {}, {}, {}
        with open(Kd_fpath) as f:
            line = next(f)
            assert line.startswith('# Target:')
            self.target = line.strip().split(': ')[1]
            line = next(f)
            assert line.startswith('# Neg Control')
            self.neg_control_target = line.strip().split(': ')[1]
            line = next(f)
            assert line.startswith('# Concentration\tImin\tImax')
            line = next(f)
            while not line.startswith('#'):
                conc, imn, imx = map(float, line.strip().split())
                self.concentrations.append(conc)
                self.Imin.append(imn)
                self.Imax.append(imx)
                line = next(f)
            assert line.startswith('# Seq')
            for line in f:
                if line.startswith('#'):
                    continue
                words = line.strip().split()
                seq = words[0]
                assert seq not in self.Kd, seq
                Kd, Kd_err, ABA, ABA_err = map(float, words[1:])
                self.Kd[seq] = Kd
                self.Kd_error[seq] = Kd_err
                self.ABA[seq] = ABA
                self.ABA_error[seq] = ABA_err
        self.neg_control_Kd = self.Kd[self.neg_control_target]
        self.log_neg_control_Kd = np.log(self.neg_control_Kd)
        self.target_ABA = self.ABA[self.target]

    def ABA_given_Kd(self, Kd):
        if Kd is None:
            return None
        return self.log_neg_control_Kd - np.log(Kd)


def fixed_delta_y_hyperbola(delta_y):
    """
    :param concentrations: array of titrant concentrations
    :param yint: Y-intercept
    :param delta_y: Total change in signal
    :param kd: Dissociation constant
    :return: array of Y values

    """
    return lambda (concentrations, yint, kd): yint + ((delta_y * concentrations) / (concentrations + kd))


def fit_hyperbola(concentrations, signals, delta_y=None):
    """
    :param concentrations: X-axis values representing concentrations in arbitrary units
    :param signals: Y-axis values representing some kind of signal. Don't normalize this. These are batched by
    concentration, and each batch can have variable numbers of observations

    :return:
        yint: the Y-intercept of the fit, often the background signal
        yint_stddev: standard deviation of the error of yint
        delta_y: total height of the fit
        delta_y_stddev: standard deviation of the error of delta_y
        kd: the dissociation constant
        kd_stddev: the standard deviation of the error of kd

    """
    if delta_y is None:
        (yint, delta_y, kd), covariance = curve_fit(hyperbola,
                                                    concentrations,
                                                    signals,
                                                    bounds=((-np.inf, 0.0, 10**-100),
                                                            (np.inf, np.inf, np.inf)))
    else:
        (yint, kd), covariance = curve_fit(fixed_delta_y_hyperbola(delta_y),
                                           concentrations,
                                           signals,
                                           bounds=((-np.inf, 10 ** -100), (np.inf, np.inf)))
    return yint, delta_y, kd


def fit_all_kds(read_name_intensities, h5_fpaths, process_count=36):
    all_concentrations = [misc.parse_concentration(fpath) for fpath in h5_fpaths]
    minimum_required_observations = max(len(all_concentrations) - 3, 5)
    for result in lomp.parallel_map(read_name_intensities.items(),
                                    _thread_fit_kd,
                                    args=(all_concentrations, minimum_required_observations),
                                    process_count=process_count):
        if result is not None:
            yield result


def _thread_fit_kd(read_name_intensities, all_concentrations, minimum_required_observations):
    read_name, intensities = read_name_intensities
    fitting_concentrations = []
    fitting_intensities = []
    for intensity, concentration in zip(intensities, all_concentrations):
        if intensity is np.nan:
            continue
        fitting_intensities.append(intensity)
        fitting_concentrations.append(concentration)
    if len(fitting_intensities) < minimum_required_observations:
        return None
    kd, yint, delta_y = fit_kd(fitting_concentrations, fitting_intensities)
    if kd is None:
        return None
    return read_name, kd, yint, delta_y


def fit_kd(all_concentrations, all_intensities):
    """ all_intensities is a list of dicts, with read_name: intensity"""
    try:
        yint, delta_y, kd = fit_hyperbola(all_concentrations, all_intensities)
    except (FloatingPointError, RuntimeError, Exception):
        return None, None, None
    else:
        return kd, yint, delta_y


def bootstrap_kd_uncertainty(all_concentrations, all_intensities):
    kds = []
    indexes = list(range(max([len(i) for i in all_intensities])))
    for i in range(BOOTSTRAP_ROUNDS):
        sample_of_indexes = np.random.choice(indexes, min(MAX_BOOTSTRAP_SAMPLE_SIZE, len(indexes)), replace=True)
        sample_of_intensities = []
        concentrations = []
        for n, concentration in enumerate(all_concentrations):
            concentration_subsample = []
            for index in sample_of_indexes:
                try:
                    intensity = all_intensities[n][index]
                except IndexError:
                    continue
                concentration_subsample.append(intensity)
            if concentration_subsample:
                sample_of_intensities.append(concentration_subsample)
                concentrations.append(concentration)
        try:
            _, _, kd = fit_hyperbola(concentrations, sample_of_intensities)
        except (FloatingPointError, RuntimeError, Exception) as e:
            continue
        else:
            kds.append(kd)
    if not kds:
        return None
    return np.std(kds)


# def build_intensity_concentration_array(sequence_intensities):
#     all_concentrations = []
#     all_intensities = []
#     for h5_filename, intensities in sorted(sequence_intensities.items(),
#                                            key=lambda hi: misc.parse_concentration(hi[0])):
#         concentration = misc.parse_concentration(h5_filename)
#         all_concentrations.append(concentration)
#         all_intensities.append(intensities)
#     return all_concentrations, all_intensities


def filter_reads_with_unusual_intensities(intensities):
    """
    Filters out intensity gradients where any of the measurements were absurdly high or low. Each intensity gradient
    is from a single cluster of DNA.

    :param intensities: a list of numpy arrays, potentially with np.nan values

    """
    bad_indexes = set()
    assert len(set([len(intensity) for intensity in intensities])) == 1, "All reads should have the same number of observations. Missing observations should be represented by np.nan"
    for index in range(len(intensities[0])):
        index_intensities = [intensity_gradient[index] for intensity_gradient in intensities if intensity_gradient[index] is not np.nan]
        if len(index_intensities) < MINIMUM_READ_COUNT:
            continue
        q1 = np.percentile(index_intensities, 25)
        q3 = np.percentile(index_intensities, 75)
        iqr = q3 - q1
        min_range, max_range = (q1 - TUKEY_CONSTANT * iqr, q3 + TUKEY_CONSTANT * iqr)
        for n, intensity_gradient in enumerate(intensities):
            if intensity_gradient[index] is not np.nan and (intensity_gradient[index] < min_range or intensity_gradient[index] > max_range):
                bad_indexes.add(n)
    return [ints for n, ints in enumerate(intensities) if n not in bad_indexes]


def assemble_read_intensities_for_fitting(intensities):
    """ Takes a list of read names that are known to be of high quality and assembles their intensities for each
    concentration into variable-sized lists. Some of these values may be np.nan and we need to exclude those (this is
    also why not all lists will be the same size). The goal here is to give something to a non-linear curve fitting
    function that doesn't need to be processed any further. """
    assembled_intensities = []
    for index in range(len(intensities[0])):
        concentration_intensities = []
        for intensity_gradient in intensities:
            intensity = intensity_gradient[index]
            if intensity is not np.nan:
                concentration_intensities.append(intensity)
        assembled_intensities.append(concentration_intensities)
    return assembled_intensities


def assemble_fitting_inputs(assembled_intensities, all_concentrations):
    """ We have to handle the case where at one concentration, no reads have any values. This will break our fitting
    code so we need to make sure we have a list of concentrations that only matches points with data. Typically the
    reason this happens is that higher concentration images just don't align at all. """
    intensities = []
    concentrations = []
    concentrations_per_observation = []
    for cluster_intensities, concentration in zip(assembled_intensities, all_concentrations):
        if cluster_intensities:
            intensities.append(cluster_intensities)
            concentrations.append(concentration)
            for _ in intensities:
                concentrations_per_observation.append(concentration)
    return concentrations, concentrations_per_observation, intensities


def main(interesting_read_names, h5_fpaths, int_scores, data_channel):
    skipped = 0
    sequence_kds = []
    for n, (sequence, read_names) in enumerate(interesting_read_names.items()):
            if len(read_names) < MINIMUM_READ_COUNT:
                skipped += 1
                continue
            scores = int_scores.score_given_read_name_in_channel
            sequence_intensities = [[scores[h5_fpath][data_channel].get(read_name, np.nan)
                                     for h5_fpath in h5_fpaths] for read_name in read_names]
            filtered_intensities = filter_reads_with_unusual_intensities(sequence_intensities)
            assembled_intensities = assemble_read_intensities_for_fitting(filtered_intensities)

            concentrations = [misc.parse_concentration(h5_fpath) for h5_fpath in h5_fpaths]
            concentrations, concentrations_per_observation, intensities = assemble_fitting_inputs(assembled_intensities, concentrations)

            kd, uncertainty, yint, delta_y = fit_kd(concentrations, intensities)
            uncertainty = bootstrap_kd_uncertainty(concentrations_per_observation, intensities)
            if kd is not None and uncertainty is not None:
                count = len(intensities)
                sequence_kds.append((sequence, kd, uncertainty, count))
