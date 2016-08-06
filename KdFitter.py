import sys
from scipy.optimize import minimize, curve_fit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import random
import misc


class KdFitter(object):
    """
    A class to fit Kd values.

        Params:
            :IntensityScores:   int_scores
            :list:              h5_fpaths
            :str:               signal_channel
            :dict:              read_names_given_seq 
            :str:               target
            :set:               perfect_target_read_names
            :str:               neg_control_target
            :set:               neg_contol_target_read_names
            :int:               min_reads
    """
    def __init__(self, 
                 int_scores,
                 h5_fpaths,
                 signal_channel,
                 read_names_given_seq, 
                 target,
                 perfect_target_read_names,
                 neg_control_target,
                 neg_contol_target_read_names, 
                 min_reads,
                ):
        self.int_scores = int_scores
        self.signal_channel = signal_channel
        self.read_names_given_seq = read_names_given_seq
        self.target = target
        self.perfect_target_read_names = perfect_target_read_names
        self.neg_control_target = neg_control_target
        self.neg_contol_target_read_names = neg_contol_target_read_names
        self.min_reads = min_reads

        self.h5_fpaths = h5_fpaths
        self.h5_fpaths.sort(key=misc.pM_concentration_given_fpath)
        self.concentrations = np.array(map(misc.pM_concentration_given_fpath, self.h5_fpaths))
        self.nM_concentrations = np.array([conc/1000.0 for conc in self.concentrations])
        self.max_sample_size = 1000

    def find_medians_given_read_names(self, read_names):
        """
        Return the median value of given reads at each concentration.
        """
        medians = []
        for h5_fpath in self.h5_fpaths:
            # Already verified this exists in __init__
            score_given_read_name = self.int_scores.score_given_read_name_in_channel[h5_fpath][self.signal_channel]
            medians.append(
                np.median(
                    [float(score_given_read_name[read_name])
                     for read_name in read_names
                     if read_name in score_given_read_name]
                )
            )
        return medians

    def find_modes_given_read_names(self, read_names):
        """
        Return the mode value of given reads at each concentration.
        """
        modes = []
        for h5_fpath in self.h5_fpaths:
            # Already verified this exists in __init__
            score_given_read_name = self.int_scores.score_given_read_name_in_channel[h5_fpath][self.signal_channel]
            modes.append(
                misc.get_mode(
                    [float(score_given_read_name[read_name])
                     for read_name in read_names
                     if read_name in score_given_read_name]
                )
            )
        return modes

    def find_stdevs_given_read_names(self, read_names):
        """
        Return the mode value of given reads at each concentration.
        """
        stdevs = []
        for h5_fpath in self.h5_fpaths:
            # Already verified this exists in __init__
            score_given_read_name = self.int_scores.score_given_read_name_in_channel[h5_fpath][self.signal_channel]
            stdevs.append(
                np.std(
                    [float(score_given_read_name[read_name])
                     for read_name in read_names
                     if read_name in score_given_read_name]
                )
            )
        return stdevs

    def get_all_conc_and_inten_given_read_names(self, read_names):
        """
        Returns all concentration/intensities pairs in read_names, return as two lists.

        Returns:
            :list: all_concentrations
            :list: all_intensities
        """
        all_concentrations, all_intensities = [], []
        for conc, h5_fpath in zip(self.concentrations, self.h5_fpaths):
            score_given_read_name = self.int_scores.score_given_read_name_in_channel[h5_fpath][self.signal_channel]
            for read_name in read_names:
                if read_name in score_given_read_name:
                    all_concentrations.append(float(conc))
                    all_intensities.append(float(score_given_read_name[read_name]))
        return all_concentrations, all_intensities

    def get_adjusted_conc_and_inten_given_read_names(self, read_names):
        """
        Returns all concentration/intensities pairs in read_names, return as two lists, with the
        intensities adjusted by Imin and Imax to run typically between 0 and 1.

        Returns:
            :list: all_concentrations
            :list: all_intensities
        """
        all_concentrations, all_intensities = [], []
        for conc, h5_fpath, Imin, Imax in zip(self.concentrations,
                                              self.h5_fpaths,
                                              self.Imin,
                                              self.Imax):
            score_given_read_name = self.int_scores.score_given_read_name_in_channel[h5_fpath][self.signal_channel]
            for read_name in read_names:
                if read_name in score_given_read_name:
                    all_concentrations.append(float(conc))
                    all_intensities.append(float(float(score_given_read_name[read_name])-Imin)/float(Imax-Imin))
        return all_concentrations, all_intensities

    def get_adjusted_median_intensities_and_mads(self, read_names):
        """
        Returns all concentration/intensities pairs in read_names, return as two lists, with the
        intensities adjusted by Imin and Imax to run typically between 0 and 1.

        Returns:
            :list: all_concentrations
            :list: all_intensities
        """
        medians, mads = [], []
        for conc, h5_fpath, Imin, Imax in zip(self.concentrations,
                                              self.h5_fpaths,
                                              self.Imin,
                                              self.Imax):
            score_given_read_name = self.int_scores.score_given_read_name_in_channel[h5_fpath][self.signal_channel]
            intensities = []
            for read_name in read_names:
                if read_name in score_given_read_name:
                    intensities.append((float(score_given_read_name[read_name])-Imin)/float(Imax-Imin))
            med, mad = misc.median_and_median_absolute_deviation(intensities)
            medians.append(med)
            mads.append(mad)
        return medians, mads

    def find_Imin_and_background_noise(self):
        """
        Find Imin, the mode of the negative control intensities at each concentration, and standard
        deviation of negative control at each point.
        """
        self.Imin = self.find_modes_given_read_names(self.neg_contol_target_read_names)
        self.Imin_given_conc = {conc: Imin for conc, Imin in zip(self.concentrations, self.Imin)}
        self.Imin_stdev = self.find_stdevs_given_read_names(self.neg_contol_target_read_names)
        self.Imin_stdev_given_conc = {conc: Imin for conc, Imin in zip(self.concentrations, self.Imin_stdev)}

    def find_Imax(self, fit_seq_list=None):
        """
        Find Imax for each concentration, given by
            Imax = max((Fit Imax on perfect reads), median(perfect read intensities))
        """
        def Iobs(x, Kd, Imax):
            Imin = self.Imin[0] #self.Imin_given_conc[x]
            return (Imax - Imin)/(1 + (float(Kd)/x)) + Imin

        read_names = self._sample_if_needed(self.perfect_target_read_names)

        all_concentrations, all_intensities = self.get_all_conc_and_inten_given_read_names(read_names)

        def Iobs_sq_error(params):
            Kd, Imax = params
            return sum((Iobs(conc, Kd, Imax) - inten)**2
                       for conc, inten in zip(all_concentrations, all_intensities))

        res = minimize(Iobs_sq_error, x0=(1, 500), bounds=((0, None), (0, None)), method='L-BFGS-B',
                       options=dict(maxiter=10000))
        if not res.success:
            raise RuntimeError('Failed fit of perfect sequences:\n{}'.format(str(res)))

        self.perfect_Kd, self.fit_Imax = res.x

        perfect_medians = self.find_medians_given_read_names(self.perfect_target_read_names)
        self.Imax = []
        for conc, Imin, med in zip(self.concentrations, self.Imin, perfect_medians):
            fit_Iobs = Iobs(conc, self.perfect_Kd, self.fit_Imax)
            if fit_Iobs < 0.9 * self.fit_Imax:
                self.Imax.append(self.fit_Imax)
            else:
                # Using the median value as real Iobs, solve for Imax
                self.Imax.append(Imin + (1 + (self.perfect_Kd/conc)) * (med - Imin))
        self.Imax_given_conc = {conc: Imax for conc, Imax in zip(self.concentrations, self.Imax)}

        if fit_seq_list:
            self.EM_Imax_fit(fit_seq_list)

    def EM_Imax_fit(self, fit_seq_list):
        def weighted_intensity_se(Imax, conc, Imin, Imin_stdev, score_given_read_name):
            se = 0
            for seq in fit_seq_list:
                I_theoretical = (Imax - Imin) / (1.0 + self.Kds[seq]/conc) + Imin
                read_names = self._sample_if_needed(self.read_names_given_seq[seq], k=200)
                for read_name in read_names:
                    if read_name in score_given_read_name:
                        se += (I_theoretical - score_given_read_name[read_name])**2
            sigma_cluster = ((self.fit_Imax - Imin) / (1.0 + self.Kds[seq]/conc)) / 2.0 + Imin_stdev
            sigma_Imax = (self.fit_Imax - Imin)/10.0
            return se/(sigma_cluster)**2 + (self.fit_Imax - Imax)**2/(sigma_Imax)**2
            
        def fit_Imax_given_Kds():
            se = 0
            for i in range(len(self.h5_fpaths)):
                h5_fpath = self.h5_fpaths[i]
                conc = self.concentrations[i]
                Imin = self.Imin[i]
                Imin_stdev = self.Imin_stdev[i]
                curr_Imax = self.Imax[i]
                score_given_read_name = self.int_scores.score_given_read_name_in_channel[h5_fpath][self.signal_channel]
                res = minimize(weighted_intensity_se,
                               x0=curr_Imax,
                               method='Powell',
                               args=(conc, Imin, Imin_stdev, score_given_read_name),
                               options=dict(maxiter=100000))
                if not res.success:
                    sys.stdout.write('X')
                self.Imax[i] = res.x
                se += res.fun
            return se


        prev_se = float('inf')
        diff_se = float('inf')
        i = 0
        while abs(diff_se) > 1:
            sys.stdout.write('.')
            sys.stdout.flush()
            self.fit_Kds_for_given_seqs(fit_seq_list, verbose=False)
            curr_se = fit_Imax_given_Kds()
            i += 1
            if i % 10 == 0:
                diff_se = curr_se - prev_se
                prev_se = curr_se
                print 'SE={} ({})'.format(curr_se, diff_se)
        

    def _sample_if_needed(self, read_names, k=None, seed=42):
        """Helper function to reduce to <= max_sample_size"""
        if k is None:
            k = self.max_sample_size
        if len(read_names) > k:
            random.seed(seed)
            return random.sample(read_names, k)
        return read_names

    def slow_fit_Kd_given_read_names(self, read_names):
        """
        Fit Kd of one sequence given the corresponding read names.
        """
        def Iobs(conc, Kd):
            return (self.Imax_given_conc[conc] - self.Imin_given_conc[conc])/(1 + float(Kd)/conc) + self.Imin_given_conc[conc]

        read_names = self._sample_if_needed(read_names)
        all_concentrations, all_intensities = self.get_all_conc_and_inten_given_read_names(read_names)
        def Iobs_sq_error(Kd):
            return sum((Iobs(conc, Kd) - inten)**2
                       for conc, inten in zip(all_concentrations, all_intensities))
        res = minimize(Iobs_sq_error, (10000), options=dict(maxiter=10000))
        return float(res.x)

    def adjusted_fit_Kd_given_read_names(self, read_names):
        """
        Fit Kd of one sequence given the corresponding read names.
        """
        def Iobs(conc, Kd):
            return 1.0/(1.0 + (Kd/conc))

        read_names = self._sample_if_needed(read_names)
        all_concentrations, all_intensities = self.get_adjusted_conc_and_inten_given_read_names(read_names)
        popt, pcov = curve_fit(Iobs, all_concentrations, all_intensities, maxfev=10000)
        return popt[0]

    def model_fit_Kd_given_read_names(self, read_names):
        """
        Fit Kd of one sequence given the corresponding read names.
        """
        def weighted_intensity_se(Kd):
            se = 0
            for h5_fpath, conc, Imin, Imax, Imin_stdev in zip(
                self.h5_fpaths,
                self.concentrations,
                self.Imin,
                self.Imax,
                self.Imin_stdev):

                score_given_read_name = self.int_scores.score_given_read_name_in_channel[h5_fpath][self.signal_channel]

                I_theoretical = (Imax - Imin) / (1.0 + Kd/conc) + Imin
                sigma_cluster = ((self.fit_Imax - Imin) / (1.0 + Kd/conc)) / 2.0 + Imin_stdev

                loc_se = sum((I_theoretical - score_given_read_name[read_name])**2
                             for read_name in read_names if read_name in score_given_read_name)
                se += loc_se / sigma_cluster
            return se

        res = minimize(weighted_intensity_se, x0=(10,), bounds=((1e-8, None),), options=dict(maxiter=10000))
        if not res.success:
            sys.stdout.write('x')
            sys.stdout.flush()
        return float(res.x)
            
    def adjusted_median_fit_Kd_given_read_names(self, read_names):
        """
        Fit Kd of one sequence given the corresponding read names.
        """
        def Iobs(conc, Kd):
            return 1.0/(1.0 + (Kd/conc))

        read_names = self._sample_if_needed(read_names)
        medians, mads = self.get_adjusted_median_intensities_and_mads(read_names)
        popt, pcov = curve_fit(Iobs, self.concentrations, medians, sigma=mads, maxfev=10000)
        return popt[0]

    def old_fit_Kd_given_read_names(self, read_names):
        """
        Fit Kd of one sequence given the corresponding read names.
        """
        Imin = self.Imin[0]
        Imax = self.fit_Imax
        def Iobs(conc, Kd):
            return (Imax - Imin)/(1 + float(Kd)/conc) + Imin

        read_names = self._sample_if_needed(read_names)
        all_concentrations, all_intensities = self.get_all_conc_and_inten_given_read_names(read_names)
        popt, pcov = curve_fit(Iobs, all_concentrations, all_intensities, maxfev=10000)
        return popt[0]

    def bootstrap_Kds_given_read_names(self,
                                       read_names,
                                       num_bootstraps=10,
                                       sample_size=None,
                                       fit_func='adjusted_fit_Kd_given_read_names'):
        """
        Calculates bootstrap values of Kd by resampling the read_names.
        """
        if sample_size is None:
            sample_size = min(len(read_names), self.max_sample_size)
        fit_func = getattr(self, fit_func)

        read_names = list(read_names)
        bs_Kds = []
        for _ in range(num_bootstraps):
            samp_read_names = np.random.choice(read_names, size=sample_size, replace=True)
            bs_Kds.append(fit_func(samp_read_names))
        return bs_Kds

    def fit_Kds_for_given_seqs(self,
                               seqs,
                               fit_errors=False,
                               fit_func='adjusted_fit_Kd_given_read_names',
                               verbose=True):
        """
        Fit Kds for given list of seqs.
        """
        fit_func = getattr(self, fit_func)
        if not hasattr(self, 'Kds'):
            self.Kds, self.Kd_errors = {}, {}
            self.dGs, self.dG_errors = {}, {}
        for i, seq in enumerate(seqs):
            read_names = self.read_names_given_seq[seq]
            if len(read_names) < self.min_reads:
                continue
            if verbose and i % 100 == 0:
                sys.stdout.write('.')
                sys.stdout.flush()
            Kd = fit_func(read_names)
            self.Kds[seq] = Kd
            self.dGs[seq] = np.log(Kd)
            if fit_errors:
                bs_Kds = self.bootstrap_Kds_given_read_names(read_names)
                self.Kd_errors[seq] = np.std(bs_Kds)
                self.dG_errors[seq] = np.std(map(np.log, bs_Kds))

    def error_analysis_and_figs(self, axes, seq, conf_pct=90,
                                fit_func='adjusted_fit_Kd_given_read_names'):
        """
        For the given sequence, performs bootstrap analysis of errors for all numbers of clusters
        from 3 to 100.
        """
        ff = getattr(self, fit_func)
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

        
    def plot_raw_intensities(self, ax, seq):
        """
        Plot raw intensities, the fit curve, the simplified fit curve, Imin, Imax, and statistics.
        """
        read_names = self._sample_if_needed(self.read_names_given_seq[seq])

        if len(read_names) > 500:
            alpha = 0.01
        else:
            alpha = 0.1

        for read_name in read_names:
            nM_concs, intensities = [], []
            for conc, h5_fpath in zip(self.concentrations, self.h5_fpaths):
                score_given_read_name = self.int_scores.score_given_read_name_in_channel[h5_fpath][self.signal_channel]
                if read_name in score_given_read_name:
                    nM_concs.append(conc/1000.0)
                    intensities.append(float(score_given_read_name[read_name]))
            ax.plot(nM_concs, intensities, 'b', alpha=alpha)
        ax.set_xscale('log')

        try:
            Kd = self.Kds[seq]
        except:
            self.fit_Kds_for_given_seqs([seq])
            Kd = self.Kds[seq]

        nM_Kd = Kd/1000.0
        Imin, Imax = min(self.Imin), self.fit_Imax
        x = np.logspace(np.log10(self.concentrations[0]), np.log10(self.concentrations[-1]), 200)
        y = [(Imax - Imin)/(1 + (float(Kd)/xx)) + Imin for xx in x]
        ax.plot(x/1000, y, 'r', label='Simplified Fit Curve', linewidth=2.5)

        y = [(self.Imax_given_conc[conc] - self.Imin_given_conc[conc])/(1 + (float(Kd)/conc)) + self.Imin_given_conc[conc] 
             for conc in self.concentrations]
        ax.plot(self.nM_concentrations, y, 'k--', label='Fit Curve', linewidth=2.5)

        ax.plot(self.nM_concentrations, self.Imin, 'ko', label='$I_{min}$')
        ax.plot(self.nM_concentrations, self.Imax, 'o', color='darkgoldenrod', label='$I_{max}$')
        legend = ax.legend()
        legend.get_frame().set_facecolor('white')

        ax.set_xlabel('Concentration (nM)', fontsize=18)
        ax.set_ylabel('Intensity', fontsize=18)
        ax.set_ylim((0, 2*(max(self.Imax))))
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        text_x = 10**((3*np.log10(xlim[0]) + np.log10(xlim[1]))/4.0)
        
        ax.text(text_x, 1.5 * Imax, 
                '$K_d = %.2f$ nM' % (nM_Kd), 
                fontsize=16, va='center', ha='center')

        inc = (ylim[1] - ylim[0]) / 3
        oom = int(np.log10(inc))
        inc -= inc % max(1, int(0.05 * 10**oom))
        ax.yaxis.set_major_locator(MultipleLocator(inc))
        
        ax.set_axis_bgcolor('white')
        ax.grid(False)

    def plot_adjusted_intensities(self, ax, seq):
        """
        Plot adjusted intensities and fit curve.
        """
        read_names = self._sample_if_needed(self.read_names_given_seq[seq])

        try:
            Kd = self.Kds[seq]
        except:
            self.fit_Kds_for_given_seqs([seq])
            Kd = self.Kds[seq]
        nM_Kd = Kd/1000.0

        if len(read_names) > 500:
            alpha = 0.01
        elif len(read_names) > 100:
            alpha = 0.05
        elif len(read_names) > 50:
            alpha = 0.1
        else:
            alpha = 0.2

        for read_name in read_names:
            nM_concs, intensities = [], []
            for i, (conc, h5_fpath) in enumerate(zip(self.concentrations, self.h5_fpaths)):
                score_given_read_name = self.int_scores.score_given_read_name_in_channel[h5_fpath][self.signal_channel]
                if read_name in score_given_read_name:
                    nM_concs.append(conc/1000.0)
                    intensities.append((float(score_given_read_name[read_name])-self.Imin[i])/(self.Imax[i]-self.Imin[i]))
            ax.plot(nM_concs, intensities, 'b', alpha=alpha)
        ax.set_xscale('log')

        Imin, Imax = min(self.Imin), self.fit_Imax
        x = np.logspace(np.log10(self.concentrations[0]), np.log10(self.concentrations[-1]), 200)
        y = [1.0/(1 + (float(Kd)/xx)) for xx in x]
        ax.plot(x/1000, y, 'r', linewidth=2.5)

        ax.set_xlabel('Concentration (nM)')
        ax.set_ylabel('Intensity')
        ax.set_ylim((-0.25, 2.5))
        ax.set_yticks(range(3))
        ax.set_yticklabels(['$I_{min}$', '$I_{max}$', '$2 I_{max}$'])
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        text_x = 10**((3*np.log10(xlim[0]) + np.log10(xlim[1]))/4.0)
        
        ax.text(text_x, 1.5,
                '$K_d = %.2f$ nM' % (nM_Kd), 
                fontsize=16, va='center', ha='center')
        
        ax.set_axis_bgcolor('white')
        ax.grid(False)

    def write_values(self, out_fpath):
        with open(out_fpath, 'w') as out:
            out.write('# Concentration\tImin\tImax\n')
            for conc, Imin, Imax in zip(self.concentrations, self.Imin, self.Imax):
                out.write('\t'.join(map(str, map(float, (conc, Imin, Imax)))) + '\n')
            out.write('# Negative Control Sequence: {}\n'.format(self.neg_control_target))
            out.write('\t'.join(['# Seq', 'Kd (pM)', 'Kd error', 'ABA (kB T)', 'ABA error']) + '\n')
            neg_control_dG = self.dGs[self.neg_control_target]
            out.write(
                '\n'.join(
                    '\t'.join(map(str, [
                        seq,
                        self.Kds[seq],
                        self.Kd_errors[seq],
                        neg_control_dG - self.dGs[seq],
                        self.dG_errors[seq]
                    ]))
                    for seq in self.Kds.keys()
                )
            )
