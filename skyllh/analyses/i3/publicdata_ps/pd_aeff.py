# -*- coding: utf-8 -*-

import numpy as np

from scipy import interpolate
from scipy import integrate

from skyllh.core.binning import (
    get_bincenters_from_binedges,
    get_bin_indices_from_lower_and_upper_binedges,
)
from skyllh.core.storage import create_FileLoader

from skyllh.analyses.i3.publicdata_ps.utils import FctSpline2D


def load_effective_area_array(pathfilenames):
    """Loads the (nbins_decnu, nbins_log10enu)-shaped 2D effective
    area array from the given data file.

    Parameters
    ----------
    pathfilename : str | list of str
        The file name of the data file.

    Returns
    -------
    aeff_decnu_log10enu : (nbins_decnu, nbins_log10enu)-shaped 2D ndarray
        The ndarray holding the effective area for each
        (dec_nu,log10(E_nu/GeV)) bin.
    decnu_binedges_lower : (nbins_decnu,)-shaped ndarray
        The ndarray holding the lower bin edges of the dec_nu axis.
    decnu_binedges_upper : (nbins_decnu,)-shaped ndarray
        The ndarray holding the upper bin edges of the dec_nu axis.
    log10_enu_binedges_lower : (nbins_log10enu,)-shaped ndarray
        The ndarray holding the lower bin edges of the log10(E_nu/GeV) axis.
    log10_enu_binedges_upper : (nbins_log10enu,)-shaped ndarray
        The ndarray holding the upper bin edges of the log10(E_nu/GeV) axis.
    """
    loader = create_FileLoader(pathfilenames=pathfilenames)
    data = loader.load_data()
    renaming_dict = {
        'log10(E_nu/GeV)_min': 'log10_enu_min',
        'log10(E_nu/GeV)_max': 'log10_enu_max',
        'Dec_nu_min[deg]':     'decnu_min',
        'Dec_nu_max[deg]':     'decnu_max',
        'A_Eff[cm^2]':         'a_eff'
    }
    data.rename_fields(renaming_dict, must_exist=True)

    # Convert the true neutrino declination from degrees to radians.
    data['decnu_min'] = np.deg2rad(data['decnu_min'])
    data['decnu_max'] = np.deg2rad(data['decnu_max'])

    # Determine the binning for energy and declination.
    log10_enu_binedges_lower = np.unique(data['log10_enu_min'])
    log10_enu_binedges_upper = np.unique(data['log10_enu_max'])
    decnu_binedges_lower = np.unique(data['decnu_min'])
    decnu_binedges_upper = np.unique(data['decnu_max'])

    if(len(log10_enu_binedges_lower) != len(log10_enu_binedges_upper)):
        raise ValueError('Cannot extract the log10(E/GeV) binning of the '
            'effective area from data file "{}". The number of lower and upper '
            'bin edges is not equal!'.format(str(loader.pathfilename_list)))
    if(len(decnu_binedges_lower) != len(decnu_binedges_upper)):
        raise ValueError('Cannot extract the dec_nu binning of the effective '
            'area from data file "{}". The number of lower and upper bin edges '
            'is not equal!'.format(str(loader.pathfilename_list)))

    nbins_log10_enu = len(log10_enu_binedges_lower)
    nbins_decnu = len(decnu_binedges_lower)

    # Construct the 2d array for the effective area.
    aeff_decnu_log10enu = np.zeros(
        (nbins_decnu, nbins_log10_enu), dtype=np.double)

    decnu_idx = np.digitize(
        0.5*(data['decnu_min'] +
             data['decnu_max']),
        decnu_binedges_lower) - 1
    log10enu_idx = np.digitize(
        0.5*(data['log10_enu_min'] +
             data['log10_enu_max']),
        log10_enu_binedges_lower) - 1

    aeff_decnu_log10enu[decnu_idx, log10enu_idx] = data['a_eff']

    return (
        aeff_decnu_log10enu,
        decnu_binedges_lower,
        decnu_binedges_upper,
        log10_enu_binedges_lower,
        log10_enu_binedges_upper
    )


class PDAeff(object):
    """This class provides a representation of the effective area provided by
    the public data.
    """
    def __init__(
            self, pathfilenames, **kwargs):
        """Creates an effective area instance by loading the effective area
        data from the given file.
        """
        super().__init__(**kwargs)

        (
            self._aeff_decnu_log10enu,
            self._decnu_binedges_lower,
            self._decnu_binedges_upper,
            self._log10_enu_binedges_lower,
            self._log10_enu_binedges_upper
        ) = load_effective_area_array(pathfilenames)

        # Note: self._aeff_decnu_log10enu is numpy 2D ndarray of shape
        # (nbins_decnu, nbins_log10enu).

        # Cut the energies where all effective areas are zero.
        m = np.sum(self._aeff_decnu_log10enu, axis=0) > 0
        self._aeff_decnu_log10enu = self._aeff_decnu_log10enu[:,m]
        self._log10_enu_binedges_lower = self._log10_enu_binedges_lower[m]
        self._log10_enu_binedges_upper = self._log10_enu_binedges_upper[m]

        self._decnu_binedges = np.concatenate(
            (self._decnu_binedges_lower,
             self._decnu_binedges_upper[-1:])
        )
        self._log10_enu_binedges = np.concatenate(
            (self._log10_enu_binedges_lower,
             self._log10_enu_binedges_upper[-1:])
        )

    @property
    def decnu_binedges(self):
        """(read-only) The bin edges of the neutrino declination axis in
        radians.
        """
        return self._decnu_binedges

    @property
    def sin_decnu_binedges(self):
        """(read-only) The sin of the bin edges of the neutrino declination
        in radians.
        """
        return np.sin(self._decnu_binedges)

    @property
    def decnu_bincenters(self):
        """(read-only) The bin center values of the neutrino declination axis in
        radians.
        """
        return get_bincenters_from_binedges(self._decnu_binedges)

    @property
    def n_decnu_bins(self):
        """(read-only) The number of bins of the neutrino declination axis.
        """
        return len(self._decnu_binedges) - 1

    @property
    def log10_enu_binedges(self):
        """(read-only) The bin edges of the log10(E_nu/GeV) neutrino energy
        axis.
        """
        return self._log10_enu_binedges

    @property
    def log10_enu_bincenters(self):
        """(read-only) The bin center values of the log10(E_nu/GeV) neutrino
        energy axis.
        """
        return get_bincenters_from_binedges(self._log10_enu_binedges)

    @property
    def n_log10_enu_bins(self):
        """(read-only) The number of bins of the log10 neutrino energy axis.
        """
        return len(self._log10_enu_binedges) - 1

    @property
    def aeff_decnu_log10enu(self):
        """(read-only) The effective area in cm^2 as (n_decnu,n_log10enu)-shaped
        2D numpy ndarray.
        """
        return self._aeff_decnu_log10enu

    def create_sin_decnu_log10_enu_spline(self):
        """Creates a FctSpline2D object representing a 2D spline of the
        effective area in sin(dec_nu)-log10(E_nu/GeV)-space.

        Returns
        -------
        spl : FctSpline2D instance
            The FctSpline2D instance representing a spline in the
            sin(dec_nu)-log10(E_nu/GeV)-space.
        """
        spl = FctSpline2D(
            self._aeff_decnu_log10enu,
            self.sin_decnu_binedges,
            self.log10_enu_binedges
        )
        return spl

    def get_aeff_for_decnu(self, decnu):
        """Retrieves the effective area as function of log10_enu.

        Parameters
        ----------
        decnu : float
            The true neutrino declination.

        Returns
        -------
        aeff : (n,)-shaped numpy ndarray
            The effective area in cm^2 for the given true neutrino declination
            as a function of log10 true neutrino energy.
        """
        decnu_idx = np.digitize(decnu, self._decnu_binedges) - 1

        aeff = self._aeff_decnu_log10enu[decnu_idx]

        return aeff

    #def get_detection_pd_for_sin_true_dec(self, sin_true_dec, true_e):
        #"""Calculates the detection probability density p(E_nu|sin_dec) in
        #unit GeV^-1 for the given true energy values.

        #Parameters
        #----------
        #sin_true_dec : float
            #The sin of the true declination.
        #true_e : (n,)-shaped 1d numpy ndarray of float
            #The values of the true energy in GeV for which the probability
            #density value should get calculated.

        #Returns
        #-------
        #det_pd : (n,)-shaped 1d numpy ndarray of float
            #The detection probability density values for the given true energy
            #value.
        #"""
        #aeff = self.get_aeff_for_sin_true_dec(sin_true_dec)

        #dE = np.diff(np.power(10, self.log_true_e_binedges))

        #det_pdf = aeff / np.sum(aeff) / dE

        #x = np.power(10, self.log_true_e_bincenters)
        #y = det_pdf
        #tck = interpolate.splrep(x, y, k=1, s=0)

        #det_pd = interpolate.splev(true_e, tck, der=0)

        #return det_pd

    #def get_detection_pd_in_log10E_for_sin_true_dec(
            #self, sin_true_dec, log10_true_e):
        #"""Calculates the detection probability density p(E_nu|sin_dec) in
        #unit log10(GeV)^-1 for the given true energy values.

        #Parameters
        #----------
        #sin_true_dec : float
            #The sin of the true declination.
        #log10_true_e : (n,)-shaped 1d numpy ndarray of float
            #The log10 values of the true energy in GeV for which the
            #probability density value should get calculated.

        #Returns
        #-------
        #det_pd : (n,)-shaped 1d numpy ndarray of float
            #The detection probability density values for the given true energy
            #value.
        #"""
        #aeff = self.get_aeff_for_sin_true_dec(sin_true_dec)

        #dlog10E = np.diff(self.log_true_e_binedges)

        #det_pdf = aeff / np.sum(aeff) / dlog10E

        #spl = interpolate.splrep(
            #self.log_true_e_bincenters, det_pdf, k=1, s=0)

        #det_pd = interpolate.splev(log10_true_e, spl, der=0)

        #return det_pd

    def get_detection_prob_for_decnu(
            self, decnu, enu_min, enu_max, enu_range_min, enu_range_max):
        """Calculates the detection probability for a given true neutrino energy
        range for a given neutrino declination.

        Parameters
        ----------
        decnu : float
            The neutrino declination in radians.
        enu_min : float
            The minimum energy in GeV.
        enu_max : float
            The maximum energy in GeV.
        enu_range_min : float
            The minimum energy in GeV of the entire energy range.
        enu_range_max : float
            The maximum energy in GeV of the entire energy range.

        Returns
        -------
        det_prob : float
            The neutrino energy detection probability.
        """
        enu_binedges = np.power(10, self.log10_enu_binedges)

        # Get the bin indices for the lower and upper energy range values.
        (lidx,) = get_bin_indices_from_lower_and_upper_binedges(
            enu_binedges[:-1],
            enu_binedges[1:],
            np.array([enu_range_min])
        )
        if enu_range_max >= enu_binedges[-1]:
            uidx = len(enu_binedges)-2
        else:
            (uidx,) = get_bin_indices_from_lower_and_upper_binedges(
                enu_binedges[:-1],
                enu_binedges[1:],
                np.array([enu_range_max])
            )
            # Note: The get_bin_indices_from_lower_and_upper_binedges function
            #       is based on the lower edges. So by definition the upper bin
            #       index is one too large.
            uidx -= 1

        aeff = self.get_aeff_for_decnu(decnu)
        aeff = aeff[lidx+1:uidx+1]
        enu_binedges = enu_binedges[lidx+1:uidx+2]

        dE = np.diff(enu_binedges)

        daeff_dE = aeff / dE

        # Create a spline representation that spans the entire enu range.
        x = np.empty((len(enu_binedges)+1,), dtype=np.double)
        x[0] = enu_binedges[0]
        x[1:-1] = get_bincenters_from_binedges(enu_binedges)
        x[-1] = enu_binedges[-1]

        y = np.empty((len(enu_binedges)+1,), dtype=np.double)
        y[0] = daeff_dE[0]
        y[1:-1] = daeff_dE
        y[-1] = daeff_dE[-1]

        spl = interpolate.splrep(
            x,
            y,
            xb=enu_range_min,
            xe=enu_range_max,
            k=1,
            s=0
        )

        def _eval_spl_func(x):
            return interpolate.splev(x, spl, der=0, ext=1)

        norm = integrate.quad(
            _eval_spl_func,
            enu_range_min,
            enu_range_max,
            limit=200,
            full_output=1
        )[0]

        integral = integrate.quad(
            _eval_spl_func,
            enu_min,
            enu_max,
            limit=200,
            full_output=1
        )[0]

        det_prob = integral / norm

        return det_prob

    #def get_aeff_integral_for_sin_true_dec(
            #self, sin_true_dec, log_true_e_min, log_true_e_max):
        #"""Calculates the integral of the effective area using the trapezoid
        #method.

        #Returns
        #-------
        #integral : float
            #The integral in unit cm^2 GeV.
        #"""
        #aeff = self.get_aeff_for_sin_true_dec(sin_true_dec)

        #integral = (
            #(np.power(10, log_true_e_max) -
             #np.power(10, log_true_e_min)) *
            #0.5 *
            #(np.interp(log_true_e_min, self.log_true_e_bincenters, aeff) +
             #np.interp(log_true_e_max, self.log_true_e_bincenters, aeff))
        #)

        #return integral

    #def get_aeff(self, sin_true_dec, log_true_e):
        #"""Retrieves the effective area for the given sin(dec_true) and
        #log(E_true) value pairs.

        #Parameters
        #----------
        #sin_true_dec : (n,)-shaped 1D ndarray
            #The sin(dec_true) values.
        #log_true_e : (n,)-shaped 1D ndarray
            #The log(E_true) values.

        #Returns
        #-------
        #aeff : (n,)-shaped 1D ndarray
            #The 1D ndarray holding the effective area values for each value
            #pair. For value pairs outside the effective area data zero is
            #returned.
        #"""
        #valid = (
            #(sin_true_dec >= self.sin_true_dec_binedges[0]) &
            #(sin_true_dec <= self.sin_true_dec_binedges[-1]) &
            #(log_true_e >= self.log_true_e_binedges[0]) &
            #(log_true_e <= self.log_true_e_binedges[-1])
        #)
        #sin_true_dec_idxs = np.digitize(
            #sin_true_dec[valid], self.sin_true_dec_binedges) - 1
        #log_true_e_idxs = np.digitize(
            #log_true_e[valid], self.log_true_e_binedges) - 1

        #aeff = np.zeros((len(valid),), dtype=np.double)
        #aeff[valid] = self.aeff_arr[sin_true_dec_idxs,log_true_e_idxs]

        #return aeff