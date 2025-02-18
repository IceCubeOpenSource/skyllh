# -*- coding: utf-8 -*-
# Authors:
#   Dr. Martin Wolf <mail@martin-wolf.org>

import numpy as np

from skyllh.core.debugging import (
    get_logger,
)
from skyllh.core.parameters import (
    ParameterModelMapper,
)
from skyllh.core.pdfratio import (
    SigSetOverBkgPDFRatio,
)
from skyllh.core.py import (
    module_class_method_name,
)


class PDSigSetOverBkgPDFRatio(
        SigSetOverBkgPDFRatio):
    def __init__(
            self,
            sig_pdf_set,
            bkg_pdf,
            cap_ratio=False,
            **kwargs):
        """Creates a PDFRatio instance for the public data.
        It takes a signal PDF set for different discrete gamma values.

        Parameters
        ----------
        sig_pdf_set : instance of PDSignalEnergyPDFSet
            The PDSignalEnergyPDFSet instance holding the set of signal energy
            PDFs.
        bkg_pdf : instance of PDDataBackgroundI3EnergyPDF
            The PDDataBackgroundI3EnergyPDF instance holding the background
            energy PDF.
        cap_ratio : bool
            Switch whether the S/B PDF ratio should get capped where no
            background is available. Default is False.
        """
        self._logger = get_logger(module_class_method_name(self, '__init__'))

        super().__init__(
            sig_pdf_set=sig_pdf_set,
            bkg_pdf=bkg_pdf,
            **kwargs)

        # Construct the instance for the fit parameter interpolation method.
        self._interpolmethod = self.interpolmethod_cls(
            func=self._get_ratio_values,
            param_grid_set=sig_pdf_set.param_grid_set)

        self.cap_ratio = cap_ratio
        if self.cap_ratio:
            self._logger.info('The energy PDF ratio will be capped!')

            # Calculate the ratio value for the phase space where no background
            # is available. We will take the p_sig percentile of the signal
            # like phase space.
            ratio_perc = 99

            # Get the log10 reco energy values where the background pdf has
            # non-zero values.
            n_logE = bkg_pdf.get_binning('log_energy').nbins
            n_sinDec = bkg_pdf.get_binning('sin_dec').nbins
            bd = bkg_pdf._hist_logE_sinDec > 0
            log10_e_bc = bkg_pdf.get_binning('log_energy').bincenters
            self.ratio_fill_value_dict = dict()
            for sig_pdf_key in sig_pdf_set.pdf_keys:
                sigpdf = sig_pdf_set[sig_pdf_key]
                sigvals = sigpdf.get_pd_by_log10_reco_e(log10_e_bc)
                sigvals = np.broadcast_to(sigvals, (n_sinDec, n_logE)).T
                r = sigvals[bd] / bkg_pdf._hist_logE_sinDec[bd]
                # Remove possible inf values.
                r = r[np.invert(np.isinf(r))]
                val = np.percentile(r[r > 1.], ratio_perc)
                self.ratio_fill_value_dict[sig_pdf_key] = val
                self._logger.info(
                    f'The cap value for the energy PDF ratio key {sig_pdf_key} '
                    f'is {val}.')

        # Create cache variables for the last ratio value and gradients in
        # order to avoid the recalculation of the ratio value when the
        # ``get_gradient`` method is called (usually after the ``get_ratio``
        # method was called).
        self._cache_tdm_trial_data_state_id = None
        self._cache_fitparams_hash = None
        self._cache_ratio = None
        self._cache_grads = None

    @property
    def cap_ratio(self):
        """Boolean switch whether to cap the ratio where no background
        information is available (True) or use the smallest possible floating
        point number greater than zero as background pdf value (False).
        """
        return self._cap_ratio

    @cap_ratio.setter
    def cap_ratio(self, b):
        self._cap_ratio = b

    def _is_cached(
            self,
            tdm,
            fitparams_hash):
        """Checks if the ratio and gradients for the given hash of local fit
        parameters are already cached.

        Parameters
        ----------
        tdm : instance of TrialDataManager
            The instance of TrialDataManager holding the trial data events.
        fitparams_hach : int
            The hash value of the local fit parameter values.

        Returns
        -------
        check : bool
            ``True`` if the ratio and gradient values are already cached,
            ``False`` otherwise.
        """
        if (self._cache_tdm_trial_data_state_id == tdm.trial_data_state_id) and\
           (self._cache_fitparams_hash == fitparams_hash) and\
           (self._cache_ratio is not None) and\
           (self._cache_grads is not None):
            return True

        return False

    def _get_hash_of_local_sig_fit_param_values(
            self,
            src_params_recarray):
        """Gets the hash of the values of the local signal fit parameters from
        the given ``src_params_recarray``.

        Parameters
        ----------
        src_params_recarray : instance of ndarray
            The (N_sources,)-shaped structured numpy ndarray holding the local
            parameter names and values of the sources.

        Returns
        -------
        hash : int
            The hash of the (N_fitparams, N_sources)-shaped tuple of tuples
            holding the values of the local signal fit parameters.
        """
        values = []
        for param_name in self.sig_param_names:
            if ParameterModelMapper.is_local_param_a_fitparam(
                    param_name, src_params_recarray):
                values.append(tuple(src_params_recarray[param_name]))

        values = tuple(values)

        return hash(values)

    def _get_ratio_values(
            self,
            tdm,
            eventdata,
            gridparams_recarray,
            n_values):
        """Select the signal PDF for the given fit parameter grid point and
        evaluates the S/B ratio for all the trial data events and sources.
        """
        n_sources = len(gridparams_recarray)

        ratio = np.empty((n_values,), dtype=np.double)

        same_pdf_for_all_sources = True
        if len(gridparams_recarray) > 1:
            for pname in gridparams_recarray.dtype.fields.keys():
                if not np.all(np.isclose(np.diff(gridparams_recarray[pname]), 0)):
                    same_pdf_for_all_sources = False
                    break
        if same_pdf_for_all_sources:
            # Special case where the grid parameter values are the same for all
            # sources for all grid parameters
            gridparams = dict(
                zip(gridparams_recarray.dtype.fields.keys(),
                    gridparams_recarray[0])
            )
            sig_pdf_key = self.sig_pdf_set.make_key(gridparams)
            sig_pdf = self.sig_pdf_set.get_pdf(sig_pdf_key)
            (ratio, sig_grads) = sig_pdf.get_pd(
                tdm=tdm,
                params_recarray=None)
        else:
            # General case, we need to loop over the sources.
            for (sidx, interpol_param_values) in enumerate(gridparams_recarray):
                m_src = np.zeros((n_sources), dtype=np.bool_)
                m_src[sidx] = True
                m_values = tdm.get_values_mask_for_source_mask(m_src)

                gridparams = dict(
                    zip(gridparams_recarray.dtype.fields.keys(),
                        interpol_param_values)
                )
                sig_pdf_key = self.sig_pdf_set.make_key(gridparams)
                sig_pdf = self.sig_pdf_set.get_pdf(sig_pdf_key)
                (sig_pd, sig_grads) = sig_pdf.get_pd(
                    tdm=tdm,
                    params_recarray=None)

                ratio[m_values] = sig_pd[m_values]

        (bkg_pd, bkg_grads) = self.bkg_pdf.get_pd(tdm=tdm)
        (bkg_pd,) = tdm.broadcast_selected_events_arrays_to_values_arrays(
            (bkg_pd,))

        m_nonzero_bkg = bkg_pd > np.finfo(np.double).resolution
        m_zero_bkg = np.invert(m_nonzero_bkg)
        if np.any(m_zero_bkg):
            ev_idxs = np.where(m_zero_bkg)[0]
            self._logger.debug(
                f'For {len(ev_idxs)} events the background probability is '
                f'zero. The event indices of these events are: {ev_idxs}')

        np.divide(
            ratio,
            bkg_pd,
            where=m_nonzero_bkg,
            out=ratio)

        if self._cap_ratio:
            ratio[m_zero_bkg] = self.ratio_fill_value_dict[sig_pdf_key]
        else:
            np.divide(
                ratio,
                np.finfo(np.double).resolution,
                where=m_zero_bkg,
                out=ratio)

        # Check for positive inf values in the ratio and set the ratio to a
        # finite number. Here we choose the maximum value of float32 to keep
        # room for additional computational operations.
        m_inf = np.isposinf(ratio)
        ratio[m_inf] = np.finfo(np.float32).max

        return ratio

    def _calculate_ratio_and_grads(
            self,
            tdm,
            src_params_recarray,
            fitparams_hash):
        """Calculates the ratio and ratio gradient values for all the trial data
        events and sources given the fit parameters using the interpolation
        method for the fit parameter. It caches the results.
        """
        (ratio, grads) = self._interpolmethod(
            tdm=tdm,
            eventdata=None,
            params_recarray=src_params_recarray)

        # Cache the ratio and gradient values.
        self._cache_fitparams_hash = fitparams_hash
        self._cache_ratio = ratio
        self._cache_grads = grads

    def get_ratio(
            self,
            tdm,
            src_params_recarray,
            tl=None):
        """Calculates the PDF ratio values for all events and sources.

        Parameters
        ----------
        tdm : instance of TrialDataManager
            The instance of TrialDataManager holding the trial data events for
            which the PDF ratio values should get calculated.
        src_params_recarray : instance of numpy structured ndarray | None
            The (N_sources,)-shaped numpy structured ndarray holding the
            parameter names and values of the sources.
            See the documentation of the
            :meth:`skyllh.core.parameters.ParameterModelMapper.create_src_params_recarray`
            method for more information.
        tl : instance of TimeLord | None
            The optional instance of TimeLord that should be used to measure
            timing information.

        Returns
        -------
        ratios : instance of ndarray
            The (N_values,)-shaped 1d numpy ndarray of float holding the PDF
            ratio value for each trial event and source.
        """
        fitparams_hash = self._get_hash_of_local_sig_fit_param_values(
            src_params_recarray)

        # Check if the ratio value is already cached.
        if self._is_cached(
                tdm=tdm,
                fitparams_hash=fitparams_hash):
            return self._cache_ratio

        self._calculate_ratio_and_grads(
            tdm=tdm,
            src_params_recarray=src_params_recarray,
            fitparams_hash=fitparams_hash)

        return self._cache_ratio

    def get_gradient(
            self,
            tdm,
            src_params_recarray,
            fitparam_id,
            tl=None):
        """Retrieves the PDF ratio gradient for the global fit parameter
        ``fitparam_id`` for each trial data event and source, given the given
        set of parameters ``src_params_recarray`` for each source.

        Parameters
        ----------
        tdm : instance of TrialDataManager
            The instance of TrialDataManager holding the trial data events for
            which the PDF ratio gradient values should get calculated.
        src_params_recarray : instance of numpy structured ndarray | None
            The (N_sources,)-shaped numpy structured ndarray holding the
            parameter names and values of the sources.
            See the documentation of the
            :meth:`skyllh.core.parameters.ParameterModelMapper.create_src_params_recarray`
            method for more information.
        fitparam_id : int
            The name of the fit parameter for which the gradient should get
            calculated.
        tl : instance of TimeLord | None
            The optional TimeLord instance that should be used to measure
            timing information.

        Returns
        -------
        grad : instance of ndarray
            The (N_values,)-shaped numpy ndarray holding the gradient values
            for all sources and trial events w.r.t. the given global fit
            parameter.
        """
        fitparams_hash = self._get_hash_of_local_sig_fit_param_values(
            src_params_recarray)

        # Calculate the gradients if they are not calculated yet.
        if not self._is_cached(
            tdm=tdm,
            fitparams_hash=fitparams_hash
        ):
            self._calculate_ratio_and_grads(
                tdm=tdm,
                src_params_recarray=src_params_recarray,
                fitparams_hash=fitparams_hash)

        tdm_n_sources = tdm.n_sources

        grad = np.zeros((tdm.get_n_values(),), dtype=np.float64)

        # Loop through the parameters of the signal PDF set and match them with
        # the global fit parameter.
        for (pidx, pname) in enumerate(
                self._sig_pdf_set.param_grid_set.params_name_list):
            if pname not in src_params_recarray.dtype.fields:
                continue
            p_gpidxs = src_params_recarray[f'{pname}:gpidx']
            src_mask = p_gpidxs == (fitparam_id + 1)
            n_sources = np.count_nonzero(src_mask)
            if n_sources == 0:
                continue
            if n_sources == tdm_n_sources:
                # This parameter applies to all sources, hence to all values,
                # and hence it's the only local parameter contributing to the
                # global parameter fitparam_id.
                return self._cache_grads[pidx]

            # The current parameter does not apply to all sources.
            # Create a values mask that matches a given source mask.
            m_values = tdm.get_values_mask_for_source_mask(src_mask)
            grad[m_values] = self._cache_grads[pidx][m_values]

        return grad
