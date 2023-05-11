# -*- coding: utf-8 -*-

"""Setup the time-dependent analysis. For now this only works on a single
dataset.
"""

import numpy as np

from skyllh.analyses.i3.publicdata_ps.backgroundpdf import (
    PDDataBackgroundI3EnergyPDF,
)
from skyllh.analyses.i3.publicdata_ps.detsigyield import (
    PDSingleParamFluxPointLikeSourceI3DetSigYieldBuilder,
)
from skyllh.analyses.i3.publicdata_ps.pdfratio import (
    PDSigSetOverBkgPDFRatio,
)
from skyllh.analyses.i3.publicdata_ps.signalpdf import (
    PDSignalEnergyPDFSet,
)
from skyllh.analyses.i3.publicdata_ps.utils import (
    # create_energy_cut_spline,
    tdm_field_func_psi,
)
from skyllh.core.analysis import (
    SingleSourceMultiDatasetLLHRatioAnalysis as Analysis,
)
from skyllh.core.backgroundpdf import (
    BackgroundTimePDF,
)
from skyllh.core.debugging import (
    get_logger,
)
from skyllh.core.event_selection import (
    SpatialBoxEventSelectionMethod,
)
from skyllh.core.expectation_maximization import (
    em_fit,
)
from skyllh.core.minimizer import (
    LBFGSMinimizerImpl,
    Minimizer,
)
from skyllh.core.minimizers.iminuit import (
    IMinuitMinimizerImpl,
)
from skyllh.core.model import (
    DetectorModel,
)
from skyllh.core.parameters import (
    Parameter,
    ParameterModelMapper,
)
from skyllh.core.pdfratio import (
    SigOverBkgPDFRatio,
)
from skyllh.core.progressbar import (
    ProgressBar,
)
from skyllh.core.random import (
    RandomStateService,
)
from skyllh.core.scrambling import (
    DataScrambler,
)
from skyllh.core.signal_generator import (
    MultiDatasetSignalGenerator,
)
from skyllh.core.signalpdf import (
    FixedBoxSignalTimePDF,
    FixedGaussianSignalTimePDF,
    RayleighPSFPointSourceSignalSpatialPDF,
)
from skyllh.core.smoothing import (
    BlockSmoothingFilter,
)
from skyllh.core.source_hypo_grouping import (
    SourceHypoGroup,
    SourceHypoGroupManager,
)
from skyllh.core.test_statistic import (
    WilksTestStatistic,
)
from skyllh.core.trialdata import (
    TrialDataManager,
)
from skyllh.core.utils.analysis import (
    pointlikesource_to_data_field_array
)
from skyllh.i3.background_generation import (
    FixedScrambledExpDataI3BkgGenMethod,
)
from skyllh.i3.backgroundpdf import (
    DataBackgroundI3SpatialPDF,
)
from skyllh.i3.livetime import (
    I3Livetime,
)
from skyllh.i3.scrambling import (
    I3SeasonalVariationTimeScramblingMethod,
)
from skyllh.physics.flux_model import (
    BoxTimeFluxProfile,
    PowerLawEnergyFluxProfile,
    SteadyPointlikeFFM,
)


def create_signal_time_pdf(
        grl,
        gauss=None,
        box=None):
    """Creates the signal time PDF, either a gaussian or a box shaped PDF.

    Parameters
    ----------
    grl : instance of numpy structured ndarray
        The structured numpy ndarray holding the good-run-list data.
    gauss : dict | None
        None or dictionary with {"mu": float, "sigma": float}.
    box : dict | None
        None or dictionary with {"start": float, "end": float}.

    Returns
    -------
    pdf : instance of PDF
        The created time PDF instance.
    """
    if (gauss is None) and (box is None):
        raise TypeError(
            'Either gauss or box have to be specified as time pdf.')

    if gauss is not None:
        pdf = FixedGaussianSignalTimePDF(
            grl=grl,
            mu=gauss['mu'],
            sigma=gauss['sigma'])
    elif box is not None:
        pdf = FixedBoxSignalTimePDF(
            grl=grl,
            start=box['start'],
            end=box['end'])

    return pdf


def change_signal_time_pdf_of_llhratio_function(
        ana,
        gauss=None,
        box=None):
    """Changes the signal time PDF of the log-likelihood ratio function.

    Parameters
    ----------
    gauss : dict | None
        None or dictionary with {"mu": float, "sigma": float}.
    box : dict | None
        None or dictionary with {"start": float, "end": float}.
    """
    grl = ana.data_list[0].grl

    time_sigpdf = create_signal_time_pdf(
        grl=grl,
        gauss=gauss,
        box=box)

    pdfratio = ana.llhratio.llhratio_list[0].pdfratio

    # pdfratio is an instance of PDFRatioProduct.
    # The first item is the PDF ratio product of the spatial and energy PDF
    # ratios. The second item is the time PDF ratio.
    pdfratio.pdfratio2.sig_pdf = time_sigpdf

    # TODO: Change detector signal yield with flare livetime in sample
    # (1 / grl_norm in pdf), rebuild the histograms if it is changed.


def get_energy_spatial_signal_over_background(
        ana,
        fitparam_values,
        tl=None):
    """Returns the signal over background ratio for
    (spatial_signal * energy_signal) / (spatial_background * energy_background).

    Parameter
    ---------
    fitparam_values : instance of ndarray
        The (N_fitparams,)-shaped numpy ndarray holding the values of the global
        fit parameters, e.g. ns and gamma.

    Returns
    -------
    ratio : 1d ndarray
        Product of spatial and energy signal over background pdfs.
    """
    tdm = ana.tdm_list[0]

    pdfratio = ana.llhratio.llhratio_list[0].pdfratio

    # pdfratio is an instance of PDFRatioProduct.
    # The first item is the PDF ratio product of the spatial and energy PDF
    # ratios. The second item is the time PDF ratio.
    pdfratio = pdfratio.pdfratio1

    src_params_recarray = ana.pmm.create_src_params_recarray(
        gflp_values=fitparam_values)

    ratio = pdfratio.get_ratio(
        tdm=tdm,
        src_params_recarray=src_params_recarray,
        tl=tl)

    return ratio


def change_fluxmodel_gamma(
        ana,
        gamma):
    """Sets the given gamma value to the flux model of the single source.

    Parameter
    ---------
    ana : instance of SingleSourceMultiDatasetLLHRatioAnalysis
        The analysis that should be used.
    gamma : float
        Spectral index for the flux model.
    """
    ana.shg_mgr.shg_list[0].fluxmodel.set_params({'gamma': gamma})
    ana.change_shg_mgr(shg_mgr=ana.shg_mgr)


def change_signal_time(
        ana,
        gauss=None,
        box=None):
    """Change the signal injection to gauss or box.

    Parameters
    ----------
    ana : instance of SingleSourceMultiDatasetLLHRatioAnalysis
        The analysis that should be used.
    gauss : dict | None
        None or dictionary {"mu": float, "sigma": float}.
    box : dict | None
        None or dictionary {"start" : float, "end" : float}.
    """
    # FIXME
    ana.sig_generator.set_flare(box=box, gauss=gauss)


def calculate_TS(
        ana,
        em_results,
        rss):
    """Calculate the best TS value for the expectation maximization gamma scan.

    Parameters
    ----------
    ana : instance of SingleSourceMultiDatasetLLHRatioAnalysis
        The analysis that should be used.
    em_results : 1d ndarray of tuples
        Gamma scan result.
    rss : instance of RandomStateService
        The instance of RandomStateService that should be used to generate
        random numbers from.

    Returns
    -------
    max_TS : float
        The maximal TS value of all maximized time hypotheses.
    best_time : tuple
        The tuple(gamma from em scan [float], best fit mean time [float],
        best fit width [float]).
    best_fitparam_values : instance of numpy ndarray
        The instance of numpy ndarray holding the fit parameter values of the
        overall best fit result.
    """
    max_TS = 0
    best_time = None
    best_fitparam_values = None
    for em_result in em_results:
        change_signal_time_pdf_of_llhratio_function(
            ana=ana,
            gauss={
                'mu': em_result['mu'],
                'sigma': em_result['sigma']})

        (log_lambda_max, fitparam_values, status) = ana.llhratio.maximize(
            rss=rss)

        TS = ana.calculate_test_statistic(
            log_lambda=log_lambda_max,
            fitparam_values=fitparam_values)

        if TS > max_TS:
            max_TS = TS
            best_time = em_result
            best_fitparam_values = fitparam_values

    return (max_TS, best_time, best_fitparam_values)


def run_gamma_scan_single_flare(
        ana,
        remove_time=None,
        gamma_min=1,
        gamma_max=5,
        n_gamma=51):
    """Run em for different gammas in the signal energy pdf

    Parameters
    ----------
    ana : instance of SingleSourceMultiDatasetLLHRatioAnalysis
        The analysis that should be used.
    remove_time : float
        Time information of event that should be removed.
    gamma_min : float
        Lower bound for gamma scan.
    gamma_max : float
        Upper bound for gamma scan.
    n_gamma : int
        Number of steps for gamma scan.

    Returns
    -------
    results : instance of numpy structured ndarray
        The numpy structured ndarray with fields

        gamma : float
            The spectral index value.
        mu : float
            The determined mean value of the gauss curve.
        sigma : float
            The determoned standard deviation of the gauss curve.
        ns_em : float
            The scaling factor of the flare.
    """
    dtype = [
        ('gamma', np.float64),
        ('mu', np.float64),
        ('sigma', np.float64),
        ('ns_em', np.float64)
    ]
    results = np.empty(n_gamma, dtype=dtype)

    time = ana._tdm_list[0].get_data('time')

    for (i, g) in enumerate(np.linspace(gamma_min, gamma_max, n_gamma)):
        ratio = get_energy_spatial_signal_over_background(ana, {"gamma": g})
        (mu, sigma, ns) = em_fit(
            time,
            ratio,
            n=1,
            tol=1.e-200,
            iter_max=500,
            weight_thresh=0,
            initial_width=5000,
            remove_x=remove_time)
        results[i] = (g, mu[0], sigma[0], ns[0])

    return results


def unblind_flare(
        ana,
        remove_time=None):
    """Run EM on unscrambled data. Similar to the original analysis, remove the alert event.

    Parameters
    ----------
    ana : instance of SingleSourceMultiDatasetLLHRatioAnalysis
        The analysis that should be used.
    remove_time : float
        Time information of event that should be removed.
        In the case of the TXS analysis: remove_time=58018.8711856

    Returns
    -------
    results :
    array with "gamma", "mu", "sigma", and scaling factor for flare "ns_em"
    """
    rss = RandomStateService(seed=1)
    ana.unblind(
        rss=rss)
    results = run_gamma_scan_single_flare(
        ana=ana,
        remove_time=remove_time)

    return results


def create_analysis(  # noqa: C901
        datasets,
        source,
        gauss=None,
        box=None,
        refplflux_Phi0=1,
        refplflux_E0=1e3,
        refplflux_gamma=2.0,
        ns_seed=100.0,
        ns_min=0.,
        ns_max=1e3,
        gamma_seed=3.0,
        gamma_min=1.,
        gamma_max=5.,
        kde_smoothing=False,
        minimizer_impl="LBFGS",
        cut_sindec=None,
        spl_smooth=None,
        cap_ratio=False,
        compress_data=False,
        keep_data_fields=None,
        evt_sel_delta_angle_deg=10,
        construct_bkg_generator=True,
        construct_sig_generator=True,
        tl=None,
        ppbar=None,
        logger_name=None):
    """Creates the Analysis instance for this particular analysis.

    Parameters:
    -----------
    datasets : list of Dataset instances
        The list of Dataset instances, which should be used in the
        analysis.
    source : PointLikeSource instance
        The PointLikeSource instance defining the point source position.
    gauss : None or dictionary with mu, sigma
        None if no Gaussian time pdf. Else dictionary with {"mu": float, "sigma": float} of Gauss
    box : None or dictionary with start, end
        None if no Box shaped time pdf. Else dictionary with {"start": float, "end": float} of box.
    refplflux_Phi0 : float
        The flux normalization to use for the reference power law flux model.
    refplflux_E0 : float
        The reference energy to use for the reference power law flux model.
    refplflux_gamma : float
        The spectral index to use for the reference power law flux model.
    ns_seed : float
        Value to seed the minimizer with for the ns fit.
    ns_min : float
        Lower bound for ns fit.
    ns_max : float
        Upper bound for ns fit.
    gamma_seed : float | None
        Value to seed the minimizer with for the gamma fit. If set to None,
        the refplflux_gamma value will be set as gamma_seed.
    gamma_min : float
        Lower bound for gamma fit.
    gamma_max : float
        Upper bound for gamma fit.
    kde_smoothing : bool
        Apply a KDE-based smoothing to the data-driven background pdf.
        Default: False.
    minimizer_impl : str | "LBFGS"
        Minimizer implementation to be used. Supported options are "LBFGS"
        (L-BFG-S minimizer used from the :mod:`scipy.optimize` module), or
        "minuit" (Minuit minimizer used by the :mod:`iminuit` module).
        Default: "LBFGS".
    cut_sindec : list of float | None
        sin(dec) values at which the energy cut in the southern sky should
        start. If None, np.sin(np.radians([-2, 0, -3, 0, 0])) is used.
    spl_smooth : list of float
        Smoothing parameters for the 1D spline for the energy cut. If None,
        [0., 0.005, 0.05, 0.2, 0.3] is used.
    cap_ratio : bool
        If set to True, the energy PDF ratio will be capped to a finite value
        where no background energy PDF information is available. This will
        ensure that an energy PDF ratio is available for high energies where
        no background is available from the experimental data.
        If kde_smoothing is set to True, cap_ratio should be set to False!
        Default is False.
    compress_data : bool
        Flag if the data should get converted from float64 into float32.
    keep_data_fields : list of str | None
        List of additional data field names that should get kept when loading
        the data.
    evt_sel_delta_angle_deg : float
        The delta angle in degrees for the event selection optimization methods.
    construct_bkg_generator : bool
        Flag if the background generator should be constructed (``True``) or not
        (``False``).
    construct_sig_generator : bool
        Flag if the signal generator should be constructed (``True``) or not
        (``False``).
    tl : TimeLord instance | None
        The TimeLord instance to use to time the creation of the analysis.
    ppbar : ProgressBar instance | None
        The instance of ProgressBar for the optional parent progress bar.
    logger_name : str | None
        The name of the logger to be used. If set to ``None``, ``__name__`` will
        be used.

    Returns
    -------
    ana : SingleSourceMultiDatasetLLHRatioAnalysis
        The Analysis instance for this analysis.
    """
    if logger_name is None:
        logger_name = __name__
    logger = get_logger(logger_name)

    if len(datasets) != 1:
        raise RuntimeError(
            'This analysis supports only analyses with only single datasets '
            'at the moment!')

    if (gauss is None) and (box is None):
        raise ValueError(
            'No time pdf specified (box or gauss)!')
    if (gauss is not None) and (box is not None):
        raise ValueError(
            'Time PDF cannot be both gaussian and box shaped. '
            'Please specify only one shape.')

    # Create the minimizer instance.
    if minimizer_impl == 'LBFGS':
        minimizer = Minimizer(LBFGSMinimizerImpl())
    elif minimizer_impl == 'minuit':
        minimizer = Minimizer(IMinuitMinimizerImpl(ftol=1e-8))
    else:
        raise NameError(
            f"Minimizer implementation `{minimizer_impl}` is not supported "
            "Please use `LBFGS` or `minuit`.")

    # Define the flux model.
    fluxmodel = SteadyPointlikeFFM(
        Phi0=refplflux_Phi0,
        energy_profile=PowerLawEnergyFluxProfile(
            E0=refplflux_E0,
            gamma=refplflux_gamma))

    # Define the fit parameter ns.
    param_ns = Parameter(
        name='ns',
        initial=ns_seed,
        valmin=ns_min,
        valmax=ns_max)

    # Define the fit parameter gamma.
    param_gamma = Parameter(
        name='gamma',
        initial=gamma_seed,
        valmin=gamma_min,
        valmax=gamma_max)

    # Define the detector signal yield builder for the IceCube detector and this
    # source and flux model.
    # The sin(dec) binning will be taken by the builder automatically from the
    # Dataset instance.
    gamma_grid = param_gamma.as_linear_grid(delta=0.1)
    detsigyield_builder =\
        PDSingleParamFluxPointLikeSourceI3DetSigYieldBuilder(
            param_grid=gamma_grid)

    # Define the signal generation method.
    sig_gen_method = None

    # Create a source hypothesis group manager with a single source hypothesis
    # group for the single source.
    shg_mgr = SourceHypoGroupManager(
        SourceHypoGroup(
            sources=source,
            fluxmodel=fluxmodel,
            detsigyield_builders=detsigyield_builder,
            sig_gen_method=sig_gen_method))

    # Define a detector model for the ns fit parameter.
    detector_model = DetectorModel('IceCube')

    # Define the parameter model mapper for the analysis, which will map global
    # parameters to local source parameters.
    pmm = ParameterModelMapper(
        models=[detector_model, source])
    pmm.def_param(param_ns, models=detector_model)
    pmm.def_param(param_gamma, models=source)

    logger.info(str(pmm))

    # Define the test statistic.
    test_statistic = WilksTestStatistic()

    # Create the Analysis instance.
    ana = Analysis(
        shg_mgr=shg_mgr,
        pmm=pmm,
        test_statistic=test_statistic,
        sig_generator_cls=MultiDatasetSignalGenerator,
    )

    # Define the event selection method for pure optimization purposes.
    # We will use the same method for all datasets.
    event_selection_method = SpatialBoxEventSelectionMethod(
        shg_mgr=shg_mgr,
        delta_angle=np.deg2rad(evt_sel_delta_angle_deg))

    # Prepare the spline parameters for the signal generator.
    if cut_sindec is None:
        cut_sindec = np.sin(np.radians([-2, 0, -3, 0, 0]))
    if spl_smooth is None:
        spl_smooth = [0., 0.005, 0.05, 0.2, 0.3]
    if len(spl_smooth) < len(datasets) or len(cut_sindec) < len(datasets):
        raise AssertionError(
            'The length of the spl_smooth and of the cut_sindec must be equal '
            f'to the length of datasets: {len(datasets)}.')

    # Add the data sets to the analysis.
    pbar = ProgressBar(len(datasets), parent=ppbar).start()
    data_list = []
    for (ds_idx, ds) in enumerate(datasets):
        data = ds.load_and_prepare_data(
            keep_fields=keep_data_fields,
            compress=compress_data,
            tl=tl)
        data_list.append(data)

        sin_dec_binning = ds.get_binning_definition('sin_dec')
        log_energy_binning = ds.get_binning_definition('log_energy')

        # Create the spatial PDF ratio instance for this dataset.
        spatial_sigpdf = RayleighPSFPointSourceSignalSpatialPDF(
            dec_range=np.arcsin(sin_dec_binning.range))
        spatial_bkgpdf = DataBackgroundI3SpatialPDF(
            data_exp=data.exp,
            sin_dec_binning=sin_dec_binning)
        spatial_pdfratio = SigOverBkgPDFRatio(
            sig_pdf=spatial_sigpdf,
            bkg_pdf=spatial_bkgpdf)

        # Create the energy PDF ratio instance for this dataset.
        energy_sigpdfset = PDSignalEnergyPDFSet(
            ds=ds,
            src_dec=source.dec,
            flux_model=fluxmodel,
            param_grid_set=gamma_grid,
            ppbar=ppbar
        )
        smoothing_filter = BlockSmoothingFilter(nbins=1)
        energy_bkgpdf = PDDataBackgroundI3EnergyPDF(
            data_exp=data.exp,
            logE_binning=log_energy_binning,
            sinDec_binning=sin_dec_binning,
            smoothing_filter=smoothing_filter,
            kde_smoothing=kde_smoothing)

        energy_pdfratio = PDSigSetOverBkgPDFRatio(
            sig_pdf_set=energy_sigpdfset,
            bkg_pdf=energy_bkgpdf,
            cap_ratio=cap_ratio)

        pdfratio = spatial_pdfratio * energy_pdfratio

        # Create the time PDF ratio instance for this dataset.
        if (gauss is not None) or (box is not None):
            livetime = I3Livetime.from_grl_data(
                grl_data=data.grl)
            time_bkgpdf = BackgroundTimePDF(
                livetime=livetime,
                time_flux_profile=BoxTimeFluxProfile.from_start_and_stop_time(
                    start=livetime.time_start,
                    stop=livetime.time_stop))
            if gauss is not None:
                time_sigpdf = FixedGaussianSignalTimePDF(
                    grl=data.grl,
                    mu=gauss['mu'],
                    sigma=gauss['sigma'])
            elif box is not None:
                time_sigpdf = FixedBoxSignalTimePDF(
                    grl=data.grl,
                    start=box['start'],
                    end=box['end'])
            time_pdfratio = SigOverBkgPDFRatio(
                sig_pdf=time_sigpdf,
                bkg_pdf=time_bkgpdf,
            )

            pdfratio = pdfratio * time_pdfratio

        # Create a trial data manager and add the required data fields.
        tdm = TrialDataManager()
        tdm.add_source_data_field(
            name='src_array',
            func=pointlikesource_to_data_field_array)
        tdm.add_data_field(
            name='psi',
            func=tdm_field_func_psi,
            dt='dec',
            is_srcevt_data=True)

        # energy_cut_spline = create_energy_cut_spline(
        #     ds,
        #     data.exp,
        #     spl_smooth[ds_idx])#
        #
        sig_generator = None
        # if construct_sig_generator is True:
        #     pass # FIXME

        ana.add_dataset(
            ds,
            data,
            pdfratio=pdfratio,
            tdm=tdm,
            event_selection_method=event_selection_method,
            sig_generator=sig_generator)

        pbar.increment()
    pbar.finish()

    ana.construct_services(
        ppbar=ppbar)

    ana.llhratio = ana.construct_llhratio(
        minimizer=minimizer,
        ppbar=ppbar)

    # Define the data scrambler with its data scrambling method, which is used
    # for background generation.

    # FIXME: Support multiple datasets for the DataScrambler.
    data_scrambler = DataScrambler(
        I3SeasonalVariationTimeScramblingMethod(
            data_list[0]))
    bkg_gen_method = FixedScrambledExpDataI3BkgGenMethod(data_scrambler)
    ana.bkg_gen_method = bkg_gen_method

    if construct_bkg_generator is True:
        ana.construct_background_generator()

    if construct_sig_generator is True:
        ana.construct_signal_generator()

    return ana
