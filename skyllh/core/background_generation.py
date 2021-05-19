# -*- coding: utf-8 -*-

import abc
import numpy as np

from skyllh.core.optimize import (
    AllEventSelectionMethod,
    EventSelectionMethod
)
from skyllh.core.py import (
    float_cast,
    func_has_n_args,
    issequenceof
)
from skyllh.core.debugging import get_logger
from skyllh.core.scrambling import DataScrambler
from skyllh.core.timing import TaskTimer
from skyllh.core.config import CFG


logger = get_logger(__name__)


class BackgroundGenerationMethod(object, metaclass=abc.ABCMeta):
    """This is the abstract base class for a detector specific background
    generation method.
    """

    def __init__(self):
        """Constructs a new background generation method instance.
        """
        super(BackgroundGenerationMethod, self).__init__()

    def change_source_hypo_group_manager(self, src_hypo_group_manager):
        """Notifies the background generation method about an updated
        SourceHypoGroupManager instance.

        Parameters
        ----------
        src_hypo_group_manager : SourceHypoGroupManager instance
            The new SourceHypoGroupManager instance.
        """
        pass

    @abc.abstractmethod
    def generate_events(self, rss, dataset, data, mean, tl=None, **kwargs):
        """This method is supposed to generate a `mean` number of background
        events for the given dataset and its data.

        Parameters
        ----------
        rss : instance of RandomStateService
            The instance of RandomStateService that should be used to generate
            random numbers from.
        dataset : instance of Dataset
            The Dataset instance describing the dataset for which background
            events should get generated.
        data : instance of DatasetData
            The DatasetData instance holding the data of the dataset for which
            background events should get generated.
        mean : float
            The mean number of background events to generate.
        tl : instance of TimeLord | None
            The optional instance of TimeLord that should be used to collect
            timing information about this method.
        **kwargs
            Additional keyword arguments, which might be required for a
            particular background generation method.

        Returns
        -------
        n_bkg : int
            The number of generated background events.
        bkg_events : DataFieldRecordArray
            The instance of DataFieldRecordArray holding the generated
            background events. The number of events in this array might be less
            than `n_bkg` if an event selection method was used for optimization
            purposes. The difference `n_bkg - len(bkg_events)` is then the
            number of pure background events in the generated background event
            sample.
        """
        pass


class MCDataSamplingBkgGenMethod(BackgroundGenerationMethod):
    """This class implements the method to generate background events from
    monte-carlo (MC) data by sampling events from the MC data set according to a
    probability value given for each event. Functions can be provided to get the
    mean number of background events and the probability of each monte-carlo
    event.
    """
    def __init__(
        self, get_event_prob_func, get_mean_func=None, unique_events=False,
        data_scrambler=None, mc_inplace_scrambling=False,
        keep_mc_data_fields=None, pre_event_selection_method=None):
        """Creates a new instance of the MCDataSamplingBkgGenMethod class.

        Parameters
        ----------
        get_event_prob_func : callable
            The function to get the background probability of each monte-carlo
            event. The call signature of this function must be
            `__call__(dataset, data, events)`,
            where `dataset` and `data` are `Dataset` and `DatasetData` instances
            of the data set for which background events needs to get generated.
            The `events` argument holds the actual set of events, for which the
            background event probabilities need to get calculated.
        get_mean_func : callable | None
            The function to get the mean number of background events.
            The call signature of this function must be
            `__call__(dataset, data, events)`,
            where `dataset` and `data` are `Dataset` and `DatasetData` instances
            of the data set for which background events needs to get generated.
            The `events` argument holds the actual set of events, for which the
            mean number of background events should get calculated.
            This argument can be `None`, which means that the mean number of
            background events to generate needs to get specified through the
            `generate_events` method. However, if an event selection method is
            provided, this argument cannot be None.
        unique_events : bool
            Flag if unique events should be drawn from the monte-carlo (True),
            or if events can be drawn several times (False). Default is False.
        data_scrambler : instance of DataScrambler | None
            If set to an instance of DataScrambler, the drawn monte-carlo
            background events will get scrambled. This can ensure more
            independent data trials. It is especially important when monte-carlo
            statistics are low.
        mc_inplace_scrambling : bool
            Flag if the scrambling of the monte-carlo data should be done
            inplace, i.e. without creating a copy of the MC data first.
            Default is False.
        keep_mc_data_fields : str | list of str | None
            The MC data field names that should be kept in order to be able to
            calculate the background events rates by the functions
            ``get_event_prob_func`` and ``get_mean_func``. All other MC fields
            will get droped due to computational efficiency reasons.
        pre_event_selection_method : instance of EventSelectionMethod | None
            If set to an instance of EventSelectionMethod, this method will
            pre-select the MC events that will be used for later background
            event generation. Using this pre-selection a large portion of the
            MC data can be reduced prior to background event generation.
        """
        super(MCDataSamplingBkgGenMethod, self).__init__()

        self.get_event_prob_func = get_event_prob_func
        self.get_mean_func = get_mean_func
        self.unique_events = unique_events
        self.data_scrambler = data_scrambler
        self.mc_inplace_scrambling = mc_inplace_scrambling
        self.keep_mc_data_field_names = keep_mc_data_fields
        self.pre_event_selection_method = pre_event_selection_method

        if((pre_event_selection_method is not None) and (get_mean_func is None)):
            raise ValueError('If an event pre-selection method is provided, a '
                             'get_mean_func needs to be provided as well!')

        # Define cache members to cache the background probabilities for each
        # monte-carlo event. The probabilities change only if the data changes.
        self._cache_data_id = None
        self._cache_mc_pre_selected = None
        self._cache_mc_event_bkg_prob = None
        self._cache_mc_event_bkg_prob_pre_selected = None
        self._cache_mean = None

    @property
    def get_event_prob_func(self):
        """The function to obtain the background probability for each
        monte-carlo event of the data set.
        """
        return self._get_event_prob_func
    @get_event_prob_func.setter
    def get_event_prob_func(self, func):
        if(not callable(func)):
            raise TypeError('The get_event_prob_func property must be a '
                'callable!')
        if(not func_has_n_args(func, 3)):
            raise TypeError('The function provided for the get_event_prob_func '
                'property must have 3 arguments!')
        self._get_event_prob_func = func

    @property
    def get_mean_func(self):
        """The function to obtain the mean number of background events for the
        data set. This can be None, which means that the mean number of
        background events to generate needs to be specified through the
        `generate_events` method.
        """
        return self._get_mean_func
    @get_mean_func.setter
    def get_mean_func(self, func):
        if(func is not None):
            if(not callable(func)):
                raise TypeError('The get_mean_func property must be a '
                    'callable!')
            if(not func_has_n_args(func, 3)):
                raise TypeError('The function provided for the get_mean_func '
                    'property must have 3 arguments!')
        self._get_mean_func = func

    @property
    def unique_events(self):
        """Flag if unique events should be drawn from the monto-carlo (True),
        or if the same event can be drawn multiple times from the monte-carlo.
        """
        return self._unique_events
    @unique_events.setter
    def unique_events(self, b):
        if(not isinstance(b, bool)):
            raise TypeError('The unique_events property must be of type bool!')
        self._unique_events = b

    @property
    def data_scrambler(self):
        """The DataScrambler instance that should be used to scramble the drawn
        monte-carlo background events to ensure more independent data trials.
        This is especially important when monte-carlo statistics are low. It is
        `None`, if no data scrambling should be used.
        """
        return self._data_scrambler
    @data_scrambler.setter
    def data_scrambler(self, scrambler):
        if(scrambler is not None):
            if(not isinstance(scrambler, DataScrambler)):
                raise TypeError('The data_scrambler property must be an instance '
                    'of DataScrambler!')
        self._data_scrambler = scrambler

    @property
    def mc_inplace_scrambling(self):
        """Flag if the scrambling of the monte-carlo data should be done
        inplace, i.e. without creating a copy of the MC data first.
        """
        return self._mc_inplace_scrambling
    @mc_inplace_scrambling.setter
    def mc_inplace_scrambling(self, b):
        if(not isinstance(b, bool)):
            raise TypeError('The mc_inplace_scrambling property must be of '
                'type bool!')
        self._mc_inplace_scrambling = b

    @property
    def keep_mc_data_field_names(self):
        """The MC data field names that should be kept in order to be able to
        calculate the background events rates by the functions
        ``get_event_prob_func`` and ``get_mean_func``. All other MC fields
        will get droped due to computational efficiency reasons.
        """
        return self._keep_mc_data_field_names
    @keep_mc_data_field_names.setter
    def keep_mc_data_field_names(self, names):
        if(names is None):
            names = []
        elif(isinstance(names, str)):
            names = [ names ]
        elif(not issequenceof(names, str)):
            raise TypeError('The keep_mc_data_field_names must be None, an '
                'instance of type str, or a sequence of objects of type str!')
        self._keep_mc_data_field_names = names

    @property
    def pre_event_selection_method(self):
        """The instance of EventSelectionMethod that pre-selects the MC events,
        which can be considered for background event generation.
        """
        return self._pre_event_selection_method
    @pre_event_selection_method.setter
    def pre_event_selection_method(self, method):
        if(method is not None):
            if(not isinstance(method, EventSelectionMethod)):
                raise TypeError('The pre_event_selection_method property must '
                    'be None, or an instance of EventSelectionMethod!')
            # If the event selection method selects all events, it's equivalent
            # to have it set to None, because then no operation has to be
            # performed.
            if(isinstance(method, AllEventSelectionMethod)):
                method = None
        self._pre_event_selection_method = method

    def change_source_hypo_group_manager(self, src_hypo_group_manager):
        """Changes the SourceHypoGroupManager instance of the
        pre-event-selection method. Also it invalides the data cache of this
        background generation method.

        Parameters
        ----------
        src_hypo_group_manager : SourceHypoGroupManager instance
            The new SourceHypoGroupManager instance.
        """
        if(self._pre_event_selection_method is not None):
            self._pre_event_selection_method.change_source_hypo_group_manager(
                src_hypo_group_manager)

        # Invalidate the data cache.
        self._cache_data_id = None

    def generate_events(
            self, rss, dataset, data, mean=None, poisson=True, tl=None):
        """Generates a `mean` number of background events for the given dataset
        and its data.

        Parameters
        ----------
        rss : instance of RandomStateService
            The instance of RandomStateService that should be used to generate
            random numbers from.
        dataset : instance of Dataset
            The Dataset instance describing the dataset for which background
            events should get generated.
        data : instance of DatasetData
            The DatasetData instance holding the data of the dataset for which
            background events should get generated.
        mean : float | None
            The mean number of background events to generate.
            Can be `None`. In that case the mean number of background events is
            obtained through the `get_mean_func` function.
        poisson : bool
            If set to True (default), the actual number of generated background
            events will be drawn from a Poisson distribution with the given mean
            value of background events.
            If set to False, the argument ``mean`` specifies the actual number
            of generated background events.
        tl : instance of TimeLord | None
            The optional instance of TimeLord that should be used to collect
            timing information about this method.

        Returns
        -------
        n_bkg : int
            The number of generated background events for the data set.
        bkg_events : instance of DataFieldRecordArray
            The instance of DataFieldRecordArray holding the generated
            background events. The number of events can be less than `n_bkg`
            if an event selection method is used.
        """
        tracing_enabled = CFG['debugging']['enable_tracing']

        # Create aliases to avoid dot-lookup.
        self__pre_event_selection_method = self._pre_event_selection_method

        # Check if the data set has changed. In that case need to get new
        # background probabilities for each monte-carlo event and a new mean
        # number of background events.
        data_id = id(data)
        if(self._cache_data_id != data_id):
            if(tracing_enabled):
                logger.debug(
                    f'DatasetData instance id of dataset "{dataset.name}" '
                    f'changed from {self._cache_data_id} to {data_id}')
            # Cache the current id of the data.
            self._cache_data_id = data_id

            # Create a copy of the MC data with all MC data fields removed,
            # except the specified MC data fields to keep for the
            # ``get_mean_func`` and ``get_event_prob_func`` functions.
            keep_field_names = list(set(
                CFG['dataset']['analysis_required_exp_field_names'] +
                data.exp_field_names +
                self._keep_mc_data_field_names
            ))
            data_mc = data.mc.copy(keep_fields=keep_field_names)

            if(self._get_mean_func is not None):
                with TaskTimer(tl, 'Calculate total MC background mean.'):
                    self._cache_mean = self._get_mean_func(
                        dataset, data, data_mc)

            with TaskTimer(tl, 'Calculate MC background event probability cache.'):
                self._cache_mc_event_bkg_prob = self._get_event_prob_func(
                    dataset, data, data_mc)

            if(self__pre_event_selection_method is not None):
                with TaskTimer(tl, 'Pre-select MC events.'):
                    (self._cache_mc_pre_selected,
                     mc_pre_selected_mask_idxs) =\
                    self__pre_event_selection_method.select_events(
                        data_mc, retidxs=True, tl=tl)
                self._cache_mc_event_bkg_prob_pre_selected = self._cache_mc_event_bkg_prob[mc_pre_selected_mask_idxs]
            else:
                self._cache_mc_pre_selected = data_mc


        if(mean is None):
            if(self._cache_mean is None):
                raise ValueError('No mean number of background events and no '
                    'get_mean_func were specified! One of the two must be '
                    'specified!')
            mean = self._cache_mean
        else:
            mean = float_cast(mean, 'The mean number of background events must '
                'be castable to type float!')

        # Draw the number of background events from a poisson distribution with
        # the given mean number of background events. This will be the number of
        # background events for this data set.
        n_bkg = (int(rss.random.poisson(mean)) if poisson else
                 int(np.round(mean, 0)))

        # Apply only event pre-selection before choosing events.
        data_mc_selected = self._cache_mc_pre_selected

        # Calculate the mean number of background events for the pre-selected
        # MC events.
        if(self__pre_event_selection_method is None):
            # No selection at all, use the total mean.
            mean_selected = mean
        else:
            with TaskTimer(tl, 'Calculate selected MC background mean.'):
                mean_selected = self._get_mean_func(
                    dataset, data, data_mc_selected)

        # Calculate the actual number of background events for the selected
        # events.
        p_binomial = mean_selected / mean
        with TaskTimer(tl, 'Get p array.'):
            if(self__pre_event_selection_method is None):
                p = self._cache_mc_event_bkg_prob
            else:
                # Pre-selection.
                p = self._cache_mc_event_bkg_prob_pre_selected / p_binomial
        n_bkg_selected = int(np.around(n_bkg * p_binomial, 0))

        # Draw the actual background events from the selected events of the
        # monto-carlo data set.
        with TaskTimer(tl, 'Draw MC background indices.'):
            bkg_event_indices = rss.random.choice(
                data_mc_selected.indices,
                size=n_bkg_selected,
                p=p,
                replace=(not self._unique_events))
        with TaskTimer(tl, 'Select MC background events from indices.'):
            bkg_events = data_mc_selected[bkg_event_indices]

        # Scramble the drawn MC events if requested.
        if(self._data_scrambler is not None):
            with TaskTimer(tl, 'Scramble MC background data.'):
                bkg_events = self._data_scrambler.scramble_data(
                    rss, bkg_events, copy=False)

        # Remove MC specific data fields from the background events record
        # array. So the result contains only experimental data fields. The list
        # of experimental data fields is defined as the unique set of the
        # required experimental data fields defined by the data set, and the
        # actual experimental data fields (in case there are additional kept
        # data fields by the user).
        with TaskTimer(tl, 'Remove MC specific data fields from MC events.'):
            exp_field_names = list(set(
                CFG['dataset']['analysis_required_exp_field_names'] +
                data.exp_field_names))
            bkg_events.tidy_up(exp_field_names)

        return (n_bkg, bkg_events)
