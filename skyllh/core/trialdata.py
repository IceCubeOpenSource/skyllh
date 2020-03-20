# -*- coding: utf-8 -*-

"""The trialdata module of SkyLLH provides a trial data manager class that
manages the data of an analysis trial. It provides also possible additional data
fields and their calculation, which are required by the particular analysis.
The rational behind this manager is to compute data fields only once, which can
then be used by different analysis objects, like PDF objects.
"""

import numpy as np

from skyllh.core.py import (
    func_has_n_args,
    issequenceof
)
from skyllh.core.storage import DataFieldRecordArray


class DataField(object):
    """This class defines a data field and its calculation that is used by an
    Analysis class instance. The calculation is defined through an external
    function.
    """

    def __init__(self, name, func, fitparam_names=None):
        """Creates a new instance of DataField that might depend on fit
        parameters.

        Parameters
        ----------
        name : str
            The name of the data field. It serves as the identifier for the
            data field.
        func : callable
            The function that calculates the values of this data field. The call
            signature must be
            `__call__(tdm, src_hypo_group_manager, fitparams)`,
            where `tdm` is the TrialDataManager instance holding the event data,
            `src_hypo_group_manager` is the SourceHypoGroupManager instance, and
            `fitparams` is the dictionary with the current fit parameter names
            and values. If the data field depends solely on source parameters,
            the call signature must be `__call__(tdm, src_hypo_group_manager)`
            instead.
        fitparam_names : sequence of str | None
            The sequence of str instances specifying the names of the fit
            parameters this data field depends on. If set to None, the data
            field does not depend on any fit parameters.
        """
        super(DataField, self).__init__()

        self.name = name
        self.func = func

        if(fitparam_names is None):
            fitparam_names = []
        if(not issequenceof(fitparam_names, str)):
            raise TypeError('The fitparam_names argument must be None or a '
                            'sequence of str instances!')
        self._fitparam_name_list = list(fitparam_names)

        # Define the list of fit parameter values for which the fit parameter
        # depend data field values have been calculated for.
        self._fitparam_value_list = [None]*len(self._fitparam_name_list)

        # Define the member variable that holds the numpy ndarray with the data
        # field values.
        self._values = None

        # Define the most efficient `calculate` method for this kind of data
        # field.
        if(func_has_n_args(self._func, 2)):
            self.calculate = self._calc_source_values
        elif(len(self._fitparam_name_list) == 0):
            self.calculate = self._calc_static_values
        else:
            self.calculate = self._calc_fitparam_dependent_values

    @property
    def name(self):
        """The name of the data field.
        """
        return self._name

    @name.setter
    def name(self, name):
        if(not isinstance(name, str)):
            raise TypeError('The name property must be an instance of str!')
        self._name = name

    @property
    def func(self):
        """The function that calculates the data field values.
        """
        return self._func

    @func.setter
    def func(self, f):
        if(not callable(f)):
            raise TypeError('The func property must be a callable object!')
        if((not func_has_n_args(f, 2)) and
           (not func_has_n_args(f, 3))):
            raise TypeError('The func property must be a function with 2 or 3 '
                            'arguments!')
        self._func = f

    @property
    def values(self):
        """(read-only) The calculated data values of the data field.
        """
        return self._values

    def _calc_source_values(
            self, tdm, src_hypo_group_manager, fitparams):
        """Calculates the data field values utilizing the defined external
        function. The data field values solely depend on source parameters.
        """
        self._values = self._func(tdm, src_hypo_group_manager)

    def _calc_static_values(
            self, tdm, src_hypo_group_manager, fitparams):
        """Calculates the data field values utilizing the defined external
        function, that are static and only depend on source parameters.

        Parameters
        ----------
        tdm : instance of TrialDataManager
            The TrialDataManager instance this data field is part of and is
            holding the event data.
        src_hypo_group_manager : instance of SourceHypoGroupManager
            The instance of SourceHypoGroupManager, which defines the groups of
            source hypotheses.
        fitparams : dict
            The dictionary holding the current fit parameter names and values.
            By definition this dictionary is empty.
        """
        self._values = self._func(tdm, src_hypo_group_manager, fitparams)

    def _calc_fitparam_dependent_values(
            self, tdm, src_hypo_group_manager, fitparams):
        """Calculate data field values utilizing the defined external
        function, that depend on fit parameter values. We check if the fit
        parameter values have changed.

        Parameters
        ----------
        tdm : instance of TrialDataManager
            The TrialDataManager instance this data field is part of and is
            holding the event data.
        src_hypo_group_manager : instance of SourceHypoGroupManager
            The instance of SourceHypoGroupManager, which defines the groups of
            source hypotheses.
        fitparams : dict
            The dictionary holding the current fit parameter names and values.
        """
        if(self._values is None):
            # It's the first time this method is called, so we need to calculate
            # the data field values for sure.
            self._values = self._func(tdm, src_hypo_group_manager, fitparams)
            return

        for (idx, fitparam_name) in enumerate(self._fitparam_name_list):
            if(fitparams[fitparam_name] != self._fitparam_value_list[idx]):
                # This current fit parameter value has changed. So we need to
                # re-calculate the data field values.
                self._fitparam_value_list = [
                    fitparams[name] for name in self._fitparam_name_list
                ]
                self._values = self._func(
                    tdm, src_hypo_group_manager, fitparams)
                break


class TrialDataManager(object):
    """The TrialDataManager class manages the event data for an analysis trial.
    It provides possible additional data fields and their calculation. New data fields can be defined via the `add_data_field` method. Whenever a new trial
    is being initialized the data fields get re-calculated. The data trial
    manager is provided to the PDF evaluation method. Hence, data fields are
    calculated only once.
    """

    def __init__(self, index_field_name=None):
        """Creates a new TrialDataManager instance.

        Parameters
        ----------
        index_field_name : str | None
            The name of the field that should be used as primary index field.
            If provided, the events will be sorted along this data field. This
            might be useful for run-time performance.
        """
        super(TrialDataManager, self).__init__()

        self.index_field_name = index_field_name

        # Define the list of data fields that depend only on the source
        # parameters.
        self._source_data_fields = []
        self._source_data_field_reg = dict()

        # Define the list of data fields that are static, i.e. don't depend on
        # any fit parameters. These fields have to be calculated only once when
        # a new evaluation data is available.
        self._static_data_fields = []
        self._static_data_field_reg = dict()

        # Define the list of data fields that depend on fit parameters. These
        # data fields have to be re-calculated whenever a fit parameter value
        # changes.
        self._fitparam_data_fields = []
        self._fitparam_data_field_reg = dict()

        # Define the member variable that will hold the raw events for which the
        # data fields get calculated.
        self._events = None

        # We store an integer number for the trial data state and increase it
        # whenever the state of the trial data changed. This way other code,
        # e.g. PDFs, can determine when the data changed and internal caches
        # must be flushed.
        self._trial_data_state_id = -1

    @property
    def index_field_name(self):
        """The name of the primary index data field. If not None, events will
        be sorted by this data field.
        """
        return self._index_field_name

    @index_field_name.setter
    def index_field_name(self, name):
        if(name is not None):
            if(not isinstance(name, str)):
                raise TypeError('The index_field_name property must be an '
                                'instance of type str!')
        self._index_field_name = name

    @property
    def events(self):
        """The DataFieldRecordArray instance holding the data events, which
        should get evaluated.
        """
        return self._events

    @events.setter
    def events(self, arr):
        if(not isinstance(arr, DataFieldRecordArray)):
            raise TypeError('The events property must be an instance of '
                            'DataFieldRecordArray!')
        self._events = arr

    @property
    def n_events(self):
        """(read-only) The number of events which should get evaluated.
        """
        return len(self._events)

    @property
    def trial_data_state_id(self):
        """(read-only) The integer ID number of the trial data. This ID number
        can be used to determine when the trial data has changed its state.
        """
        return self._trial_data_state_id

    def __contains__(self, name):
        """Checks if the given data field is defined in this data field manager.

        Parameters
        ----------
        name : str
            The name of the data field.

        Returns
        -------
        check : bool
            True if the data field is defined in this data field manager,
            False otherwise.
        """
        if((name in self._source_data_field_reg) or
           (name in self._static_data_field_reg) or
           (name in self._fitparam_data_field_reg)):
            return True

        return False

    def change_source_hypo_group_manager(self, src_hypo_group_manager):
        """Recalculate the source data fields.

        Parameters
        ----------
        src_hypo_group_manager : instance of SourceHypoGroupManager
            The SourceHypoGroupManager manager that defines the groups of
            source hypotheses.
        """
        self.calculate_source_data_fields(src_hypo_group_manager)

    def initialize_for_new_trial(self, src_hypo_group_manager, events):
        """Initializes the trial data manager for a new trial. It sets the raw
        events and calculates the static data fields.

        Parameters
        ----------
        src_hypo_group_manager : SourceHypoGroupManager instance
            The instance of SourceHypoGroupManager that defines the source
            hypothesis groups.
        events : DataFieldRecordArray instance
            The DataFieldRecordArray instance holding the raw events.
        """
        # Sort the events by the index field, if a field was provided.
        if(self._index_field_name is not None):
            events.sort_by_field(self._index_field_name)

        # Set the events property, so that the calculation functions of the data
        # fields can access them.
        self.events = events

        # Now calculate all the static data fields.
        self.calculate_static_data_fields(src_hypo_group_manager)

        # Increment the trial data state ID.
        self._trial_data_state_id += 1

    def add_source_data_field(self, name, func):
        """Adds a new data field to the manager. The data field must depend
        solely on source parameters.

        Parameters
        ----------
        name : str
            The name of the data field. It serves as the identifier for the
            data field.
        func : callable
            The function that calculates the data field values. The call
            signature must be
            `__call__(tdm, src_hypo_group_manager, fitparams)`, where
            `tdm` is the TrialDataManager instance holding the event data,
            `src_hypo_group_manager` is the SourceHypoGroupManager instance,
            and `fitparams` is an unused interface argument.
        """
        if(name in self):
            raise KeyError('The data field "%s" is already defined!' % (name))

        data_field = DataField(name, func)

        self._source_data_fields.append(data_field)
        self._source_data_field_reg[name] = data_field

    def add_data_field(self, name, func, fitparam_names=None):
        """Adds a new data field to the manager.

        Parameters
        ----------
        name : str
            The name of the data field. It serves as the identifier for the
            data field.
        func : callable
            The function that calculates the data field values. The call
            signature must be
            `__call__(tdm, src_hypo_group_manager, fitparams)`, where
            `tdm` is the TrialDataManager instance holding the event data,
            `src_hypo_group_manager` is the SourceHypoGroupManager instance,
            and `fitparams` is the dictionary with the current fit parameter
            names and values.
        fitparam_names : sequence of str | None
            The sequence of str instances specifying the names of the fit
            parameters this data field depends on. If set to None, it means that
            the data field does not depend on any fit parameters.
        """
        if(name in self):
            raise KeyError('The data field "%s" is already defined!' % (name))

        data_field = DataField(name, func, fitparam_names)

        if(fitparam_names is not None):
            self._fitparam_data_fields.append(data_field)
            self._fitparam_data_field_reg[name] = data_field
        else:
            self._static_data_fields.append(data_field)
            self._static_data_field_reg[name] = data_field

    def calculate_source_data_fields(self, src_hypo_group_manager):
        """Calculates the data values of the data fields that solely depend on
        source parameters.

        Parameters
        ----------
        src_hypo_group_manager : instance of SourceHypoGroupManager
            The instance of SourceHypoGroupManager, which defines the groups of
            source hypotheses.
        """
        if(len(self._source_data_fields) == 0):
            return

        fitparams = None
        for data_field in self._source_data_fields:
            data_field.calculate(self, src_hypo_group_manager, fitparams)

        self._trial_data_state_id += 1

    def calculate_static_data_fields(self, src_hypo_group_manager):
        """Calculates the data values of the data fields that do not depend on
        any source or fit parameters.

        Parameters
        ----------
        src_hypo_group_manager : instance of SourceHypoGroupManager
            The instance of SourceHypoGroupManager, which defines the groups of
            source hypotheses.
        """
        if(len(self._static_data_fields) == 0):
            return

        fitparams = dict()
        for data_field in self._static_data_fields:
            data_field.calculate(self, src_hypo_group_manager, fitparams)

        self._trial_data_state_id += 1

    def calculate_fitparam_data_fields(self, src_hypo_group_manager, fitparams):
        """Calculates the data values of the data fields that depend on fit
        parameter values.

        Parameters
        ----------
        src_hypo_group_manager : instance of SourceHypoGroupManager
            The instance of SourceHypoGroupManager, which defines the groups of
            source hypotheses.
        fitparams : dict
            The dictionary holding the fit parameter names and values.
        """
        if(len(self._fitparam_data_fields) == 0):
            return

        for data_field in self._fitparam_data_fields:
            data_field.calculate(self, src_hypo_group_manager, fitparams)

        self._trial_data_state_id += 1

    def get_data(self, name):
        """Gets the data for the given data field name. The data is stored
        either in the raw events record ndarray or in one of the additional
        defined data fields. Data from the raw events record ndarray is
        prefered.

        Parameters
        ----------
        name : str
            The name of the data field for which to retrieve the data.

        Returns
        -------
        data : numpy ndarray
            The data of the requested data field.

        Raises
        ------
        KeyError
            If the given data field is not defined.
        """
        if(name in self._events.field_name_list):
            return self._events[name]
        if(name in self._source_data_field_reg):
            return self._source_data_field_reg[name].values
        if(name in self._static_data_field_reg):
            return self._static_data_field_reg[name].values
        if(name in self._fitparam_data_field_reg):
            return self._fitparam_data_field_reg[name].values

        raise KeyError('The data field "%s" is not defined!' % (name))
