# -*- coding: utf-8 -*-

from collections import (
    defaultdict,
)
import numpy as np

from skyllh.core.dataset import (
    Dataset,
    DatasetData,
)
from skyllh.core.progressbar import (
    ProgressBar,
)
from skyllh.core.py import (
    classname,
    issequenceof,
)
from skyllh.core.source_hypo_grouping import (
    SourceHypoGroupManager,
)


class DetSigYieldService(
        object):
    """This class provides a service to build and hold detector signal yield
    instances for multiple datasets and source hypothesis groups.
    """

    def __init__(
            self,
            shg_mgr,
            dataset_list,
            data_list,
            ppbar=None,
            **kwargs):
        """Creates a new DetSigYieldService instance.
        """
        super().__init__(
            **kwargs)

        self._set_shg_mgr(shg_mgr)

        self.dataset_list = dataset_list
        self.data_list = data_list

        self._arr = self.construct_detsigyield_array(
            ppbar=ppbar)

    @property
    def shg_mgr(self):
        """(read-only) The instance of SourceHypoGroupManager providing the list
        of source hypothesis groups.
        """
        return self._shg_mgr

    @property
    def dataset_list(self):
        """The list of instance of Dataset for which the detector signal yields
        should be built.
        """
        return self._dataset_list

    @dataset_list.setter
    def dataset_list(self, datasets):
        if not issequenceof(datasets, Dataset):
            raise TypeError(
                'The dataset_list property must be a sequence of Dataset '
                'instances! '
                f'Its current type is {classname(datasets)}!')
        self._dataset_list = list(datasets)

    @property
    def data_list(self):
        """The list of instance of DatasetData for which the detector signal
        yields should be built.
        """
        return self._data_list

    @data_list.setter
    def data_list(self, datas):
        if not issequenceof(datas, DatasetData):
            raise TypeError(
                'The data_list property must be a sequence of DatasetData '
                'instances! '
                f'Its current type is {classname(datas)}!')
        self._data_list = list(datas)

    @property
    def arr(self):
        """(read-only) The (N_datasets, N_source_hypo_groups)-shaped numpy
        ndarray of object holding the constructed DetSigYield instances.
        """
        return self._arr

    @property
    def n_datasets(self):
        """(read-only) The number of datasets this service was created for.
        """
        return self._arr.shape[0]

    @property
    def n_shgs(self):
        """(read-only) The number of source hypothesis groups this service was
        created for.
        """
        return self._arr.shape[1]

    def _set_shg_mgr(self, mgr):
        """Sets the internal member variable to the given instance of
        SourceHypoGroupManager.
        """
        if not isinstance(mgr, SourceHypoGroupManager):
            raise TypeError(
                'The shg_mgr argument must be an instance of '
                'SourceHypoGroupManager! '
                f'Its current type is {classname(mgr)}!')

        self._shg_mgr = mgr

    def change_shg_mgr(
            self,
            shg_mgr,
            ppbar=None,
    ):
        """Changes the instance of SourceHypoGroupManager of this service. This
        will also rebuild the detector signal yields.
        """
        self._set_shg_mgr(shg_mgr)

        self._arr = self.construct_detsigyield_array(
            ppbar=ppbar)

    def construct_detsigyield_array(
            self,
            ppbar=None,
    ):
        """Creates a (N_datasets, N_source_hypo_groups)-shaped numpy ndarray of
        object holding the constructed DetSigYield instances.

        Parameters
        ----------
        ppbar : instance of ProgressBar | None
            The instance of ProgressBar of the optional parent progress bar.

        Returns
        -------
        detsigyield_arr : instance of numpy ndarray
            The (N_datasets, N_source_hypo_groups)-shaped numpy ndarray of
            object holding the constructed DetSigYield instances.
        """
        n_datasets = len(self._dataset_list)

        detsigyield_arr = np.empty(
            (n_datasets,
             self._shg_mgr.n_src_hypo_groups),
            dtype=object
        )

        pbar = ProgressBar(
            self._shg_mgr.n_src_hypo_groups * n_datasets,
            parent=ppbar).start()
        for (g, shg) in enumerate(self._shg_mgr.shg_list):
            fluxmodel = shg.fluxmodel
            detsigyield_builder_list = shg.detsigyield_builder_list

            if (len(detsigyield_builder_list) != 1) and\
               (len(detsigyield_builder_list) != n_datasets):
                raise ValueError(
                    'The number of detector signal yield builders '
                    f'({len(detsigyield_builder_list)}) is not 1 and does not '
                    f'match the number of datasets ({n_datasets})!')

            for (j, (dataset, data)) in enumerate(zip(self._dataset_list,
                                                      self._data_list)):
                if len(detsigyield_builder_list) == 1:
                    # Only one detsigyield builder was defined,
                    # so we use it for all datasets.
                    detsigyield_builder = detsigyield_builder_list[0]
                else:
                    detsigyield_builder = detsigyield_builder_list[j]

                detsigyield = detsigyield_builder.construct_detsigyield(
                    dataset=dataset,
                    data=data,
                    fluxmodel=fluxmodel,
                    livetime=data.livetime,
                    ppbar=pbar)
                detsigyield_arr[j, g] = detsigyield

                pbar.increment()
        pbar.finish()

        return detsigyield_arr


class SrcDetSigYieldWeightsService(
        object):
    r"""This class provides a service for the source detector signal yield
    weights, which are the product of the source weights with the detector
    signal yield, denoted with :math:`a_{j,k}(\vec{p}_{\mathrm{s}_k})` in the
    math formalism documentation.

    .. math::

        a_{j,k}(\vec{p}_{\mathrm{s}_k}) = W_k
            \mathcal{Y}_{\mathrm{s}_{j,k}}(\vec{p}_{\mathrm{s}_k})

    The service has a method to calculate the weights and a method to retrieve
    the weights. The weights are stored internally.
    """

    @staticmethod
    def create_src_recarray_list_list(
            detsigyield_service,
    ):
        """Creates a list of numpy record ndarrays, one for each source
        hypothesis group suited for evaluating the detector signal yield
        instance of that source hypothesis group.

        Parameters
        ----------
        detsigyield_service : instance of DetSigYieldService
            The instance of DetSigYieldService providing the
            (N_datasets, N_source_hypo_groups)-shaped 2D ndarray of
            DetSigYield instances, one for each dataset and source hypothesis
            group combination.

        Returns
        -------
        src_recarray_list_list : list of list of numpy record ndarrays
            The (N_datasets,N_source_hypo_groups)-shaped list of list of the
            source numpy record ndarrays, one for each dataset and source
            hypothesis group combination, which is needed for
            evaluating a particular detector signal yield instance.
        """
        n_datasets = detsigyield_service.n_datasets
        n_shgs = detsigyield_service.n_shgs
        shg_list = detsigyield_service.shg_mgr.shg_list

        src_recarray_list_list = []
        for ds_idx in range(n_datasets):
            src_recarray_list = []
            for shg_idx in range(n_shgs):
                shg = shg_list[shg_idx]
                src_recarray_list.append(
                    detsigyield_service.arr[ds_idx][shg_idx].sources_to_recarray(
                        shg.source_list))

            src_recarray_list_list.append(src_recarray_list)

        return src_recarray_list_list

    @staticmethod
    def create_src_weight_array_list(
            shg_mgr,
    ):
        """Creates a list of numpy 1D ndarrays holding the source weights, one
        for each source hypothesis group.

        Parameters
        ----------
        shg_mgr : instance of SourceHypoGroupManager
            The instance of SourceHypoGroupManager defining the source
            hypothesis groups with their sources.

        Returns
        -------
        src_weight_array_list : list of numpy 1D ndarrays
            The list of 1D numpy ndarrays holding the source weights, one for
            each source hypothesis group.
        """
        src_weight_array_list = [
            np.array([src.weight for src in shg.source_list])
            for shg in shg_mgr.shg_list
        ]
        return src_weight_array_list

    def __init__(
            self,
            detsigyield_service,
            **kwargs,
    ):
        """Creates a new SrcDetSigYieldWeightsService instance.

        Parameters
        ----------
        detsigyield_service : instance of DetSigYieldService
            The instance of DetSigYieldService providing the
            (N_datasets, N_source_hypo_groups)-shaped array of DetSigYield
            instances, one instance for each combination of dataset and source
            hypothesis group.
        """
        super().__init__(
            **kwargs)

        self.detsigyield_service = detsigyield_service

        # Create the list of list of source record arrays for each combination
        # of dataset and source hypothesis group.
        self._src_recarray_list_list = type(self).create_src_recarray_list_list(
            detsigyield_service=self._detsigyield_service)

        # Create the list of 1D ndarrays holding the source weights for each
        # source hypothesis group.
        self._src_weight_array_list = type(self).create_src_weight_array_list(
            shg_mgr=self._detsigyield_service.shg_mgr)

        self._a_jk = None
        self._a_jk_grads = None

    @property
    def shg_mgr(self):
        """(read-only) The instance of SourceHypoGroupManager defining the
        source hypothesis groups.
        """
        return self._detsigyield_service.shg_mgr

    @property
    def detsigyield_service(self):
        """The instance of DetSigYieldService providing the
        (N_datasets, N_source_hypo_groups)-shaped array of DetSigYield
        instances.
        """
        return self._detsigyield_service

    @detsigyield_service.setter
    def detsigyield_service(self, service):
        if not isinstance(service, DetSigYieldService):
            raise TypeError(
                'The detsigyield_service property must be an instance of '
                'DetSigYieldService! '
                f'Its current type is {classname(service)}!')
        self._detsigyield_service = service

    @property
    def detsigyield_arr(self):
        """(read-only) The (N_datasets, N_source_hypo_groups)-shaped 2D numpy
        ndarray holding the DetSigYield instances for each source hypothesis
        group.
        """
        return self._detsigyield_service.arr

    @property
    def n_datasets(self):
        """(read-only) The number of datasets this service was created for.
        """
        return self._detsigyield_service.n_datasets

    @property
    def n_shgs(self):
        """(read-only) The number of source hypothesis groups this service was
        created for.
        """
        return self._detsigyield_service.n_shgs

    @property
    def src_recarray_list_list(self):
        """(read-only) The (N_datasets,N_source_hypo_groups)-shaped list of list
        of the source numpy record ndarrays, one for each dataset and source
        hypothesis group combination, which is needed for evaluating a
        particular detector signal yield instance.
        """
        return self._src_recarray_list_list

    def change_shg_mgr(
            self,
            shg_mgr,
    ):
        """Re-creates the internal source numpy record arrays needed for the
        detector signal yield calculation.

        Parameters
        ----------
        shg_mgr : instance of SourceHypoGroupManager
            The new SourceHypoGroupManager instance.
        """
        if id(shg_mgr) != id(self._detsigyield_service.shg_mgr):
            raise ValueError(
                'The provides instance of SourceHypoGroupManager does not '
                'match the instance of the detector signal yield service!')

        self._src_recarray_list_list = type(self).create_src_recarray_list_list(
            detsigyield_service=self._detsigyield_service)

        self._src_weight_array_list = type(self).create_src_weight_array_list(
            shg_mgr=self._detsigyield_service.shg_mgr)

    def calculate(
            self,
            src_params_recarray):
        """Calculates the source detector signal yield weights for each source
        and their derivative w.r.t. each global floating parameter. The result
        is stored internally as:

            a_jk : instance of ndarray
                The (N_datasets,N_sources)-shaped numpy ndarray holding the
                source detector signal yield weight for each combination of
                dataset and source.
            a_jk_grads : dict
                The dictionary holding the (N_datasets,N_sources)-shaped numpy
                ndarray with the derivatives w.r.t. the global fit parameter
                the SrcDetSigYieldWeightsService depend on. The dictionary's key
                is the index of the global fit parameter.

        Parameters
        ----------
        src_params_recarray : instance of numpy record ndarray
            The numpy record ndarray of length N_sources holding the local
            source parameters. See the documentation of
            :meth:`skyllh.core.parameters.ParameterModelMapper.create_src_params_recarray`
            for more information about this record array.
        """
        n_datasets = self.n_datasets

        shg_mgr = self._detsigyield_service.shg_mgr

        self._a_jk = np.empty(
            (n_datasets, shg_mgr.n_sources,),
            dtype=np.double)

        self._a_jk_grads = defaultdict(
            lambda: np.zeros(
                (n_datasets, shg_mgr.n_sources),
                dtype=np.double))

        sidx = 0
        for (shg_idx, (shg, src_weights)) in enumerate(zip(
                shg_mgr.shg_list, self._src_weight_array_list)):

            shg_n_src = shg.n_sources

            shg_src_slice = slice(sidx, sidx+shg_n_src)

            shg_src_params_recarray = src_params_recarray[shg_src_slice]

            for ds_idx in range(n_datasets):
                detsigyield = self._detsigyield_service.arr[ds_idx, shg_idx]
                src_recarray = self._src_recarray_list_list[ds_idx][shg_idx]

                (Yg, Yg_grads) = detsigyield(
                    src_recarray=src_recarray,
                    src_params_recarray=shg_src_params_recarray)

                self._a_jk[ds_idx][shg_src_slice] = src_weights * Yg

                for gpidx in Yg_grads.keys():
                    self._a_jk_grads[gpidx][ds_idx, shg_src_slice] =\
                        src_weights * Yg_grads[gpidx]

            sidx += shg_n_src

    def get_weights(self):
        """Returns the source detector signal yield weights and their
        derivatives w.r.t. the global fit parameters.

        Returns
        -------
        a_jk : instance of ndarray
            The (N_datasets, N_sources)-shaped numpy ndarray holding the
            source detector signal yield weight for each combination of
            dataset and source.
        a_jk_grads : dict
            The dictionary holding the (N_datasets, N_sources)-shaped numpy
            ndarray with the derivatives w.r.t. the global fit parameter
            the SrcDetSigYieldWeightsService depend on. The dictionary's key
            is the index of the global fit parameter.
        """
        return (self._a_jk, self._a_jk_grads)
