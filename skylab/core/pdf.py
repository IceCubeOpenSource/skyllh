# -*- coding: utf-8 -*-

import abc

from skylab.core.py import typename
from skylab.core.parameters import ParameterGridSet, make_params_hash

class PDF(object):
    """The abstract base class for all probability distribution functions (PDF)
    models.
    All PDF model classes must be derived from this class.
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, *args, **kwargs):
        # Make sure that multiple inheritance can be used. This super call will
        # invoke the __init__ method of a possible second inheritance.
        super(PDF, self).__init__(*args, **kwargs)

    @abc.abstractmethod
    def assert_is_valid_for_exp_data(self, data_exp):
        """This abstract method is supposed to check if this PDF is valid for
        all the given experimental data. This means, it needs to check if there
        is a PDF value for each data event that will be used in the likelihood
        evaluation. This is just a seatbelt. The method must raise a ValueError
        if the PDF is not valid for the given experimental data.

        Parameters
        ----------
        data_exp : numpy record ndarray
            The array holding the experimental data to check for.

        Errors
        ------
        ValueError
            If the PDF is not valid for the given experimental data.
        """
        pass

    @abc.abstractmethod
    def get_prob(self, events, params):
        """This abstract method is supposed to calculate the probability for
        the specified events given the specified parameters.

        Parameters
        ----------
        events : numpy record ndarray
            The numpy record ndarray holding the data events for which the
            probability should be calculated for. What data fields are required
            is defined by the derived PDF class and depends on the application.
        params : dict
            The dictionary containing the parameters needed for the calculation.
            What parameters these are is defined by derived PDF class and
            depends on the application.
        """
        pass


class SpatialPDF(PDF):
    """This is the abstract base class for a spatial PDF model.
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, *args, **kwargs):
        super(SpatialPDF, self).__init__(*args, **kwargs)


class EnergyPDF(PDF):
    """This is the abstract base class for an energy PDF model.
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, *args, **kwargs):
        super(EnergyPDF, self).__init__(*args, **kwargs)


class TimePDF(PDF):
    """This is the abstract base class for a time PDF model.
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, *args, **kwargs):
        super(TimePDF, self).__init__(*args, **kwargs)


class PDFSet(object):
    """This class describes a set of PDF objects which are related to each other
    via different values of a set of fit parameters. A signal PDF usually
    consists of multiple same-kind PDFs for different signal fit parameters.
    In general background PDFs could have fit parameters, too.

    This class has the ``fitparams_gid_set`` property holding the set of fit
    parameter grids. Also it holds a dictionary with the PDFs for the different
    sets of fit parameter values. The type of the PDF objects is defined through
    the ``pdf_type`` property. PDF objects of type ``pdf_type`` can be added
    via the ``add_pdf`` method and retrieved via the ``get_pdf`` method.
    """
    def __init__(self, pdf_type, *args, **kwargs):
        """Constructor method. Gets called when the an instance of a class is
        created which derives from this PDFSet class.

        Parameters
        ----------
        pdf_type : type
            The PDF class that can be added to the set.
        """
        if(not issubclass(pdf_type, PDF)):
            raise TypeError('The pdf_type argument must be a subclass of PDF!')
        self._pdf_type = pdf_type
        self.fitparams_grid_set = ParameterGridSet()
        self._params_hash_pdf_dict = dict()

    @property
    def pdf_type(self):
        """(read-only) The PDF type which can be added to the PDF set.
        """
        return self._pdf_type

    @property
    def fitparams_grid_set(self):
        """The ParameterGridSet object defining the value grids of the different
        fit parameters.
        """
        return self._fitparams_grid_set
    @fitparams_grid_set.setter
    def fitparams_grid_set(self, obj):
        if(not isinstance(obj, ParameterGridSet)):
            raise TypeError('The fitparams_gid_set property must be an object of type ParameterGridSet!')
        self._fitparams_grid_set = obj

    @property
    def gridfitparams_list(self):
        """(read-only) The list of dictionaries of all the fit parameter
        permutations on the grid.
        """
        return self.fitparams_grid_set.parameter_permutation_dict_list

    @property
    def pdf_keys(self):
        """(read-only) The list of stored PDF object keys.
        """
        return self._params_hash_pdf_dict.keys()

    def items(self):
        """Returns the list of 2-element tuples for the PDF stored in this
        PDFSet object.
        """
        return self._params_hash_pdf_dict.items()

    def add_pdf(self, pdf, params):
        """Adds the given PDF object for the given parameters to the internal
        registry.

        Parameters
        ----------
        pdf : pdf_type
            The object derived from ``pdf_type`` that should be added.
        params : dict
            The dictionary with the PDF parameters, which identify the PDF
            object.

        Errors
        ------
        KeyError
            If the given PDF was already added for the given set of parameters.
        """
        if(not isinstance(pdf, self.pdf_type)):
            raise TypeError('The pdf argument must be an instance of %s!'%(typename(self.pdf_type)))
        if(not isinstance(params, dict)):
            raise TypeError('The params argument must be of type dict!')

        params_hash = make_params_hash(params)
        if(params_hash in self._params_hash_pdf_dict):
            raise KeyError('The signal PDF with parameters %s was already added!'%(str(params)))
        self._params_hash_pdf_dict[params_hash] = signalpdf

    def get_pdf(self, params):
        """Retrieves the PDF object for the given set of signal
        parameters.

        Parameters
        ----------
        params : dict | int
            The dictionary with the parameters for which the PDF object should
            get retrieved. If an integer is given, it is assumed to be the PDF
            key.

        Returns
        -------
        pdf : pdf_type
            The pdf_type object for the given parameters.

        Errors
        ------
        KeyError
            If no PDF object was created for the given set of parameters.
        """
        if(isinstance(params, int)):
            params_hash = params
        elif(isinstance(params, dict)):
            params_hash = make_params_hash(params)
        else:
            raise TypeError('The params argument must be of type dict or int!')

        if(params_hash not in self._params_hash_pdf_dict):
            raise KeyError('No PDF was created for the parameter set "%s"!'%(str(params)))

        pdf = self._params_hash_pdf_dict[params_hash]
        return pdf


class IsBackgroundPDF(object):
    """This is a classifier class that can be used by other classes to indicate
    that the class describes a background PDF. This is useful for type checking.
    """
    def __init__(self, *args, **kwargs):
        """Constructor method. Gets called when the an instance of a class is
        created which derives from this IsBackgroundPDF class.
        """
        super(IsBackgroundPDF, self).__init__(*args, **kwargs)


class IsSignalPDF(object):
    """This is a classifier class that can be used by other classes to indicate
    that the class describes a signal PDF.
    """
    def __init__(self, *args, **kwargs):
        """Constructor method. Gets called when the an instance of a class is
        created which derives from this IsSignalPDF class.
        """
        super(IsSignalPDF, self).__init__(*args, **kwargs)

