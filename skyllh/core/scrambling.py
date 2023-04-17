# -*- coding: utf-8 -*-

import abc

import numpy as np

from skyllh.core.times import TimeGenerator


class DataScramblingMethod(object, metaclass=abc.ABCMeta):
    """Base class (type) for implementing a data scrambling method.
    """

    def __init__(self):
        super(DataScramblingMethod, self).__init__()

    @abc.abstractmethod
    def scramble(self, rss, data):
        """The scramble method implements the actual scrambling of the given
        data, which is method dependent. The scrambling must be performed
        in-place, i.e. it alters the data inside the given data array.

        Parameters
        ----------
        rss : RandomStateService
            The random state service providing the random number
            generator (RNG).
        data : instance of DataFieldRecordArray
            The DataFieldRecordArray containing the to be scrambled data.

        Returns
        -------
        data : DataFieldRecordArray
            The given DataFieldRecordArray holding the scrambled data.
        """
        pass


class UniformRAScramblingMethod(DataScramblingMethod):
    """The UniformRAScramblingMethod method performs right-ascention scrambling
    uniformly within a given RA range. By default it's (0, 2\pi).

    Note: This alters only the ``ra`` values of the data!
    """
    def __init__(self, ra_range=None):
        """Initializes a new RAScramblingMethod instance.

        Parameters
        ----------
        ra_range : tuple | None
            The two-element tuple holding the range in radians within the RA
            values should get drawn from. If set to None, the default (0, 2\pi)
            will be used.
        """
        super(UniformRAScramblingMethod, self).__init__()

        self.ra_range = ra_range

    @property
    def ra_range(self):
        """The two-element tuple holding the range within the RA values
        should get drawn from.
        """
        return self._ra_range
    @ra_range.setter
    def ra_range(self, ra_range):
        if(ra_range is None):
            ra_range = (0, 2*np.pi)
        if(not isinstance(ra_range, tuple)):
            raise TypeError('The ra_range property must be a tuple!')
        if(len(ra_range) != 2):
            raise ValueError('The ra_range tuple must contain 2 elements!')
        self._ra_range = ra_range

    def scramble(self, rss, data):
        """Scrambles the given data uniformly in right-ascention.

        Parameters
        ----------
        rss : RandomStateService
            The random state service providing the random number
            generator (RNG).
        data : instance of DataFieldRecordArray
            The DataFieldRecordArray instance containing the to be scrambled
            data.

        Returns
        -------
        data : DataFieldRecordArray
            The given DataFieldRecordArray holding the scrambled data.
        """
        dt = data['ra'].dtype
        data['ra'] = rss.random.uniform(
            *self.ra_range, size=len(data)).astype(dt)
        return data


class TimeScramblingMethod(DataScramblingMethod):
    """The TimeScramblingMethod class provides a data scrambling method to
    perform data coordinate scrambling based on a generated time. It draws a
    random time from a time generator and transforms the horizontal (local)
    coordinates into equatorial coordinates using a specified transformation
    function.
    """
    def __init__(self, timegen, hor_to_equ_transform):
        """Initializes a new time scramling method instance.

        Parameters
        ----------
        timegen : TimeGenerator
            The time generator that should be used to generate random MJD times.
        hor_to_equ_transform : callable
            The transformation function to transform coordinates from the
            horizontal system into the equatorial system.

            The call signature must be:

                __call__(azi, zen, mjd)

            The return signature must be: (ra, dec)

        """
        super(TimeScramblingMethod, self).__init__()

        self.timegen = timegen
        self.hor_to_equ_transform = hor_to_equ_transform

    @property
    def timegen(self):
        """The TimeGenerator instance that should be used to generate random MJD
        times.
        """
        return self._timegen
    @timegen.setter
    def timegen(self, timegen):
        if(not isinstance(timegen, TimeGenerator)):
            raise TypeError('The timegen property must be an instance of TimeGenerator!')
        self._timegen = timegen

    @property
    def hor_to_equ_transform(self):
        """The transformation function to transform coordinates from the
        horizontal system into the equatorial system.
        """
        return self._hor_to_equ_transform
    @hor_to_equ_transform.setter
    def hor_to_equ_transform(self, transform):
        if(not callable(transform)):
            raise TypeError('The hor_to_equ_transform property must be a callable object!')
        self._hor_to_equ_transform = transform

    def scramble(self, rss, data):
        """Scrambles the given data based on random MJD times, which are
        generated from a TimeGenerator instance. The event's right-ascention and
        declination coordinates are calculated via a horizontal-to-equatorial
        coordinate transformation and the generated MJD time of the event.

        Parameters
        ----------
        rss : RandomStateService
            The random state service providing the random number
            generator (RNG).
        data : instance of DataFieldRecordArray
            The DataFieldRecordArray instance containing the to be scrambled
            data.

        Returns
        -------
        data : DataFieldRecordArray
            The given DataFieldRecordArray holding the scrambled data.
        """
        mjds = self.timegen.generate_times(rss, len(data))
        data['time'] = mjds
        (data['ra'], data['dec']) = self.hor_to_equ_transform(
            data['azi'], data['zen'], mjds)
        return data


class DataScrambler(object):
    def __init__(self, method):
        """Creates a data scrambler instance with a given defined scrambling
        method.

        Parameters
        ----------
        method : DataScramblingMethod
            The instance of DataScramblingMethod that defines the method of
            the data scrambling.
        """
        self.method = method

    @property
    def method(self):
        """The underlaying scrambling method that should be used to scramble
        the data. This must be an instance of the DataScramblingMethod class.
        """
        return self._method
    @method.setter
    def method(self, method):
        if(not isinstance(method, DataScramblingMethod)):
            raise TypeError('The data scrambling method must be an instance '
                'of DataScramblingMethod!')
        self._method = method

    def scramble_data(self, rss, data, copy=False):
        """Scrambles the given data by calling the scramble method of the
        scrambling method class, that was configured for the data scrambler.
        If the ``inplace_scrambling`` property is set to False, a copy of the
        data is created before the scrambling is performed.

        Parameters
        ----------
        rss : RandomStateService
            The random state service providing the random number generator
            (RNG).
        data : instance of DataFieldRecordArray
            The DataFieldRecordArray instance holding the data, which should get
            scrambled.
        copy : bool
            Flag if a copy of the given data should be made before scrambling
            the data. The default is False.

        Returns
        -------
        data : DataFieldRecordArray
            The given DataFieldRecordArray instance with the scrambled data.
            If the ``inplace_scrambling`` property is set to True, this output
            array is the same array as the input array, otherwise it's a new
            array.
        """
        if(copy):
            data = data.copy()

        data = self._method.scramble(rss, data)

        return data


class TimeDepScrambling(DataScramblingMethod):
    """The TimeScramblingMethod class provides a data scrambling method to
        perform data coordinate scrambling based on a generated time. It draws a
        random time from a time generator and transforms the horizontal (local)
        coordinates into equatorial coordinates using a specified transformation
        function.
        """
    def __init__(self, events):
            """Initializes a new time scrambling instance.

            Parameters
            ----------
            timegen : TimeGenerator
                The time generator that should be used to generate random MJD times.
            hor_to_equ_transform : callable
                The transformation function to transform coordinates from the
                horizontal system into the equatorial system.

                The call signature must be:

                    __call__(azi, zen, mjd)

                The return signature must be: (ra, dec)

            """
            super(TimeDepScrambling, self).__init__()

            # get the times and the events for this for correct weighting
            self.weights = []
            for start, stop in zip(events.grl['start'], events.grl['stop']):
                mask_grl = (events.exp["time"] > start) & (events.exp['time'] < stop)
                # how many events in each run (relative to all events)
                self.weights.append(len(events.exp[mask_grl]) / len(events.exp['time']))
            
            # renormalize
            weight_sum = sum(self.weights)
            self.weights = [x / weight_sum for x in self.weights]
            
            self.grl = events.grl
            

    # this function could also be in utils, but I'm not sure in which utils.
    def azimuth_ra_converter(self, angles_in, mjd):
        """Rotate angles_in (right ascension / azimuth) according to the time mjd.
        The result is (azimuth / right ascension) since the formula is symmetric.
        This assumes the rotation can be approximated so that the axis is at Pole,
        which neglects the exact position of IceCube and all astronomical effects.
        
        Parameters
        ----------
        angles_in : angle (in rad), azimuth or zenith
        mjd : time in MJD

        Returns
        -------
        angle : (in rad), zenith or azimuth
        
        """

        # constants
        sidereal_length = 0.997269566  # sidereal day = length * solar day
        sidereal_offset = 2.54199002505
        sidereal_day_residuals = ((mjd / sidereal_length) % 1)
        angles_out = sidereal_offset + 2 * np.pi * sidereal_day_residuals - angles_in
        angles_out = np.mod(angles_out, 2 * np.pi)

        return angles_out
    

    def scramble(self, rss, data):
        """Scrambles the given data based on random MJD times, which are
        generated from a TimeGenerator instance. The event's right-ascension and
        declination coordinates are calculated via a horizontal-to-equatorial
        coordinate transformation and the generated MJD time of the event.

        Parameters
        ----------
        rss : RandomStateService
            The random state service providing the random number
            generator (RNG).
        data : instance of DataFieldRecordArray
            The DataFieldRecordArray instance containing the to be scrambled
            data.

        Returns
        -------
        data : DataFieldRecordArray
            The given DataFieldRecordArray holding the scrambled data.
        """

        # from which runs to draw the events, weighted with number of events in each run
        random_runs = rss.random.choice(self.grl['start'].size, size=len(data["time"]), p=self.weights)

        # draw the random times
        times = rss.random.uniform(
            self.grl['start'][random_runs],
            self.grl['stop'][random_runs])

        # get the correct right ascension 
        data['time'] = times
        data['ra'] = self.azimuth_ra_converter(data['azi'], times)

        return data
