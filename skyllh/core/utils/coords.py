# -*- coding: utf-8 -*-

import numpy as np

from astropy import (
    units,
)
from astropy.coordinates import (
    AltAz,
    ICRS,
    SkyCoord,
)
from astropy.time import (
    Time,
)

from skyllh.core.model import (
    DetectorModel,
)

def get_equ_to_hor_transform(detector_model):
    """Returns the equ_to_hor_transform function for the given detector model,
    i.e. detector location.

    Parameters
    ----------
    detector_model : instance of DetectorModel
        The instance of DetectorModel defining the location of the detector.

    Returns
    -------
    equ_to_hor_transform : callable
        The transform function to transform equatorial coordinates into
        horizontal coordinates. The function has the following call signature:

            __call__(ra, dec, mjd)

    """
    if not isinstance(detector_model, DetectorModel):
        raise TypeError(
            'The detector_model argument must be an instance of DetectorModel!')

    location = detector_model.location

    def equ_to_hor_transform(src_ra, src_dec, mjd):
        """Transforms the given horizontal coordinates into equatorial
        coordinates using astropy.

        Parameters
        ----------
        src_ra : instance of numpy.ndarray
            The (N,)-shaped numpy.ndarray holding the right-ascensions in
            radians.
        src_dec : instance of numpy.ndarray
            The (N,)-shaped numpy.ndarray holding the declinations angles in
            radians.
        mjd : instance of numpy.ndarray
            The (N,)-shaped numpy.ndarray holding the MJD times.
        location : instance of astropy.coordinates.EarthLocation
            The location on Earth of these horizontal coordinates.

        Returns
        -------
        azi : instance of numpy.ndarray
            The (N,)-shaped numpy.ndarray holding the azimuth angles in
            radians.
        alt : instance of numpy.ndarray
            The (N,)-shaped numpy.ndarray holding the altitude angles in
            radians.
        """
        equ = SkyCoord(
            ra=src_ra,
            dec=src_dec,
            frame="icrs",
            unit="rad"
        )

        hor = equ.tranform_to(
            AltAz(obstime=Time(mjd, format='mjd'),
                  location=location)
        )
        return (hor.az.to(units.rad).value, hor.alt.to(units.rad).value)

    return equ_to_hor_transform

def get_hor_to_equ_transform(detector_model):
    """Returns the hor_to_equ_transform function for the given detector model,
    i.e. detector location.

    Parameters
    ----------
    detector_model : instance of DetectorModel
        The instance of DetectorModel defining the location of the detector.

    Returns
    -------
    hor_to_equ_transform : callable
        The transform function to transform horizontal coordinates into
        equatorial coordinates. The function has the following call signature:

            __call__(azi, alt, mjd)

    """
    if not isinstance(detector_model, DetectorModel):
        raise TypeError(
            'The detector_model argument must be an instance of DetectorModel!')

    location = detector_model.location

    def hor_to_equ_transform(azi, alt, mjd):
        """Transforms the given horizontal coordinates into equatorial
        coordinates using astropy.

        Parameters
        ----------
        azi : instance of numpy.ndarray
            The (N,)-shaped numpy.ndarray holding the azimuth angles in
            radians.
        alt : instance of numpy.ndarray
            The (N,)-shaped numpy.ndarray holding the altitude angles in
            radians.
        mjd : instance of numpy.ndarray
            The (N,)-shaped numpy.ndarray holding the MJD times.
        location : instance of astropy.coordinates.EarthLocation
            The location on Earth of these horizontal coordinates.

        Returns
        -------
        ra : instance of numpy.ndarray
            The (N,)-shaped numpy.ndarray holding the right-ascensions in
            radians.
        dec : instance of numpy.ndarray
            The (N,)-shaped numpy.ndarray holding the declinations in
            radians.
        """
        hor = AltAz(
            az=azi*units.rad,
            alt=alt*units.rad,
            obstime=Time(mjd, format='mjd'),
            location=location)

        equ = hor.transform_to(ICRS())

        return (equ.ra.to(units.rad).value, equ.dec.to(units.rad).value)

    return hor_to_equ_transform


def rotate_spherical_vector(ra1, dec1, ra2, dec2, ra3, dec3):
    """Calculates the rotation matrix R to rotate the spherical vector
    (ra1,dec1) onto the direction (ra2,dec2), and performs this rotation on the
    spherical vector (ra3,dec3).

    In practice (ra1,dec1) refers to the true location of a MC event,
    (ra2,dec2) the true location of the signal source, and (ra3,dec3) the
    reconstructed location of the MC event, which should get rotated according
    to the rotation of the two true directions.
    """
    # Make sure, the inputs are 1D arrays.
    ra1 = np.atleast_1d(ra1)
    dec1 = np.atleast_1d(dec1)
    ra2 = np.atleast_1d(ra2)
    dec2 = np.atleast_1d(dec2)
    ra3 = np.atleast_1d(ra3)
    dec3 = np.atleast_1d(dec3)

    assert (
        len(ra1) == len(dec1) ==
        len(ra2) == len(dec2) ==
        len(ra3) == len(dec3)
    ), 'All input argument arrays must be of the same length!'

    N_event = len(ra1)

    # Calculate the space angle alpha between vector 1 and vector 2, and
    # correct for possible rounding errors.
    cos_alpha = (np.cos(ra2 - ra1) * np.cos(dec1) * np.cos(dec2)
                 + np.sin(dec1) * np.sin(dec2))
    cos_alpha[cos_alpha > 1] = 1
    cos_alpha[cos_alpha < -1] = -1
    alpha = np.arccos(cos_alpha)

    # Define the three 3D-vectors in spherical coordinate system. Each vector
    # is a (N_event,3)-shaped 2D array.
    vec1 = np.vstack([np.cos(ra1) * np.cos(dec1),
                      np.sin(ra1) * np.cos(dec1),
                      np.sin(dec1)]).T
    vec2 = np.vstack([np.cos(ra2) * np.cos(dec2),
                      np.sin(ra2) * np.cos(dec2),
                      np.sin(dec2)]).T
    vec3 = np.vstack([np.cos(ra3) * np.cos(dec3),
                      np.sin(ra3) * np.cos(dec3),
                      np.sin(dec3)]).T

    # Calculate the normalized rotation axis vector, nrot. nrot is a
    # (N_event,3)-shaped ndarray.
    nrot = np.cross(vec1, vec2)
    norm = np.sqrt(np.sum(nrot**2, axis=1))
    nrot[norm > 0] /= norm[np.newaxis, norm > 0].T

    # Define the diagonal 3D unit matrix.
    one = np.diagflat(np.ones(3))

    # Calculate the rotation matrix R_i for each event i and perform the
    # rotation on vector 3 for each event.
    vec = np.empty((N_event, 3), dtype=np.float64)

    sin_alpha = np.sin(alpha)
    twopi = 2*np.pi
    # Remap functions to avoid Python's (.)-resolution millions of times.
    (np_outer, np_dot, np_roll, np_diag, np_T) = (
     np.outer, np.dot, np.roll, np.diag, np.transpose)
    for i in range(N_event):
        cos_alpha_i = cos_alpha[i]
        nrot_i = nrot[i]
        nrotTnrot_i = np_outer(nrot_i, nrot_i)

        # Calculate cross product matrix, nrotx_i:
        # A[ij] = x_i * y_j - y_i * x_j
        skv = np_roll(np_roll(np_diag(nrot_i), shift=1, axis=1), shift=-1, axis=0)
        nrotx_i = skv - np_T(skv)

        # Calculate rotation matrix, R_i.
        R_i = ((1. - cos_alpha_i) * nrotTnrot_i
               + one*cos_alpha_i
               + sin_alpha[i] * nrotx_i)
        vec[i] = np_dot(R_i, np_T(vec3[i]))

    ra = np.arctan2(vec[:, 1], vec[:, 0])
    ra += np.where(ra < 0., twopi, 0.)
    dec = np.arcsin(vec[:, 2])

    return (ra, dec)


def rotate_signal_events_on_sphere(
        src_ra,
        src_dec,
        evt_true_ra,
        evt_true_dec,
        evt_reco_ra,
        evt_reco_dec,
):
    """Rotate signal events on a sphere to a given source position preserving
    position angle and separation (great circle distance) between the event's
    true and reco directions.

    Parameters
    ----------
    src_ra : instance of numpy.ndarray
        The (N_events,)-shaped 1D numpy.ndarray holding the true right-ascension
        of the source.
    src_dec : instance of numpy.ndarray
        The (N_events,)-shaped 1D numpy.ndarray holding the true declination
        of the source.
    evt_true_ra : instance of numpy.ndarray
        The (N_events,)-shaped 1D numpy.ndarray holding the true right-ascension
        of the MC event.
    evt_true_dec : instance of numpy.ndarray
        The (N_events,)-shaped 1D numpy.ndarray holding the true declination of
        the MC event.
    evt_reco_ra : instance of numpy.ndarray
        The (N_events,)-shaped 1D numpy.ndarray holding the reconstructed
        right-ascension of the MC event.
    evt_reco_dec : instance of numpy.ndarray
        The (N_events,)-shaped 1D numpy.ndarray holding the reconstructed
        declination of the MC event.

    Returns
    -------
    rot_evt_reco_ra : instance of numpy.ndarray
        The (N_events,)-shaped 1D numpy.ndarray holding the rotated
        reconstructed event right-ascension.
    rot_evt_reco_dec : instance of numpy.ndarray
        The (N_events,)-shaped 1D numpy.ndarray holding the rotated
        reconstructed event declination.
    """
    assert (
        len(src_ra) == len(src_dec) ==
        len(evt_true_ra) == len(evt_true_dec) ==
        len(evt_reco_ra) == len(evt_reco_dec)
    ), 'All input argument arrays must be of the same length!'

    v_source = SkyCoord(src_ra, src_dec, frame="icrs", unit="rad")
    v_evt_true = SkyCoord(evt_true_ra, evt_true_dec, frame="icrs", unit="rad")
    v_evt_reco = SkyCoord(evt_reco_ra, evt_reco_dec, frame="icrs", unit="rad")

    position_angle = v_evt_true.position_angle(v_evt_reco)
    separation = v_evt_true.separation(v_evt_reco)

    v_rotated = v_source.directional_offset_by(position_angle, separation)
    (rot_evt_reco_ra, rot_evt_reco_dec) = (v_rotated.ra.rad, v_rotated.dec.rad)

    return (rot_evt_reco_ra, rot_evt_reco_dec)


def angular_separation(ra1, dec1, ra2, dec2, psi_floor=None):
    """Calculates the angular separation on the sphere between two vectors on
    the sphere.

    Parameters
    ----------
    ra1 : instance of numpy.ndarray
        The (N_events,)-shaped numpy.ndarray holding the right-ascension or
        longitude coordinate of the first vector in radians.
    dec1 : instance of numpy.ndarray
        The (N_events,)-shaped numpy.ndarray holding declination or latitude
        coordinate of the first vector in radians.
    ra2 : instance of numpy.ndarray
        The (N_events,)-shaped numpy.ndarray holding the right-ascension or
        longitude coordinate of the second vector in radians.
    dec2 : instance of numpy.ndarray
        The (N_events,)-shaped numpy.ndarray holding declination coordinate of
        the second vector in radians.
    psi_floor : float | None
        If not ``None``, specifies the floor value of psi.

    Returns
    -------
    psi : instance of numpy.ndarray
        The (N_events,)-shaped numpy.ndarray holding the calculated angular
        separation value of each event.
    """
    delta_ra = np.abs(ra1 - ra2)
    delta_dec = np.abs(dec1 - dec2)

    x = np.sin(delta_dec / 2.)**2. +\
        np.cos(dec1) * np.cos(dec2) * np.sin(delta_ra / 2.)**2.

    # Handle possible floating precision errors.
    x[x < 0.] = 0.
    x[x > 1.] = 1.

    psi = 2. * np.arcsin(np.sqrt(x))

    if psi_floor is not None:
        psi = np.where(psi < psi_floor, psi_floor, psi)

    return psi
