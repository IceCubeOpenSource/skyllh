# -*- coding: utf-8 -*-

"""Plotting module to plot IceCube specific background PDF objects.
"""
import numpy as np

from matplotlib.axes import Axes
from matplotlib.colors import LogNorm

from skylab.core.py import classname
from skylab.i3.backgroundpdf import I3BackgroundSpatialPDF

class I3BackgroundSpatialPDFPlotter(object):
    """Plotter class to plot an I3BackgroundSpatialPDF object.
    """
    def __init__(self, pdf):
        """Creates a new plotter object for plotting an I3BackgroundSpatialPDF
        object.

        Parameters
        ----------
        pdf : I3BackgroundSpatialPDF
            The PDF object to plot.
        """
        self.pdf = pdf

    @property
    def pdf(self):
        """The PDF object to plot.
        """
        return self._pdf
    @pdf.setter
    def pdf(self, pdf):
        if(not isinstance(pdf, I3BackgroundSpatialPDF)):
            raise TypeError('The pdf property must be an object of instance I3BackgroundSpatialPDF!')
        self._pdf = pdf

    def plot(self, axes):
        """Plots the spatial PDF. It uses the sin(dec) binning of the PDF to
        propperly represent the resolution of the PDF in the drawing.

        Parameters
        ----------
        axes : mpl.axes.Axes
            The matplotlib Axes object on which the PDF should get drawn to.
        """
        if(not isinstance(axes, Axes)):
            raise TypeError('The axes argument must be an instance of matplotlib.axes.Axes!')

        # By construction the I3BackgroundSpatialPDF does not depend on
        # right-ascention. Hence, we only need a single bin for the
        # right-ascention.
        sin_dec_binning = self.pdf.get_binning('sin_dec')
        pdfprobs = np.zeros((1, sin_dec_binning.nbins))

        sin_dec_points = sin_dec_binning.bincenters
        events = np.zeros((pdfprobs.size,),
                          dtype=[('sin_dec', np.float)])
        for (i, sin_dec) in enumerate(sin_dec_points):
            events['sin_dec'][i] = sin_dec

        event_probs = self.pdf.get_prob(events)

        for i in range(len(events)):
            pdfprobs[0,i] = event_probs[i]

        ra_axis = self.pdf.axes.get_axis('ra')
        (left, right, bottom, top) = (ra_axis.vmin, ra_axis.vmax,
                                      sin_dec_binning.lower_edge, sin_dec_binning.upper_edge)
        axes.imshow(pdfprobs.T, extent=(left, right, bottom, top), origin='lower',
                    norm=LogNorm(), interpolation='none')
        axes.set_xlabel('ra')
        axes.set_ylabel('sin_dec')
        axes.set_title(classname(self.pdf))
