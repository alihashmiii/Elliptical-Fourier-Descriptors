#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A Python implementation of the method described in [#a]_ and [#b]_ for
calculating Fourier coefficients for characterizing
closed contours.
References
----------
.. [#a] F. P. Kuhl and C. R. Giardina, “Elliptic Fourier Features of a
   Closed Contour," Computer Vision, Graphics and Image Processing,
   Vol. 18, pp. 236-258, 1982.
.. [#b] Oivind Due Trier, Anil K. Jain and Torfinn Taxt, “Feature Extraction
   Methods for Character Recognition - A Survey”, Pattern Recognition
   Vol. 29, No.4, pp. 641-662, 1996
Created by hbldh <henrik.blidh@nedomkull.com> on 2016-01-30.
- Edited by Ali Hashmi on 2017-02-03
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from ellipse import plot_ellipse
from contourexample import contours

import numpy as np

try:
    _range = xrange
except NameError:
    _range = range


def elliptic_fourier_descriptors(contour, order=10, normalize = False):
    """Calculate elliptical Fourier descriptors for a contour.
    :param numpy.ndarray contour: A contour array of size ``[M x 2]``.
    :param int order: The order of Fourier coefficients to calculate.
    :param bool normalize: If the coefficients should be normalized;
        see references for details.
    :return: A ``[order x 4]`` array of Fourier coefficients.
    :rtype: :py:class:`numpy.ndarray`
    """
    dxy = np.diff(contour, axis=0)
    dt = np.sqrt((dxy ** 2).sum(axis=1))
    t = np.concatenate([([0., ]), np.cumsum(dt)])
    T = t[-1]

    phi = (2 * np.pi * t) / T

    coeffs = np.zeros((order, 4))
    for n in _range(1, order + 1):
        const = T / (2 * n * n * np.pi * np.pi)
        phi_n = phi * n
        d_cos_phi_n = np.cos(phi_n[1:]) - np.cos(phi_n[:-1])
        d_sin_phi_n = np.sin(phi_n[1:]) - np.sin(phi_n[:-1])
        a_n = const * np.sum((dxy[:, 0] / dt) * d_cos_phi_n)
        b_n = const * np.sum((dxy[:, 0] / dt) * d_sin_phi_n)
        c_n = const * np.sum((dxy[:, 1] / dt) * d_cos_phi_n)
        d_n = const * np.sum((dxy[:, 1] / dt) * d_sin_phi_n)
        coeffs[n - 1, :] = a_n, b_n, c_n, d_n

    if normalize:
        coeffs = normalize_efd(coeffs)

    return coeffs


def normalize_efd(coeffs, size_invariant = True):
    """Normalizes an array of Fourier coefficients.
    See [#a]_ and [#b]_ for details.
    :param numpy.ndarray coeffs: A ``[n x 4]`` Fourier coefficient array.
    :param bool size_invariant: If size invariance normalizing should be done as well.
        Default is ``True``.
    :return: The normalized ``[n x 4]`` Fourier coefficient array.
    :rtype: :py:class:`numpy.ndarray`
    """
    # Make the coefficients have a zero phase shift from
    # the first major axis. Theta_1 is that shift angle.
    theta_1 = 0.5 * np.arctan2(
        2 * ((coeffs[0, 0] * coeffs[0, 1]) + (coeffs[0, 2] * coeffs[0, 3])),
        ((coeffs[0, 0] ** 2) - (coeffs[0, 1] ** 2) + (coeffs[0, 2] ** 2) - (coeffs[0, 3] ** 2)))
    # Rotate all coefficients by theta_1.
    for n in _range(1, coeffs.shape[0] + 1):
        coeffs[n - 1, :] = np.dot(
            np.array([[coeffs[n - 1, 0], coeffs[n - 1, 1]],
                      [coeffs[n - 1, 2], coeffs[n - 1, 3]]]),
            np.array([[np.cos(n * theta_1), -np.sin(n * theta_1)],
                      [np.sin(n * theta_1), np.cos(n * theta_1)]])).flatten()

    # Make the coefficients rotation invariant by rotating so that
    # the semi-major axis is parallel to the x-axis.
    psi_1 = np.arctan2(coeffs[0, 2], coeffs[0, 0])
    psi_rotation_matrix = np.array([[np.cos(psi_1), np.sin(psi_1)],
                                    [-np.sin(psi_1), np.cos(psi_1)]])
    # Rotate all coefficients by -psi_1.
    for n in _range(1, coeffs.shape[0] + 1):
        coeffs[n - 1, :] = psi_rotation_matrix.dot(
            np.array([[coeffs[n - 1, 0], coeffs[n - 1, 1]],
                      [coeffs[n - 1, 2], coeffs[n - 1, 3]]])).flatten()

    if size_invariant:
        # Obtain size-invariance by normalizing.
        coeffs /= np.abs(coeffs[0, 0])

    return coeffs

def Length_Major_axis(coeffs):
    """returns the length of the Main elliptical axis. First locates the semi-major axix`
    """
    # Make the coefficients have a zero phase shift from
    # the first major axis. Theta_1 is that shift angle.
    theta_1 = 0.5 * np.arctan2(
        2 * ((coeffs[0, 0] * coeffs[0, 1]) + (coeffs[0, 2] * coeffs[0, 3])),
        ((coeffs[0, 0] ** 2) - (coeffs[0, 1] ** 2) + (coeffs[0, 2] ** 2) - (coeffs[0, 3] ** 2)))

    coeffs_N=np.array(coeffs)
    # Rotate all coefficients by theta_1.
    for n in _range(1, coeffs.shape[0] + 1):
        coeffs_N[n - 1, :] = np.dot(
            np.array([[coeffs_N[n - 1, 0], coeffs_N[n - 1, 1]],
                      [coeffs_N[n - 1, 2], coeffs_N[n - 1, 3]]]),
            np.array([[np.cos(n * theta_1), -np.sin(n * theta_1)],
                      [np.sin(n * theta_1), np.cos(n * theta_1)]])).flatten()

    # Make the coefficients rotation invariant by rotating so that
    # the semi-major axis is parallel to the x-axis.
    psi_1 = np.arctan2(coeffs_N[0, 2], coeffs_N[0, 0])
    psi_rotation_matrix = np.array([[np.cos(psi_1), np.sin(psi_1)],
                                    [-np.sin(psi_1), np.cos(psi_1)]])
    # Rotate all coefficients by -psi_1.
    for n in _range(1, coeffs.shape[0] + 1):
        coeffs_N[n - 1, :] = psi_rotation_matrix.dot(
            np.array([[coeffs_N[n - 1, 0], coeffs_N[n - 1, 1]],
                      [coeffs_N[n - 1, 2], coeffs_N[n - 1, 3]]])).flatten()

    Length_MAxis=np.sqrt(coeffs_N[0, 0]**2+coeffs_N[0, 2]**2)

    return Length_MAxis


def calculate_dc_coefficients(contour):
    """Calculate the :math:`A_0` and :math:`C_0` coefficients of the elliptic Fourier series.
    :param numpy.ndarray contour: A contour array of size ``[M x 2]``.
    :return: The :math:`A_0` and :math:`C_0` coefficients.
    :rtype: tuple
    """
    contour = np.array(contour);
    dxy = np.diff(contour, axis=0)
    dt = np.sqrt((dxy ** 2).sum(axis=1))
    t = np.concatenate([([0., ]), np.cumsum(dt)])
    T = t[-1]

    xi = np.cumsum(dxy[:, 0]) - (dxy[:, 0] / dt) * t[1:]
    A0 = (1 / T) * np.sum(((dxy[:, 0] / (2 * dt)) * np.diff(t ** 2)) + xi * dt)
    delta = np.cumsum(dxy[:, 1]) - (dxy[:, 1] / dt) * t[1:]
    C0 = (1 / T) * np.sum(((dxy[:, 1] / (2 * dt)) * np.diff(t ** 2)) + delta * dt)

    # A0 and CO relate to the first point of the contour array as origin.
    # Adding those values to the coefficients to make them relate to true origin.
    return contour[0, 0] + A0, contour[0, 1] + C0

def Length_Major_Minor_axis(coeffs):
    """returns the length and width of the Main elliptical axis. First locates the semi-major axix`
    """
    # Make the coefficients have a zero phase shift from
    # the first major axis. Theta_1 is that shift angle.

    theta_1 = 0.5 * np.arctan2(
        2 * ((coeffs[0, 0] * coeffs[0, 1]) + (coeffs[0, 2] * coeffs[0, 3])),
        ((coeffs[0, 0] ** 2) - (coeffs[0, 1] ** 2) + (coeffs[0, 2] ** 2) - (coeffs[0, 3] ** 2)))

    coeffs_N=np.array(coeffs)
    # Rotate all coefficients by theta_1.
    for n in _range(1, coeffs.shape[0] + 1):
        coeffs_N[n - 1, :] = np.dot(
            np.array([[coeffs_N[n - 1, 0], coeffs_N[n - 1, 1]],
                      [coeffs_N[n - 1, 2], coeffs_N[n - 1, 3]]]),
            np.array([[np.cos(n * theta_1), -np.sin(n * theta_1)],
                      [np.sin(n * theta_1), np.cos(n * theta_1)]])).flatten()

    # Make the coefficients rotation invariant by rotating so that
    # the semi-major axis is parallel to the x-axis.
    psi_1 = np.arctan2(coeffs_N[0, 2], coeffs_N[0, 0])
    psi_rotation_matrix = np.array([[np.cos(psi_1), np.sin(psi_1)],
                                    [-np.sin(psi_1), np.cos(psi_1)]])
    # Rotate all coefficients by -psi_1.
    for n in _range(1, coeffs.shape[0] + 1):
        coeffs_N[n - 1, :] = psi_rotation_matrix.dot(
            np.array([[coeffs_N[n - 1, 0], coeffs_N[n - 1, 1]],
                      [coeffs_N[n - 1, 2], coeffs_N[n - 1, 3]]])).flatten()

    Length_MAxis=np.sqrt(coeffs_N[0, 0]**2+coeffs_N[0, 2]**2)
    Length_MinAxis=np.sqrt(coeffs_N[0, 1]**2+coeffs_N[0, 3]**2)

    return Length_MAxis, Length_MinAxis


def DescriptorProperties(coeffs):
    # same as fourierHarmonicAnalysis
    """returns the lengths and orientations of the all ellipses."""

    save = [];
    coeffs = np.array(coeffs)
    coeffs_N = coeffs
    # Rotate all coefficients by theta.
    for n in _range(1, coeffs.shape[0] + 1):
        theta = 0.5 * np.arctan2(
            2 * ((coeffs[n-1, 0] * coeffs[n-1, 1]) + (coeffs[n-1, 2] * coeffs[n-1, 3])),
            ((coeffs[n-1, 0] ** 2) - (coeffs[n-1, 1] ** 2) + (coeffs[n-1, 2] ** 2) - (coeffs[n-1, 3] ** 2)))

        coeffs_N[n - 1, :] = np.dot(
            np.array([[np.cos(theta), np.sin(theta)],
                      [-np.sin(theta), np.cos(theta)]]),
                      np.array([[coeffs_N[n - 1, 0], coeffs_N[n - 1, 2]],
                                [coeffs_N[n - 1, 1], coeffs_N[n - 1, 3]]])).flatten()


        Length_MAxis=np.sqrt(coeffs_N[n-1, 0]**2+coeffs_N[n-1, 1]**2)
        Length_MinAxis=np.sqrt(coeffs_N[n-1, 2]**2+coeffs_N[n-1, 3]**2)
        phi_1 = np.arctan2(coeffs_N[n-1, 1],coeffs_N[n-1, 0])*180/np.pi
        phi_2 = np.arctan2(coeffs_N[n-1, 3],coeffs_N[n-1, 2])*180/np.pi
        save.append([Length_MAxis,Length_MinAxis,phi_1,phi_2])
    return np.swapaxes(np.array(save),1,1)


def plot_efd(coeffs, ell_coeff=None, locus=(0., 0.), image=None, contour=None, ellipse=None, n=300):
    """Plot a ``[2 x (N / 2)]`` grid of successive truncations of the series.
    .. note::
        Requires `matplotlib <http://matplotlib.org/>`_!
    :param numpy.ndarray coeffs: ``[N x 4]`` Fourier coefficient array.
    :param list, tuple or numpy.ndarray locus:
        The :math:`A_0` and :math:`C_0` elliptic locus in [#a]_ and [#b]_.
    :param int n: Number of points to use for plotting of Fourier series.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Cannot plot: matplotlib was not installed.")
        return

    N = coeffs.shape[0]
    N_half = int(np.ceil(N / 2))
    n_rows = 2

    t = np.linspace(0, 1.0, n)
    xt = np.ones((n,)) * locus[0]
    yt = np.ones((n,)) * locus[1]
    contour = np.array(contour);
    for n in _range(coeffs.shape[0]):
        xt += (coeffs[n, 0] * np.cos(2 * (n + 1) * np.pi * t)) + \
              (coeffs[n, 1] * np.sin(2 * (n + 1) * np.pi * t))
        yt += (coeffs[n, 2] * np.cos(2 * (n + 1) * np.pi * t)) + \
              (coeffs[n, 3] * np.sin(2 * (n + 1) * np.pi * t))
        ax = plt.subplot2grid((n_rows, N_half), (n // N_half, n % N_half))
        ax.set_title(str(n + 1))
        if contour is not None:
            ax.plot(contour[:, 1], contour[:, 0], 'c--', linewidth=2)
        ax.plot(yt, xt, 'r', linewidth=2)
        if image is not None:
            ax.imshow(image, plt.cm.gray)
        if ellipse is not None:

            semimaj=ell_coeff[1][0]/2
            semimin=ell_coeff[1][1]/2
            phi=ell_coeff[2]/180*np.pi
            theta = np.linspace(0,2*np.pi,1e4)
            r = semimaj*semimin / np.sqrt((semimin*np.cos(theta))**2 + (semimaj*np.sin(theta))**2)
            x = r*np.cos(theta)
            y = r*np.sin(theta)
            data = np.array([x,y])
            R = np.array([[np.cos(phi),-np.sin(phi)],[np.sin(phi),np.cos(phi)]])
            data = np.dot(R,data)
            data[0] += ell_coeff[0][1]
            data[1] += ell_coeff[0][0]

            ax.plot(data[0],data[1],color='b',linestyle='-')

    plt.show()

def fourierHarmonicAnalysis(coeff):
    components = [];
    for lst in (coeff):
        theta = 0.5 * np.arctan2(
            2*(lst[0] * lst[1] + lst[2] * lst[3]),
            (lst[0]**2 - lst[1]**2 + lst[2]**2 - lst[3]**2))
        coeffMat = np.array([[lst[0],lst[1]],[lst[2],lst[3]]])
        rotationMat = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]);
        componentMat = np.dot(coeffMat,rotationMat).flatten();
        a,b,c,d = componentMat
        psi = np.arctan2(c,a)*180/np.pi;
        phi = np.arctan2(d,b)*180/np.pi;
        L = np.sqrt(a**2 + c**2);
        W = np.sqrt(b**2 + d**2);
        components.append([L,W,psi,phi])
    return np.swapaxes(np.array(components),1,1)


lst = elliptic_fourier_descriptors(contours, order=10, normalize=False);
print("<EFA coefficients> ","\n", lst)
print("<fourierHarmonicAnalysis>","\n",fourierHarmonicAnalysis(lst));
plot_efd(lst,contour=contours,locus = calculate_dc_coefficients(contours),n=300)
