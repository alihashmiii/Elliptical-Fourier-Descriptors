#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A Python implementation of the method to calculate Fourier coefficients for characterizing
closed contours.
original EFA code in python was created by hbldh <henrik.blidh@nedomkull.com> on 2016-01-30
based on [#a] and [#b].
- Ali Hashmi transformed the EFA code into Classes on 2017-02-03 and introduced code for LOCO-EFA
on 2017-08-08 (from [#c])
References
----------
.. [#a] F. P. Kuhl and C. R. Giardina, “Elliptic Fourier Features of a
   Closed Contour," Computer Vision, Graphics and Image Processing,
   Vol. 18, pp. 236-258, 1982.
.. [#b] Oivind Due Trier, Anil K. Jain and Torfinn Taxt, “Feature Extraction
   Methods for Character Recognition - A Survey”, Pattern Recognition
   Vol. 29, No.4, pp. 641-662, 1996
 .. [#c] Yara E. Sanchez-Corrales, Matthew Hartley, Jop van Rooij,
Athanasius F. M. Maree, Veronica A. Grieneisen - Morphometrics of complex cell shapes: Lobe
Contribution Elliptic Fourier Analysis (LOCO-EFA)
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from ellipse import plot_ellipse
from contourexample import contours
import numpy as np
import csv

try:
    _range = xrange
except NameError:
    _range = range

# to convert a csv example to coordinates
filename = "C:\\Users\\Ali Hashmi\\Downloads\\DEVELOP_2017_156778\\cell_outline.csv"
with open(filename,'rU') as openfile:
    contours = [list(map(int,rec)) for rec in csv.reader(openfile)]


# class definitions
class EFAcoefficients:
    def __init__(self,contour,order):
        self.contour = contour
        self.order = order

    def elliptic_fourier_descriptors(self, normalize = True):

        """Calculate elliptical Fourier descriptors for a contour.
        :param numpy.ndarray contour: A contour array of size ``[M x 2]``.
        :param int order: The order of Fourier coefficients to calculate.
        :param bool normalize: If the coefficients should be normalized;
        see references for details.
        :return: A ``[order x 4]`` array of Fourier coefficients.
        :rtype: :py:class:`numpy.ndarray`
        """
        dxy = np.diff(self.contour, axis=0)
        dt = np.sqrt((dxy ** 2).sum(axis=1))
        t = np.concatenate([([0., ]), np.cumsum(dt)])
        T = t[-1]

        phi = (2 * np.pi * t) / T

        coeffs = np.zeros((self.order, 4))
        for n in _range(1, self.order + 1):
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
            coeffs = self.normalize_efd(coeffs)

        return coeffs

    @staticmethod
    def normalize_efd(coeffs, size_invariant = True):
        """Normalizes an array of Fourier coefficients.
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
            coeffs[n - 1, :] = np.dot(np.array([[coeffs[n - 1, 0], coeffs[n - 1, 1]],
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

    def calculate_dc_coefficients(self):
        """Calculate the :math:`A_0` and :math:`C_0` coefficients of the elliptic Fourier series.
        :param numpy.ndarray contour: A contour array of size ``[M x 2]``.
        :return: The :math:`A_0` and :math:`C_0` coefficients.
        :rtype: tuple
        """
        contour = np.array(self.contour)
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

class ShapeDescriptor(EFAcoefficients):
    # initialize
    def __init__(self,contour,order,normalized = False):
        EFAcoefficients.__init__(self,contour,order)
        self.coefficients = self.elliptic_fourier_descriptors(normalize=normalized)
        self.dcOffset = self.calculate_dc_coefficients()

    def printEFAcoeffs(self):
        """ prints the EFA coefficients """
        print(self.coefficients)

    def printdcOffset(self):
        """ prints the starting offsets """
        print(self.dcOffset)

    def printHarmonicinRSpace(self):
        """ changing EFA coefficients to lengths and respective angles """
        coeff = self.coefficients
        components = []
        for lst in (coeff):
            theta = 0.5 * np.arctan2(
                2*(lst[0] * lst[1] + lst[2] * lst[3]),
                (lst[0]**2 - lst[1]**2 + lst[2]**2 - lst[3]**2))
            coeffMat = np.array([[lst[0],lst[1]],[lst[2],lst[3]]])
            rotationMat = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
            componentMat = np.dot(coeffMat,rotationMat).flatten()
            a,b,c,d = componentMat
            psi = np.arctan2(c,a)*180/np.pi
            phi = np.arctan2(d,b)*180/np.pi
            L = np.sqrt(a**2 + c**2)
            W = np.sqrt(b**2 + d**2)
            components.append([L,W,psi,phi])
        print(np.swapaxes(np.array(components),1,1))

    def plot_efd(self, ell_coeff = None, image = None, ellipse = None, n = 300):
        """Plot a ``[2 x (N / 2)]`` grid of successive truncations of the series.
        .. note::
            Requires `matplotlib <http://matplotlib.org/>`_!
        :param numpy.ndarray coeffs: ``[N x 4]`` Fourier coefficient array.
        :param list, tuple or numpy.ndarray locus:
            The :math:`A_0` and :math:`C_0` elliptic locus in [#a]_ and [#b]_.
        :param int n: Number of points to use for plotting of Fourier series.
        """
        contour = self.contour
        coeffs = self.coefficients
        locus = self.dcOffset
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
        contour = np.array(contour)
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

    def LobeContribution(self,DEBUG=False):
        """ yields lobe contributions for holomorphic and nonholomorphic shapes """
        modes = self.coefficients

        starmat,lcL = [],[]
        lambdaminusvec,lambdaplusvec,zetaplusvec,zetaminusvec = ([] for i in range(4))
        lczetaplus,lclambdaplus,lczetaminus,lclambdaminus = \
        map(lambda x: [0]*len(modes), range(4))
        lcaplus,lcbplus,lccplus,lcdplus = ([0]*len(modes) for i in range(4))
        lcaminus,lcbminus,lccminus,lcdminus = ([0]*len(modes) for i in range(4))

        MAXERROR = 1e-13

        tau = 0.5*np.arctan2(2.0*(modes[0][0]*modes[0][1]+modes[0][2]*modes[0][-1]),
        modes[0][0]**2 + modes[0][2]**2 - modes[0][1]**2 - modes[0][-1]**2)
        alphaprime = modes[0][0]*np.cos(tau) + modes[0][1]*np.sin(tau)
        gammaprime = modes[0][2]*np.cos(tau) + modes[0][-1]*np.sin(tau)
        rho = np.arctan2(gammaprime,alphaprime)

        if rho < 0:
            tau += np.pi

        for j in range(len(modes)):
            alphastar = (modes[j][0])*np.cos((j+1)*tau) + (modes[j][1])*np.sin((j+1)*tau)
            betastar = - modes[j][0]*np.sin((j+1)*tau) + modes[j][1]*np.cos((j+1)*tau)
            gammastar = modes[j][2]*np.cos((j+1)*tau) + modes[j][-1]*np.sin((j+1)*tau)
            deltastar = - modes[j][2]*np.sin((j+1)*tau) + modes[j][-1]*np.cos((j+1)*tau)
            starmat.append([alphastar,betastar,gammastar,deltastar])

        r = starmat[0][0]*starmat[0][-1] - starmat[0][1]*starmat[0][-2]

        if r<0:
            for i,val in enumerate(starmat):
                starmat[i][1] *= -1
                starmat[i][-1] *= -1
        starmat = np.array(starmat)
        offx,offy = self.dcOffset
        a,b,c,d = (starmat[:,i] for i in [0,1,2,-1])

        if DEBUG:
            print("modified EFA coefficients\n=================================")
            [print("mode:",i+1,(a[i],b[i],c[i],d[i])) for i in range(len(modes))]
            print("\n")
        if DEBUG:
            print("Lambda matrices:\n==========================================")
        for i in range(len(starmat)):
            phi = 0.5*np.arctan2(2.0*(a[i]*b[i]+c[i]*d[i]),
        a[i]**2 + c[i]**2 - b[i]**2 - d[i]**2)

            aprime = a[i]*np.cos(phi) + b[i]*np.sin(phi)
            bprime = - a[i]*np.sin(phi) + b[i]*np.cos(phi)
            cprime = c[i]*np.cos(phi) + d[i]*np.sin(phi)
            dprime = -c[i]*np.sin(phi) + d[i]*np.cos(phi)

            theta = np.arctan2(cprime,aprime)

            lambda1 = aprime*np.cos(theta) + cprime*np.sin(theta)
            lambda12 = bprime*np.cos(theta) + dprime*np.sin(theta)
            lambda21 = - aprime*np.sin(theta) + cprime*np.cos(theta)
            lambda2 = - bprime*np.sin(theta) + dprime*np.cos(theta)

            if(DEBUG or abs(lambda12)>MAXERROR or abs(lambda21)>MAXERROR or
            lambda1<0 or lambda1 < abs(lambda2)):
                if abs(lambda12)>MAXERROR or abs(lambda21) > MAXERROR:
                    print("off-diagonal lambda matrix unequal to zero:\n")
                if lambda1 < 0:
                    print("Warning: lambda1 negative:\n")
                if lambda1 < abs(lambda2):
                    print("Warning: lambda 1 <|lambda2|:\n")
                print("mode: {0}\n(lambda1,lambda12):{1} {2} \n(lambda21,lambda2):{3} {4}\n".format(i+1,
                lambda1,lambda12,lambda21,lambda2))

            lambdaplus = (lambda1 + lambda2)/2.0
            lambdaminus = (lambda1 - lambda2)/2.0
            zetaplus = theta - phi
            zetaminus = - theta - phi

            zetaplusvec.append(zetaplus)
            lambdaplusvec.append(lambdaplus)
            lambdaminusvec.append(lambdaminus)
            zetaminusvec.append(zetaminus)

        lclambdapluszeroM = lambdaplusvec[1]
        lczetapluszeroM = zetaplusvec[1]

        lclambdaplus[0] = lambdaplusvec[0]
        lczetaplus[0] = zetaplusvec[0]

        for i in range(1,len(modes)-1):
            lclambdaplus[i]=lambdaplusvec[i+1]
            lczetaplus[i] = zetaplusvec[i+1]

        for i in range(1,len(modes)):
            lclambdaminus[i] = lambdaminusvec[i-1]
            lczetaminus[i] = zetaminusvec[i-1]

        if True:
            print("Ln quadruplets:\n===========================================")
            print("lc-EFA mode 0", (lclambdapluszeroM,lczetapluszeroM),"\n")
            for i in range(len(modes)):
                print("lc-EFA mode {0}\n(lambdaplus: {1} lambdaminus:{2})\
                \n(zetaplus: {3},zetaminus: {4})\n".format(i+1,
                lclambdaplus[i],lclambdaminus[i],lczetaplus[i],lczetaminus[i]))

        for i,val in enumerate(modes):
            lcaplus[i] = lclambdaplus[i]*np.cos(lczetaplus[i])
            lcbplus[i]= -lclambdaplus[i]*np.sin(lczetaplus[i])
            lccplus[i] = lclambdaplus[i]*np.sin(lczetaplus[i])
            lcdplus[i] = lclambdaplus[i]*np.cos(lczetaplus[i])
            lcaminus[i] = lclambdaminus[i]*np.cos(lczetaminus[i])
            lcbminus[i] = -lclambdaminus[i]*np.sin(lczetaminus[i])
            lccminus[i] = -lclambdaminus[i]*np.sin(lczetaminus[i])
            lcdminus[i] = -lclambdaminus[i]*np.cos(lczetaminus[i])

        if DEBUG:
            print("LC coefficients:\n===========================================")
            for i in range(len(modes)):
                print("\nmode: {0} A+\n{1} {2}\n{3} {4}".format(i+1,lcaplus[i],lcbplus[i],lccplus[i],lcdplus[i]))
                print("mode: {0} A-\n{1} {2}\n{3} {4}".format(i+1,lcaminus[i],lcbminus[i],
                lccminus[i],lcdminus[i]))

        for i in range(len(modes)):
            lcL.append(np.sqrt(lclambdaplus[i]**2 + lclambdaminus[i]**2
            + 2*lclambdaplus[i]*lclambdaminus[i]*np.cos(lczetaplus[i] -
            lczetaminus[i]-2*lczetaplus[0])))

        if True:
            print("\nLn scalar:\n===============================================")
            for i in range(len(modes)):
                print("LC-EFA mode {0}\tLn={1}".format(i+1,lcL[i]))

examplecontour = ShapeDescriptor(contours,10,normalized = False)

print("EFA coefficients:")
examplecontour.printEFAcoeffs()
print("\n")
print("---- [length_major, length_minor, angle_1, angle_2] ----")
examplecontour.printHarmonicinRSpace()
print("\n")
examplecontour.plot_efd()
print("LobeContribution:")
examplecontour.LobeContribution()
