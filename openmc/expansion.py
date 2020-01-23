import numpy as np
import openmc
from collections.abc import Iterable
from scipy.integrate import quad


def legendre_from_expcoef(coef, domain=(-1, 1)):
    """Return a Legendre series object based on expansion coefficients.

    Given a list of coefficients from FET tally and a array of down, return
    the numpy Legendre object.

    Parameters
    ----------
    coef : Iterable of float
        A list of coefficients of each term in Legendre polynomials
    domain : (2,) List of float
        Domain of the Legendre polynomial

    Returns
    -------
    numpy.polynomial.Legendre
        A numpy Legendre series class

    """

    n = np.arange(len(coef))
    c = (2*n + 1) * np.asarray(coef) / (domain[1] - domain[0])
    return np.polynomial.Legendre(c, domain)


class Expansion(object):
    """Abstract Polynomial Class for creating polynomials.
    """
    def __init__(self, coef):
        self.coef = np.asarray(coef)


class ZernikeRadial(Expansion):
    """Create radial only Zernike polynomials given coefficients and domain.

    The radial only Zernike polynomials are defined as in
    :class:`ZernikeRadialFilter`.

    Parameters
    ----------
    coef : Iterable of float
        A list of coefficients of each term in radial only Zernike polynomials
    radius : float
        Domain of Zernike polynomials to be applied on. Default is 1.
    r : float
        Position to be evaluated, normalized on radius [0,1]

    Attributes
    ----------
    order : int
        The maximum (even) order of Zernike polynomials.
    radius : float
        Domain of Zernike polynomials to be applied on. Default is 1.
    norm_coef : iterable of float
        The list of coefficients of each term in the polynomials after
        normailization.

    """
    def __init__(self, coef, radius=1):
        super().__init__(coef)
        self._order = 2 * (len(self.coef) - 1)
        self.radius = radius
        norm_vec = (2 * np.arange(len(self.coef)) + 1) / (np.pi * radius**2)
        self._norm_coef = norm_vec * self.coef

    @property
    def order(self):
        return self._order

    def __call__(self, r):
        import openmc.lib as lib
        if isinstance(r, Iterable):
            return [np.sum(self._norm_coef * lib.calc_zn_rad(self.order, r_i))
                    for r_i in r]
        else:
            return np.sum(self._norm_coef * lib.calc_zn_rad(self.order, r))


class Exponential(Expansion):
    """Create radial exponential expansion given coefficients and domain.

    The radial exponetial series are defined as in 
    :class:`ExponentialFilter`.

    Parameters
    ----------
    coef : Iterable of float
        A list of raw coefficients of each term in Exponential series
    radius : float
        Domain of Exponential polynomials to be applied on. Default is 1.
    exponent : float
        User-specified exponent to for a unique class of basis set
    r : float
        Position to be evaluated, normalized on radius [0,1]

    Attributes
    ----------
    order : int
        The maximum (even) order of expoential polynomials.
    radius : float
        Domain of Zernike polynomials to be applied on. Default is 1.
    norm_coef : iterable of float
        The list of coefficients of each term in the polynomials after
        normailization.

    """
    def __init__(self, coef, exponent, radius=1):
        super().__init__(coef)
        self._order = len(coef) - 1
        self._exponent = exponent
        self.radius = radius
        inn_prod_mat = np.empty([len(coef),len(coef)])
        in_mat = np.empty([len(coef),len(coef)])
        for i in range(len(coef)):
            for j in range(len(coef)):
                intergrand = lambda r: np.exp((i+j)*r**exponent)*(2*np.pi*r)
                y, err = quad(intergrand, 0, 1) 
                inn_prod_mat[i,j] = y
        self._norm_coef = np.linalg.solve(inn_prod_mat,coef)   

    @property
    def order(self):
        return self._order

    @property
    def exponent(self):
        return self._exponent

    @property
    def norm_coef(self):
        return self._norm_coef

    def __call__(self, r):
        import openmc.lib as lib
        if isinstance(r, Iterable):
            return [np.sum(self._norm_coef * lib.calc_exp(self.order, 
                self.exponent, r_i))/(self.radius**2) for r_i in r]

        else:
            return np.sum(self._norm_coef * lib.calc_exp(self.order, 
                self.exponent, r))/(self.radius**2)


