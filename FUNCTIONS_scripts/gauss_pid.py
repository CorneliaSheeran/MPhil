"""
Quick implementation of a few PID functions for bivariate Gaussian
distributions.

Key references:

    Williams PL, Beer RD (2010). Nonnegative decomposition of multivariate
    information. arXiv:1004.2515.

    Barrett AB (2015). Exploration of synergistic and redundant information
    sharing in static and dynamical Gaussian systems. Physical Review E.

    Ince RA (2017). Measuring multivariate redundant information with pointwise
    common change in surprisal. Entropy.

    James RG, Emenheiser J, Crutchfield JP (2018). Unique information via
    dependency constraints. Journal of Physics A: Mathematical and Theoretical.

    Kay JW, Ince RA (2018). Exact partial information decompositions for
    Gaussian systems based on dependency constraints. Entropy.

Pedro Mediano, Dec 2020
"""
import numpy as np
from scipy.stats import multivariate_normal as mvn
from collections import namedtuple

PID = namedtuple('PID', 'Red UnX UnY Syn')

__all__ = ['PID', 'ccs', 'mmi', 'dep']


def _gauss_mi(S, idx):
    """
    Quick-and-dirty implementation of the mutual information between two sets
    of variables in a multivariate Gaussian distribution.

    Parameters
    ----------
    S : np.ndarray
        Covariance matrix. Must be square and positive definite.
    idx: iterable of iterables
        Iterable with two elements, each containing the index of variables to compute MI

    Returns
    -------
    mi : float
    """
    assert len(idx) == 2, 'Index list must contain exactly two items'
    allidx = [i for l in idx for i in l]
    Hx  = 0.5*np.linalg.slogdet(S[np.ix_(idx[0], idx[0])])[1]
    Hy  = 0.5*np.linalg.slogdet(S[np.ix_(idx[1], idx[1])])[1]
    Hxy = 0.5*np.linalg.slogdet(S[np.ix_(allidx, allidx)])[1]
    mi = Hx + Hy - Hxy

    return mi


def _fill_pid(S, idx, red=None, unx=None, uny=None, syn=None):
    """
    Given a covariance matrix and a PID atom, computes mutual information and
    solves the PID linear system of equations.

    Parameters
    ----------
    S : np.ndarray
        Covariance matrix. Must be square and positive definite.
    idx : list of lists
        List with three items, each containing the index of both sources and
        the target in the covariance matrix.
    red, unx, uny, syn : float
        Partial information atoms. Exactly one has to be given.

    Returns
    -------
    pid : named tuple
        Named tuple with information decomposition results
    """
    notnone = [a is not None for a in [red, unx, uny, syn]]
    assert np.sum(notnone) == 1, 'Exactly one of red, unx, uny, and syn has to be given'

    # Pre-compute MIs, which are needed anyway
    Ix  = _gauss_mi(S, [idx[0], idx[2]])
    Iy  = _gauss_mi(S, [idx[1], idx[2]])
    Ixy = _gauss_mi(S, [idx[0] + idx[1], idx[2]])

    # Assign remaining atoms depending on which atom was given
    if red is not None:
        unx = Ix - red
        uny = Iy - red
        syn = Ixy - unx - uny - red

    elif unx is not None:
        red = Ix - unx
        uny = Iy - red
        syn = Ixy - unx - uny - red

    elif uny is not None:
        red = Iy - uny
        unx = Ix - red
        syn = Ixy - unx - uny - red

    elif syn is not None:
        red = Ix + Iy - Ixy + syn
        unx = Ix - red
        uny = Iy - red

    # Assemble PID named tuple and return
    return PID(red, unx, uny, syn)


def _prepare_corr(X, Y, Z):
    """
    Compute joint correlation matrix from input arrays, reshaping from 1D to 2D
    as necessary.

    Parameters
    ---------
    X, Y, Z : np.ndarray of floats
        Input data. For data with T samples: if 2D, each array must have
        shape[1] == T, and if 1D, len == T

    Returns
    -------
    S : np.ndarray of floats
        Correlation matrix
    idx : list of lists
        List with three items, each containing the index of both sources and
        the target in the covariance matrix.
    """
    if X.ndim > 2 or Y.ndim > 2 or Z.ndim > 2:
        raise ValueError('Every input array must be of dimension 2 at most')

    ensure_vertical = lambda M: M if M.ndim == 2 else M[:,np.newaxis]
    vX, vY, vZ = ensure_vertical(X), ensure_vertical(Y), ensure_vertical(Z)

    S = np.corrcoef(np.hstack([vX, vY, vZ]).T)
    dimX, dimY, dimZ = vX.shape[1], vY.shape[1], vZ.shape[1]
    idx = [list(range(dimX)), 
           list(range(dimX, dimX + dimY)),
           list(range(dimX + dimY, dimX + dimY + dimZ))]

    return S, idx


def mmi(X, Y, Z):
    """
    Compute all PID atoms from sources X,Y to target Z, assuming they are all
    jointly Gaussian distributed and using Barrett's Minimum Mutual Information
    (MMI) PID.

    Parameters
    ----------
    X, Y, Z : np.ndarray of floats
        Input data. For data with T samples: if 2D, each array must have
        shape[1] == T, and if 1D, len == T

    Returns
    -------
    pid : named tuple
        Named tuple with information decomposition results
    """
    S, idx = _prepare_corr(X, Y, Z)
    red = np.min([_gauss_mi(S, [idx[0], idx[2]]), _gauss_mi(S, [idx[1], idx[2]])])
    return _fill_pid(S, idx, red=red)


def ccs(X, Y, Z, nb_samples=100000):
    """
    Compute all PID atoms from sources X,Y to target Z, assuming they are all
    jointly Gaussian distributed and using Ince's Common Change in Surprisal
    (CCS) PID.

    Parameters
    ----------
    X, Y, Z : np.ndarray of floats
        Input data. For data with T samples: if 2D, each array must have
        shape[1] == T, and if 1D, len == T
    nb_samples : int
        Number of Monte Carlo samples to use to approximate CCS (default: 1e5)

    Returns
    -------
    pid : named tuple
        Named tuple with information decomposition results
    """

    S, idx = _prepare_corr(X, Y, Z)

    sample = mvn.rvs(cov=S, size=nb_samples)

    h   = lambda x, i: -mvn.logpdf(sample[:,i], cov=S[np.ix_(i,i)])
    mi  = lambda x, src, tgt: h(x, src) + h(x, tgt) - h(x, src + tgt);

    # Local MI and co-I for each MC sample
    ix  = mi(sample, idx[0], idx[2])
    iy  = mi(sample, idx[1], idx[2])
    ixy = mi(sample, idx[0] + idx[1], idx[2])
    coinfo = ix + iy - ixy

    # Note: we define the coinfo to be positive when redundancy > synergy, so
    # the sign is changed wrt table 3 in Ince2017.

    signs = np.sign(np.vstack([ix, iy, ixy, coinfo]))

    # Samples for which all local quantities have the same sign
    p = np.sum(np.abs(signs - signs[0]), 0) == 0

    # Average coinfo in the relevant samples to compute redundancy and return
    red = coinfo[p].mean()

    return _fill_pid(S, idx, red=red)


def dep(X, Y, Z):
    """
    Compute all PID atoms from sources X,Y to target Z, assuming they are all
    jointly Gaussian distributed and using James' dependency contraints (dep)
    PID.

    Parameters
    ----------
    X, Y, Z : np.ndarray of floats
        Input data. For data with T samples: if 2D, each array must have
        shape[1] == T, and if 1D, len == T

    Returns
    -------
    pid : named tuple
        Named tuple with information decomposition results
    """
    S, idx = _prepare_corr(X, Y, Z)
    dimX, dimY, dimZ = tuple(map(len, idx))

    # Note: code based on:
    # https://github.com/robince/partial-info-decomp/blob/master/calc_pi_Idep_mvn.m

    # Extract blockwise correlation components
    P = S[np.ix_(idx[0], idx[1])]
    Q = S[np.ix_(idx[0], idx[2])]
    R = S[np.ix_(idx[1], idx[2])]

    # Standard mutual informations
    Ix  = _gauss_mi(S, [idx[0], idx[2]])
    Iy  = _gauss_mi(S, [idx[1], idx[2]])
    Ixy = _gauss_mi(S, [idx[0] + idx[1], idx[2]])

    # Necessary quantities
    halflogdet = lambda V: 0.5*np.linalg.slogdet(V)[1]
    inum = halflogdet(np.eye(dimY) - R@(Q.T)@Q@(R.T))
    iden = halflogdet(np.eye(dimZ) - (Q.T)@Q) + halflogdet(np.eye(dimZ) - (R.T)@R)
    knum = halflogdet(np.eye(dimY) - (P.T)@P)
    kden = halflogdet(S)

    # Dependency lattice edges
    b = Ix;
    i = inum - iden - Iy;
    k = knum - kden - Iy;

    # Take the minimum of all dependency lattice edges and fill PID
    unx = np.min([b, i, k])

    return _fill_pid(S, idx, unx=unx)

