import numpy as np
from phd import stats as stats

"""
Title:
    thermo.py
Last update:
2021-01-22
Author(s):
    Manuel Razo-Mejia
Purpose:
    This module deals with the definition of functions derived from the 
    thermodynamic models of gene regulation.
"""

def pact_log(iptg, ea, ei, epsilon=4.5, n=2):
    '''
    Returns the probability of a repressor being active as described
    by the MWC model.

    Parameter
    ---------
    iptg : array-like.
        Concentrations of inducer on which to evaluate the function.
        All values must be positive.
    ea, ei : float.
        Minus log of the dissociation constants of the active and the
        inactive states respectively.
    epsilon : float.
        Positive log of the energy difference between the active and the
        inactive state.
    n : int
        Number of inducer binding sites.
    Returns
    -------
    pact : float.
        probability of a repressor of being in the active state.
        Active state is defined as the state that can bind to the DNA.

    Raises
    ------
    ValueError
        Thrown if any value of iptg concentration are negative.
    '''
    # Ensure that all values of iptg are positive.
    if (iptg < 0).any():
        raise ValueError('iptg array cannot have negative values.')

    pact = (1 + iptg * np.exp(ea))**n / ((1 + iptg * np.exp(ea))**n +
                                         np.exp(-epsilon) *
                                         (1 + iptg * np.exp(ei))**n)

    return pact



def fold_change_log(
    iptg, ea, ei, epsilon, R, epsilon_r, n=2, quaternary_state=2, nonspec_sites=4.6e6
):
    """
    Returns the gene expression fold change according to the
    thermodynamic model with the extension that takes into account the
    effect of the inducer.

    Parameter
    ---------
    iptg : array-like.
        Concentrations of inducer on which to evaluate the function
    ea, ei : float.
        Minus log of the dissociation constants of the active and the
        inactive states respectively
    epsilon : float.
        Energy difference between the active and the inactive state
    R : array-like.
        Repressor copy number for each of the strains. The length of
        this array should be equal to the iptg array. If only one value
        of the repressor is given it is asssume that all the data points
        should be evaluated with the same repressor copy number
    epsilon_r : array-like
        Repressor binding energy. The length of this array
        should be equal to the iptg array. If only one value of the
        binding energy is given it is asssume that all the data points
        should be evaluated with the same repressor copy number
    quaternary_state: int
        Prefactor in front of R in fold-change. Default is 2
        indicating that there are two functional heads per repressor molecule.
        This value must not be zero.
    nonspec_sites : int
        Number of nonspecific binding sites in the system.
        This value must be greater than 0.

    Returns
    -------
    fold_change : float.
        Gene expression fold change as dictated by the thermodynamic model.

    Raises
    ------
    ValueError
        Thrown if any entry of the IPTG vector, number of repressors,
        quaternary prefactor, or number of nonspecific binding sites is
        negative. This is also thrown if the quaternary
        state  or number of nonspecific binding sites is 0.


    """
    # Ensure that IPTG values and R is positive.
    if (iptg < 0).any() or (R < 0).any():
        raise ValueError("iptg and R must be positive.")
    if (quaternary_state <= 0) or (nonspec_sites <= 0):
        raise ValueError(
            "quaternary_state  and nonspec_sites must be greater\
        than zero."
        )

    return (
        1
        + quaternary_state
        * R
        / nonspec_sites
        * pact_log(iptg, ea, ei, epsilon, n)
        * (1 + np.exp(-epsilon))
        * np.exp(-epsilon_r)
    ) ** -1


def pact(IPTG, K_A, K_I, e_AI):
    """
    Computes the probability that a repressor is active
    Parameters
    ----------
    IPTG : array-like
        Array of IPTG concentrations in uM
    K_A : float
        Dissociation constant for active repressor
    K_I : float
        Dissociation constant for inactive repressor
    e_AI : float
        Energetic difference between the active and inactive state
    Returns
    -------
    probability that repressor is active
    """
    pact = (1 + IPTG * 1 / K_A) ** 2 / (
        ((1 + IPTG * 1 / K_A)) ** 2 + np.exp(-e_AI) * (1 + IPTG * 1 / K_I) ** 2
    )
    return pact


def fold_change(IPTG, K_A, K_I, e_AI, R, Op):
    """
    Computes fold-change for simple repression
    Parameters
    ----------
    IPTG : array-like
        Array of IPTG concentrations in uM
    K_A : float
        Dissociation constant for active repressor
    K_I : float
        Dissociation constant for inactive repressor
    e_AI : float
        Energetic difference between the active and inactive state
    R : float
        Number of repressors per cell
    Op : float
        Operator binding energy
    Returns
    -------
    probability that repressor is active
    """
    return 1 / (1 + R / 5e6 * pact(IPTG, K_A, K_I, e_AI) * np.exp(-Op))

def leakiness(num_rep, ep_r, ep_ai, n_ns=4.6e6):
    """
    Leakiness value of the simple-repression motif
    """
    pact = 1 / (1 + np.exp(-ep_ai))
    return (1 + pact * (num_rep / n_ns) * np.exp(-ep_r)) ** -1


def saturation(num_rep, ep_r, ep_ai, ka_ki, n_sites=2, n_ns=4.6e6):
    """
    Maximum fold-change of the simple-repression motif
    """
    pact = 1 / (1 + np.exp(-ep_ai) * ka_ki ** n_sites)
    return (1 + pact * (num_rep / n_ns) * np.exp(-ep_r)) ** -1


def saturation_cred_region(
    num_rep,
    ep_r,
    ep_ai,
    ka_flatchain,
    ki_flatchain,
    n_sites=2,
    n_ns=4.6e6,
    mass_frac=0.95,
):
    """
    Credible region of the maximum fold-change of the simple repression motif
    """
    pact = 1 / (1 + np.exp(-ep_ai) * (ka_flatchain / ki_flatchain) ** n_sites)
    cred_region = np.zeros([2, len(num_rep)])
    for i, R in enumerate(num_rep):
        fc = (1 + pact * (R / n_ns) * np.exp(-ep_r)) ** -1
        cred_region[:, i] = stats.hpd(fc, mass_frac)
    return cred_region


def dyn_range(num_rep, ep_r, ka_ki, ep_ai=4.5, n_sites=2, n_ns=4.6e6):
    """
    Dynamic range of the simple-repression motif
    """
    pact_leak = 1 / (1 + np.exp(-ep_ai))
    pact_sat = 1 / (1 + np.exp(-ep_ai) * (ka_ki) ** n_sites)
    leak = (1 + pact_leak * (num_rep / n_ns) * np.exp(-ep_r)) ** -1
    sat = (1 + pact_sat * (num_rep / n_ns) * np.exp(-ep_r)) ** -1
    return sat - leak

def dyn_cred_region(
    num_rep, ka_flatchain, ki_flatchain, ep_r, mass_frac=0.95, epsilon=4.5
):
    """
    Credible region of the dynamic range of the simple repression motif
    """
    cred_region = np.zeros([2, len(num_rep)])
    ka_ki = ka_flatchain / ki_flatchain
    for i, R in enumerate(num_rep):
        drng = dyn_range(R, ep_r, ka_ki, ep_ai=epsilon)
        cred_region[:, i] = stats.hpd(drng, mass_frac)
    return cred_region


def EC50(K_A, K_I, e_AI, R, Op):
    """
    Computes the concentration at which half of the repressors are in the active state
    Parameters
    ----------
    K_A : float
        Dissociation constant for active repressor
    K_I : float
        Dissociation constant for inactive repressor
    e_AI : float
        Energetic difference between the active and inactive state

    Returns
    -------
    Concentration at which half of repressors are active (EC50)
    """
    t = (
        1
        + (R / 4.6e6) * np.exp(-Op)
        + (K_A / K_I) ** 2 * (2 * np.exp(-e_AI) + 1 + (R / 4.6e6) * np.exp(-Op))
    )
    b = (
        2 * (1 + (R / 4.6e6) * np.exp(-Op))
        + np.exp(-e_AI)
        + (K_A / K_I) ** 2 * np.exp(-e_AI)
    )
    return K_A * ((K_A / K_I - 1) / (K_A / K_I - (t / b) ** (1 / 2)) - 1)


def ec50_cred_region(num_rep, Op, e_AI, K_A, K_I, mass_frac=0.95):
    """
    Computes the credible region of the concentration at which half of the
    repressors are in the active state
    Parameters
    ----------
    K_A : float
        Dissociation constant for active repressor
    K_I : float
        Dissociation constant for inactive repressor
    e_AI : float
        Energetic difference between the active and inactive state

    Returns
    -------
    Concentration at which half of repressors are active (EC50)
    """
    cred_region = np.zeros([2, len(num_rep)])
    for i, R in enumerate(num_rep):
        t = (
            1
            + (R / 4.6e6) * np.exp(-Op)
            + (K_A / K_I) ** 2 * (2 * np.exp(-e_AI) + 1 + (R / 4.6e6) * np.exp(-Op))
        )
        b = (
            2 * (1 + (R / 4.6e6) * np.exp(-Op))
            + np.exp(-e_AI)
            + (K_A / K_I) ** 2 * np.exp(-e_AI)
        )
        ec50_rng = K_A * ((K_A / K_I - 1) / (K_A / K_I - (t / b) ** (1 / 2)) - 1)
        cred_region[:, i] = stats.hpd(ec50_rng, mass_frac)
    return cred_region


def effective_Hill(K_A, K_I, e_AI, R, Op):
    """
    Computes the effective Hill coefficient
    Parameters
    ----------
    K_A : float
        Dissociation constant for active repressor
    K_I : float
        Dissociation constant for inactive repressor
    e_AI : float
        Energetic difference between the active and inactive state
    Returns
    -------
    effective Hill coefficient
    """
    c = EC50(K_A, K_I, e_AI, R, Op)
    return (
        2
        / (
            fold_change(c, K_A, K_I, e_AI, R, Op)
            - fold_change(0, K_A, K_I, e_AI, R, Op)
        )
        * (
            -((fold_change(c, K_A, K_I, e_AI, R, Op)) ** 2)
            * R
            / 5e6
            * np.exp(-Op)
            * 2
            * c
            * np.exp(-e_AI)
            * (
                1 / K_A * (1 + c / K_A) * (1 + c / K_I) ** 2
                - 1 / K_I * (1 + c / K_A) ** 2 * (1 + c / K_I)
            )
            / ((1 + c / K_A) ** 2 + np.exp(-e_AI) * (1 + c / K_I) ** 2) ** 2
        )
    )


def effective_hill_cred(num_rep, Op, e_AI, K_A, K_I, mass_frac=0.95):
    """
    Credible region for the effective hill coefficient
    Parameters
    ----------
    K_A : float
        Dissociation constant for active repressor
    K_I : float
        Dissociation constant for inactive repressor
    e_AI : float
        Energetic difference between the active and inactive state
    Returns
    -------
    effective Hill coefficient
    """
    cred_region = np.zeros([2, len(num_rep)])
    for i, R in enumerate(num_rep):
        # Compute the EC50
        c = EC50(K_A, K_I, e_AI, R, Op)
        # Compute the hill
        e_hill = (
            2
            / (
                fold_change(c, K_A, K_I, e_AI, R, Op)
                - fold_change(0, K_A, K_I, e_AI, R, Op)
            )
            * (
                -((fold_change(c, K_A, K_I, e_AI, R, Op)) ** 2)
                * R
                / 5e6
                * np.exp(-Op)
                * 2
                * c
                * np.exp(-e_AI)
                * (
                    1 / K_A * (1 + c / K_A) * (1 + c / K_I) ** 2
                    - 1 / K_I * (1 + c / K_A) ** 2 * (1 + c / K_I)
                )
                / ((1 + c / K_A) ** 2 + np.exp(-e_AI) * (1 + c / K_I) ** 2) ** 2
            )
        )
        cred_region[:, i] = stats.hpd(e_hill, mass_frac)

    return cred_region

def mcmc_cred_region(iptg, flatchain, R, epsilon_r, mass_frac=0.95, epsilon=4.5):
    """
    PROJECT: MWC Induction
    This function takes every element in the MCMC flatchain and computes
    the fold-change for each iptg concentration returning at the end the
    indicated mass_frac fraction of the fold change.

    Parameters
    ----------
    iptg : array-like.
        iptg concentrations on which evaluate the fold change
    flatchain : array-like.
        MCMC traces for the two MWC parameteres.
        flatchain[:,0] = ea flat-chain
        flatchain[:,1] = ei flat-chain
    R : float.
        Mean repressor copy number.
    epsilon_r : float.
        Repressor binding energy.
    mass_frac : float with 0 < mass_frac <= 1
        The fraction of the probability to be included in
        the HPD.  For example, `massfrac` = 0.95 gives a
        95% HPD.
    epsilon : float.
        Energy difference between active and inactive state.

    Returns
    -------
    cred_region : array-like
        array of 2 x len(iptg) with the upper and the lower fold-change HPD
        bound for each iptg concentration
    """
    # initialize the array to save the credible region
    cred_region = np.zeros([2, len(iptg)])

    # loop through iptg concentrations, compute all the fold changes and
    # save the HPD for each concentration
    for i, c in enumerate(iptg):
        fc = fold_change_log(c, flatchain[:, 0], flatchain[:, 1], epsilon, R, epsilon_r)
        cred_region[:, i] = stats.hpd(fc, mass_frac)

    return cred_region