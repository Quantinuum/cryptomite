"""
The Raz extractor [Raz2005]_ takes two inputs of length
'n_1, n_2', such that 'n_1/2 > n_2'. This implementation
is based on the efficient construction described in [Fore2025]_,
which requires a known irreducible trinomial for the field
GF_2^{n_1/2}.
"""
from __future__ import annotations

from math import ceil, floor, log2
from typing import cast

from cryptomite._cryptomite import NTT
from cryptomite.utils import BitsT, log_2

__all__ = ['Raz']


class Raz:
    """
    The Raz extractor [Raz2005]_ based on the efficient
    construction described in [Fore2025]_.
    """
    trinomial_s = {3: 1,
                   7: 1,
                   15: 1,
                   31: 3,
                   63: 1,
                   127: 7,
                   255: 52,
                   521: 32,
                   1279: 216,
                   2281: 715,
                   3217: 67,
                   4423: 271,
                   23209: 1530,
                   44497: 8575,
                   110503: 25230,
                   132049: 7000,
                   756839: 279695,
                   859433: 170340,
                   3021377: 361604,    # Hodgkin
                   6972593: 3037958,   # Bibury
                   24036583: 8412642,  # Judy-anne
                   25964951: 880890,   # t25a
                   30402457: 2162059,  # Florence
                   32582657: 5110722,  # Priscilla
                   42643801: 55981,    # t42a
                   43112609: 3569337,  # t43a
                   74207281: 9156813,  # t44a
                   }

    def __init__(self, n1: int, m: int, trinomial=None):
        """
        Initialize a Raz extractor.

        Parameters
        ----------
        n1: int
            The length of the first input (in bits).
            **Note:** GF_2^{n1/2} must have a known irreducible trinomial.
        m : int
            The length of the extractor output (in bits).
        trinomial : int
            An optional parameter which defines an irreducible trinomial over
            GF_2^{n1/2}, i.e., x^{n1/2} + x^{trinomial} + 1 is an irreducible
            trinomial. **Note:** This does not test whether the trinomial is
            irreducible, because it is computationally expensive to do so.
            Please ensure if you provide a value that it is correct, or use
            one of the known values.
        """
        assert (m <= n1/2)
        self.n = int(n1/2)
        self.m = m
        self.logp = log_2(self.n)+1
        self.pad_amount = (1 << self.logp) - self.n
        if trinomial is None:
            if self.n not in self.trinomial_s:
                raise ValueError('GF(2^(n1/2)) must have a known irreducible trinomial.')  # noqa: E501
            self.s = self.trinomial_s[self.n]
        else:
            self.s = trinomial
        self.ntt = NTT(self.logp)

    def __poly_reduce(self, x: BitsT):
        r = self.n
        s = self.s
        for i in range(self.n-1, -1, -1):
            red = x[r+i] % 2
            x[r+i] = 0
            x[s+i] = (x[s+i] % 2) ^ red
            x[i] = (x[i] % 2) ^ red

    def __gf_add(self, x: BitsT, y: BitsT) -> BitsT:
        return map(lambda b1, b2: b1 ^ b2, x, y)

    def __gf_add_one(self, input: BitsT) -> BitsT:
        x = input.copy()
        x[0] = x[0] ^ 1
        return x

    def gf_mul(self, x: BitsT, y: BitsT) -> BitsT:
        conv_output = cast(BitsT, [e % 2 for e in
                                   self.ntt.conv_and_reduce(x, y,
                                                            self.n, self.s)])
        return conv_output

    def extract(self, input1: BitsT, input2: BitsT) -> BitsT:
        """
        Perform randomness extraction.

        Parameters
        ----------
        input1 : list of bits (0s and 1s)
            The first input, consisting of n_1 bits.
        input2 : list of bits (0s and 1s)
            The second input, consisting of n_2 < n_1/2 bits.

        Returns
        -------
        list of bits (0s and 1s)
            The extractor output bits, of length m.
        """
        assert len(input1) >= 2*self.n
        assert 0 < len(input2) <= self.n

        x1, x2 = input1[0:self.n], input1[self.n:]
        y = input2 + [0] * (self.n - len(input2))

        pad_amount = self.pad_amount
        x1 = x1 + [0]*pad_amount
        x2 = x2 + [0]*pad_amount
        y = y + [0]*pad_amount

        cur_delta = self.gf_mul(y, x1)
        product = self.__gf_add_one(cur_delta)  # GF 1 times lambda_0
        cur_delta = self.gf_mul(cur_delta, cur_delta)  # delta^2

        # First loop is starting from product = product * lambda_1,
        # i.e. (delta^(2^1) + 1)), followed by cur_delta = delta^4
        # up to product * (lambda_{logp-1}, i.e. delta^{2^(logp-1))
        # and cur_delta = delta^(2^logp)
        # Note: we could unwrap the final multiplication at self.logp-1
        # to save the unnecessary calculation of delta^(2^logp),
        # but it happens in parallel anyway and makes the code harder to read.
        for _i in range(1, self.logp):
            # product times lambda_1 ... lambda_{l-1}
            (product, cur_delta) = self.ntt.raz_iteration(product, cur_delta,
                                                          self.n, self.s)
        conv_output = self.gf_mul(product, x2)
        return cast(BitsT, [conv_output[i] for i in range(self.m)])

# ------- UTILITY FUNCTIONS -------


def log2_error_raz(n_1: int,
                   k_1: float,
                   k_2: float,
                   m: int,
                   l: int,
                   p: int):
    """
    Compute an upper bound on the logarithm base 2
    of the error for the efficient weak version of
    Raz's extractor presented in [Fore2025]_.

    Parameters
    ----------
    n_1 : int
        The length of the first input (in bits).
    k_1 : float
        The min-entropy of the first input.
    k_2 : float
        The min-entropy of the second input.
    m : int
        The length of the extractor output (in bits).
    l : int
        The logarithm base 2 of p'.
        **Note:** the efficient construction requires p'
        to be a power of 2, i.e. l is an integer.
    p : int
        A free parameter that is an even integer that
        satisfies p <= 2^l / m.

    Returns
    -------
    float :
        An upper bound on the logarithm base 2 of the
        extractor error.
    """
    log_gamma_bound = (n_1 - k_1)/p + max((l - n_1/2 + 1)/p,
                                          log2(p) - k_2/2) + 1
    return log_gamma_bound + m/2


def opt_error_raz(n_1: int,
                  k_1: float,
                  n_2: int,
                  k_2: float,
                  m: int,
                  max_tests_basic=1,
                  max_tests_detailed=1000,
                  detailed_opt=False, verbose=False):
    """
    Compute an upper bound on the logarithm base 2
    of the error for the efficient weak version of
    Raz's extractor presented in [Fore2025]_.

    Parameters
    ----------
    n_1 : int
        The length of the first input (in bits).
    k_1 : float
        The min-entropy of the first input.
    n_2 : int
        The length of the second input (in bits).
    k_2 : float
        The min-entropy of the second input.
    m : int
        The length of the extractor output (in bits).
    max_tests_basic : int
        The maximum number of interations for the
        basic parameter optimisation method, i.e. when
        detailed_opt is set to False (default: 0).
    max_tests_detailed : int
        The maximum number of interations for the
        intense parameter optimisation method, i.e.
        when detailed_opt is set to True (default: 1000).
    detailed_opt : bool
        Flag to indicate the intensity of the
        optimisation performed (default: False).
    verbose : bool
        If True, prints all parameters found (default: False).

    Returns
    -------
    float :
        An upper bound on the logarithm base 2 of the
        extractor error.
    """
    # Ensure input parameters meet required constraints.
    assert 0 < n_2 <= n_1 / 2
    assert 0 < k_1 < n_1 and 0 < k_2 < n_2
    assert 0 < m <= n_1 / 2
    assert n_1 % 2 == 0
    assert max_tests_basic > 0
    assert max_tests_detailed > 0

    # Compute maximum possible l value based on input parameters.
    l_max = int(n_2 + floor(log2(n_1 / 2)))

    # Cap exponent to avoid overflow in 2^x computations.
    max_pow_for_overflow = 32

    # Initialize variables to track the best (minimal)
    # log2 error and corresponding parameters.
    min_log2_error, best_l, best_p = 0, 'Not found', 'Not found'

    # Estimate a good initial value for l based on m and (n_1 - k_1).
    l_use = max(floor(log2(m * (n_1 - k_1))), 1)

    # Define the range of l values to explore around l_use.
    max_plus = min(floor(l_max - l_use),
                   (max_tests_basic - 1) // 2)
    max_minus = min(floor(l_use - ceil(log2(m)) - 1),
                    (max_tests_basic - 1) // 2)
    ls = [i for i in range(l_use - max_minus, l_use + max_plus + 1)]

    # Coarse search: iterate over candidate l values.
    for current_l in ls:
        # Compute the maximum number of p values to try for the current l.
        p_half_max = int((2**(current_l - log2(m)))//2)
        # Generate candidate p values: even integers starting from 2.
        p_values = [2 * phalf + 2 for phalf in range(p_half_max)]
        for current_p in p_values:
            # Evaluate the log2 of the error for current parameters.
            eps = log2_error_raz(n_1, k_1, k_2, m, current_l, current_p)
            # Update best found parameters if error improves.
            if eps < min_log2_error:
                min_log2_error, best_l, best_p = eps, current_l, current_p

    # If detailed optimisation is enabled, perform a more exhaustive search.
    if detailed_opt:
        # Generate a list of l values to try with finer granularity.
        num_values = min(l_max, max_tests_detailed)
        step_size = (l_max - (ceil(log2(m)) + 1)) / (num_values - 1)
        ls = [
            int(round((ceil(log2(m)) + 1) + i * step_size)
                ) for i in range(num_values)]
        total = len(ls)
        # Define progress milestones for verbose output.
        milestones = {int(total * p / 100) for p in range(10, 101, 10)}
        for i, current_l in enumerate(ls):
            # Prevent overflow when computing p values.
            if current_l - log2(m) - 1 < max_pow_for_overflow:
                p_half_max = int(2**(current_l - log2(m) - 1))
            else:
                p_half_max = int(floor(2 ** max_pow_for_overflow))
            # Sample p_half values uniformly for detailed testing.
            num_values = min(p_half_max, max_tests_detailed)
            step_size = (p_half_max - 1) / num_values
            phalf_values = [
                int(round(1 + i * step_size)
                    ) for i in range(num_values + 1)]

            for current_phalf in phalf_values:
                current_p = 2*current_phalf
                eps = log2_error_raz(n_1, k_1, k_2, m, current_l, current_p)
                # Update best found parameters if error improves.
                if eps < min_log2_error:
                    min_log2_error, best_l, best_p = eps, current_l, current_p
            # Print progress if enabled and at a milestone.
            if verbose and i in milestones:
                percent = int((i / total) * 100)
                print(f'[{percent}% Completed] ({i}/{total})')

        if verbose:
            print(f'[100% Completed] ({total}/{total})')
            print('Performing final refinement...')

        # Final refinement around best_l after detailed search.
        if best_l != 'Not found':
            # Compute the range of l values to explore around best_l.
            max_plus = min(floor(l_max - best_l), max_tests_detailed // 2)
            max_minus = min(floor(best_l - ceil(log2(m)) - 1),
                            max_tests_detailed // 2)
            # Generate candidate l values around best_l.
            ls = [i for i in range(best_l - max_minus, best_l + max_plus + 1)]
            # Iterate over candidate l values for final refinement.
            for current_l in ls:
                # Prevent overflow when computing p values.
                if current_l - log2(m) - 1 < max_pow_for_overflow:
                    p_half_max = int(2**(current_l - log2(m) - 1))
                else:
                    p_half_max = int(floor(2 ** max_pow_for_overflow))
                # Sample p_half values uniformly for detailed testing.
                num_values = min(p_half_max, max_tests_detailed)
                step_size = (p_half_max - 1) / num_values
                phalf_values = [
                    int(round(1 + i * step_size)
                        ) for i in range(num_values + 1)]
                # Iterate over candidate p values.
                for current_phalf in phalf_values:
                    current_p = 2*current_phalf
                    eps = log2_error_raz(n_1,
                                         k_1,
                                         k_2,
                                         m,
                                         current_l,
                                         current_p)
                    if eps < min_log2_error:
                        min_log2_error = eps
                        best_l, best_p = current_l, current_p
    if verbose:
        print(f'log2 of minimum error found: {min_log2_error}')
        print(f'Corresponding l value (p prime = 2^l): {best_l}')
        print(f'Corresponding p value: {best_p}')
    # Return the minimal log2 error found.
    return min_log2_error


def calc_raz_out(n_1: int,
                 k_1: float,
                 n_2: int,
                 k_2: float,
                 log2_error_tol: float,
                 m_init=1,
                 max_tests_basic=1,
                 max_tests_detailed=1000,
                 detailed_opt=False,
                 verbose=False):
    """
    Compute the maximum output length for the efficient
    weak version of Raz's extractor presented in
    [Fore2025]_ that satisfies a given error tolerance.

    Parameters
    ----------
    n_1 : int
        The length of the first input (in bits).
    k_1 : float
        The min-entropy of the first input.
    n_2 : int
        The length of the second input (in bits).
    k_2 : float
        The min-entropy of the second input.
    log2_error_tol : float
        The logarithm base 2 of the acceptable extractor error.
        Must be negative, as the extractor error is 2^log2_error.
    m_init : int
        The initial value for the output length (default: 1).
    max_tests_basic : int
        The maximum number of interations for the
        basic parameter optimisation method, i.e. when
        detailed_opt is set to False (default: 0).
    max_tests_detailed : int
        The maximum number of interations for the
        intense parameter optimisation method, i.e.
        when detailed_opt is set to True (default: 1000).
    detailed_opt : bool
        Flag to indicate the intensity of the
        optimisation performed (default: False).
    verbose : bool
        If True, prints all parameters found (default: False).

    Returns
    -------
    int :
        The maximum output length that satisfies the
        error tolerance.
    """
    max_m = 0
    initial_tests = 100
    steps = min(floor(k_2) - m_init, initial_tests - 1)
    step_size = (floor(k_2) - m_init) / steps
    ms = [int(round(m_init + i * step_size)) for i in range(steps + 1)]

    for m in ms:
        if opt_error_raz(n_1, k_1, n_2, k_2, m,
                         max_tests_basic=max_tests_basic,
                         max_tests_detailed=max_tests_detailed,
                         detailed_opt=detailed_opt,
                         verbose=False) <= log2_error_tol:
            max_m = m
        else:
            break
    if detailed_opt:
        max_m += 1
        while opt_error_raz(n_1, k_1, n_2, k_2, max_m,
                            max_tests_basic=max_tests_basic,
                            max_tests_detailed=max_tests_detailed,
                            detailed_opt=detailed_opt,
                            verbose=False) <= log2_error_tol:
            max_m += 1
        max_m -= 1
    if verbose:
        print(f'Maximum output length found: {max_m}')
        opt_error_raz(n_1, k_1, n_2, k_2, max_m,
                      max_tests_basic=max_tests_basic,
                      max_tests_detailed=max_tests_detailed,
                      detailed_opt=detailed_opt,
                      verbose=True)
    if max_m > n_1 / 2:
        max_m = n_1 / 2
    return max_m

# ------ MAIN FUNCTION -------


@staticmethod
def from_params(
        n_1: int,
        k_1: float,
        n_2: int,
        k_2: float,
        log2_error: float,
        detailed_opt=False,
        verbose: bool = False) -> Raz:
    """
    Generate a weak version of the efficient Raz
    extractor from [Fore2025]_ with valid
    parameters based on input constraints.

    Parameters
    ----------
    n_1 : int
        The length of the first input (in bits).
    k_1 : float
        The min-entropy of the first input.
    n_2 : int
        The length of the second input (in bits).
    k_2 : float
        The min-entropy of the second input.
    log2_error : float
        The logarithm (base 2) of the acceptable extractor error.
        Must be negative, as the extractor error is 2^log2_error.
    detailed_opt : bool
        Flag to indicate the intensity of the optimisation
        performed (default: False).
    verbose : bool
        If True, prints the parameters used for
        extraction (default: True).

    Returns
    -------
    Raz
        A configured Raz extractor.

    Raises
    ------
    ValueError
        If the output length is non-positive.
    """
    m = calc_raz_out(n_1, k_1, n_2, k_2, log2_error,
                     detailed_opt=detailed_opt, verbose=verbose)
    if m <= 0:
        raise ValueError('Output length must be positive.')
    if verbose:
        print(
            f'--- New Raz Extractor Parameters ---\n'
            f'Input Length 1 (n_1): {n_1}, '
            f'Min Entropy of Input 1 (k_1): {k_1}, '
            f'Input Length 2 (n_2): {n_2}, '
            f'Min Entropy of Input 2 (k_2): {k_2}, '
            f'Output Length (m): {m}, '
            f'Extraction Error (log2_error): {log2_error}. '
        )
        print("""Adjust the length of the input
                and (weak) seed accordingly.""")
    return Raz(n_1, m)
