"""
This is a module.
"""
from __future__ import annotations

from math import sqrt
from typing import Literal, Sequence

# For the CompactBoolList
from typing import overload

from cryptomite._cryptomite import BigNTT, NTT

__all__ = ['is_prime', 'prime_facto', 'previous_prime', 'next_prime',
           'closest_prime', 'previous_na_set', 'next_na_set',
           'closest_na_set', 'suggest_extractor', 'von_neumann',

           'CompactBoolList'
           ]


BitT = Literal[0, 1]
BitsT = Sequence[BitT]

class CompactBoolList:
    """
    A compact list of booleans using bit-level storage, 
    designed to interface with functions expecting list[bool].
    """

    ############ Python native methods ##############
    def __init__(self, iterable: BitsT|Sequence[bool]|Sequence|None= None):
        self._data = bytearray()
        self._length = 0
        self._index: int = 0

        if iterable is not None:
            for item in iterable:
                self.append(item)

    def __len__(self):
        return self._length

    def _get_byte_and_bit_index(self, index: int) -> tuple[int, int]:
        """ Helper functions for getter """ 
        
        # Handle negative indexing for convenience
        if index < 0:
            index += self._length
        
        byte_index = index // 8
        bit_index = index % 8
        return byte_index, bit_index

    def _set_byte_and_bit_index(self, index: int, value: BitT|bool):
        """ Helper functions for setter """ 
        
        byte_idx, bit_idx = self._get_byte_and_bit_index(index)
        if value:
            self._data[byte_idx] |= (1 << bit_idx)
        else:
            self._data[byte_idx] &= ~(1 << bit_idx)

    ### Type check consistency
    @overload
    def __getitem__(self, index: int) -> bool: ...

    @overload
    def __getitem__(self, index: slice) -> CompactBoolList: ...

    # Actual  function signature
    def __getitem__(self, index: int|slice) -> bool|CompactBoolList: 
        # --- Slice Handling  ---
        if isinstance(index, slice):
            start, stop, step = index.indices(self._length)
            return self.__class__([self[i] for i in range(start, stop, step)])

        # --- Integer Handling ---
        # Boundary check for index
        if not (-self._length <= index < self._length):
             raise IndexError("Index out of range")
        
        byte_idx, bit_idx = self._get_byte_and_bit_index(index)
        mask = 1 << bit_idx
        return (self._data[byte_idx] & mask) != 0    

    def __setitem__(self, index: int, value: BitT|bool):
        # Boundary check for index
        if not (-self._length <= index < self._length):
             raise IndexError("Index out of range")

        self._set_byte_and_bit_index(index, value)

    def __iter__(self):
        return self

    def __next__(self) -> bool:
        if self._index < self._length:
            value = self.__getitem__(self._index)
            self._index += 1
            return value
        else:
            self._index = 0
            raise StopIteration

    def append(self, value: BitT|bool):
        if self._length % 8 == 0:
            self._data.append(0)
        
        self._set_byte_and_bit_index(self._length, value)            
        self._length += 1

    def __repr__(self) -> str:
        list_representation = list(self) 
        return f"{self.__class__.__name__}({list_representation!r})"

    ############## Custom methods ##############
    def to_list(self) -> list[bool]:
        """Converts to a standard list[bool] for functions expecting it."""
        return list(self)

def log_2(n: int) -> int:
    """
    Take the ceiling of the base 2 logarithm of an integer.

    Parameters
    ----------
    n : int
        The input integer.

    Returns
    -------
    int
        The ceiling of the base 2 logarithm of n.
    """
    x = 0
    while n > 0:
        n >>= 1
        x += 1
    return x


def conv(l: int, source1: Sequence[int], source2: Sequence[int]) -> list[int]:
    """
    Perform a cyclic convolution of size 2^l.

    Parameters
    ----------
    l : int
        The base 2 logarithm of the size of the convolution.
    source1: list of int
        The first input vector.
    source2: list of int
        The second input vector.

    Returns
    -------
    list[int]
        The output of the convolution.
    """
    L = 1 << l
    assert len(source1) == len(source2) == L
    ntt = BigNTT(l) if l > 30 else NTT(l)
    # ntt_source1 = ntt.ntt(source1, False)
    # ntt_source2 = ntt.ntt(source2, False)
    # mul_source = ntt.mul_vec(ntt_source1, ntt_source2)
    # conv_output = ntt.ntt(mul_source, True)
    return ntt.conv(source1, source2)


def is_prime(n: int) -> bool:
    """
    Checks whether an integer is prime.

    Parameters
    ----------
    n : int
        The integer to check for primality.

    Returns
    -------
    bool
        Whether n is prime.
    """
    for i in range(2, round(sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True


def prime_facto(n: int) -> tuple[list[int], list[int]]:
    """
    Factorizes an integer into its prime factors and their powers.
    It is required for the later function: na_set.

    Parameters
    ----------
    n : int
        The integer to check.

    Returns
    -------
    list
        A list representing the prime factorization of n.
    """
    factors = []
    i = 2
    while i <= round(sqrt(n)) + 1:
        if n % i == 0:
            factors.append(i)
            n = n // i
        else:
            i = i + 1
    if n != 1:
        factors.append(n)
    factors2 = [factors[0]]
    powers = [1]
    p = 0
    # this will put 24 = 2 * 2 * 2 * 3 in the form
    # factors=(2,3) and powers=(3,1) for (2^3) * (3^1)
    for i in range(1, len(factors)):
        if factors[i] == factors[i - 1]:
            powers[p] = powers[p] + 1
        if factors[i] != factors[i - 1]:
            powers.append(1)
            p = p + 1
            factors2.append(factors[i])
    return factors2, powers


def previous_prime(k: int) -> int:
    """
    Finds the largest prime number less than or equal to the given input.

    Parameters
    ----------
    k : int
        The upper limit for checking prime numbers.

    Returns
    -------
    int
        The largest prime number less than or equal to k.
    """
    k -= 1
    if k % 2 != 0:
        k = k - 1
    stop = False
    while not stop:
        stop = True
        while not is_prime(k + 1):
            k = k - 2
    return k + 1


def next_prime(k: int) -> int:
    """
    Finds the smallest prime number greater than or equal to the given input.

    Parameters
    ----------
    k : int
        The lower limit for checking prime numbers.

    Returns
    -------
    int
        The smallest prime number greater than or equal to k.
    """
    k -= 1
    if k % 2 != 0:
        k = k + 1
    stop = False
    while not stop:
        stop = True
        while not is_prime(k + 1):
            k = k + 2
    return k + 1


def closest_prime(k: int) -> int:
    """
    Finds the closest prime number to the given input.

    k : int
        The input value to find the closest prime number to.

    Returns
    -------
    int
        The closest prime number to k.
    """
    next_p = next_prime(k)
    previous_p = previous_prime(k)
    if next_p - k >= k - previous_p:
        out = previous_p
    else:
        out = next_p
    return out


def previous_na_set(k: int) -> int:
    """
    Finds the largest prime number with primitive root 2
    less than or equal to the given input.

    Parameters
    ----------
    k : int
        The upper limit for checking prime numbers with
        primitive root 2.

    Returns
    -------
    int
        The largest prime number with primitive root 2
        less than or equal to k.
    """
    k -= 1
    if k % 2 != 0:
        k = k - 1
    stop = False
    while not stop:
        stop = True
        while not is_prime(k + 1):
            k = k - 2
        primes, _ = prime_facto(k)
        for prime in primes:
            if pow(2, k // prime, k + 1) == 1:
                stop = False
                k = k - 2
                break
    return k + 1


def next_na_set(k: int) -> int:
    """
    Finds the smallest prime number with primitive root 2
    greater than or equal to the given input.

    Parameters
    ----------
    k : int
        The lower limit for checking prime numbers with
        primitive root 2.

    Returns
    -------
    int
        The smallest prime number with primitive root 2
        greater than or equal to k.
    """
    k -= 1
    if k % 2 != 0:
        k = k + 1
    stop = False
    while not stop:
        stop = True
        while not is_prime(k + 1):
            k = k + 2
        primes, _ = prime_facto(k)
        for prime in primes:
            if pow(2, k // prime, k + 1) == 1:
                stop = False
                k = k + 2
                break
    return k + 1


def closest_na_set(k: int) -> int:
    """
    Finds the closest prime number to the given input.

    k : int
        The input value to find the closest prime
        with primitive root 2 to.

    Returns
    -------
    int
        The closest prime number with primitive root 2
        to k.
    """
    next_p = next_na_set(k)
    previous_p = previous_na_set(k)
    if next_p - k >= k - previous_p:
        out = previous_p
    else:
        out = next_p
    return out


def suggest_extractor(n_1: int, exchangeable_sequence: bool,
                      efficiency_required: bool) -> str:
    """
    Suggests the best extractor for a user, based on Fig.2
    from the technical paper.

    Parameters
    ----------
    n_1 : int
        The length of the first input (in bits).
    exchangeable_sequence : bool
        Boolean input indicating whether the source forms an
        exchangeable sequence.
    efficiency_required : bool
        Boolean input indicating whether the user requires efficient
        extraction.

    Returns
    -------
    string
        The suggested extractor.
    """
    if exchangeable_sequence:
        out = 'Von Neumann'
    else:
        if n_1 <= 10**6 or efficiency_required:
            out = 'Circulant'
        else:
            out = 'Trevisan'
    return out


def von_neumann(bits: BitsT) -> BitsT:
    """
    Perform extraction using Von-Neumann extractor.
    """
    return [x for x, y in zip(bits[::2], bits[1::2]) if x != y]
