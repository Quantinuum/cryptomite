from cryptomite.raz import Raz

import pytest


@pytest.mark.parametrize('n', [3, 7, 15])
@pytest.mark.parametrize('input_length', ['one', 'short', 'full'])
def test_raz_zero_extension(n, input_length, monkeypatch):
    extractor = Raz(2 * n, n)
    length = {'one': 1, 'short': n - 1, 'full': n}[input_length]
    input1 = [i % 2 for i in range(2 * n)]
    input2 = [(i + 1) % 2 for i in range(length)]
    original_input2 = input2.copy()
    multiply = extractor.gf_mul

    def checked_multiply(x, y):
        # Reject undersized buffers before the native NTT reads past them.
        assert len(x) == len(y) == 1 << extractor.logp
        return multiply(x, y)

    monkeypatch.setattr(extractor, 'gf_mul', checked_multiply)
    expected = extractor.extract(input1, input2 + [0] * (n - length))

    assert extractor.extract(input1, input2) == expected
    assert input2 == original_input2
