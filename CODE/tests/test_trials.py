"""The Monte-Carlo statistics of calc_TP_M1E2.

The two shape statistics are Excel's, because that is what the author
checks the sampling against; they are pinned here to values computed by
hand from the same definitions.
"""

import math

import pytest

from toolbox import trials


# --------------------------------------------------------------------------
# the deviates
# --------------------------------------------------------------------------

def test_one_seed_gives_one_stream():
    a = trials.Randomizer('sn8')
    b = trials.Randomizer('sn8')
    assert [a.normal(1.0, 0.1) for _ in range(20)] == \
           [b.normal(1.0, 0.1) for _ in range(20)]


def test_two_seeds_do_not():
    a = trials.Randomizer('one')
    b = trials.Randomizer('two')
    assert [a.normal(0.0, 1.0) for _ in range(5)] != \
           [b.normal(0.0, 1.0) for _ in range(5)]


def test_a_zero_sigma_leaves_the_value_alone():
    rng = trials.Randomizer(1)
    assert rng.normal(3.5, 0.0) == 3.5


def test_the_deviates_have_the_mean_and_width_asked_for():
    rng = trials.Randomizer('wide')
    got = [rng.normal(10.0, 2.0) for _ in range(20000)]
    mean = sum(got) / len(got)
    sigma = math.sqrt(sum((x - mean) ** 2 for x in got) / len(got))
    assert mean == pytest.approx(10.0, abs=0.05)
    assert sigma == pytest.approx(2.0, abs=0.05)


# --------------------------------------------------------------------------
# the running sums
# --------------------------------------------------------------------------

def test_an_accumulator_is_the_population_mean_and_sigma():
    a = trials.Accumulator()
    for x in (1.0, 2.0, 3.0, 4.0):
        a.add(x)
    assert a.n == 4
    assert a.mean == 2.5
    assert a.sigma == pytest.approx(math.sqrt(1.25))


def test_an_empty_accumulator_is_zero_rather_than_an_error():
    a = trials.Accumulator()
    assert (a.n, a.mean, a.sigma) == (0, 0.0, 0.0)


def test_rounding_error_never_makes_a_negative_variance():
    a = trials.Accumulator()
    for _ in range(1000):
        a.add(1e8)
    assert a.sigma == 0.0


# --------------------------------------------------------------------------
# the two statistics the Perl prints
# --------------------------------------------------------------------------

def test_stdev_is_the_rms_deviation_from_the_reference():
    # Spread about the data's own mean, plus the offset to the reference,
    # in quadrature (`calc_TP_M1E2.pl:1076`).
    data = [1.0, 3.0]
    assert trials.stdev(2.0, data) == pytest.approx(1.0)     # no offset
    assert trials.stdev(4.0, data) == pytest.approx(math.sqrt(1.0 + 4.0))


def test_stdev_of_nothing_is_nothing():
    assert trials.stdev(1.0, []) is None


def test_kurtosis_of_a_flat_sample_is_the_excel_value():
    # For the six values -2,-1,0,0,1,2 about a mean of 0: the population
    # variance is 10/6, so sum((x-m)^2/v)^2 = 2*(2.4^2) + 2*(0.6^2) = 12.24,
    # and K = (12.24*6*7/5 - 3*25) / (4*3).
    data = [-2.0, -1.0, 0.0, 0.0, 1.0, 2.0]
    expected = (12.24 * 6 * 7 / 5 - 3 * 25) / 12
    assert trials.kurtosis(0.0, data) == pytest.approx(expected)


def test_kurtosis_needs_four_points_and_some_spread():
    assert trials.kurtosis(0.0, [1.0, 2.0, 3.0]) is None
    assert trials.kurtosis(1.0, [1.0] * 10) is None


# --------------------------------------------------------------------------
# the histograms
# --------------------------------------------------------------------------

@pytest.mark.parametrize('value, key', [
    (0.04, 0), (0.05, 1), (0.14, 1), (0.15, 2),
    (-0.04, 0), (-0.05, -1), (-0.15, -2), (1.0, 10),
])
def test_a_key_rounds_away_from_zero_keeping_its_sign(value, key):
    assert trials.rounded_key(10.0 * value) == key


def test_a_histogram_counts_by_tenths_of_a_sigma():
    h = trials.Histogram('x')
    for value in (0.0, 0.02, 0.11, -0.11):
        h.add(value)
    assert h.rows() == [(-0.1, 1), (0.0, 2), (0.1, 1)]
    assert h.values == []


def test_a_histogram_keeps_the_values_when_asked():
    h = trials.Histogram('x')
    h.add(0.3, keep=True)
    h.add(0.4)
    assert h.values == [0.3]
    assert bool(h)
