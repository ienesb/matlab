"""Deterministic, read-only checks of the manuscript's bound identities.

Run: python -B matlab/crb/test/check_bound_identities.py
Dependency: NumPy. No simulations, input files, or output files are used.

Exact finite-alphabet enumeration checks the mean, conjugate covariance,
complementary covariance, and real-score covariance for Gray-QPSK/BSC pairs;
the single-target no-ICI covariance reduction in both TF and DD; and the
population mean-projection loss. A deterministic positive-definite FIM checks
unknown-gain Schur-complement and inverse-block attenuation scaling.

All matrix comparisons use rtol=2e-11 and atol=2e-10. These tolerances include
roundoff from summing all 65,536 states in the four-symbol check. Data and bit
flips are fully enumerated, not Monte Carlo sampled. Random-looking operators
are fixed by seed 9273 and do not introduce sampling uncertainty.
"""

import numpy as np


RTOL = 2e-11
ATOL = 2e-10


def check(actual, expected, name):
    np.testing.assert_allclose(
        actual, expected, rtol=RTOL, atol=ATOL, err_msg=name
    )


def enumerate_pairs(n, p):
    """All 4**n data vectors crossed with all 2**(2*n) bit-flip vectors."""
    codes = np.arange(1 << (2 * n), dtype=np.int64)
    bits = (codes[:, None] >> np.arange(2 * n)) & 1
    s = ((1 - 2 * bits[:, 0::2]) + 1j * (1 - 2 * bits[:, 1::2])) / np.sqrt(2)
    detected_bits = bits[:, None, :] ^ bits[None, :, :]
    sh = (
        (1 - 2 * detected_bits[:, :, 0::2])
        + 1j * (1 - 2 * detected_bits[:, :, 1::2])
    ) / np.sqrt(2)
    counts = bits.sum(axis=1)
    flip_weights = p**counts * (1 - p) ** (2 * n - counts)
    weights = np.broadcast_to(flip_weights / len(bits), sh.shape[:2]).ravel()
    true = np.broadcast_to(s[:, None, :], sh.shape).reshape(-1, n)
    check(weights.sum(), 1.0, "enumeration probabilities")
    return true, sh.reshape(-1, n), weights


def weighted_cov(values, weights, complementary=False):
    centered = values - weights @ values
    right = centered if complementary else centered.conj()
    return np.einsum("t,tk,tl->kl", weights, centered, right)


def covariance_formula(kernels, p):
    size = len(kernels)
    conjugate = np.empty((size, size), dtype=complex)
    complementary = np.empty_like(conjugate)
    for k in range(size):
        for l in range(size):
            left, right = kernels[k], kernels[l]
            dl, dr = np.diag(left), np.diag(right)
            conjugate[k, l] = (
                4 * p * np.trace(left @ right.conj().T)
                - 4 * p**2 * np.sum(dl * dr.conj())
            )
            complementary[k, l] = 4 * p**2 * (
                np.trace(left @ right) - np.sum(dl * dr)
            )
    return conjugate, complementary


def check_general_covariance(rng, sigma2):
    kernels = rng.normal(size=(3, 2, 2)) + 1j * rng.normal(size=(3, 2, 2))
    for p in (0.0, 0.03, 0.17, 0.43):
        true, detected, weights = enumerate_pairs(2, p)
        error = true - detected
        bilinear = np.einsum("ti,kij,tj->tk", error.conj(), kernels, detected)
        check(
            weights @ bilinear,
            -2 * p * np.trace(kernels, axis1=1, axis2=2),
            "BSC bilinear mean",
        )
        conjugate, complementary = covariance_formula(kernels, p)
        check(weighted_cov(bilinear, weights), conjugate, "conjugate covariance")
        check(
            weighted_cov(bilinear, weights, complementary=True),
            complementary,
            "complementary covariance",
        )
        score = 2 / sigma2 * bilinear.real
        expected = 2 / sigma2**2 * (conjugate + complementary).real
        check(weighted_cov(score, weights), expected, "real-score covariance")
    print("PASS: general BSC moments, 256 states at each of four error rates")


def check_single_target(sigma2):
    p, alpha, ts = 0.17, 0.8 + 0.3j, 1.2
    true, detected, weights = enumerate_pairs(4, p)
    n, m = np.tile(np.arange(2), 2), np.repeat(np.arange(2), 2)
    phase = np.diag(np.exp(-2j * np.pi * n * 0.23 + 2j * np.pi * m * ts * 0.15))
    channel = alpha * phase
    derivatives = np.array([
        phase, 1j * phase,
        -2j * np.pi * alpha * phase @ np.diag(n),
        2j * np.pi * alpha * ts * phase @ np.diag(m),
    ])
    kernels = np.array([channel.conj().T @ derivative for derivative in derivatives])
    fim = 2 / sigma2 * np.einsum("kij,lij->kl", derivatives.conj(), derivatives).real
    expected = 4 * p * (1 - p) * abs(alpha)**2 / sigma2 * fim
    f2 = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    unitary = np.kron(f2, f2)
    dd_kernels = np.array([unitary.conj().T @ k @ unitary for k in kernels])
    for name, current in (("OFDM", kernels), ("OTFS", dd_kernels)):
        bilinear = np.einsum("ti,kij,tj->tk", (true - detected).conj(), current, detected)
        score = 2 / sigma2 * bilinear.real
        check(weighted_cov(score, weights), expected, name + " single-target reduction")
    print("PASS: single-target no-ICI covariance reduction, 65,536 states per basis")


def check_attenuation(rng, sigma2):
    gradient = rng.normal(size=(18, 8)) + 1j * rng.normal(size=(18, 8))
    fim = 2 / sigma2 * np.real(gradient.conj().T @ gradient)
    c, gain_count = 0.66, 4
    scaling = np.diag([1.0] * gain_count + [c] * 4)
    attenuated = scaling @ fim @ scaling

    def schur(matrix):
        a = matrix[:gain_count, :gain_count]
        b = matrix[:gain_count, gain_count:]
        return matrix[gain_count:, gain_count:] - b.T @ np.linalg.solve(a, b)

    check(schur(attenuated), c**2 * schur(fim), "nuisance-gain Schur scaling")
    check(
        np.linalg.inv(attenuated)[gain_count:, gain_count:],
        np.linalg.inv(fim)[gain_count:, gain_count:] / c**2,
        "kinematic block of full inverse",
    )
    h0 = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    h1 = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    true, detected, weights = enumerate_pairs(2, (1 - c) / 2)
    residual = true @ h0.T - detected @ h1.T
    loss = weights @ np.sum(abs(residual)**2, axis=1)
    expected = np.linalg.norm(h1 - c * h0)**2 + (1 - c**2) * np.linalg.norm(h0)**2
    check(loss, expected, "population mean-projection loss")
    print("PASS: population projection and unknown-gain attenuation identities")


if __name__ == "__main__":
    rng = np.random.default_rng(9273)
    check_general_covariance(rng, sigma2=0.37)
    check_single_target(sigma2=0.37)
    check_attenuation(rng, sigma2=0.37)
    print(f"All identities passed (rtol={RTOL:g}, atol={ATOL:g}).")
