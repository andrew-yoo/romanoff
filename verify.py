from __future__ import annotations

import argparse
import hashlib
import json
from bisect import bisect_right
from collections import Counter
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
from itertools import product
from math import lcm, prod
from pathlib import Path
import sys

BASE_PRIMES = [3, 5, 7, 11, 13, 17, 19, 29, 31, 37, 41, 61, 73]

ADDED_CLUSTERS = {
    11: [23, 89],
    13: [8191],
    23: [47],
    29: [233, 1103, 2089],
    37: [223],
    43: [431, 9719],
    47: [2351, 4513],
    53: [6361],
    73: [439],
    79: [2687],
    83: [167],
    131: [263],
    179: [359, 1433],
    191: [383],
    233: [1399],
    239: [479, 1913],
    251: [503],
}


ALL_PRIMES = sorted(BASE_PRIMES + [q for qs in ADDED_CLUSTERS.values() for q in qs])

G2026_UPPER = Fraction(490249407811155, 10**15)


@dataclass(frozen=True)
class CertificateResult:
    number_of_primes: int
    M: int
    phi_M: int
    T: int
    seed_histogram_sha256: str
    seed_histogram_entries: int
    added_multiplier_entries: int
    rational_upper_sum_numerator: int
    rational_upper_sum_denominator: int
    rational_upper_sum_decimal: str
    g2026_upper_bound_decimal: str
    improvement_decimal: str


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_seed_histogram(path: Path) -> dict[int, int]:
    hist: dict[int, int] = {}
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 2:
            raise ValueError(f"{path}:{lineno}: expected two columns, got {raw!r}")
        nu = int(parts[0])
        count = int(parts[1])
        if not (0 <= nu <= 2520):
            raise ValueError(f"{path}:{lineno}: nu out of range: {nu}")
        if count <= 0:
            raise ValueError(f"{path}:{lineno}: nonpositive count: {count}")
        if nu in hist:
            raise ValueError(f"{path}:{lineno}: duplicate nu value: {nu}")
        hist[nu] = count
    if not hist:
        raise ValueError(f"{path}: empty histogram")
    return hist


def ord_mod_2(q: int) -> int:
    x = 2 % q
    d = 1
    while x != 1:
        x = (x * 2) % q
        d += 1
    return d


def cluster_polynomial(order: int, primes: list[int]) -> dict[int, int]:
    d = order
    dist: dict[int, int] = {0: 1}
    for q in primes:
        actual = ord_mod_2(q)
        require(actual == d, f"order mismatch for q={q}: got {actual}, expected {d}")
        new: Counter[int] = Counter()
        for s, count in dist.items():
            non_power = q - d
            if non_power:
                new[s] += count * non_power
            if s:
                new[s] += count * s
            if s < d:
                new[s + 1] += count * (d - s)
        dist = dict(new)
    poly = {d - s: count for s, count in dist.items()}
    require(sum(poly.values()) == prod(primes), "cluster polynomial total weight mismatch")
    return poly


def brute_force_cluster_polynomial(order: int, primes: list[int]) -> dict[int, int]:
    power_residues = []
    for q in primes:
        require(ord_mod_2(q) == order, f"order mismatch in brute force cluster for q={q}")
        residues = []
        x = 1 % q
        for _ in range(order):
            residues.append(x)
            x = (x * 2) % q
        power_residues.append(residues)
    counts: Counter[int] = Counter()
    for residue_tuple in product(*(range(q) for q in primes)):
        h = 0
        for r in range(order):
            if all(residue_tuple[i] != power_residues[i][r] for i in range(len(primes))):
                h += 1
        counts[h] += 1
    return dict(counts)


def multiplier_distribution() -> Counter[int]:
    mult: Counter[int] = Counter({1: 1})
    for order, primes in ADDED_CLUSTERS.items():
        poly = cluster_polynomial(order, primes)
        new: Counter[int] = Counter()
        for m, weight in mult.items():
            for h, c in poly.items():
                new[m * h] += weight * c
        mult = new
    return mult


def log2_lower_bound(terms: int = 40) -> Fraction:
    return sum(Fraction(2, (2*j + 1) * 3**(2*j + 1)) for j in range(terms))


def decimal(frac: Fraction, prec: int = 100) -> Decimal:
    getcontext().prec = prec
    return Decimal(frac.numerator) / Decimal(frac.denominator)


def run_self_tests() -> None:
    for order, primes in (
        (11, [23, 89]),
        (23, [47]),
        (83, [167]),
    ):
        expected = brute_force_cluster_polynomial(order, primes)
        actual = cluster_polynomial(order, primes)
        require(actual == expected, f"cluster polynomial self-test failed for order {order}")

    lb10 = log2_lower_bound(10)
    lb20 = log2_lower_bound(20)
    lb40 = log2_lower_bound(40)
    require(lb10 < lb20 < lb40, "log2 lower bound should increase with more terms")
    require(Fraction(69, 100) < lb40 < Fraction(70, 100), "log2 lower bound sanity check failed")


def verify_certificate(seed_path: Path) -> CertificateResult:
    seed_digest = hashlib.sha256(seed_path.read_bytes()).hexdigest()
    base_histogram_13_prime = load_seed_histogram(seed_path)
    require(set(base_histogram_13_prime) <= set(range(2521)), "seed histogram contains nu > 2520")
    require(
        sum(base_histogram_13_prime.values()) == prod(BASE_PRIMES),
        "seed histogram total does not match base prime product",
    )

    M = prod(ALL_PRIMES)
    phi_M = prod(q - 1 for q in ALL_PRIMES)
    T = 1
    for q in ALL_PRIMES:
        T = lcm(T, ord_mod_2(q))

    mult = multiplier_distribution()
    added_prime_product = prod(q for qs in ADDED_CLUSTERS.values() for q in qs)
    require(sum(mult.values()) == added_prime_product, "multiplier distribution total mismatch")
    require(prod(BASE_PRIMES) * added_prime_product == M, "prime product mismatch")

    items = sorted(mult.items())
    factors = [h for h, _ in items]
    prefix_weight = [0]
    prefix_weighted_factor = [0]
    for h, weight in items:
        prefix_weight.append(prefix_weight[-1] + weight)
        prefix_weighted_factor.append(prefix_weighted_factor[-1] + weight * h)
    total_weight = prefix_weight[-1]

    log2_lb = log2_lower_bound(40)
    trivial_term = Fraction(1, 2 * M)
    prime_progression_scale = Fraction(1, T * phi_M) / log2_lb

    S_upper = Fraction(0, 1)
    for nu, delta in base_histogram_13_prime.items():
        if nu == 0:
            continue
        threshold = trivial_term / (prime_progression_scale * nu)
        threshold_floor = threshold.numerator // threshold.denominator
        idx = bisect_right(factors, threshold_floor)
        small_weighted_factor = prefix_weighted_factor[idx]
        large_weight = total_weight - prefix_weight[idx]
        S_upper += delta * (
            prime_progression_scale * nu * small_weighted_factor
            + trivial_term * large_weight
        )

    require(S_upper < G2026_UPPER, "39-prime bound does not improve on G2026 (36-prime) bound!")

    return CertificateResult(
        number_of_primes=len(ALL_PRIMES),
        M=M,
        phi_M=phi_M,
        T=T,
        seed_histogram_sha256=seed_digest,
        seed_histogram_entries=len(base_histogram_13_prime),
        added_multiplier_entries=len(mult),
        rational_upper_sum_numerator=S_upper.numerator,
        rational_upper_sum_denominator=S_upper.denominator,
        rational_upper_sum_decimal=str(decimal(S_upper)),
        g2026_upper_bound_decimal=str(decimal(G2026_UPPER, 30)),
        improvement_decimal=str(decimal(G2026_UPPER - S_upper, 30)),
    )


def render_result(result: CertificateResult) -> str:
    return "\n".join(
        [
            "39-prime Romanoff upper-density certificate",
            f"number of primes = {result.number_of_primes}",
            f"M = {result.M}",
            f"phi(M) = {result.phi_M}",
            f"T = {result.T}",
            f"seed histogram sha256 = {result.seed_histogram_sha256}",
            f"seed histogram entries = {result.seed_histogram_entries}",
            f"added multiplier entries = {result.added_multiplier_entries}",
            f"rational upper sum = {result.rational_upper_sum_decimal}",
            f"proved: upper density < {result.rational_upper_sum_decimal[:17]}",
            f"improvement over G2026 bound = {result.improvement_decimal}",
        ]
    )


def write_json_summary(result: CertificateResult, path: Path) -> None:
    payload = asdict(result)
    payload["rational_upper_sum"] = {
        "numerator": result.rational_upper_sum_numerator,
        "denominator": result.rational_upper_sum_denominator,
        "decimal": result.rational_upper_sum_decimal,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Verify the C45 Romanoff upper-density certificate with 39 primes.")
    parser.add_argument(
        "seed_histogram",
        nargs="?",
        default="seed_histogram_13_prime.txt",
        type=Path,
        help="Path to the generated 13-prime seed histogram.",
    )
    parser.add_argument(
        "--json-out",
        type=Path,
        default=None,
        help="Optional path for a machine-readable certificate summary.",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run exact helper self-tests before verifying the certificate.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.self_test:
        run_self_tests()
        print("Python verifier self-tests passed", file=sys.stderr)

    result = verify_certificate(args.seed_histogram)
    print(render_result(result))
    if args.json_out is not None:
        write_json_summary(result, args.json_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
