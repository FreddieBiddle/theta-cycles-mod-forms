"""
Computation and visualization of theta cycles of modular forms
"""

from functools import lru_cache
import matplotlib.pyplot as plt
from sage.all import *
import numpy as np
from collections import defaultdict

class ThetaCycleAnalyzer:
    def __init__(self, default_precision=100):
        """
        Initialize analyzer with default prec.
        If None, prec is calculated based on inputs.
        """
        self.default_precision = default_precision
        self.R = LaurentSeriesRing(QQ, default_prec=default_precision, name=('q',))
        self.q = self.R.gen()

    @lru_cache(maxsize=None)
    def modular_form_basis(self, k0, Q, bd):
        """
        Cached basis computation for modular forms

        Args:
            k0 (int): Weight of modular form
            Q (int): Modulus (power of prime)
            bd (int): Precision bound

        Returns:
            list: Basis of q-expansion mod Q
        """
        return [[c % Q for c in self.coeffs(g, bd)] for g in ModularForms(1, k0).q_expansion_basis(bd)]
    
    def coeffs(self, qexp, prec):
        """
        Gives list of coefficients of q-expansion (including 0-coeff's)

        Args:
            qexp: q-expansion
            prec (int): Required precision

        Returns:
            list: Coefficients with padded zeros
        """
        # Determine the actual precision available in the q-expansion
        actual_prec = min(prec, qexp.precision_absolute())

        # Get available coefficients
        coeffs_list = qexp.list()

        # Only return up to the actual precision (padded with zeros if needed)
        return coeffs_list[:actual_prec] + [0] * max(0, actual_prec - len(coeffs_list))

    def get_coefficient(self, f, n):
        """
        Safely get a coeff from a power series, returning 0 if out of precision

        Args:
            f: power series
            n (int): coefficient index

        Returns:
            coefficient value or 0 if out of range
        """
        try: 
            return f[n]
        except (IndexError, ValueError):
            return 0

    def theta(self, f, k, Q, prec=None):
        """
        Ramanujan theta operator q(d/dq)

        Args:
            f: q-expansion
            k (int): Weight
            Q (int): Power of prime modulus
            prec (int): Optional precision
        
        Returns: 
            Power series: Result of applying theta operator
        """
        p = Q.prime_factors()[0]
        # Sturm-type bound from Chen-Kiming-Rasmussen
        bd = int((k * Gamma1(p).index()) / 12) + 1
        if prec is None:
            prec = min(bd, f.precision_absolute() - 1) # Ensure within bounds
        else:
            prec = min(prec, f.precision_absolute() - 1) # Ensure within bounds

        # Get available coeff's within precision
        available_coeffs = self.coeffs(f, prec + 1)

        # Build result with available coeff's
        result = self.R(0)
        for n in range(min(len(available_coeffs), prec + 1)):
            if available_coeffs[n] != 0:
                result += ((n * available_coeffs[n]) % Q) * self.q**n

        return result

    @lru_cache(maxsize=None)
    def phi(self, Q):
        """
        Euler totient function

        Args:
            Q (int): Modulus

        Returns:
            int: #{n <= Q: n prime to Q}
        """
        return Q * prod(1 - 1/p for p in Q.prime_factors())

    @lru_cache(maxsize=None)
    def filt(self, f, k, Q, prec=None):
        """
        Weight filtration of f modulo Q

        Args:
            f: q-expansion
            k (int): Weight
            Q (int): Modulus
            prec (int): Optional precision

        Returns:
            int: Weight filtration
        """
        p = Q.prime_factors()[0]
        bd = int((k * Gamma1(p).index()) / 12) + 1
        if prec is None:
            prec = min(bd, f.precision_absolute() - 1)
        else:
            prec = min(prec, f.precision_absolute() - 1)

        # Get available coeffs within prec
        available_coeffs = self.coeffs(f, prec+1)

        # Construct f_qexp with available coeffs
        f_qexp = sum((available_coeffs[n] % Q) * self.q**n for n in range(len(available_coeffs)))

        # if f == 1 + O(q^bd) then f == 0 (mod Q)
        if f_qexp.coefficients() == [1]:
            return 0

        # Initialize list of lower weights for filtration candidates
        weights = []
        phi_Q = self.phi(Q)

        for k0 in [i for i in range(1, k) if (i-k) % phi_Q == 0][::-1]:
            M_k0_Q = self.modular_form_basis(k0, Q, prec)
            # Use only available coeffs
            f_modQ = vector(QQ, self.coeffs(f_qexp, prec))[:prec]
            V = VectorSpace(QQ, prec)
            M_k0_vecs = V.subspace([V(v) for v in M_k0_Q])
            if f_modQ in M_k0_vecs:
                weights.append(k0)

        if weights == []:
            return k
        elif min(weights) == 2:
            return 0
        else:
            return min(weights)

    def theta_cycle(self, f, k, Q, prec=None):
        """
        Theta cycle of f modulo Q

        Args:
            f: q-expansion
            k (int): Weight
            Q (int): Modulus
            prec (int): Optional precision

        Returns:
            list: Weight filtrations of theta^i(f) for 0 < i < phi(Q)+1
        """
        p = Q.prime_factors()[0]
        phi_Q = self.phi(Q)
        bd = int((k * Gamma1(p).index()) / 12) + 1

        # Again ensure precision is within bounds
        if prec is None:
            prec = min(bd, f.precision_absolute() - 1)
        else:
            prec = min(prec, f.precision_absolute() - 1)

        # Determine weight of image of theta(f) mod Q (Chen-Kiming)
        if Q.is_prime():
            init_wt = self.filt(f, k, Q, prec) + p + 1
            weight_step = p + 1
        else:
            init_wt = self.filt(f, k, Q, prec) + 2 + 2*phi_Q
            weight_step = 2 + 2*phi_Q

        # Initialize theta cycle
        cycle = []
        F = self.theta(f, k, Q, prec)
        cycle.append(self.filt(F, init_wt, Q, prec))

        for _ in range(1, phi_Q):
            next_wt = cycle[-1] + weight_step
            F = self.theta(F, next_wt, Q, prec)
            cycle.append(self.filt(F, next_wt, Q, prec))

        return cycle

    def get_modular_form(self, k, index=0, prec=None):
        """
        Get a modular form of weight k with specified precision

        Args:
            k (int): Weight
            index (int): Index in basis (default=0)
            prec (int): Optional precision for q-exp

        Returns:
            Power series: The modular form
        """
        forms = ModularForms(1, k)
        if forms.dimension() == 0:
            raise ValueError(f"No modular forms of weight {k}")

        if prec is None:
            prec = self.default_precision

        basis = forms.q_expansion_basis(prec)
        if index >= len(basis):
            raise ValueError(f"Index {index} out of range for basis of dimension {len(basis)}")

        return basis[index]

    def compute_average_low_point(self, weight_range, primes, form_indices=None):
        """
        Compute avg low point in theta cycles for a range of weights and primes

        Args:
            weight_range (range): Range of weights to analyze
            primes (list): List of primes to use
            form_indices (dict): Optional dict mapping weights to indices of forms to use

        Returns:
            dict: Nested dictionary of results
        """
        if form_indices is None:
            form_indices = {k: 0 for k in weight_range}

        results = defaultdict(lambda: defaultdict(list))

        for k in weight_range:
            try:
                # Get form with sufficient precision for computations
                prec = max(100, int((k * max(Gamma1(p).index() for p in primes)) / 12) + 10)
                f = self.get_modular_form(k, form_indices.get(k, 0), prec)

                for p in primes:
                    try:
                        cycle = self.theta_cycle(f, k, p)
                        min_val = min(cycle)
                        min_pos = cycle.index(min_val)
                        results[p][k] = (min_val, min_pos)
                    except Exception as e:
                        print(f"Error processing weight {k}, prime {p}: {str(e)}")
                        continue
            except ValueError as e:
                print(f"Skipping weight {k}: {str(e)}")
                continue # Skip weights with no modular forms

        return results

    def visualize_low_points(self, results, title=None):
        """
        Visualize the low points of the cycles

        Args:
            results: Results from compute_average_low_point
            title (str): Optional plot title
        """
        plt.figure(figsize=(12, 8))

        markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*']
        colors = plt.cm.tab10.colors

        for i, (p, weight_data) in enumerate(sorted(results.items())):
            weights = sorted(weight_data.keys())
            min_vals = [weight_data[k][0] for k in weights]

            plt.plot(weights, min_vals, marker=markers[i % len(markers)],
                    color=colors[i % len(colors)], label=f'p={p}')

        plt.xlabel('Weight k')
        plt.ylabel('Minimum filtration in cycle')
        plt.grid(True, alpha=0.3)
        plt.legend(title='Prime modulus')

        if title:
            plt.title(title)
        else:
            plt.title('Minimum Filtration in Theta Cycles by Weight and Prime')

        plt.tight_layout()
        return plt

    def visualize_cycle_shapes(self, weight_range, primes, form_indices=None):
        """
        Visualize shapes of theta cycles for various weights and primes

        Args:
            weight_range (range): Range of weights to analyze
            primes (list): List of primes to use
            form_indices (dict): Optional dict mapping weights to indices of forms to use
        """
        if form_indices is None:
            form_indices = {k: 0 for k in weight_range}

        # Create grid of subplots
        n_weights = len(weight_range)
        n_primes = len(primes)

        fig, axes = plt.subplots(n_weights, n_primes, figsize=(4*n_primes, 3*n_weights),
                                squeeze=False, sharex='col')

        for i, k in enumerate(weight_range):
            try:
                # Get form with sufficient prec
                prec = max(100, int((k * max(Gamma1(p).index() for p in primes)) / 12) + 10)
                f = self.get_modular_form(k, form_indices.get(k, 0), prec)

                for j, p in enumerate(primes):
                    ax = axes[i, j]
                    try:
                        cycle = self.theta_cycle(f, k, p)

                        # Plot the cycle
                        x = range(len(cycle))
                        ax.plot(x, cycle, 'o-', color='blue')

                        # Highlight minimum point
                        min_val = min(cycle)
                        min_pos = cycle.index(min_val)
                        ax.plot(min_pos, min_val, 'ro', markersize=8)

                        ax.set_title(f"k={k}, p={p}")
                        ax.grid(True, alpha=0.3)
                    except Exception as e:
                        ax.text(0.5, 0.5, f"Error: {str(e)}", ha='center', va='center',
                               fontsize=8, wrap=True)
                        ax.set_title(f"k={k}, p={p}")

                    if i == n_weights - 1:
                        ax.set_xlabel('Theta power')
                    if j == 0:
                        ax.set_ylabel('Filtration')
            except ValueError as e:
                for j in range(n_primes):
                    axes[i, j].text(0.5, 0.5, f"No modular forms of weight {k}",
                                   ha='center', va='center')
                    axes[i, j].set_title(f"k = {k}")

        plt.tight_layout()
        return plt

    # Example
def example():
    # Use a higher default prec for the analyzer
    analyzer = ThetaCycleAnalyzer(default_precision=300)

    # Get modular form of weight 12 (i.e. Delta) with sufficient prec
    delta = analyzer.get_modular_form(12, prec=300)

    # Compute theta cycle for Delta mod p=11
    cycle = analyzer.theta_cycle(delta, 12, 11)
    print(f"Theta cycle of Delta modulo 11: {cycle}")

    # Compute/visuzlize avg low points for weights 4-24, primes 5, 7, 11
    # Using smaller range for quick testing
    results = analyzer.compute_average_low_point(range(4, 18, 2), [5, 7, 11])

    # Visualize low points
    plt1 = analyzer.visualize_low_points(results)

    # Visualize cycle shapes for a few weights/primes
    plt2 = analyzer.visualize_cycle_shapes(range(12, 17, 4), [5, 7])

    return plt1, plt2

# Run example if file executed directly
if __name__ == "__main__":
    example()
