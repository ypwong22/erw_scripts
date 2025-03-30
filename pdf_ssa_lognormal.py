"""
I tested exponential, gamma, and lognormal distributions. I also tested setting x0 to different starting values (1e-6, 1e-4, 1e-2). 

Lognormal distribution gives the most stable result.
"""

import numpy as np
from scipy.special import gamma


def pdf_expo(x, lam):
    return lam * np.exp(-lam*x)


def pdf_gam(x, theta):
    k = 0.5
    return (x**(k - 1) * np.exp(-x/theta)) / (theta**k * gamma(k))


def pdf_log(x, mu):
    sigma = 2 * np.log(mu + 1)
    return (1 / (x * sigma * np.sqrt(2 * np.pi))) * np.exp(-((np.log(x) - mu) ** 2) / (2 * sigma ** 2))

def ssa(x):
    if x < 1e-12:
        return 0
    else:
        return 69.18 * (x **(-1.24))
        # return 4 * np.pi * (x / 1e4)**2 * (1e4 * x) ** 0.33


forc_gra = 107

# lognormal distribution parameter
mu = np.log(forc_gra)

# gamma distribution parameter
theta = 4.394

# numerical integration setup
n_discret = 1000
x0 = 1
x1 = 16000

#dx_log = (np.log(x1+1) - np.log(x0+1)) / n_discret
#x_range = np.array([np.exp(np.log(x0+1) + ii*dx_log) for ii in range(n_discret)])
dx = (x1-x0) / n_discret
x_range = [x0 + ii*dx for ii in range(n_discret)]

# integration procedure using Simpson's rule
num = ssa(x0) * pdf_log(x0, mu)
prob = pdf_log(x0, mu)

for ii in range(1, n_discret-1):
    if np.mod(ii, 2) == 1:
        num += 4 * ssa(x_range[ii]) * pdf_log(x_range[ii], mu)
        prob += 4 * pdf_log(x_range[ii], mu)
    else:
        num += 2 * ssa(x_range[ii]) * pdf_log(x_range[ii], mu)
        prob += 2 * pdf_log(x_range[ii], mu)

num += ssa(x1) * pdf_log(x1, mu)
prob += ssa(x1) * pdf_log(x1, mu)

num *= (dx/3)
prob *= (dx/3)

# estimated SSA over the distribution
ssa_mean = num / prob

# print the outcomes
print(ssa_mean, num, prob)
print(ssa(forc_gra))