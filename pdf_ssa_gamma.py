import numpy as np
from scipy.special import gamma


def gamma_pdf(x, theta):
    """
    Gamma PDF with shape k=0.5, and scale theta.
    """
    k = 0.5
    return(x**(k - 1) * np.exp(-x/theta)) / (theta**k * gamma(k))


def u_pdf(u, theta):
    """ Gamma PDF with k = 0.5 transformed by x = u**2
        pdf(x)dx = pdf(u^2) * 2u du
        note gamma(0.5) = np.sqrt(pi)
    """
    return 2 / np.sqrt(theta) / np.sqrt(np.pi) * np.exp(-u**2 / theta)


def ssa(x):
    """ Area of a sphere adjusted by roughtness, divided by density and volume of the sphere """
    # return np.where(x<1e-8, 0, 4 * np.pi * (x / 1e4)**2 * (1e4 * x) ** 0.33)
    # return np.where(x<1e-8, 0, 20.89 * (x**(-0.67)))
    return np.where(x<1e-8, 0, 69.18 * (x **(-1.24)))


def expected_ssa(theta, gra, N=10000):
    """
    Compute the expected value of ssa(x) over the gamma PDF using a change of variable.
    
    Parameters
    ----------
    theta : float
        Scale parameter of the gamma PDF.
    x_max : float, optional
        Upper limit for the integration in x-space.
    N : int, optional
        Number of subintervals in u-space (should be even for Simpson's rule).

    Returns
    -------
    float
        Expected value E[ssa(x)].
    """
    # Set up u-space limits: x = u^2 so u goes from 0 to sqrt(x_max)

    # note: if we start from a small u_min, e.g. 1e-3, probability can integrate to one, but we 
    #       will get exploding grain sizes. The model is therefore sensitive to this function
    u_min = np.sqrt(0.01 * gra) # start from 1
    u_max = np.sqrt(100 * gra)

    # Ensure N is even
    if N % 2 == 1:
        N += 1

    # Create grid in u-space
    u = np.linspace(u_min, u_max, N+1)
    du = (u_max - u_min) / N

    # Compute x from u
    x = u**2

    # Transform the integrals: dx = 2u du
    # Numerator: ssa(x)*pdf(x)*dx = ssa(u^2)*pdf(u^2)*2u du
    num_integrand = ssa(x) * u_pdf(u, theta)
    # Denominator: pdf(x)*dx = pdf(u^2)*2u du
    den_integrand = u_pdf(u, theta)

    # Simpson's rule integration function
    def simpson(y):
        return (du / 3.0) * (y[0] + 4.0 * np.sum(y[1:-1:2]) + 2.0 * np.sum(y[2:-1:2]) + y[-1])

    num_int = simpson(num_integrand)
    den_int = simpson(den_integrand)

    return num_int / den_int

if __name__ == "__main__":
    theta_value = 4.394
    gra = 107
    expected_value = expected_ssa(theta_value, gra, N=1000)
    print("Expected value of ssa(x) over the gamma PDF:", expected_value)
