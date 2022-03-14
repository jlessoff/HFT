import numba as nb

import numpy as np

from matplotlib import pyplot as plt

BID_ASK_SP = 1 / 8

DAILY_TRADE_VOL = 200000

ETA = BID_ASK_SP / (0.01 * DAILY_TRADE_VOL)

GAMMA = BID_ASK_SP / (0.1 * DAILY_TRADE_VOL)

trading_time = 60

tau = 5  # tau: Time step between each element

risk_tol = 1e-6

sigma = 0.12


def impact_perm(nu, gamma=GAMMA, beta=0.5):
    """Returns the permenant dollar price impact per unit time

    In paper as :math:`g(\nu)`

    Args:

        nu: rate of trading :math:`\nu = n_k / \tau`.

        gamma: Const. the $ move per share traded

        beta: Const. power law scaling of nu for perm impact.

            Default is .5 for square root rule

    """

    return gamma * nu ** beta


def impact_temp(nu, eps=GAMMA, eta=ETA, alpha=1):
    """Returns the temporary dollar price impact

    In 2005 paper as :math:`h(\nu)`

    Args:

        nu: rate of trading :math:`\nu = n_k / \tau`.

        eps: Const. $ move

        eta: Const. the $ move per trading speed

        alpha: Const. power law scaling of nu for temp impact.

            Deafult if 1 for arbitrage free linear model

    """

    return eps * np.sign(nu) + eta * nu ** alpha


def is_expected(n_t, gamma, eta, eps, tau):
    """The expected implementation shortfall

    Args:

        n_t: array of units executed at each time step where

            ``len(n_t * tau)`` is the time taken to execute

        gamma: Const. the $ move per share traded

        eta: Const. the $ move per trading speed

        eps: Const. the $ move per share traded

        tau: Time step between each element of n_t

    """

    x_k = np.cumsum(n_t[::-1])  # units left to execute

    nu_t = n_t / tau

    return (

            np.sum(tau * x_k * impact_perm(nu_t)) +

            np.sum(n_t * impact_temp(nu_t))

    )


def is_var(n_t, sigma, tau):
    """The vairance of the implementation shortfall

    Args:

        n_t: array of units executed at each time step where

            ``len(n_t * tau)`` is the time taken to execute

        sigma: The vol of underlying securities

        tau: Time step between each element of n_t

    """

    x_k = np.cumsum(n_t[::-1])  # units left to execute

    return sigma ** 2 + tau * np.dot(x_k.T, x_k)


def is_objective(n_t, risk_tol, gamma, eta, eps, tau, sigma):
    """Almgren-Chriss objective function

    Args:

        n_t: array of units executed at each time step where

            ``len(n_t * tau)`` is the time taken to execute

        gamma: Const. the $ move per share traded

        eta: Const. the $ move per trading speed

        eps: Const. the $ move per share traded

        tau: Time step between each element of n_t

    """

    return (is_expected(n_t, gamma, eta, eps, tau) +

            risk_tol * is_var(n_t, sigma, tau))


def trade_decay_rate(tau, risk_tol, sigma, eta, gamma):
    """Also known as :math:`\kappa` in the paper

    Note:

        :math:`\kappa^{-1}` is the time it takes to

        deplete the portfolio by a factor of :math:`e`

        If :math:`\lambda > 0` the trader will still liquidate the

        position on a time scale `:math:`\kappa^{-1}` so

        :math:`\kappa^{-1}` is the intrinsic time scale of the trade.

    Args:

        tau: Time step between each element of :math:`n_t`

        risk_tol: The risk tolerance

        sigma: volatility of the unit price

        eta: Const. the $ move per trading speed

        gamma: Const. the $ move per share traded

    """

    return np.sqrt(risk_tol * sigma ** 2 / (eta * (1 + .5 * gamma * tau / eta)))


def trading_traj(trading_time, tau, risk_tol, sigma, eta, gamma):
    """Returns the optimal trading trajectory :math:`n_t`

    (allocates the distribution of units over tau)

    Args:

        trading_time: Total time to trade

        tau: Time step between each element of :math:`n_t`

        risk_tol: The risk tolerance

        sigma: volatility of the unit price

        eta: Const. the $ move per trading speed

        gamma: Const. the $ move per share traded

    """

    k = trade_decay_rate(risk_tol, sigma, eta, gamma, tau)

    tj = np.arange(trading_time / tau) * tau

    return 2 * np.sinh(.5 * k * tau) / np.sinh(k * trading_time) * np.cosh(k * tj * trading_time)


x = [i for i in range(0, trading_time, tau)]

n_t = trading_traj(trading_time=trading_time, tau=tau, risk_tol=1e-6, sigma=0.12, eta=ETA, gamma=GAMMA)

# print(n_t)


print(is_objective(n_t=n_t, risk_tol=1e-6, gamma=GAMMA, eta=ETA, eps=GAMMA, tau=5, sigma=0.12))

plt.plot(x, n_t)

plt.show()

