from SALib.sample import saltelli
from SALib.analyze import sobol
from matplotlib import pyplot as plt
from SALib.test_functions import Ishigami
import numpy as np
# def perm_impact(q0,q,t,d,v):
#     one=(abs(q0-(q*t)))/(((q**(t/3))))
#     two=(v*t)*(d*t)
#     perm_impact=one*two
#     return perm_impact
#
# print(perm_impact(500,20,300,))


def block_price(q,s,k,psi,gamma,theta,V,n):
    """S -   price
       q- quantity
       k-perm market impact factor -higher impact, higher premium
       psi -proportional fees and width of bid ask sprad -higher fee, higher premium
       gamma-risk adversion -premium is increasing function of gamma..
       theta- volatility- premim increasing function of vol
       psi- measurement of convex behavior
       V - liquid - more liquidity, lower premium
       n-higher execution cost, higher premium

    """
    #mtm value
    mtm_val=q*s

    #permanent market impact
    perm_impact=(k/2)*(q**2)

    #linear execution costs
    one=psi*q
    #nonlinear execution costs
    two=(n**(1/(1+psi)))/(psi**(psi/(psi+1)))
    three=((1+psi)**2)/(1+3*psi)
    four=((gamma*theta**2)/(2*V))**(psi/(1+psi))
    five=q**((1+(3*psi))/(1+psi))
    exec_cost=((two)*(three)*(four)*(five))

    blockprice=mtm_val-perm_impact-one-exec_cost
    return blockprice

problem = {'num_vars': 6,
           'names': ['k','psi','gamma','theta','V','n'],
           'bounds': [[0,1],
                      [0, 1],
                      [0, 1],
                      [0.0,1],
                      [0,60000000],
                      [0,40000000]]
           }

vals = saltelli.sample(problem, 2^6)
# print(vals)
# Y=block_price(100,45,0.0001,0.75,0.0000051,0.21,3000000,20000000)
y = np.array([block_price(1000,45,*params) for params in vals])
print(y)
sobol_indices = sobol.analyze(problem, y, print_to_console=True)
print(sobol_indices)
print(sobol_indices['ST'])
print(sobol_indices['S1'])

sobol_indices.plot()
plt.show()





def block_price_constraint(q,s,k,psi,gamma,theta,V,n):
    """S -   price
       q- quantity
       k-perm market impact factor -higher impact, higher premium
       psi -proportional fees and width of bid ask sprad -higher fee, higher premium
       gamma-risk adversion -premium is increasing function of gamma..
       theta- volatility- premim increasing function of vol
       psi- measurement of convex behavior
       V - liquid - more liquidity, lower premium
       n-higher execution cost, higher premium

    """
    #mtm value
    mtm_val=q*s

    #permanent market impact
    perm_impact=(k/2)*(q**2)

    #linear execution costs
    lin_ex=psi*q
    one=(n**(1/(1+theta)))
    two= ((gamma*theta**2)/(6*theta*V))**(theta/1+theta)
    three=q**((3*theta+1)/(1+theta))
    exec_cost=one*two*three

    #price risk
    one=(theta*n)**(1/(1+theta))
    two= ((gamma*theta**2)/(6*theta*V))**(theta/1+theta)
    three=q**((3*theta+1)/(1+theta))
    price_risk= one*two*three
    blockprice=mtm_val-perm_impact-lin_ex-exec_cost-price_risk
    return blockprice


Y=block_price(100,45,0.0001,0.75,0.0000051,0.21,3000000,20000000)
print(Y)
Y=block_price_constraint(100,45,0.0001,0.75,0.0000051,0.21,3000000,20000000)
print(Y)


problem = {'num_vars': 6,
           'names': ['k','psi','gamma','theta','V','n'],
           'bounds': [[0,1],
                      [0, 1],
                      [0, 1],
                      [0,1],
                      [0,60000000],
                      [0,40000000]]
           }

vals = saltelli.sample(problem, 2^6)
# print(vals)
# Y=block_price(100,45,0.0001,0.75,0.0000051,0.21,3000000,20000000)
y = np.array([block_price_constraint(1000,45,*params) for params in vals])
print(y)
sobol_indices = sobol.analyze(problem, y, print_to_console=True)
print(sobol_indices)
print(sobol_indices['S1'])

sobol_indices.plot()
plt.show()




#
# sobol_indices = [sobol.analyze(problem, Y,calc_second_order=False) for Y in y.T]
# S1s = np.array([s['S1'] for s in sobol_indices])
#
# fig = plt.figure(figsize=(10, 6), constrained_layout=True)
# gs = fig.add_gridspec(2, 2)
#
# ax0 = fig.add_subplot(gs[:, 0])
# ax1 = fig.add_subplot(gs[0, 1])
# ax2 = fig.add_subplot(gs[1, 1])
#
# for i, ax in enumerate([ax1, ax2]):
#     ax.plot(x, S1s[:, i],
#             label=r'S1$_\mathregular{{{}}}$'.format(problem["names"][i]),
#             color='black')
#     ax.set_xlabel("x")
#     ax.set_ylabel("First-order Sobol index")
#
#     ax.set_ylim(0, 1.04)
#
#     ax.yaxis.set_label_position("right")
#     ax.yaxis.tick_right()
#
#     ax.legend(loc='upper right')
#
# ax0.plot(x, np.mean(y, axis=0), label="Mean", color='black')
#
# # in percent
# prediction_interval = 95
#
# ax0.fill_between(x,
#                  np.percentile(y, 50 - prediction_interval/2., axis=0),
#                  np.percentile(y, 50 + prediction_interval/2., axis=0),
#                  alpha=0.5, color='black',
#                  label=f"{prediction_interval} % prediction interval")
#
# ax0.set_xlabel("x")
# ax0.set_ylabel("y")
# ax0.legend(title=r"$y=a+b\cdot x^2$",
#            loc='upper center')._legend_box.align = "left"
#
# plt.show()
# analyse
# for i in range(len(vals)):
#   Y[i][:] = sp.odeint(iron,[vals[i][6],vals[i][7],vals[i][8]],t,\
#     args=(vals[i][0],vals[i][1],vals[i][2],vals[i][3],vals[i][4],vals[i][5],I))[len(XYZ)-1]
#








#
#
#
# x = np.linspace(-1, 1, 100)
# y = np.array([block_price(x, *params) for params in param_values])
# sobol_indices = [sobol.analyze(problem, Y) for Y in y.T]
# print(sobol_indices)
# S1s = np.array([s['S1'] for s in sobol_indices])
#
# fig = plt.figure(figsize=(10, 6), constrained_layout=True)
# gs = fig.add_gridspec(2, 2)
#
# ax0 = fig.add_subplot(gs[:, 0])
# ax1 = fig.add_subplot(gs[0, 1])
# ax2 = fig.add_subplot(gs[1, 1])
#
# for i, ax in enumerate([ax1, ax2]):
#     ax.plot(x, S1s[:, i],
#             label=r'S1$_\mathregular{{{}}}$'.format(problem["names"][i]),
#             color='black')
#     ax.set_xlabel("x")
#     ax.set_ylabel("First-order Sobol index")
#
#     ax.set_ylim(0, 1.04)
#
#     ax.yaxis.set_label_position("right")
#     ax.yaxis.tick_right()
#
#     ax.legend(loc='upper right')
#
# ax0.plot(x, np.mean(y, axis=0), label="Mean", color='black')
#
# # in percent
# prediction_interval = 95
#
# ax0.fill_between(x,
#                  np.percentile(y, 50 - prediction_interval/2., axis=0),
#                  np.percentile(y, 50 + prediction_interval/2., axis=0),
#                  alpha=0.5, color='black',
#                  label=f"{prediction_interval} % prediction interval")
#
# ax0.set_xlabel("x")
# ax0.set_ylabel("y")
# ax0.legend(title=r"$y=a+b\cdot x^2$",
#            loc='upper center')._legend_box.align = "left"
#
# plt.show()
# Y=block_price(100,45,0.0001,0.75,0.0000051,0.21,3000000,20000000)
#
# # Si = sobol.analyze(problem, Y, print_to_console=True)
#
# # print(block_price(100,45,0.0001,0.75,0.0000051,0.21,3000000,20000000))
#
#
# # x = np.linspace(-1, 1, 100)
# # y = np.array([block_price(x, *params) for params in param_values])
# # sobol_indices = [sobol.analyze(problem, Y) for Y in y.T]
#
# # def execution_cost(q0, psi,n, gamma, theta, V):
# #     one=psi*q0-(n**(1/(1+theta)))
# #     two= ((gamma*theta**2)/(6*theta*V))**(theta/1+theta)
# #     three=q0**((3*theta+1)/(1+theta))
# #     execution_cost=one-two*three
# #     return execution_cost
# #
# # def price_(q0, psi,n, gamma, theta, V):
# #     one=(-theta*n)**(1/(1+theta))
# #     two= ((gamma*theta**2)/(6*theta*V))**(theta/1+theta)
# #     three=q0**((3*theta+1)/(1+theta))
#
#
#
# # def perm_market_impact(theta, delta, W, vt, q0, q,T):
# #
# #     for t in T:
# #
# #         one= (theta*delta*W)
# #         two=
# #
# #
# #
# #
