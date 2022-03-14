import numba as nb
import statistics
from SALib.sample import saltelli
from SALib.analyze import sobol


import numpy as np
import matplotlib.pyplot as plt
def trading_traj_discrete(q0, trading_time,tau, alpha):
    front=q0*(np.sinh(alpha*(trading_time-tau))/np.sinh(alpha*trading_time))
    return front
y=[]
for i in range(299,301):
    curve=(trading_traj_discrete(50000,300,i,0.005))
    y.append(curve)
print(y)
plt.plot(y)
plt.show()


def trading_traj_discrete_nonhom(xk, trading_time, tau, k, tauk):
    front_nonhom=xk*(np.sinh(k*(trading_time-tau)))/(np.sinh(k*(trading_time-tauk)))
    return front_nonhom


y=[]
for i in range(0,300):
    curve=(trading_traj_discrete_nonhom(50000,300,i,.005,0))
    y.append(curve)
print(y)
plt.plot(y)
plt.show()


def trading_traj_continuous(q0, trading_time,tau, sigma, eta, gamma):
    """tau- time
       q0- quantity
       theta- volatility- volatility increases, the trader wants to go faster
       eta-scale parameter for execution costs paid by trader, when increases, will trade slower
       gamma-risk adversion - as gamma incrases, trader moves faster.  faster in beginning when gamma large

    """
    #linear volume
    vol=((2.5*0.000001)*tau)+6.5*0.0001
    #u shaped volume
    vol2=((2.5*0.000001)*tau)**2+((6.5*0.001)*tau)+2.5*0.0001
    #v shaped volume
    vol3= 0.00001*abs(0.001-tau)+6.5*0.0001
    #volume iwth jump at time t
    vol4= 0.0009+vol if  200 <= tau  else vol



    qt= (gamma * (sigma**2) * vol)/(2*eta)
    qt2= (gamma * (sigma**2) * vol2)/(2*eta)
    qt3= (gamma * (sigma**2) * vol3)/(2*eta)
    qt4= (gamma * (sigma**2) * vol4)/(2*eta)


    front=q0*((np.sinh((np.sqrt(qt))*(trading_time-tau)))/((np.sinh(np.sqrt(qt)*trading_time))))
    front2=q0*((np.sinh((np.sqrt(qt2))*(trading_time-tau)))/((np.sinh(np.sqrt(qt2)*trading_time))))
    front3=q0*((np.sinh((np.sqrt(qt3))*(trading_time-tau)))/((np.sinh(np.sqrt(qt3)*trading_time))))
    front4=q0*((np.sinh((np.sqrt(qt4))*(trading_time-tau)))/((np.sinh(np.sqrt(qt4)*trading_time))))

    return  front,front2, front3, front4





y=[]
for i in range(0,301):
    curve=(trading_traj_continuous(50000,300,i,0.5,0.00005,0.000008))
    y.append(curve)
print(y)
plt.plot(y)

plt.show()

def trading_traj_cont_mult(zj0, kj, trading_time, tau):
    front_nonhom=zj0*(np.sinh(kj*(trading_time-tau)))/(np.sinh(kj*(trading_time)))
    return front_nonhom

y=[]
for i in range(0,600):
    curve=(trading_traj_cont_mult(50000,0.05,600,i))
    y.append(curve)
print(y)
plt.plot(y)
plt.show()



problem = {'num_vars': 3,
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