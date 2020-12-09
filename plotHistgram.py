#https://jb102.blogspot.com/2017/10/22-histogram.html
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from scipy.stats import uniform
from scipy.stats import poisson
from scipy.stats import beta
def plot_uniform(min,max):
    f1=open('uniform.csv',"r")
    x_uniform=f1.read().split('\n')[:-1]
    x_uniform=[float(s) for s in x_uniform]
    X1=np.arange(min-1,max+1,0.1)
    Y1=uniform.pdf(X1,min,max)
    plt.hist(x_uniform,bins=30,density=True)
    plt.plot(X1,Y1,color='r')
    plt.show()

def plot_norm(mu,sigma):
    f2=open('norm.csv',"r")
    x_norm=f2.read().split('\n')[:-1]
    x_norm=[float(s) for s in x_norm]
    X2=np.arange(mu-4*sigma,mu+4*sigma,0.1)
    Y2=norm.pdf(X2,mu,sigma)
    plt.hist(x_norm,bins=30,density=True)
    plt.plot(X2,Y2,color='r')
    plt.show()

def plot_beta(a,b):
    f3=open('beta.csv',"r")
    x_beta=f3.read().split('\n')[:-1]
    x_beta=[float(s) for s in x_beta]
    X3=np.arange(0,1,0.01)
    Y3=beta.pdf(X3,a,b)
    plt.hist(x_beta,bins=30,density=True)
    plt.plot(X3,Y3,color='r')
    plt.show()

def plot_poisson(lamb):
    f4=open('poisson.csv',"r")
    x_poisson=f4.read().split('\n')[:-1]
    x_poisson=[float(s) for s in x_poisson]
    X4=np.arange(0,lamb*4,1)
    Y4=poisson.pmf(X4,lamb)
    weights=np.ones_like(x_poisson)/len(x_poisson)
    plt.hist(x_poisson,bins=30,weights=weights)
    plt.plot(X4,Y4,color='r')
    plt.show()

if __name__ == '__main__':
    plot_uniform(0,10)
    plot_norm(0,5)
    plot_beta(2,4)
    plot_poisson(5)



