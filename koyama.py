#!/usr/bin/python
# This script calculates the intramolecular distribution using the Koyama distribution.
# You must specify the number of monomers per chain, the bond length, and the stiffness.
# Created by J. McCarty 
from math import *
from numpy import *
# specify parameters
Nch=100 # number of monomers per chain
l=1.54 # bond length
q= 0.785 # stiffness parameter
p=(3*(q**2.0)-1.0)/2.0
print p

# other parameters 
DR=0.015 # r-spacing
NM=2**14  # number of points
DK=pi/(NM*DR) # k-spacing
# initialize arrays
kval = arange(NM)
kval = (kval+1)*DK
omega = ones(NM)
sum1 = zeros(NM)
print 'starting summation'

for ia in range(Nch):
    i=ia+1
    for ja in range(Nch):
        j=ja+1
        n=abs(i-j)
        d1=n**2*((1+q)/(1-q))**2.0
        d2=-n*(1+((2*q)/((1-q)**3.0))*(6+5*q+3*q**2.0)-((1+q)/(1-q))**2.0*((4*p)/(1-p)))
        d3=((2*q)/((1-q)**4.0))*(4+11*q+12*q**2.0)-((4*p)/(1-p))*(1+(8*q)/((1-q)**3.0)+(((1+q)/(1-q))**2.0)*(p/(1-p)))
	nb=float(n)
        d4=-1.0*(q**nb)*((8*q)/((1-q)**3.0))*(n*(1+3*q)+(1+2*q+3*q**2.0)/(1-q)-((2*p)/((q-p)**2.0))*(n*(1-q)*(q-p)+2*q**2.0-p*q-p))
        d5=-1.0*(6*(q)**(2.0*nb+2))/((1-q)**4.0)
        d6=(p**nb)*((4.0/(1-p))*(1+(8*q)/((1-q)**3.0)-(((1+q)/(1-q))**2.0)*(1-p/(1-p)))-((16*(q**2.0))/(((1-q)**3.0)*(q-p)**2.0))*(q+q**2.0-2*p))
        d=d1+d2+d3+d4+d5+d6
        Dij=(2.0/3.0)*d
        if n==0:
            rsqij=0.0
	    omega=ones(NM)
	else:       
            rsqij=l**2.0*n*((1+q)/(1-q)-(2*q*(1-q**nb))/(n*(1-q)**2.0))
            rfourthij=rsqij**2.0+Dij*(l**4.0)
	    CC=((1.0/2.0)*(5-3*(rfourthij/(rsqij)**2.0))**(1.0/2.0))        
            BB=(CC*rsqij)**(1.0/2.0)
	    AA=((rsqij*(1-CC)/6)**(1.0/2.0))
            omega = (sin(BB*kval)/(BB*kval))*exp(-1.0*AA**2.0*kval**2.0)
	sum1=sum1+omega 
    print 'step',i,'of',Nch
    
avomeg=sum1/(float(Nch))
dataout = column_stack((kval,avomeg))
nm=str(Nch)
name='w'+nm+'koyama.dat'
savetxt(name,dataout) 
 
print 'done'



