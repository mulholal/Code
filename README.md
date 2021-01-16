# Coding Assignments
# Aim:
# Numerical differentiation using Python
# backward, forward, central method
# computation of relative errors
print("Alexandra Mulholland 17336557")
import numpy as np
import matplotlib.pyplot as plt


#Section 1
#Defining the function to be differentiated 
def function(x):
    """ define function to be differentiated"""
    f=np.cos(x)
    return f
#Defining the analytic second derivative
def analytic(x):
    """analytic second derivative"""
    f= -1*np.cos(x)
    return f
#Defining the central difference equation 
#Does not suffer from 'cancellation catastrophe'- No subtraction of 2 similar sized numbers
def central(x,h):
    f = (function(x+h)-2*function(x)+function(x-h))/(h*h)
    return f

#Section 2
#absolute error of central difference approximation
def error_abs(x,h):
    f= abs(analytic(x)-central(x,h))
    return f 

#Section 3
#Relative error for the central difference approximation
def error_rel(x,h):
    f=error_abs(x,h)/abs(analytic(x))
    return f

#x ranging from 0 to 4pi in steps of 0.0001
x=np.arange(0,4*np.pi,0.0001)

plt.figure(1)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Plot of f(x)=cos(x)')
plt.plot(x,function(x))
plt.grid(True)
plt.show()

#analytical plot for the second derivative
plt.figure(2)
plt.xlabel('x')
plt.ylabel('second derivative f"(x)')
plt.title('Analytical result for the second derivative')
plt.plot(x,analytic(x))
plt.grid(True)
plt.show()

#fixing h as pi/10 and plotting second derivative using central difference method
plt.figure(3)
plt.xlabel('x')
plt.ylabel('second derivative')
plt.title('Central difference result for second derivative for h=pi/10')
plt.plot(x,central(x,np.pi/10.0))
plt.grid(True)
plt.show()

#precision of central difference approximation for x=pi
print("Precision test:",central(np.pi,np.pi/10))

#Graph of comparison of the analytical and numerical results
plt.figure(4)
plt.xlabel('x')
plt.ylabel('f"(x)')
plt.title('Comparing the Analytical and Central Difference results')
plt.plot(x,analytic(x),label='analytic')
plt.plot(x,central(x,np.pi/10.0),label='Central Difference (h=pi/10)') #size step of h arbitrary
plt.legend()
plt.show()

#Precision of analytical for x=pi
print("Precision:",analytic(np.pi))
#This gives a result of 1 and so proves the accuracy of the analytical result.

#absolute error for the central difference approximation
#where h=pi/10
plt.figure(5)
plt.xlabel('x')
plt.ylabel('absolute error')
plt.title('Plot of absolute error for Central Difference approximation')
plt.plot(x,error_abs(x,0.1*np.pi))
plt.grid(True)
plt.show()

#Can see that max error occurs at npi
#min error occurs and npi/2

h=np.arange(np.pi/50000,np.pi/10000,np.pi/1000000000)
#absolute error for varying h for x=pi, as this was a point of maximum absolute error
#Produce zoom in non log graph 
plt.figure(6)
plt.title('Absolute error of Central difference approximation,varying h for x=$\pi$',y=1.04)
plt.xlabel('step size, h')
plt.ylabel('absolute error')
plt.plot(h,error_abs(np.pi,h))
plt.grid(True)
plt.show()

#Plot of log relative error for log varying step size h, for x=pi
plt.figure(7)
plt.title('Plot of relative error of Central Difference method, varying h for x=$\pi$',y=1.04)
plt.xlabel('step size, h')
plt.ylabel('relative error')
plt.plot(h,error_rel(np.pi,h),color='b')
plt.grid(True)
plt.show()

h=np.arange(np.pi/50000,np.pi/10000,np.pi/1000000000)
#as h gets smaller, 2 terms become more alike- increasingly small numbers
#v small number minus v small number, will lose significant figures
def subtractive(x,h):
    f = ((function(x+h)-function(x))-(function(x)-function(x-h)))/h**2
    return f
#absolute error of central subtractive
def subtracterror_abs(x,h):
    f=abs(analytic(x)-subtractive(x,h))
    return f
#relative error of central subtractive
def subtracterror_rel(x,h):
	f = abs((analytic(x)-subtractive(x,h))/analytic(x))
	return f

#Plot of subtractive cancellation absolute error 
#for x=pi
plt.figure(8)
plt.title('Absolute error vs step size for subtractive central difference method',y=1.04)
plt.xlabel('step size h')
plt.ylabel('absolute error')
plt.plot(h,subtracterror_abs(np.pi,h))
plt.grid(True)
plt.show()

#Plot of subtractive cancellation relative error 
#for x=pi
plt.figure(9)
plt.title('Relative error vs step size for subtractive central difference method',y=1.04)
plt.xlabel('step size h')
plt.ylabel('Relative error')
plt.plot(h,subtracterror_rel(np.pi,h),color='g')
plt.grid(True)
plt.show()

plt.figure(10)
plt.xlabel('h-step size')
plt.ylabel('relative error')
plt.title('Comparison of relative error for original and subtractive central difference methods',y=1.04)#
plt.plot(h,error_rel(np.pi,h), label='original',color='b')
plt.plot(h,subtracterror_rel(np.pi,h), label='subtractive',color='g') #for x=pi and log-log plot
plt.legend()
plt.show()
