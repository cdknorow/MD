import numpy as np
# Integrate using simpsons rule
# int_a^b f(x) dx ~ (b-a)/6 * (f(a) + 4 * f((a+b)/2) + f(b))
def simpson_integration(f,x):
    A = [0]
    m = [0]
    for i in range(0,x.shape[0]-2,2):
        A.append((x[i+2]-x[i])/6. * (f[i]+4*f[i+1]+f[i+2])+A[-1])
        m.append(x[i+2])
    return A, m


def f(x):
	return x**2

x = np.array([0,1,2,3,4,5,6,8])
y1 = f(x)
I1 = simpson_integration(y1,x)
print I1