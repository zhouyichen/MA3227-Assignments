import numpy as np
from numpy import linalg as LA

def fixed_point(x0, g, TOL, N):
	x = x0
	for i in range(N):
		x1 = g(x)
		print x1
		print 'residue', abs(x1 - x)
		if abs(x1 - x) < TOL:
			
			return x1
		x = x1

def g(x):
	return np.pi + np.sin(x/2)/2

def dg(x):
	return np.cos(x/2)/4

TOL = 10 ** (-9)

x0 = 0.0
N = 4
# fixed_point(x0, g, TOL, N)

def app_equal(x, y):
	return abs(x - y) < 10 ** (-8)

def f(x):
	return np.pi + np.sin(x/2)/2 - x

def df(x):
	return np.cos(x/2)/4 - 1

def newton(p0, f, df, TOL, N):
	print 'Newton method'
	for i in range(N):
		if app_equal(df(p0), 0):
			print 'Zero derivative'
			return;
		print 'f=', f(p0), 'f\'=', df(p0) 
		p = p0 - f(p0)/df(p0)
		print p
		print 'residue', abs(p - p0)
		if (abs(p - p0) < TOL):
			return p
		p0 = p

# newton(x0, f, df, TOL, N)

def inf_norm(x, y):
	error = max(abs(x[i] - y[i]) for i in range(len(x)))
	print 'error', error
	return error

def fixed_point_system(x0, G, TOL, N):
	x = x0
	for i in range(N):
		x1 = G(x)
		print i + 1
		print x1
		if inf_norm(x1, x) < TOL:
			print'converged'
			return x1
		x = x1

def G(x):
	a, b = x
	x1 = (a ** 2 + b ** 2 + 8) / 10
	x2 = (a * (b ** 2) + a + 8) / 10
	return [x1, x2]

x0 = [0.0, 0.0]
TOL = 10 ** (-3)

fixed_point_system(x0, G, TOL, 10)
