import numpy as np
from numpy import linalg as LA

n = 4
A = np.zeros([n, n])
b = np.zeros([n, 1])
h = 1.0 / (n + 1)
exact_y = np.zeros([n, 1])

def q1(i):
	return 1/15.0

def q2(i):
	return 2/15.0

def q3(i):
	return 2/15.0

def q4(i):
	return (3 * i**2 - 3 * i + 1) / 15.0

def q5(i):
	return (-2 * i**2 + 4 * i /3.0 - 1/3.0) / 125.0

def q6(i):
	return (-2 * i**2 - 4 * i /3.0 - 1/3.0) / 125.0

def ai(i):
	return q4(i) + q4(i+1) + q2(i) + q3(i)
	# print 2 * (i**2 + 1)/5.0
	# return 4 * (3 * i **2 +2)/15.0

def aip1(i):
	return -q4(i+1) + q1(i)

def aim1(i):
	return -q4(i) + q1(i-1)

def bi(i):
	return q5(i) + q6(i)

def y_exact(x):
	return x**2 - x

# def phi(i):
# 	def helper(x):

# 	return helper


for i in range(n):
	A[i,i] = ai(i+1)
	if (i > 0):
		A[i,i-1] = aim1(i+1)
	if (i < n-1):
		A[i,i+1] = aip1(i+1)
	b[i] = bi(i+1)
	exact_y[i] = y_exact(h * (i+1))
print A
print b
c = LA.solve(A, b)

print c
print exact_y

print exact_y-c

