import numpy as np

def jacobi(mat, b):
	upper = -np.triu(mat , 1)
	lower = -np.tril(mat , -1)
	d_inv = np.linalg.inv(np.diag(np.diag(mat)))

	T = d_inv.dot(lower + upper)
	c = d_inv.dot(b)
	print 'T'
	print T
	print 'c: ', c
	return (T, c)

def gauss_seidel(mat, b):
	upper = -np.triu(mat , 1)
	lower = -np.tril(mat , -1)
	inv = np.linalg.inv(np.diag(np.diag(mat)) - lower)

	T = inv.dot(upper)
	c = inv.dot(b)
	print 'T'
	print T
	print 'c: ', c
	return (T, c)

def iterate(method, mat, b, number):
	T, c = method(mat, b)
	x = np.zeros(T.shape[0])
	for i in range(number):
		x = T.dot(x) + c
		print x
	return x

a1 = np.array([[-2, 1, 1/2.0], [1, -2, -1/2.0], [0, 1, 2]]);
b1 = np.array([4, -4, 0])

a2 = np.array([
	[10, 5, 0, 0],
	[5, 10, -4, 0],
	[0, -4, 8, -1],
	[0, 0, -1, 5]
	]);

b2 = np.array([6, 25, -11, -11])

# iterate(gauss_seidel, a2, b2, 2)
# print np.linalg.solve(a1, b1)

def radius(mat):
	values, vectors = np.linalg.eig(mat)
	print values
	magnitude = np.absolute(values)
	print (magnitude)
	return max(magnitude)

a3 = np.array([
	[2, -1, 1],
	[2, 2, 2],
	[-1, -1, 2]
	])
b3 = np.array([-1, 4, -5])

# print radius(gauss_seidel(a3, b3)[0])