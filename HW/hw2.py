import numpy as np

def SOR_gauss_seidel(mat, b, w):
	upper = -np.triu(mat , 1)
	lower = -np.tril(mat , -1)
	diag_mat = np.diag(np.diag(mat))
	inv = np.linalg.inv(diag_mat - w * lower)

	T = inv.dot((1 - w) * diag_mat + w * upper)
	print inv
	c = w * inv.dot(b)
	print 'T'
	print T
	print 'c: ', c
	return (T, c)


def SOR(method, mat, b, x0,  w, number):
	T, c = method(mat, b, w)
	x = x0
	for i in range(number):
		x = T.dot(x) + c
		print x
	return x


def conjugate_gradient(A, b, x0, number):
	r = b - A.dot(x0)
	v = r
	x = x0
	for i in range(number):
		
		previous_r_norm_square = np.transpose(r).dot(r)[0][0]
		print 'r k-1: ', previous_r_norm_square
		A_v = A.dot(v)
		print 'Av:'
		print A_v
		print 'v dot Av:'
		print np.transpose(v).dot(A_v)[0][0]
		t = previous_r_norm_square / np.transpose(v).dot(A_v)[0][0]

		print 't:' , t
		x = x + t * v
		print "x"
		print x
		r = r - t * A_v
		print "r"
		print r
		s = np.transpose(r).dot(r)[0][0] / previous_r_norm_square
		print 'r', np.transpose(r).dot(r)[0][0]
		print 's', s
		v = r + s * v
		print "v"
		print v
	return x


a1 = np.array([[10.0, -1.0, 0.0], [-1.0, 10.0, -2.0], [0.0, -2.0, 10.0]]);
b1 = np.array([[9], [7], [6]])
x01 = np.zeros([3,1])
w = 1.1

# SOR(SOR_gauss_seidel, a1, b1, x0, w, 2)

# a1 = np.array([[4.0, 3.0, 0.0], [3.0, 4.0, -1.0], [0.0, -1.0, 4.0]]);
# b1 = np.array([24, 30, -24])
# x0 = np.ones(a1.shape[0])
# w = 1.25

# SOR(SOR_gauss_seidel, a1, b1, x0, w, 7)


a2 = np.array([
	[10, 5, 0, 0],
	[5, 10, -4, 0],
	[0, -4, 8, -1],
	[0, 0, -1, 5]
	]);

b2 = np.array([[6], [25], [-11], [-11]])

x02 = np.zeros([4,1])
w = 1.1

# SOR(SOR_gauss_seidel, a2, b2, x0, w, 2)

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


A = np.array([
	[0.2, 0.1, 1.0, 1.0, 0.0],
	[0.1, 4.0, -1.0, 1.0, -1.0],
	[1.0, -1.0, 60.0, 0.0 ,-2.0],
	[1.0, 1.0, 0.0, 8.0, 4.0],
	[0.0, -1.0, -2.0, 4.0, 700.0]
])

b = np.array([[1.0],[2.0],[3.0],[4.0],[5.0]])

x0 = np.zeros([3, 1])
# conjugate_gradient(a1, b1, x01, 2)
conjugate_gradient(a2, b2, x02, 2)
