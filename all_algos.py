import numpy as np
from numpy import linalg as LA
from math import *

def iterate(method, mat, b, x0,  w, number):
	T, c = method(mat, b)
	x = x0
	for i in range(number):
		x = T.dot(x) + c
		print x
	return x

def radius(mat):
	values, vectors = LA.eig(mat)
	print values
	magnitude = np.absolute(values)
	print (magnitude)
	return max(magnitude)

def app_equal(x, y):
	return abs(x - y) < 10 ** (-8)

def inf_norm_diff(x, y):
	error = max(abs(x[i] - y[i]) for i in range(len(x)))
	print 'error', error
	return error

def jacobi(mat, b, x0, number):
	def get_Tc(mat, b):
		upper = -np.triu(mat , 1)
		lower = -np.tril(mat , -1)
		d_inv = LA.inv(np.diag(np.diag(mat)))

		T = d_inv.dot(lower + upper)
		c = d_inv.dot(b)
		print 'T'
		print T
		print 'c: ', c
		return (T, c)
	return iterate(get_Tc, mat, b, x0, 1, number)

def gauss_seidel(mat, b, x0, number):
	def get_Tc(mat, b):
		upper = -np.triu(mat , 1)
		lower = -np.tril(mat , -1)
		inv = LA.inv(np.diag(np.diag(mat)) - lower)

		T = inv.dot(upper)
		c = inv.dot(b)
		print 'T'
		print T
		print 'c: ', c
		return (T, c)
	return iterate(get_Tc, mat, b, x0, 1, number)

def SOR_gauss_seidel(mat, b, x0, w, number):
	def get_Tc(mat, b, w):
		upper = -np.triu(mat , 1)
		lower = -np.tril(mat , -1)
		diag_mat = np.diag(np.diag(mat))
		inv = LA.inv(diag_mat - w * lower)

		T = inv.dot((1 - w) * diag_mat + w * upper)
		print inv
		c = w * inv.dot(b)
		print 'T'
		print T
		print 'c: ', c
		return (T, c)
	return  iterate(get_Tc, mat, b, x0, w, number)

def conjugate_gradient(A, b, x0, number):
	r = b - A.dot(x0)
	v = r
	x = x0
	for i in range(number):
		print 'previous r:'
		print r
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

def fixed_point(x0, g, TOL, N):
	x = x0
	for i in range(N):
		x1 = g(x)
		print x1
		print 'residue', abs(x1 - x)
		if abs(x1 - x) < TOL:
			
			return x1
		x = x1

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

def fixed_point_system(x0, G, TOL, N):
	x = x0
	for i in range(N):
		x1 = G(x)
		print i + 1
		print x1
		if inf_norm_diff(x1, x) < TOL:
			print'converged'
			return x1
		x = x1

class LinearIVP:
	def __init__(self, p, q, r, a, b, ya, yap):
		self.p = p
		self.q = q
		self.r = r
		self.a = a
		self.b = b
		self.ya = ya
		self.yap = yap

	def modified_euler_method(self, n):
		h = (self.b - self.a)/n
		u = self.yap
		y = self.ya
		result = [y]
		for i in range(n):
			x = self.a + i * h
			x1 = x + h
			k_y1 = h * u
			print 'k_y1', k_y1
			k_u1 = h * (self.p(x) * u + self.q(x) * y + self.r(x))
			print 'k_u1', k_u1
			k_y2 = h * (u + k_u1)
			print 'k_y2', k_y2
			k_u2 = h * (self.p(x1) * (u + k_u1) + self.q(x1) * (y + k_y1) + self.r(x1))
			print 'k_u2', k_u2
			y += (k_y1 + k_y2)/2.0
			print 'y', y 
			result.append(y)
			u += (k_u1 + k_u2)/2.0
			print 'u', u
		return result

class NonlinearIVP:
	def __init__(self, f, a, b, ya, yap):
		self.f = f
		self.a = a
		self.b = b
		self.ya = ya
		self.yap = yap

	def euler_method(self, n, h):
		x = self.a
		y = self.ya
		v = self.yap
		y_list = [y]
		for i in range(n):
			sy = v
			sv = self.f(x, y, v)
			y = y + h * sy
			v = v + h * sv
			x = x + h
			print 'x =', x, 'y =', y, 'v =', v
			y_list.append(y)
		return y_list

class LinearBVPShooting:
	def __init__(self, p, q, r, a, b, ya, yb):
		self.p = p
		self.q = q
		self.r = r
		self.a = a
		self.b = b
		self.ya = ya
		self.yb = yb

	def ivp_1(self):
		print 'IVP 1:______________________________'
		return Ivp(self.p, self.q, self.r, self.a, self.b, self.ya, 0)

	def ivp_2(self):
		print 'IVP 2:______________________________'
		return Ivp(self.p, self.q, lambda x:0, self.a, self.b, 0, 1)

	def linear_shooting(self, n):
		y1 = self.ivp_1().modified_euler_method(n)
		y2 = self.ivp_2().modified_euler_method(n)
		beta = self.yb
		y1_b = y1[-1]
		y2_b = y2[-1]
		result = [(y1[i] + (beta - y1_b) * y2[i] / y2_b)for i in range(n+1)]
		return result

	def exact_solutions(self, exact_y, n):
		result = [self.ya]
		h = (self.b - self.a)/n
		x = self.a
		for i in range(n-1):
			x += h 
			result.append(exact_y(x))
		result.append(self.yb)
		return result

	def error(self, exact_y, n):
		approximate = self.linear_shooting(n)
		print approximate
		exact = self.exact_solutions(exact_y, n)
		print exact
		result = [abs(approximate[i] - exact[i]) for i in range(n + 1)]
		return result

class NonlinearBVPShooting:
	def __init__(self, f, fy, fv, a, b, ya, yb):
		self.f = f
		self.fy = fy
		self.fv = fv
		self.a = a
		self.b = b
		self.ya = ya
		self.yb = yb

	def newton(self, method, n, TOL, M, *t0):
		h = (self.b - self.a)/n
		k = 1
		if t0:
			t = t0[0]
		else:
			t = (self.yb - self.ya) / (self.b - self.a)
		while k <= M:
			print 't', t
			y_list, z = method(self, n, h, t)
			error = y_list[-1] - self.yb

			if abs(error) < TOL:
				print k, t, y_list
				return y_list
			t = t - error / z
			k += 1
		return y_list

	def euler_method(self, n, h, t):
		x = self.a
		y = self.ya
		v = t
		z = 0
		w = 1
		y_list = [y]
		for i in range(n):
			x1 = x + h
			sy = v
			sv = f(x, y, v)
			sz = w
			sw = fy(x, y, v) * z + fv(x, y, v) * w
			y = y + h * sy
			v = v + h * sv
			z = z + h * sz
			w = w + h * sw
			print 'y =', y, 'v =', v, 'z =', z, 'w =', w
			x = x1
			y_list.append(y)
		return [y_list, z]

	def modified_euler_method(self, n, h, t):
		x = self.a
		y = self.ya
		v = t
		z = 0
		w = 1
		y_list = [y]
		for i in range(n):
			x1 = x + h
			sy = v
			sv = f(x, y, v)
			sz = w
			sw = fy(x, y, v) * z + fv(x, y, v) * w
			y0 = y + h * sy
			v0 = v + h * sv
			z0 = z + h * sz
			w0 = w + h * sw
			print 'y* =', y0, 'v* =', v0, 'z* =', z0, 'w* =', w0
			v += (sv + f(x1, y0, v0)) * h / 2
			y += (sy + v0) * h / 2
			w += (sw + fy(x1, y0, v0) * z0 + fv(x1, y0, v0) * w0) * h / 2
			z += (sz + w0) * h / 2
			print 'y =', y, 'v =', v, 'z =', z, 'w =', w
			x = x1
			y_list.append(y)
		return [y_list, z]

	def exact_solutions(self, exact_y, n):
		result = [self.ya]
		h = (self.b - self.a)/n
		x = self.a
		for i in range(n-1):
			x += h 
			result.append(exact_y(x))
		result.append(self.yb)
		return result

	def error(self, exact_y, method, n, TOL, M, *t0):
		approximate = self.newton(method, n, TOL, M, *t0)
		print approximate
		exact = self.exact_solutions(exact_y, n)
		print exact
		result = [abs(approximate[i] - exact[i]) for i in range(n + 1)]
		return result


def linear_BVP_finite_difference(p, q, r, a, b, ya, yb, n):
	h = (b - a) / (n + 1.0)
	A = np.zeros([n, n])
	B = np.zeros([n, 1])
	for i in range(n):
		x = a + h * (i + 1)
		B[i] = - h * h * r(x)
		A[i,i] = 2 + h * h * q(x)
		aim1 = 1 + h * p(x) / 2
		a1p1 = 1 - h * p(x) / 2
		if (i > 0):
			A[i,i-1] = - aim1
		else:
			B[i] += aim1 * ya
		if (i < n-1):
			A[i,i+1] = - a1p1
		else:
			B[i] += a1p1 * yb
	print 'A:'
	print A
	print 'b:'
	print B.transpose()
	y = LA.solve(A, B)
	print 'y:'
	print y.transpose()
	return y

def nonlinear_BVP_finite_difference(f, fy, fv, a, b, ya, yb, y0, n, TOL, M):
	h = (b - a) / (n + 1.0)
	J = np.zeros([n, n])
	d = np.zeros([n, 1])
	for num in range(M):
		for i in range(n):
			x = a + h * (i + 1)
			yi = y0[i + 1]
			yim1 = y0[i]
			yip1 = y0[i + 2]
			t = (yip1 - yim1) / 2 * h
			d[i] = - (2 * yi - yim1 - yip1 + h * h * f(x, yi, t))
			J[i,i] = 2 + h * h * fy(x, yi, t)
			if (i > 0):
				J[i,i-1] = - (1 + h * fv(x, yi, t) / 2)
			if (i < n - 1):
				J[i,i+1] = - (1 - h * fv(x, yi, t) / 2)
		print 'J:'
		print J
		print 'd:'
		print d.transpose()
		v = LA.solve(J, d)
		print 'v:'
		print v.transpose()
		for i in range(1, n + 1):
			y0[i] = y0[i] + v[i - 1][0]
		if abs(max(v)) < TOL:
			print 'Converged'
			print 'y =', y0
			break
	return y0
