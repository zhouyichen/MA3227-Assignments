from math import *

class Ivp:
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

class BVP:
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

# Q1:
# q1 = BVP(lambda x:1, lambda x:2, lambda x:cos(x), 0, pi/2, -0.3, -0.1)
# q1_exact_y = lambda x: -(sin(x) + 3 * cos(x))/10
# n = 4
# q1_error = q1.error(q1_exact_y, n)
# print q1_error


class NonelinearBVP:
	def __init__(self, f, fy, fv, a, b, ya, yb):
		self.f = f
		self.fy = fy
		self.fv = fv
		self.a = a
		self.b = b
		self.ya = ya
		self.yb = yb

	def newton(self, n, TOL, M, *t0):
		h = (self.b - self.a)/n
		k = 1
		if t0:
			t = t0[0]
		else:
			t = (self.yb - self.ya) / (self.b - self.a)
		while k <= M:
			print 't', t
			y_list, z = self.modified_euler_method(n, h, t)
			error = y_list[-1] - self.yb

			if abs(error) < TOL:
				print k, t, y_list
				return y_list
			t = t - error / z
			k += 1
		return y_list

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
			sw = fv(x, y, v) * z + fv(x, y, v) * w
			y0 = y + h * sy
			v0 = v + h * sv
			z0 = z + h * sz
			w0 = w + h * sw
			print 'y* =', y0, 'v* =', v0, 'z* =', z0, 'w* =', w0
			v += (sv + f(x1, y0, v0)) * h / 2
			y += (sy + v) * h / 2
			w += (sw + fy(x1, y, v) * z0 + fv(x1, y, v) * w0) * h / 2
			z += (sz + w) * h / 2
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

	def error(self, exact_y, n, TOL, M, *t0):
		approximate = self.newton(n, TOL, M, *t0)
		print approximate
		exact = self.exact_solutions(exact_y, n)
		print exact
		result = [abs(approximate[i] - exact[i]) for i in range(n + 1)]
		return result
'''
def f(x, y, yp):
	return (32 + 2 * x**3 - y * yp) / 8
def fy(x, y, yp):
	return -yp/8
def fv(x, y, yp):
	return -y/8

test = NonelinearBVP(f, fy, fv, 1.0, 3.0, 17.0, 43/3.0)
test.newton(20, 10**(-5), 100)
'''
def f(x, y, yp):
	return -(yp ** 2) - y + log(x)
def fy(x, y, yp):
	return -1
def fv(x, y, yp):
	return -2 * yp

test = NonelinearBVP(f, fy, fv, 1.0, 2.0, 0.0, log(2))
error = test.error(log, 2, 10**(-5), 3, log(2))
print error

