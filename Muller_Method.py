import numpy as np
import math
import cmath


def evalPoly(n, a, t):
	res = 0.0
	
	temp_a = a[n]
	temp_b = temp_a
	temp_c = temp_b

	temp = n - 1

	while temp >= 1:
		temp_b = a[temp] + temp_b * t
		temp_c = temp_b + temp_c * t
		temp -= 1

	temp_b = a[temp] + temp_b * t

	return temp_b, temp_c

def muller(init_value, a):
	x1 = init_value[0]
	x2 = init_value[1]
	x3 = init_value[2]
	
	n = len(a) - 1

	y1, yy = evalPoly(n, a, x1)
	y2, yy = evalPoly(n, a, x2)
	y3, yy = evalPoly(n, a, x3)

	error_limit = 1e-6

	c1 = (y2 - y1)/(x2 - x1)
	c2 = (y3 - y2)/(x3 - x2)
	d1 = (c2 - c1)/(x3 - x1)

	s = c2 + d1*(x3 - x2)

	# print(s + np.sign(s)*(s**2 - 4*y3*d1))

	x4 = x3 - 2.0*y3/cmath.sqrt(s + np.sign(s)*(s**2 - 4.0*y3*d1))

	temp_value, yy = evalPoly(n, a, x4)

	error = np.absolute(temp_value)

	while error > error_limit:

		x1 = x2
		x2 = x3
		x3 = x4 

		y1, yy = evalPoly(n, a, x1)
		y2, yy = evalPoly(n, a, x2)
		y3, yy = evalPoly(n, a, x3)

		c1 = (y2 - y1)/(x2 - x1)
		c2 = (y3 - y2)/(x3 - x2)
		d1 = (c2 - c1)/(x3 - x1)

		s = c2 + d1*(x3 - x2)

		x4 = x3 - 2.0*y3/(s + np.sign(s)*(s**2 - 4.0*y3*d1)**(0.5))

		temp_value, yy = evalPoly(n, a, x4)

		error = np.absolute(temp_value)

	return x4




def linearRoot(a):

	return -a[0]/a[1]

def quardraticRoots(a):
	res = []

	res.append((-a[1] + (a[1]**2.0 - 4.0*a[0]*a[2]))/(2.0*a[2]))
	res.append((-a[1] - (a[1]**2.0 - 4.0*a[0]*a[2]))/(2.0*a[2]))

	return res


# def cubicRoots(a):

# 	res = []

# 	p = a[2]/a[3]
# 	q = a[1]/a[3]
# 	r = a[0]/a[3]

# 	a = (3.0*q - p**2)/3.0
# 	b = (2.0*p**3 - 9.0*p*q + 27.0*r)/27.0

# 	A = (-b/2 + (b**2/4.0 + a**3/27.0)**2)**(1/3)
# 	B = (-b/2.0 - (b**2/4.0 + a**3/27.0)**2)**(1/3)

# 	y1 = A + B
# 	y2 = -0.5*(A + B) + (3**(0.5)*(A - B)/2.0)*1.0j
# 	y3 = -0.5*(A + B) - (3**(0.5)*(A - B)/2.0)*1.0j

# 	res.append(y1 + p/3.0)
# 	res.append(y2 + p/3.0)
# 	res.append(y3 + p/3.0)

# 	return res

# def quarticRoots(a):
# 	res = []

# 	p = a[3]/a[4]
# 	q = a[2]/a[4]
# 	r = a[1]/a[4]
# 	s = a[0]/a[4]

# 	a_p = q - 3.0*p**2/8
# 	b_p = r + p**3.0/8.0 - p*q/2.0
# 	c_p = s - 3.0*p**4.0/256.0 + p**2*q/16.0 - p*r/4.0

# 	cubic_para = [4.0*q*s - r**2 - p**2*s, p*r - 4.0*s, -q, 1.0]

# 	cubic_roots = cubicRoots(cubic_para)

# 	z = 0.0

# 	for i in range(3):
# 		if math.fabs(cubic_roots[i].imag) < 1e-6:
# 			z = cubic_roots[i].real

# 	R = cmath.sqrt(p**2/4.0 - q + np.real(z))


# 	D = 0.0 + 0.0j

# 	E = 0.0 + 0.0j

# 	if np.absolute(R) < 1e-6:
# 		D = cmath.sqrt(3.0*p**2/4.0 - 2.0*q + 2.0*(z**2 - 4.0*s))
# 		E = cmath.sqrt(3.0*p**2/4.0 - 2.0*q - 2.0*(z**2 - 4.0*s))
# 	else:
# 		D = cmath.sqrt(3.0*p**2/4.0 - R**2 - 2.0*q + (4.0*p*q - 8.0*r - p**3)/(4.0*R))
# 		E = cmath.sqrt(3.0*p**2/4.0 - R**2 - 2.0*q - (4.0*p*q - 8.0*r - p**3)/(4.0*R))

# 	res.append(-p/4.0 + (R + D)/2.0)
# 	res.append(-p/4.0 + (R - D)/2.0)
# 	res.append(-p/4.0 - (R - E)/2.0)
# 	res.append(-p/4.0 - (R + E)/2.0)

# 	# print(R, D, E)

# 	return res


def acosTheta(a, b):

	if b > 0:

		return cmath.acos(-cmath.sqrt((b**2 / 4.0)/(-a**3/27.0)))

	else:
		return cmath.acos(cmath.sqrt((b**2 / 4.0)/(-a**3/27.0)))

def cubicRoots(a):

	res = []

	p = a[2]/a[3]
	q = a[1]/a[3]
	r = a[0]/a[3]

	# print(p, q, r)

	a_p = (3.0*q - p**2.0)/3.0
	b_p = (2.0*p**3.0 - 9.0*p*q + 27.0*r)/27.0

	# print(a_p, b_p)

	check = b_p**2.0/4.0 + a_p**3.0/27.0
	# print(check)

	A = 0.0 + 0.0j

	B = 0.0 + 0.0j

	y1 = 0.0 + 0.0j
	y2 = 0.0 + 0.0j
	y3 = 0.0 + 0.0j

	if check > 0:

		temp_A = -b_p/2.0 + math.sqrt(check)
		temp_B =  -b_p/2.0 - math.sqrt(check)

		if temp_A > 0:
			A = temp_A**(1/3.0)
		else:
			A = -(-temp_A)**(1/3.0)

		if temp_B > 0:
			B = temp_B**(1/3.0)
		else:
			B = -(-temp_B)**(1/3.0)

		y1 = A + B
		y2 = -0.5*(A + B) + (3**(0.5)*(A - B)/2.0)*1.0j
		y3 = -0.5*(A + B) - (3**(0.5)*(A - B)/2.0)*1.0j

	elif check == 0:

		if b > 0:
			y1 = -2 * cmath.sqrt(-a_p/3.0)
			y2 = cmath.sqrt(-a_p/3.0)
			y3 = cmath.sqrt(-a_p/3.0)
		elif b < 0:

			y1 = 2 * cmath.sqrt(-a_p/3.0)
			y2 = -cmath.sqrt(-a_p/3.0)
			y3 = -cmath.sqrt(-a_p/3.0)

		else:
			y1 = 0
			y2 = 0
			y3 = 0
	else:

		y1 = 2*cmath.sqrt(-a_p/3.0) * cmath.cos(acosTheta(a_p,b_p)/3.0 + 2*cmath.pi/3.0)
		y2 = 2*cmath.sqrt(-a_p/3.0) * cmath.cos(acosTheta(a_p,b_p)/3.0 + 2*2*cmath.pi/3.0)
		y3 = 2*cmath.sqrt(-a_p/3.0) * cmath.cos(acosTheta(a_p,b_p)/3.0 + 2*3*cmath.pi/3.0)

	res.append(y1 - p/3.0)
	res.append(y2 - p/3.0)
	res.append(y3 - p/3.0)

	return res

def quarticRoots(a):
	res = []

	p = a[3]/a[4]
	q = a[2]/a[4]
	r = a[1]/a[4]
	s = a[0]/a[4]

	# print(p, q, r, s)

	a_p = q - 3.0*p**2.0/8.0
	b_p = r + p**3.0/8.0 - p*q/2.0
	c_p = s - 3.0*p**4.0/256.0 + p**2.0*q/16.0 - p*r/4.0

	cubic_para = [4.0*q*s - r**2.0 - p**2.0*s, p*r - 4.0*s, -q, 1.0]

	# print('test')
	# print(cubic_para)
	cubic_roots = cubicRoots(cubic_para)
	# print(cubic_roots)

	z = 0.0

	for i in range(3):
		if math.fabs(cubic_roots[i].imag) < 1e-6:
			z = cubic_roots[i].real
			break


	# print(p**2/4 - q + np.real(z))
	R = cmath.sqrt(p**2.0/4.0 - q + z)

	# print(np.absolute(R))

	D = 0.0 + 0.0j

	E = 0.0 + 0.0j

	if np.absolute(R) < 1e-6:
		D = cmath.sqrt(3.0*p**2.0/4.0 - 2.0*q + 2.0*cmath.sqrt(z**2.0 - 4.0*s))
		E = cmath.sqrt(3.0*p**2.0/4.0 - 2.0*q - 2.0*cmath.sqrt(z**2.0 - 4.0*s))
	else:
		D = cmath.sqrt(3.0*p**2.0/4.0 - R**2.0 - 2.0*q + (4.0*p*q - 8.0*r - p**3.0)/(4.0*R))
		E = cmath.sqrt(3.0*p**2.0/4.0 - R**2.0 - 2.0*q - (4.0*p*q - 8.0*r - p**3.0)/(4.0*R))

	res.append(-p/4.0 + (R + D)/2.0)
	res.append(-p/4.0 + (R - D)/2.0)
	res.append(-p/4.0 - (R - E)/2.0)
	res.append(-p/4.0 - (R + E)/2.0)

	# print(R, D, E)

	return res


def deflation(a, r):

	l = len(a) - 1

	res = [0.0 + 0.0j] * l

	res[l-1] = a[l]

	for i in range(1,l):
		res[l - i - 1] = res[l - i]*r + a[l - i]

	return res

def deflation_two(a, r):

	l = len(a) - 2

	res = [0.0 + 0.0j] * l

	res[l - 1] = a[l - 1 + 2]
	res[l - 2] = a[l - 2 + 2] - r[1]*res[l - 2 + 1]

	for i  in range(3, l+1):
		res[l - i] = a[l - i + 2] - r[1]*res[l - i + 1] - r[0]*res[l - i + 2]

	return res


# test = [0, 1, 2, 2, 1]

# test_r = [1,1,1]

# print(deflation_two(test, test_r))


def newton(a, f):

	f = 0.0

	th = 1e-6

	count = 0 

	n = len(a) - 1

	# print "Iteration: ", count 
	# print "Estimated x is ",f
	# print "f(x) is ", evalFunc(f)

	while math.fabs(evalPoly(n,a,f)[0]) > th:

		f = f - evalPoly(n,a,f)[0]/evalPoly(n,a,f)[1]

		count += 1

		# print "Iteration: ", count 
		# print "Estimated x is ",f
		# print "f(x) is ", evalFunc(f)


	# f = a - evalFunc(a)/evalFuncDe(a)
	# print "Iteration: ", count + 1 
	# print "Estimated x is ", f
	# print "f(x) is ", evalFunc(f)


	return f


def mullerRoots(a, n):

	if n == 1:
		return linearRoot(a)

	if n == 2:
		return quardraticRoots(a)

	if n == 3:
		return cubicRoots(a)

	if n == 4:
		return quarticRoots(a)


	l = 0
	p = a
	deg = n

	init_v = [0.0, 1.0, 2.0]

	temp_root = 0.0 + 0.0j

	res = []

	while deg >= 5:
		temp_root = muller(init_v, a)
		# print(temp_root)
		res.append(temp_root)

		if math.fabs(np.imag(temp_root)) < 1e-6:

			a = deflation(a, temp_root)
			a = [np.real(temp) for temp in a]
			l += 1
			deg -= 1

			while np.absolute(evalPoly(deg, a, temp_root)[0]) < 1e-6:
				print(temp_root)
				res.append(temp_root)
				a = deflation(a, temp_root)
				l += 1
				deg -= 1
		else:

			re = np.real(temp_root)
			im = np.imag(temp_root)

			a = deflation_two(a, [re**2 + im**2, -2*re, 1])
			a = [np.real(temp) for temp in a]
			l += 2
			deg -= 2

			res.append(re - im*1.0j)
			while np.absolute(evalPoly(deg, a, temp_root)) < 1e-6:
				res.append(temp_root)
				res.append(re - im*1.0j)

				a = deflation_two(a, [re**2 + im**2, -2*re, 1])

				l += 2
				deg -= 2

	# print(deg)
	# print(a)

	# print(res)
	
	if deg == 1:
		r0 = linearRoot(a)
		res.append(r0)

	if deg == 2:
		r0, r1 = quardraticRoots(a)

		res.append(r0)
		res.append(r1)

	if deg == 3:
		r0, r1, r2 = cubicRoots(a)
		res.append(r0)
		res.append(r1)
		res.append(r2)

	if deg == 4:
		r0, r1, r2, r3 = quarticRoots(a)
		res.append(r0)
		res.append(r1)
		res.append(r2)
		res.append(r3)

	# for i in range(len(res)):
	# 	res[i] = newton(p, res[i])

	return res


# # p = [-6.8, 10.8, -10.8, 7.4, -3.7, 1]

init_v = [0, 1, 2]
# # res = muller(init_v, p)

# # print(res)
# # print(math.fabs(evalPoly(5, p, res)))

temp_p = [-6.8, 10.8, -10.8, 7.4, -3.7, 1]

temp_q = [0.00787276, -0.180591, -0.360995, 9.15636, -25.7634, 14.6196, 10.1887, -8.35979, -0.843121, 1]


temp_test = [0.31972438038439854, -7.256533885173501, -2.320987501231513, 2.7340233785189527, 1]


# print(quarticRoots(temp_test))

print(mullerRoots(temp_p, 5))

print(mullerRoots(temp_q, 9))

# print(mullerRoots(temp_q, 9))
# print(evalPoly(4, temp_test, -2.7684766402232786))
# res = deflation(temp_p, 1.7)
