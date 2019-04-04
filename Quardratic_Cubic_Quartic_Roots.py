import numpy as np
import math
import cmath

def linearRoot(a):

	return -a[0]/a[1]

def quardraticRoots(a):
	res = []

	res.append((-a[1] + (a[1]**2 - 4*a[0]*a[2]))/(2*a[2]))
	res.append((-a[1] - (a[1]**2 - 4*a[0]*a[2]))/(2*a[2]))

	return res

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

a = [4.0, 87.0, -23.0, 110.0]

b = [-3400.0, 0.0, -7.0, 1.34, 43.0]

test = [1,1,-2,1]
test1 = [-1,0,-1,1]


# print(cubicRoots(a))

# print(quarticRoots(b))


print(cubicRoots(test1))

