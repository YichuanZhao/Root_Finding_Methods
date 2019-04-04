import numpy as np

import math

from decimal import Decimal

def evalFunc(x):

	return math.exp(-x) - math.sin(2*x)


def evalFuncDe(x):

	return -math.exp(-x) - 2*math.cos(2*x)

def bisection():

	th = 5e-5

	l = 0.0
	r = 1.0
	res = 0.0

	count = 0

	res = (l + r)/2

	print "Iteration: ", count 
	print "Estimated bracket [a, b] is [" + str(l) + ", " + str(r) + "]"
	print "f(a) is", evalFunc(l)
	print "f(b) is", evalFunc(r)
	print "(a+b)/2 is", res
	print "f((a+b)/2 is ", evalFunc(res)

	while math.fabs(evalFunc(res)) > th:

		if evalFunc(l) * evalFunc(res) <= 0:
			l = l
			r = res
		else:
			l = res
			r = r

		res = (l+r)/2

		count += 1

		print "Iteration: ", count 
		print "Estimated bracket [a, b] is [" + str(l) + ", " + str(r) + "]"
		print "f(a) is", evalFunc(l)
		print "f(b) is", evalFunc(r)
		print "(a+b)/2 is", res
		print "f((a+b)/2) is ", evalFunc(res)


	# res = (l+r)/2
	# print "Iteration: ", count + 1 
	# print "Estimated bracket [a, b] is [" + str(l) + ", " + str(r) + "]"
	# print "f(a) is", evalFunc(l)
	# print "f(b) is", evalFunc(r)
	# print "(a+b)/2 is", res
	# print "f((a+b)/2 is ", evalFunc(res)

	return res


def modifiedRegFalsi():

	th = 5e-5

	l = 0.0
	r = 1.0
	res_p = 0.0

	count = 0

	f_a = evalFunc(l)
	f_b = evalFunc(r)

	res = (f_b*l - f_a*r)/(f_b - f_a)

	print "Iteration", count
	print "Estimated bracket [a, b] is [" + str(l) + ", " + str(r) + "]"
	print "F is", f_a
	print "G is", f_b
	print "w is ",res
	print "f(w) is ", evalFunc(res)


	while math.fabs(evalFunc(res)) > th:

		if evalFunc(l) * evalFunc(res) <= 0:
			l = l
			r = res
			f_b = evalFunc(r)

			if evalFunc(res_p) * evalFunc(res) > 0:
				f_a = f_a/2
		else:
			l = res
			r = r
			f_a = evalFunc(l)

			if evalFunc(res_p) * evalFunc(res) > 0:
				f_b = f_b/2

		
		res_p = res
		res = (f_b*l - f_a*r)/(f_b - f_a)
		
		count += 1

		print "Iteration: ", count 
		print "Estimated bracket [a, b] is [" + str(l) + ", " + str(r) + "]"
		print "F is", f_a
		print "G is", f_b
		print "w is ",res
		print "f(w) is ", evalFunc(res)

	
	# res = (f_b*l - f_a*r)/(f_b - f_a)
	# print "Estimated bracket [a, b] is [" + str(l) + ", " + str(r) + "]"
	# print "F is", f_a
	# print "G is", f_b
	# print "w is ",res
	# print "f(w) is ", evalFunc(res)


	return res

def secant():

	
	th = 5e-5

	a = 0.0

	count = 0

	print "Iteration: ", count 
	print "Estimated x is ", a
	print "f(x) is ", evalFunc(a)

	count += 1

	b = 1.0 
	print "Iteration: ", count 
	print "Estimated x is ",b
	print "f(x) is ", evalFunc(b)

	f = b

	while math.fabs(evalFunc(f)) > th:

		f = b - evalFunc(b)*(b - a)/(evalFunc(b) - evalFunc(a))

		a = b

		b = f
		
		count += 1

		print "Iteration: ", count 
		print "Estimated x is ",f
		print "f(x) is ", evalFunc(f)

	
	# f = f = b - evalFunc(b)*(b - a)/(evalFunc(b) - evalFunc(a))
	# print "Iteration: ", count + 1 
	# print "Estimated x is ", f
	# print "f(x) is ", evalFunc(f)


	return f


def newton():

	f = 0.0

	th = 5e-5

	count = 0 

	print "Iteration: ", count 
	print "Estimated x is ",f
	print "f(x) is ", evalFunc(f)

	while math.fabs(evalFunc(f)) > th:

		f = f - evalFunc(f)/evalFuncDe(f)

		count += 1

		print "Iteration: ", count 
		print "Estimated x is ",f
		print "f(x) is ", evalFunc(f)


	# f = a - evalFunc(a)/evalFuncDe(a)
	# print "Iteration: ", count + 1 
	# print "Estimated x is ", f
	# print "f(x) is ", evalFunc(f)


	return f

print("Bisection:")
bisection()
print("\n")
print("Modified Regula falsi:")
modifiedRegFalsi()
print("\n")

print("Secant:")
secant()
print("\n")
print("Newton:")
newton()