def p_list(n,k):
	from math import factorial
	a = {} ; dmtr = 0.;
	for i in range(k+1):
		nmtr = combination(n,k-i)*(k-i)**i
		dmtr += combination(n,k-i)*(k-i)**i
		a[i] = nmtr
	for item in a.items():
		a[item[0]] = item[1] / dmtr
	return(a.items())

def poisson(u,k):
	from math import exp, factorial
	p = exp(-u)*(u**k)/factorial(k)
	return(p)

def binomial(n,k,p):
	p = combination(n,k)*(p**k)*((1-p)**(n-k))
	return(p)

def combination(n,k):
	from math import factorial
	c = factorial(n)/(factorial(k)*factorial(n-k))
	return(c)

def poisson_distr(u, K):
	a = []
	for k in range(K+1):
		p = poisson(u,k)
		a.append(p)
	return(a)





