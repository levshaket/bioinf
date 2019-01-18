#Assigns random expression values based on low ('l'), medium ('m'), or high ('h') settings. 
#Binary: [0.,1.];   Discrete:[0.,0.25,0.5,0.75,1.];  Uniform: uniform distribution
def express(mode,setting):
	import random
#	random.seed()
	if mode == 'uniform':
		if setting == 'l':
			return(random.uniform(0,0.5))
		if setting == 'm':
			return(random.uniform(0.25,0.75))
		if setting == 'h':
			return(random.uniform(0.5,1))
		if setting == 'o':
			return(random.random())
	if mode == 'discrete':
		if setting == 'l':
			return(random.choice([0.,0.25,0.5]))
		if setting == 'm':
			return(random.choice([0.25,0.5,0.75]))
		if setting == 'h': 
			return(random.choice([0.5,0.75,1.]))
		if setting == 'o':
			return(random.choice([0.,0.25,0.5,0.75,1.]))
	if mode == 'binary':
		if setting == 'l':
			return(0.)
		if setting == 'm':
			return(random.choice([0.,1.]))
		if setting == 'h': 
			return(1.)
		if setting == 'o':
			return(random.choice([0.,1.]))
