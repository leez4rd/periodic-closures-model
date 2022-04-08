# binary search for unstable equilibrium 

# for a given set of parameters,
# we iteratively run the model with different initial conditions until finding the unstable equilibrium

# start at coral high and coral low, bisect
# does trajectory from middle go up or down (what is the sign of dC/dt?)

# if it goes up, choose upper interval, if it goes down, choose lower one 
# new bounds of interval are C_middle and C_low or C_high
# we repeat this process recursively ...

def find_unstable_eq(C_HI, C_LO):
	C_M = .5*(C_HI + C_LO)
	#solve ode for C_M = initial coral over 1 timestep... or just plug C_M in to dCdt formula
	coral_derivative = # whatever we get
	if coral_derivative > 0:
		find_unstable_eq(C_M, C_HI)
	elif coral_derivative < 0:
		find_unstable_eq(C_M, C_LO)
	else:
		return C_M 


