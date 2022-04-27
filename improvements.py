# idk


# find a way to avoid global variables here 

def initialize_patch_model(n, frac_nomove):
	frac_dispersed = (1-frac_nomove)*(1/(n)) #fraction of fish that disperse to other patches symmetrically
	#transition matrix for dispersal: element [i,j] of kP describes influx of P from j to i
	#doesn't include fish staying in patch (just subtracts the ones that leave) but math should still be consistent 
	global kP, P, C, M, dPs, dCs, dMs
	kP = np.empty((n,n))
	for i in range(n):
		for j in range(n):
			kP[i][j] = frac_dispersed
			if i == j:
				kP[i][j] = -frac_dispersed*(n - 1)

	P_influx = np.empty(n)
	P = np.empty(n) 
	C = np.empty(n) 
	M = np.empty(n)
	dPs = np.empty(n)
	dCs = np.empty(n)
	dMs = np.empty(n)

	if frac_nomove == 0:
		print(kP)

	#concatenate initial condition arrays 
	global X1
	X1 = [P0]*n + [C0L]*n + [M0H]*n
	global X2
	X2 = [P0]*n + [C0H]*n + [M0L]*n

# rass briggs model is only one with different state variables, so just have a separate init function 
def RB_initialize_patch_model(n, frac_nomove):
	frac_dispersed = (1-frac_nomove)*(1/(n)) #fraction of fish that disperse to other patches symmetrically
	#transition matrix for dispersal: element [i,j] of kP describes influx of P from j to i
	#doesn't include fish staying in patch (just subtracts the ones that leave) but math should still be consistent 
	global kP, P, C, Mi, Mv, dPs, dCs, dMis, dMvs
	kP = np.empty((n,n))
	for i in range(n):
		for j in range(n):
			kP[i][j] = frac_dispersed
			if i == j:
				kP[i][j] = -frac_dispersed*(n - 1)

	P_influx = np.empty(n)
	P = np.empty(n) 
	C = np.empty(n) 
	Mi = np.empty(n)
	Mv = np.empty(n)
	dPs = np.empty(n)
	dCs = np.empty(n)
	dMis = np.empty(n)
	dMvs = np.empty(n)

	if frac_nomove == 0:
		print(kP)

	#concatenate initial condition arrays -- amend to match state vec
	global X1
	X1 = [P0]*n + [C0L]*n + [M0H]*n 
	global X2
	X2 = [P0]*n + [C0H]*n + [M0L]*n


def patch_system(model, X, t, closure_length, f, m, n, poaching): 
	
	P_influx = [0]*n

	for i in range(n):
		for j in range(n):
			P_influx[i] += (kP[i][j]) *P[j]  
	

	# for this to work, we need the actual equations inside the for loop
	# but the if statements ought to come before... maybe use boolean flags? 
	# so if RB true, we just use that function in the iterating part 
	# and all that happens in the if statements is initialization of the patch model? 

	if model == 'Rass Briggs':
		P, C, Mi, Mv = X.reshape(4, n) 
		#concatenate into 1D vector to pass to next step
		return np.concatenate((dPs, dCs, dMis, dMvs), axis=0)

	elif model == 'Blackwood Mumby':
		P,C,M = X.reshape(3, n) 
		#concatenate into 1D vector to pass to next step
		return np.concatenate((dPs, dCs, dMs), axis=0)

	elif model == 'van de Leemput fbCP':
		P,C,M = X.reshape(3, n) 
		dPs[i] = P_influx[i]+ s*P[i]*(1 - (P[i] / K(sigma,C[i]))) - fishing(P[i], f)*P[i] *(square_signal(t, closure_length, i, m, n, poaching))
		dCs[i] = (i_C + r*C[i])*(1-M[i]-C[i])*(1-alpha*M[i]) - d*C[i]
		dMs[i] = (i_M+gamma*M[i])*(1-M[i]-C[i])-g*M[i]*P[i]/(g*eta*M[i]+1) 
		#concatenate into 1D vector to pass to next step
		return np.concatenate((dPs, dCs, dMs), axis=0)

	elif model == 'van de Leemput fbMP':
		P,C,M = X.reshape(3, n) 
		dPs[i] = P_influx[i]+ s*P[i]*(1 - (P[i] / K(sigma,C[i]))) - fishing(P[i], f)*P[i] *(square_signal(t, closure_length, i, m, n, poaching))
		dCs[i] = (i_C + r*C[i])*(1-M[i]-C[i])*(1-alpha*M[i]) - d*C[i]
		dMs[i] = (i_M+gamma*M[i])*(1-M[i]-C[i])-g*M[i]*P[i]/(g*eta*M[i]+1) 
		#concatenate into 1D vector to pass to next step
		return np.concatenate((dPs, dCs, dMs), axis=0)

	elif model == 'van de Leemput fbCM':

		# where to set feedback parameters ? here or earlier on ? 
		P,C,M = X.reshape(3, n) 	
		dPs[i] = P_influx[i]+ s*P[i]*(1 - (P[i] / K(sigma,C[i]))) - fishing(P[i], f)*P[i] *(square_signal(t, closure_length, i, m, n, poaching))
		dCs[i] = (i_C + r*C[i])*(1-M[i]-C[i])*(1-alpha*M[i]) - d*C[i]
		dMs[i] = (i_M+gamma*M[i])*(1-M[i]-C[i])-g*M[i]*P[i]/(g*eta*M[i]+1) 
		#concatenate into 1D vector to pass to next step
		return np.concatenate((dPs, dCs, dMs), axis=0)

	elif model == 'van de Leemput': # all feedbacks active 
		P,C,M = X.reshape(3, n) 
		dPs[i] = P_influx[i]+ s*P[i]*(1 - (P[i] / K(sigma,C[i]))) - fishing(P[i], f)*P[i] *(square_signal(t, closure_length, i, m, n, poaching))
		dCs[i] = (i_C + r*C[i])*(1-M[i]-C[i])*(1-alpha*M[i]) - d*C[i]
		dMs[i] = (i_M+gamma*M[i])*(1-M[i]-C[i])-g*M[i]*P[i]/(g*eta*M[i]+1) 

		#concatenate into 1D vector to pass to next step
		return np.concatenate((dPs, dCs, dMs), axis=0)


	else:
		print("Bad input, defaulting to Blackwood-Mumby!")
		P,C,M = X.reshape(3, n) 
		#concatenate into 1D vector to pass to next step
		return np.concatenate((dPs, dCs, dMs), axis=0)

		#set solution_vector to what each of these are returning and then put all of the above inside functions?
# return solution_vector

# to modularize this, we could just pass the iterator and the necessary arrays themselves to a function
# which will just compute the dP dC dM vector and return it 