# how do we solve the global variable problem ? 




class Model:
	# variable arguments based on model type? or just use string to specify parameters...
	def __init__(self, model_type, n, frac_nomove): 

		if model_type == 'RB':
			RB_initialize_patch_model(n, frac_nomove)
		elif model_type == 'BM':
			initialize_patch_model(n, frac_nomove)

		elif model_type == 'vdL_PC':
			initialize_patch_model(n, frac_nomove)

		elif model_type == 'vdL_MP':
			initialize_patch_model(n, frac_nomove)
		elif model_type == 'vdL_MC':
			initialize_patch_model(n, frac_nomove)
		elif model_type == 'vdL': # all feedbacks active 
			initialize_patch_model(n, frac_nomove)
		else:
			print("Bad input, defaulting to Blackwood-Mumby!")


		# parameters -- based on model type 
		# bundle initialization based on model type as well 

	def initialize_patch_model(n, frac_nomove):
		
		frac_dispersed = (1-frac_nomove)*(1/(n)) # fraction of fish that disperse to other patches symmetrically

		# transition matrix for dispersal: element [i,j] of kP describes influx of P from j to i
		kP = np.empty((n,n))
		for i in range(n):
			for j in range(n):
				kP[i][j] = frac_dispersed
				if i == j:
					kP[i][j] = -frac_dispersed*(n - 1)

		
		global kP, P, C, M, dPs, dCs, dMs


		P_influx = np.empty(n)
		P = np.empty(n) 
		C = np.empty(n) 
		M = np.empty(n)
		dPs = np.empty(n)
		dCs = np.empty(n)
		dMs = np.empty(n)

		#concatenate initial condition arrays 
		global X1
		X1 = [P0]*n + [C0L]*n + [M0H]*n
		global X2
		X2 = [P0]*n + [C0H]*n + [M0L]*n

	# rass briggs model is only one with different state variables, so just have a separate init function 
	def RB_initialize_patch_model(n, frac_nomove):

		

		global kP, P, C, Mi, Mv, dPs, dCs, dMis, dMvs

		frac_dispersed = (1-frac_nomove)*(1/(n)) # fraction of fish that disperse to other patches symmetrically
		# transition matrix for dispersal: element [i,j] of kP describes influx of P from j to i
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

	# this defines the system conditionally based on user input -- this is the function we pass to odeint method
	def patch_system(model_type):
		P_influx = [0]*n

		for i in range(n):
			for j in range(n):
				P_influx[i] += (kP[i][j]) *P[j]  
	
			if model_type == 'RB':
				return rass_briggs(X, i)
			elif model_type == 'BM':
				return blackwood(X, i)

			elif model_type == 'vdL_PC':
				return leemput(X, i)

			elif model_type == 'vdL_MP':
			# set global parameters 
				return leemput(X, i)
				
			elif model_type == 'vdL_MC':
				return leemput(X, i)
				
			elif model_type == 'vdL': # all feedbacks active 
				return leemput(X, i)
				
			else:
				print("Bad input, defaulting to Blackwood-Mumby!")
				return blackwood(X, i)

	# returns the model run for a certain set of parameters 
	def run_model(IC_set, t, closure_length, f, m, n, poaching):
		sol = odeint(patch_system, IC_set, t, args = (closure_length, f/(1-m/n), m, n, poaching))
		return sol 

	# make setting for show versus save 
	def time_series():
		if (coral_high):
			IC_set = X2
		else:
			IC_set = X1
		plt.figure()
		sol = odeint(self.patch_system(model_type), IC_set, t, args = (closure_length, f/(1-m/n), m, n, poaching))
		print(sol[:, n:2*n])
		patch_num = [x for x in range(1,n+1)]
		for i in range(len(patch_num)):
			plt.plot(sol[:, n+i],  label = 'patch % d'% patch_num[i])
		plt.xlabel('time')
		plt.ylabel('abundance')
		plt.title('spatial model')
		plt.legend(loc=0)
		txt="parameters" + "\nfishing when open: " + str(f/(1-m/n)) + "\npercent closure: " + str(m/n) +"\nclosure time: " + str(closure_length)
		plt.figtext(0.3, .31, txt, wrap=True, fontsize=8)
		plt.show()

	# make a flag for fast or slow version 
	def coral_recovery_map():
		patches = 20

	#heatmap of period vs m -- a bit messy right now 
	coral_array_HI =  np.zeros(shape=(patches,patches))
	coral_array_LO =  np.zeros(shape=(patches,patches))
	coral_array_AVG =  np.zeros(shape=(patches,patches))
	period_array = np.empty(patches)
	m_array = np.empty(patches)
	fishin = 0.35
	n = 20

	frac_nomove = 1
	for period in range(1,patches+1):
		print(period)
		print(patches)
		for m in range(patches):

			initialize_patch_model(n, frac_nomove)
			print(n, m, sep = ' ')
			displacer = 1/(1-m/float(n))
			final_coral_LO = odeint(patch_system, X1, t, args = (period, displacer*float(fishin), m, n, 0))# full_output = 1)
			#graph_sol(period*4, displacer*float(fishin), m, n, 1, False, 0)
			#graph_sol(5, 0.45, 4, 10, 1, False, 0) -- INCONSISTENT
			#avg1 = 0
			avg2 = 0
			for year in range(999- (999 % (n*period)) - 2*(n*period), 999 - (999 % (n*period))):
				#avg1 += final_coral_HI[year][n]
				avg2 += final_coral_LO[year][n]
			#avg1 = avg1 / (2*(period) + 1)
			avg2 = avg2 / (2*(period*n) + 1)
			print(avg2)
			print("________________________")
			#coral_array_HI[period-1][m] = avg1
			coral_array_LO[period-1][m] = avg2
			#coral_array_AVG[period-1][m] = 0.5 * (avg1 + avg2)
			period_array[period-1] = period
			m_array[m] = m
			show_labels = False

	plt.title('heatmap', fontsize = 20)
	sb.heatmap(coral_array_LO, vmin = 0.0, vmax = 1.0, cmap = "viridis", annot = show_labels) #YlGnBu for original color scheme
	plt.ylabel('period', fontsize = 10) # x-axis label with fontsize 15
	plt.xlabel('number of closures', fontsize = 10) # y-axis label with fontsize 15
	plt.yticks(rotation=0)
	plt.show()

	def bistable_zone():
		return None
	def find_unstable_equilibrium():
		return None
	def scenario_plot():
		return None 

	
	def K(sigma, C):
		return (1-sigma)+sigma*C

	def square_signal(t, closure_length, region, m, n, poaching):
		start = int((t % (n*closure_length))/closure_length)
		if start+m-1 >= n:
			end = (start + m - 1)%n

		else:
			end = (start + m - 1)
		if region >= start and region <= end:
			return poaching
		elif start + m - 1 >= n and (region >= start or region <= end):
			return poaching
		else:
			#determine whether to bake displacement into signal
			return (1-(m/n)*poaching)#/(1-(m/n))
		
	#this signal function is not quite working yet 
	def sigmoid_signal(t, period, p):
		if period == 0:
			return 0
		else:
			return 1.0 / (1 + math.exp(-(t % period - p * period)))

	def fishing(parrotfish, f):
		return f/(1+math.exp(-steepness*(parrotfish-shift)))




	def leemput(X, i):
		# check input 
		P,C,M = X.reshape(3, n) 
		dPs[i] = P_influx[i]+ s*P[i]*(1 - (P[i] / K(sigma,C[i]))) - fishing(P[i], f)*P[i] *(square_signal(t, closure_length, i, m, n, poaching))
		dCs[i] = (i_C + r*C[i])*(1-M[i]-C[i])*(1-alpha*M[i]) - d*C[i]
		
		dMs[i] = (i_M+gamma*M[i])*(1-M[i]-C[i])-g*M[i]*P[i]/(g*eta*M[i]+1)

		#concatenate into 1D vector to pass to next step
		return np.concatenate((dPs, dCs, dMs), axis=0)

	def rass_briggs(X, i):
		# check input
		return None 


	'''
	first idea 

	def leemput(P, C, M, which_var):
		if which_var == 'P':
			# return parrotfish equation
		elif which_var == 'C':
			# return coral equation 
		elif which_var == 'M':
			# return algae equation 
		else:
			print("error")
			exit(2)

	def rassbriggs(P, C, M, which_var):
		if which_var == 'P':
			# return parrotfish equation
		elif which_var == 'C':
			# return coral equation 
		elif which_var == 'Mv':
			# return vuln algae equation 
		elif which_var == 'Mi':
			# return invuln algae equation 
		else:
			print("error")
			exit(2)

	def blackwood(P, C, M, which_var):
		if which_var == 'P':
			# return parrotfish equation
		elif which_var == 'C':
			# return coral equation 
		elif which_var == 'M':
			# return algae equation 
		else:
			print("error")
			exit(2)
	'''