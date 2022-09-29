
import numpy as np
from scipy.integrate import odeint  
import matplotlib.pyplot as plt 
import math
import seaborn as sb
import multiprocessing as mp 
from joblib import Parallel, delayed
from matplotlib.colors import ListedColormap


class Model:
	
	# variable arguments based on model type? 
	# need a function to tinker with individual parameters, and a function to list them 
	def __init__(self, model_type, n, frac_nomove, mgmt_strat = 'periodic'): 
		self.model_type = model_type
		self.n = n 
		self.frac_nomove = frac_nomove
		self.mgmt_strat = mgmt_strat 

		self.f = 0
		self.closure_length = 0
		self.m = 0 
		self.poaching = 0


		'''
		#time -- maybe move this down to run model for dynamic stuff, same with ICs
		self.yrs = yrs #total amount of time
		self.t = np.linspace(0, yrs, yrs) #timestep array -- same number of timesteps as years 

		#initial conditions
		self.ICs = ICs

		initialize_patch_model()
		'''

	def initialize_patch_model(self, P0, C0L, C0H, M0L, M0H, M0iL = None, M0iH = None):
		# do variables defined in this function only exist in the scope of this function if it is called by constructor?
		

		if self.model_type == 'RB':
			self.RB_initialize_patch_model(P0, C0L, C0H, M0L, M0H, M0iL, M0iH)
			return 

		frac_dispersed = (1-self.frac_nomove)*(1/(self.n)) # fraction of fish that disperse to other patches symmetrically
		print("FRAC NO MOVe: ", self.frac_nomove)
		# transition matrix for dispersal: element [i,j] of kP describes influx of P from j to i
		kP = np.empty((self.n,self.n)) 
		for i in range(self.n):
			for j in range(self.n):
				kP[i][j] = frac_dispersed
				if i == j:
					kP[i][j] = -frac_dispersed*(self.n - 1)
		print(kP)
		setattr(self,'kP', kP)

		# P_influx = np.empty(n) -- I think this is unnecessary 
		setattr(self, 'P', np.zeros(self.n))
		setattr(self, 'C' , np.empty(self.n))
		setattr(self,'M', np.empty(self.n))
		setattr(self,'dPs', np.empty(self.n))
		setattr(self,'dCs', np.empty(self.n))
		setattr(self,'dMs', np.empty(self.n))
		# setattr(self, 'state_vec', [self.P,self.C,self.M])
		# concatenate initial condition arrays 
		# need to define baseline values somewhere
		# should these be attributes ?
		setattr(self, 'X1', [P0]*self.n + [C0L]*self.n + [M0H]*self.n)
		setattr(self, 'X2', [P0]*self.n + [C0H]*self.n + [M0L]*self.n)  
		
	# rass briggs model is only one with different state variables, so just have a separate init function 
	# better way might be to give initialize... an optional argument 
	# or dynamic number of params 
	def RB_initialize_patch_model(self, P0, C0L, C0H, M0vL, M0vH, M0iL, M0iH):

		frac_dispersed = (1-self.frac_nomove)*(1/(self.n)) # fraction of fish that disperse to other patches symmetrically
		# transition matrix for dispersal: element [i,j] of kP describes influx of P from j to i
		kP = np.empty((self.n,self.n))
		for i in range(self.n):
			for j in range(self.n):
				kP[i][j] = frac_dispersed
				if i == j:
					kP[i][j] = -frac_dispersed*(self.n - 1)
		setattr(self,'kP', kP)

		# P_influx = np.empty(n) -- I think this is unnecessary 
		setattr(self, 'P', np.empty(self.n))
		setattr(self, 'C' , np.empty(self.n))
		setattr(self,'Mi', np.empty(self.n))
		setattr(self,'Mv', np.empty(self.n))
		setattr(self,'dPs', np.empty(self.n))
		setattr(self,'dCs', np.empty(self.n))
		setattr(self,'dMis', np.empty(self.n))
		setattr(self,'dMvs', np.empty(self.n))
		# setattr(self, 'state_vec', [P,C,Mv,Mi])
		# instead of global variables, can we pass these directly to patch_system? 
		# concatenate initial condition arrays -- amend to match state vec
		# X1 = [P0]*n + [C0L]*n + [M0H]*n 
		# X2 = [P0]*n + [C0H]*n + [M0L]*n
		setattr(self, 'X1', [P0]*self.n + [C0L]*self.n + [M0vH]*self.n + [M0iH]*self.n)
		setattr(self, 'X2', [P0]*self.n + [C0H]*self.n + [M0vL]*self.n + [M0iL]*self.n)


	# returns the model run for a certain set of parameters 
	# maybe create a subclass for a model run to separate plotting jobs ? 
	def run_model(self, IC_set, t):
		# can we not pass self as an argument here ? 
		sol = odeint(patch_system, IC_set, t, args = (self, ))
		return sol 

	# make setting for show versus save (make these optional arguments) 
	# maybe add second way to use it by plugging in array from run_model 
	def time_series(self, IC_set, t, save, show):

		fig, ax = plt.subplots()
		sol = odeint(patch_system, IC_set, t, args = (self, ))
		patch_num = [x for x in range(1, self.n+1)]
		recov_time = self.get_coral_recovery_time(t)
		if recov_time == -1:
			print("ERROR")
			recov_time = 1


		'''
		sol_scaled = sol
		for j in range(len(sol)):
			sol_scaled[j, :] = sol[int(j / recov_time), :]
		'''
		scaling_array = np.asarray([1 / recov_time]*len(t))

		time_scaled = np.multiply(scaling_array, t)
		for i in range(self.n):

			ax.plot(time_scaled, sol[:, self.n+i],  label = 'patch % d'% (int(i) + 1))
		
		# locs, labels = plt.xticks()

		# print("locs:", locs)

		# newlabels = [int(int(loc) / int(recov_time)) for loc in locs[0:]]

		# print("newlb", newlabels)
		# ax.set_xticks([tick for tick in range(0, 6)])
		# ax.set_xticklabels(newlabels)
		ax.set_xlabel('Time (scaled to coral recovery time)')
		
		ax.set_ylabel('Coral cover (total fraction of area)')
		ax.set_title('Coral cover in each patch over time')
		plt.legend(loc=0)
		ax.set_xlim([0, len(t)/recov_time])
		ax.set_ylim([0, 1])
		# txt="parameters" + "\nfishing when open: " + str(self.f/(1-self.m/self.n)) + "\npercent closure: " + str(self.m/self.n) +"\nclosure time: " + str(self.closure_length)
		# plt.figtext(0.7, .31, txt, wrap=True, fontsize=8)
		if save:
			plt.savefig('TS' + str(self.model_type) + str(self.f) + '_' + str(self.m) +'out_of_'+ str(self.n) + 'dispersal=' + str(self.frac_nomove)+ '.jpg')
		if show:
			plt.show()
			plt.close()


	def coral_recovery_map(self, t, fishing_intensity, IC_set = None, filename = None):

		fig, (ax1,ax2) = plt.subplots(nrows=2, sharex=True, figsize=(10,10), gridspec_kw={'height_ratios': [1, 1.2*self.n]})


		P0, C0L, C0H, M0L, M0H, M0vH, M0vL, M0iH, M0iL = 0.1, 0.04, 0.4, 0.04, 0.4, 0.04, 0.4, 0.04, 0.4


		# slow version  
		# IC_set = self.X2
		MAX_TIME = len(t) # last year in the model run 

		coral_array =  np.zeros(shape=(self.n, self.n)) # array of long-term coral averages
		# CL_array = np.empty(int(0.75*self.n+1)) # array of closure lenghts 
		# m_array = np.empty(int(self.n / 2))  # array of number of closures 
		closure_lengths = np.empty(self.n)
		ms = np.empty(self.n)
		
		def fill_heatmap_cell(self, t, fishing_intensity, IC_set, closure_length, m):
			period = self.n*closure_length
			# set management parameters for this run 
			self.set_mgmt_params(closure_length, fishing_intensity, m, self.poaching)

			# solve ODE system 
			sol = odeint(patch_system, IC_set, t, args = (self, ))
			# average over coral cover of last two full rotations for a single patch (assumes symmetry, may fix that)
			avg = 0

			# dealing with very large period / millennium ratio ?
			if period / MAX_TIME < 0.5:
				for year in range(MAX_TIME - MAX_TIME % period - period, 
					MAX_TIME - MAX_TIME % period):
					avg += sol[year][self.n]
				avg = avg / ((period) + 1)
			else:
				avg = sol[MAX_TIME - 1][self.n]

			return avg


		solutions_pool = mp.Pool() # defaults to num_cpu


		final_coral_covers = Parallel(n_jobs = 4)(delayed(fill_heatmap_cell)(self, t, fishing_intensity, IC_set, closure_length, m) for m in range(self.n) for closure_length in range(1, self.n + 1))
		
		print(final_coral_covers)
		print("FINAL CORAL COVERS ABOVE")
		coral_array = (np.asarray(final_coral_covers)).reshape((self.n, self.n))

		recov_time = self.get_coral_recovery_time(t)
		if recov_time < 0:
			print("oops")
			print(recov_time)
			quit()
		ax2.set_title('Coral cover under periodic closures', fontsize = 15)
		f = lambda y: str(float(self.n*y / recov_time))[0:3]
		new_labels = [f(y) for y in range(1, self.n+1)]

		g = lambda x: str(x / self.n)[0:4]
		
		new_x_labels = [g(x) for x in range(len(ms))]
		# new_x_labels = [0.1, 0.2, 0.303, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
		# print(coral_array)
		ax2 = sb.heatmap(coral_array, vmin = 0.0, vmax = 0.8, cmap = "viridis", xticklabels = new_x_labels, yticklabels = new_labels, cbar = True) #YlGnBu for original color scheme
		ax2.invert_yaxis()
		ax2.set_ylabel('Cycle period (in terms of coral recovery time)', fontsize = 10) # x-axis label with fontsize 15
		ax2.set_xlabel('Fraction of seascape closed', fontsize = 10) # y-axis label with fontsize 15
		# ax.yticks(ax.get_yticks(), ax.get_yticks() * 3)
		# ax.set_yticks(rotation=0)
		# plt.show()45


		# ms, closure_lengths = np.meshgrid(ms, closure_lengths)

		# ax = plt.axes(projection='3d')
		# ax.plot_surface(ms, closure_lengths, coral_array,
		#                 cmap='viridis', edgecolor='none')
		

		""" Plot isoclines """ 
		ps  = np.linspace(0,self.n,100)
		isocline1 = lambda x: 10/(x+0.0000001)
		y  = np.asarray([isocline1(val) for val in ps])
		ax2.plot(ps, y, 'r--')
		isocline2 = lambda x: 35/(x+0.0000001)
		y  = np.asarray([isocline2(val) for val in ps])
		ax2.plot(ps, y, 'r--')

		isocline3 = lambda x: 55/(x+0.0000001)
		y  = np.asarray([isocline3(val) for val in ps])
		ax2.plot(ps, y, 'r--')

		isocline4 = lambda x: 89/(x+0.0000001)
		y  = np.asarray([isocline4(val) for val in ps])
		ax2.plot(ps, y, 'r--')

		isocline5 = lambda x: 125/(x+0.0000001)
		y  = np.asarray([isocline5(val) for val in ps])
		ax2.plot(ps, y, 'r--')


		# ax2.plot([self.n / 2 + 1]*100, y, 'r-')
		
		""" MPA color bar for comparison -- really hacky right now bc it involves creating a Model object within current one """

		mako = ListedColormap(sb.color_palette('viridis').as_hex())

		
		""" Make MPA colorbar for comparisons """
		# hacky but works 
		# need to simulate MPA for second colorbar to juxtapose with heatmap...
		# creating another model object within this one seems unnecessarily funky 
		
		z = Model(self.model_type, self.n, 1, mgmt_strat = 'MPA')
		mpa_corals = np.empty(self.n)
		# set initial conditions 
		z.initialize_patch_model(P0, C0L, C0H, M0L, M0H, M0iL, M0iH)
		z.load_parameters() # do this inside initializer
		for i in range(self.n):
			z.set_mgmt_params(500, fishing_intensity, i, 0)
			MPAsol = z.run_model(IC_set, t)
			print(MPAsol)
			extent = [0, self.n, 0, 1]
			coral_f = 0
			for j in range(self.n):
				coral_f += MPAsol[len(t) - 1][j+self.n]
			coral_f = coral_f / self.n # patch average
			mpa_corals[i] = coral_f
			print(mpa_corals)
		
		

		ax1.imshow(mpa_corals[np.newaxis][:], cmap=mako, aspect=.10, extent=extent, label='MPA coral cover', vmin=0, vmax=1)
		ps  = np.linspace(0,self.n,10)
		
		# ax1.plot(ps, mpa_corals)
		ax1.set_title('Coral cover under MPA')
			# ax1.set_yticks([])
		ax1.set_xlim(extent[0], extent[1])

		position = ax1.get_position()
		position = ax2.get_position()
		print(position.x0)
		# ax1.set_position([0.125, 0.17, 0.75, 1.1])
		# ax1.set_position([position.x0, position.y0 + position.y1, position.x1, 1.1])


		if filename == None:
			plt.show()
		else:
			plt.savefig(str(filename) + '.jpg')
			plt.close()


		# name = 'heatmap_viridis_coral_vs_fishing_' + str(self.model_type) + str(fishing_intensity) + '_' + str(self.n) + '.jpg'
		# plt.savefig(name)
		# plt.close()

	def bistable_zone(self, t, filename = None):
		""" plot final coral cover for different values of fishing effort for two sets of initial conditions """ 
		final_coral_high = np.empty(20)
		final_coral_low = np.empty(20)


		fishing_range = np.linspace(0, 1.5, 20)

		for i, f in enumerate(fishing_range):

			# set management parameters 
			self.set_mgmt_params(0, f, 0, self.poaching)

			# make high start solution
			high_sol = odeint(patch_system, self.X2, t, args = (self, ))

			# make low start solution 
			low_sol = odeint(patch_system, self.X1, t, args = (self, ))

			# note: this only works without periodic oscillations, which this plot assumes are not present 
			yrs = len(t)
			final_coral_high[int(i)] = high_sol[yrs - 1][self.n]

			final_coral_low[int(i)] = low_sol[yrs - 1][self.n]

		plt.plot(fishing_range, final_coral_low, label = 'coral starts low', color = 'blue')
		plt.plot(fishing_range, final_coral_high, label = 'coral starts high' , color = 'green')
		plt.xlabel('fishing effort')
		plt.ylabel('long-term coral cover')

		if filename == None:
			plt.show()
		else:
			plt.savefig(str(filename) + '.jpg')
			plt.close()

	def coral_vs_CL(self, t, fishing_intensity, IC_set = None):

		durations = [i for i in range(1,100)]

		for m in range(0, self.n):
			coral_covers = []
			for duration in durations:

				self.set_mgmt_params(duration, fishing_intensity, m, self.poaching)

				sol = odeint(patch_system, IC_set, t, args = (self, ))

				avg = 0
				# average over last cycle for one patch
				for year in range(len(t) - duration*self.n, len(t)):
					avg += sol[year][self.n]
				avg = avg / (duration*self.n)

				# avg = avg / self.n

				coral_covers.append(avg)

			plt.title('coral cover versus closure duration')
			plt.xlabel('closure duration in years')
			plt.ylabel('coral cover as fraction of seascape')
			plt.ylim([0,1])

			plt.plot(durations, coral_covers, label = 'coral when {} patches are closed'.format(m))
		plt.legend(loc = 0)
		plt.show()



	def find_unstable_equilibrium(self, t, lowC = 0.1, highC = 0.7, recursion_depth = 0):
		""" binary search for unstable equilibrium """ 


		midpoint = (lowC + highC) / 2
		print(midpoint)


		IC_mid = [0.1]*self.n + [midpoint]*self.n + [0.04]*self.n # verify this 

		# run model starting from midpoint -- could reduce runtime by making t smaller 
		mid_sol = odeint(patch_system, IC_mid, t, args = (self, ))

		if recursion_depth > 10:
			print("Close enough....")
			return midpoint

		# if coral cover grows from the midpoint, the equilibrium is above it
		if mid_sol[len(t) - 1][1] - midpoint > 0:
			print("going up...")
			new_recursion_depth = recursion_depth + 1
			return self.find_unstable_equilibrium(t, lowC = midpoint, recursion_depth = new_recursion_depth)
		# if coral cover declines from the midpoint, the equilibrium is below it 
		elif mid_sol[len(t) - 1][1] - midpoint < 0: 
			print("going down...")
			new_recursion_depth = recursion_depth + 1
			return self.find_unstable_equilibrium(t, highC = midpoint, recursion_depth = new_recursion_depth)
		else: # unstable equilibrium found!
			return midpoint 

	def scenario_plot(self, t, fishing_intensity, IC_set, filename = None):
		P0, C0L, C0H, M0L, M0H, M0vH, M0vL, M0iH, M0iL = 0.1, 0.04, 0.4, 0.04, 0.4, 0.04, 0.4, 0.04, 0.4

		crt = self.get_coral_recovery_time(t)
		if crt == -1:
			print("coral recovery time too high")
			quit()
		print("crt: ", crt)
		final_coral = np.empty(self.n)
		ms = np.empty(self.n)
		multipliers = np.asarray([0.1, 0.5, 1, 2, 4])
		crts = np.asarray([crt]*5)
		periods = np.multiply(crts, multipliers)# [0.1*crt, 0.5*crt, 1*crt, 2*crt, 4*crt]  # parametrize in terms of coral growth time? 
			
		
		# there is a cooler way to do colors than this 
		color_sequence = {periods[0]: '#1f77b4', periods[1]: '#aec7e8', periods[2]: '#ff7f0e', periods[3]:'#ffbb78', periods[4]:'#2ca02c'}

		MAX_TIME = len(t)


		for i, period in enumerate(periods):
			
			for m in range(self.n):


				# set management parameters for this run 
				self.set_mgmt_params(period / self.n, fishing_intensity, m, self.poaching)
				self.set_mgmt_params(period / self.n, fishing_intensity, m, self.poaching)

				# solve ODE system 
				sol = odeint(patch_system, IC_set, t, args = (self, ))
				avg = 0

				'''
				parallel version?
				def get_total(sol, N):
					avg +

				solutions_pool = mp.Pool() # defaults to num_cpu


				final_coral_covers = Parallel(n_jobs = 4)(delayed(get_average)()(
					for j in range(self.n) for year in range(MAX_TIME - MAX_TIME % int(period+1) 
						- int(period+1), MAX_TIME - MAX_TIME % int(period+1))

				print(final_coral_covers)
				'''

				# truncate to last full period, then average over that period 
				'''
				for year in range(MAX_TIME - MAX_TIME % int(period+1) - int(period+1), MAX_TIME - MAX_TIME % int(period+1)):
				# do we really need to do all the careful averaging stuff to get the general shape right? 
					# avg over last thirty years, all patches 
					print(year)
					for j in range(self.n):
						print(period)
						avg += sol[year][self.n + j] # only looking at one patch here
						print("year zero: ", sol[0][self.n + j])
						print("year 10: ", sol[10][self.n + j])
						print("year 50: ", sol[50][self.n + j])
						print("year 100: ", sol[100][self.n + j])
						print("running total of coral covers: ", avg)
					avg = avg / self.n
					print("averaged over patches: ", avg)
				avg = avg / period
				print("averaged over period: ", avg)
				'''

				# for j in range(self.n):
				#	avg += sol[len(t) - 1][self.n + j] 
				# avg = avg / self.n

				# plt.plot(sol)
				# plt.show()

				# crude average 
				print(sol[MAX_TIME - 1])
				avg = sol[MAX_TIME - 1][self.n]

				print("Average ",  m, "is ", avg)
				# plt.plot(t, sol[:, self.n])
				# plt.show()

					
				final_coral[m] = avg
				ms[m] = m

			# plot result for this period
			
			
			plt.plot(ms / self.n, final_coral, label = 'period = %s' % str(multipliers[i]), color = color_sequence[period])
		
		
		# plt.xlim([0, 1])
		plt.ylim([0, 1])
		plt.xlabel('Fraction closed')
		plt.ylabel('Final coral Cover')
		plt.title('Final coral state across closure scenarios - ' + str(self.model_type) + '\nfishing = ' + str(fishing_intensity))
		# plt.show()




		# make MPA line 
		z = Model(self.model_type, self.n, 1, mgmt_strat = 'MPA')
		mpa_corals = np.empty(self.n)
		# set initial conditions 
		z.initialize_patch_model(P0, C0L, C0H, M0L, M0H, M0iL, M0iH)
		z.load_parameters() # do this inside initializer
		
		# loop over frac closed 
		for i in range(self.n): 
			z.set_mgmt_params(500, fishing_intensity, i, 0)
			MPAsol = z.run_model(IC_set, t) 
			
			# loop over patches 
			total = 0 
			for j in range(self.n):
				coral_f = 0
				coral_f += MPAsol[len(t) - 1][j+self.n]
				total += coral_f
			total = total / self.n
			mpa_corals[i] = total

		arr = np.linspace(0, 1, self.n)
		plt.plot(arr, mpa_corals, label = 'MPA', color = 'black')
		plt.xlim([0, 0.66])
		if IC_set == self.X1:
			plt.legend(loc=1)
		else:
			plt.legend(loc=3)
	
		# plt.show()
		
		if filename == None:
			plt.savefig('newest_scenario_plot' + str(self.model_type) + '_' + str(fishing_intensity) + '_' + str(IC_set[self.n]) + '_' + str(self.poaching) + '_' + str(self.frac_nomove)+'.jpg')
		else:
			plt.savefig(str(filename) + '.jpg')

		plt.close()


	def get_coral_recovery_time(self, t):

		P0, C0L, C0H, M0L, M0H, M0vH, M0vL, M0iH, M0iL = 0.1, 0.04, 0.4, 0.04, 0.4, 0.04, 0.4, 0.04, 0.4


		z = Model(self.model_type, self.n, 1, mgmt_strat = 'MPA')
		# set initial conditions 
		z.initialize_patch_model(P0, C0L, C0H, M0L, M0H, M0iL, M0iH)
		z.load_parameters() # do this inside initializer

		z.set_mgmt_params(500, 0, 0, 0)

		MPAsol = z.run_model(z.X1, t)
		coralsol = MPAsol[:, self.n]
		list_of_labels = ['fish']*self.n + ['coral']*self.n + ['algae']*self.n
		# plt.plot(t, coralsol, label = 'coral')


		# plt.show()
		# np.take_along_axis(a, ai, axis=1)
		# plt.legend(loc=0)

		if self.model_type == 'RB':
			high_coral_eq = 0.5
		elif self.model_type == 'BM':
			high_coral_eq = 0.54
		else: 
			high_coral_eq = 0.7

		coral_recovery_time = -1
		for i, state in enumerate(coralsol):
			if state > high_coral_eq:
				coral_recovery_time = i
				break

		return coral_recovery_time + 10


	# modify this to take custom feedback parameters, or maybe custom anything? 
	# this doesn't work for some reason
	def load_parameters(self):
		if self.model_type == 'vdL':
			params = {
			"r": 0.303,
			"i_C" : 0.05,
			"i_M" : 0.05,
			"ext_C" : 0.0001,
			"ext_P" : 0.0001,
			"gamma" : 0.8,
			"d" : 0.1,
			"g" : 1,
			"s" : 1,
			"sigma" : .5, #strength of coral-herbivore feedback
			"eta" : 2, #strength of algae-herbivore feedback
			"alpha" : 0.5, #strength of algae-coral feedback 
			"P0" : 0.1,
			"C_HI" : .4,
			"M_LO" : .04,
			"C_LO" : .04,
			"M_HI" : .4
			}

		elif self.model_type == 'vdL_MC':
			params = {
			"r": 0.303,
			"i_C" : 0.05,
			"i_M" : 0.05,
			"ext_C" : 0.0001,
			"ext_P" : 0.0001,
			"gamma" : 0.8,
			"d" : 0.1,
			"g" : 1,
			"s" : 1,
			"sigma" : 0, #strength of coral-herbivore feedback
			"eta" : 0, #strength of algae-herbivore feedback
			"alpha" : 0.5 #strength of algae-coral feedback 
			}

		elif self.model_type == 'vdL_MP':
			params = {
			"r": 0.303,
			"i_C" : 0.05,
			"i_M" : 0.05,
			"ext_C" : 0.0001,
			"ext_P" : 0.0001,
			"gamma" : 0.8,
			"d" : 0.1,
			"g" : 1,
			"s" : 1,
			"sigma" : 0, #strength of coral-herbivore feedback
			"eta" : 2, #strength of algae-herbivore feedback
			"alpha" : 0 #strength of algae-coral feedback 
			}

		elif self.model_type == 'vdL_PC':
			params = {
			"r": 0.303,
			"i_C" : 0.05,
			"i_M" : 0.05,
			"ext_C" : 0.0001,
			"ext_P" : 0.0001,
			"gamma" : 0.8,
			"d" : 0.1,
			"g" : 1,
			"s" : 1,
			"sigma" : .5, #strength of coral-herbivore feedback
			"eta" : 0, #strength of algae-herbivore feedback
			"alpha" : 0, #strength of algae-coral feedback 
			}

		elif self.model_type == 'BM':
			params = {
			"gamma" : 0.8,
			"beta" : 1,
			"alpha" : 1,
			"s" : 0.49,
			"r" : 1,
			"d" : 0.44,
			"a" : 0.1,
			"i_C" : 0.05,
			"i_M" : 0.05
			}

		elif self.model_type == 'RB':
			params = {
			"phiC" : 0.001, #open recruitment of coral
			"phiM" : 0.0001, #open recruitment of macroalgae

			"rM" : 0.5, #production of vulnerable macroalgae from invulnerable stage"

			#growth rates
			"gTC" : 0.1, #combined rate of growth and local recruitment of corals over free space 
			"gTV" : 0.2, #growth rate of vulnerable macroalgae over free space 
			"gTI" : 0.4, #growth rate of invulnerable macroalgae over free space 
			"rH" : 0.49,

			"gamma" : 0.4, #growth of macroalgae over coral vs free space
			"omega" : 2, #maturation rate of macroalgae from vulnerable to invulnerable class "

			#death rates
			"dC" : 0.05, #death rate of coral 
			"dI" : 0.4, #death rate of invulnerable macroalgae
			"dV" : 0.58, #death rate of vulnerable macroalgae per unit biomass of herbivores "

			"K" : 20, 
			"Graze" : 0.58
			}

		for name, val in params.items():
			# initialize a new variable within class with that parameter's name and value 
			# there is probably a less hacky way to do this, but i aint got time to find it
			# exec(f"{name} = {val}", globals())
			setattr(self, name, val)

	# management parameter setter 
	def set_mgmt_params(self, closure_length, f, m, poaching):
		self.closure_length = closure_length
		self.f = f
		self.m = m
		self.poaching = poaching 



# puttin patch system outside the class so we don't have to pass self as first parameter 
# this way we can pass this function directly to odeint, and just call it from within the Model class 
# don't necessarily need outside class, just can't modify object...
# how could we bundle all of this into a library (ie import periodic_closures_modeltools) ? 


# BUG: only the first patch is being simulated? 
def patch_system(X, t, system_model):
	P_influx = [0]*system_model.n

	for i in range(system_model.n):

		# add influx at end AFTER the parrotfish have been initialized so it isn't populated with random values
		for j in range(system_model.n):
			P_influx[i] += (system_model.kP[i][j]) * X[j]  # X[j] up to n is just the parrotfish pops

		# this could be structured more nicely
		if system_model.model_type == 'RB':
			
			results = rass_briggs(X, t, i, system_model, P_influx)
			system_model.dPs[i] = results[0]
			system_model.dCs[i] = results[1]
			system_model.dMvs[i] = results[2]
			system_model.dMis[i] = results[3]

		elif system_model.model_type == 'BM':
			
			results = blackwood(X, t, i, system_model, P_influx)
			system_model.dPs[i] = results[0]
			system_model.dCs[i] = results[1]
			system_model.dMs[i] = results[2]

		elif system_model.model_type == 'vdL_PC':
			
			results = leemput(X, t, i, system_model, P_influx)
			system_model.dPs[i] = results[0]
			system_model.dCs[i] = results[1]
			system_model.dMs[i] = results[2]

		elif system_model.model_type == 'vdL_MP':
			results = leemput(X, t, i, system_model, P_influx)
			system_model.dPs[i] = results[0]
			system_model.dCs[i] = results[1]
			system_model.dMs[i] = results[2]

		elif system_model.model_type == 'vdL_MC':
			
			results = leemput(X, t, i, system_model, P_influx)
			system_model.dPs[i] = results[0]
			system_model.dCs[i] = results[1]
			system_model.dMs[i] = results[2]

		elif system_model.model_type == 'vdL': # all feedbacks active 

			results = leemput(X, t, i, system_model, P_influx)
			# print(leemput(X, t, i, system_model, P_influx))
			system_model.dPs[i] = results[0]
			system_model.dCs[i] = results[1]
			system_model.dMs[i] = results[2]
			
		else:
			print("Bad input, defaulting to Blackwood-Mumby!")
			results = blackwood(X, t, i, system_model, P_influx)
			system_model.dPs[i] = results[0]
			system_model.dCs[i] = results[1]
			system_model.dMs[i] = results[2]

		
	if system_model.model_type == 'RB':
		return np.concatenate((system_model.dPs, system_model.dCs, system_model.dMvs, system_model.dMis), axis = 0)
	else:
		return np.concatenate((system_model.dPs, system_model.dCs, system_model.dMs), axis = 0)

def rass_briggs(X, t, i, system_model, P_influx):


	P, C, Mv, Mi = X.reshape(4, system_model.n)
	T = 1 - C[i] - Mv[i] - Mi[i]
	
	# P,C,M = X.reshape(4, system_model.n) # will reshaping work since we are passing arrays of length n? 
	# T = 1 - C[i] - Mv[i] - Mi[i]
	# dC = P_influx[i]+ system_model.s*P[i]*(1 - (P[i] / K(system_model.sigma,C[i]))) - fishing(P[i], system_model.f)*P[i] *(square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching))
	dP = P_influx[i]+system_model.rH*P[i]*(1-P[i]/system_model.K) - system_model.f/(1-system_model.m/system_model.n)*P[i] *(square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching, system_model.mgmt_strat))
	# print("SQUARE SIGNAL IN PATCH ", i)
	# print(square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching))
	# print(P_influx)
	# print(square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching))


	dC = (system_model.phiC*T) + system_model.gTC*T*C[i] - system_model.gamma*system_model.gTI*Mi[i]*C[i] - system_model.dC*C[i]

	dMv = system_model.phiM*T + system_model.rM*T*Mi[i] + system_model.gTV*T*Mv[i] - system_model.dV*Mv[i] - P[i]*Mv[i]*system_model.Graze - system_model.omega * Mv[i]
	dMi = system_model.omega*Mv[i] + system_model.gTI*T*Mi[i] + system_model.gamma*system_model.gTI*Mi[i]*C[i] - system_model.dI*Mi[i]
	# print(P, C, Mv, Mi)
	# print(system_model.f)#/(1-system_model.m/system_model.n))# *P[i] *(square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching)))
	return [dP, dC, dMv, dMi]

	
	# dP = P_influx[i] + system_model.rH*P[i]*(1-P[i]/system_model.K) - system_model.f/(1-system_model.m/system_model.n)*P[i] *(square_signal(t, system_model.closure_length, i, sysem_model.m, system_model.n, system_model.poaching))
# 	dC = (system_model.phiC*T) + gTC*T*C - gamma*gTI*Mi*C - dC*C
# 	dMv = phiM*T + rM*T*Mi + gTV*T*Mv - dV*Mv - P*Mv*Graze - omega * Mv
# 	# conceptual question: why is that second term multiplied by M not Mv? 
# 	dMi = omega*Mv + gTI*T*Mi + gamma*gTI*Mi*C - dI*Mi
	# return [dP, dC, dMv, dMi]


	# check input
	# return None 

def K(sigma, C):
		return (1-sigma)+sigma*C
def BMK(C):
	return 1 - 0.5*C

def square_signal(t, closure_length, region, m, n, poaching, mgmt_strat = 'periodic'):

	if mgmt_strat == 'periodic':

		if closure_length != 0: 
			start = int((t % (n*closure_length))/closure_length)
		else:
			start = 0

		if start+m-1 >= n:
			end = (start + m - 1)%n

		else:
			end = (start + m - 1)
		# print("START: ", start)
		# print("END: ", end)
		# print("REGION NUM ", region)


		if region >= start and region <= end:
			return poaching
		elif start + m - 1 >= n and (region >= start or region <= end):
			return poaching
		else:
			# print("We are here.")
			# print("all inequalities are false... why?")
			#determine whether to bake displacement into signal
			return (1-(m/n)*poaching)#/(1-(m/n))
	elif mgmt_strat == 'MPA':
		if m == 0:
			return 1 # if we close nothing, signal does not modify fishing intensity
		if m == n:
			return poaching # if we close everything, only poaching remains  
		if region <= m:
			return poaching  # closed region 
		else: 
			return (1 - (m / n) * poaching) # open region 

	
#this signal function is not quite working yet 
def sigmoid_signal(t, period, p):
	if period == 0:
		return 0
	else:
		return 1.0 / (1 + math.exp(-(t % period - p * period)))

# fishing density dependence 
def fishing(parrotfish, f):
	steepness = 25
	shift = 0.2
	return f/(1+math.exp(-steepness*(parrotfish-shift)))

# is self keyword necessary here ? or can we call these functions regardless?
# is it necessary to send X as a parameter if it is embedded in the class? 
# should all of these be class methods which can change self.X? 
def blackwood(X, t, i, system_model, P_influx):
	
	P, C, M = X.reshape(3, system_model.n)

	# dP = s*P[i]*(1 - (P[i] / (beta*system_model.K(C[i])))) - system_model.fishing(P[i], f)*P[i]*system_model.square_signal(t, closure_length, i, m, n, poaching)
	dP = system_model.s*P[i]*(1 - (P[i] / (system_model.beta*BMK(C[i])))) - system_model.f/(1-system_model.m/system_model.n)*P[i]*square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching, system_model.mgmt_strat)
	dC = system_model.r*(1-M[i]-C[i])*C[i]-system_model.d*C[i] - system_model.a*M[i]*C[i] + 0.01*system_model.i_C*(1-M[i]-C[i])
	# need to define g(P) before this model is used 
	dM = system_model.a*M[i]*C[i] - system_model.alpha*P[i]/system_model.beta*M[i] *(1/(1-C[i])) + system_model.gamma*M[i]*(1-M[i]-C[i]) + 0.01*system_model.i_M*(1-M[i]-C[i])

	'''
	if dM < 0:
		dM = 0.01
	if dP < 0:
		dP = 0.01
	'''
	# return np.concatenate((dPs, dCs, dMs), axis=0)
	return [dP, dC, dM]
	
	'''
	P,C,M = X.reshape(3, system_model.n) # will reshaping work since we are passing arrays of length n? 
	# dC = P_influx[i]+ system_model.s*P[i]*(1 - (P[i] / K(system_model.sigma,C[i]))) - fishing(P[i], system_model.f)*P[i] *(square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching))
	dP = P_influx[i]+ system_model.s*P[i]*(1 - (P[i] / K(system_model.sigma,C[i]))) - system_model.f/(1-system_model.m/system_model.n)*P[i] *(square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching))
	# need separtate K function for blackwood
	# print(square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching))
	dC = (system_model.r*C[i])*(1-M[i]-C[i]) - system_model.a*M[i]*C[i] - system_model.d*C[i] # recruitment needed?
	
	dM = (system_model.gamma*M[i])*(1-M[i]-C[i])-system_model.g(P[i])*(1/(1-C[i]))*M[i]


	return [dP, dC, dM]
	'''

# will need to pass [self.P, self.C, self.M] to this for it to work 
def leemput(X, t, i, system_model, P_influx): # COPY THIS FORMAT FOR OTHER MODELS 
	# check input 
	P,C,M = X.reshape(3, system_model.n) # will reshaping work since we are passing arrays of length n? 
	# dC = P_influx[i]+ system_model.s*P[i]*(1 - (P[i] / K(system_model.sigma,C[i]))) - fishing(P[i], system_model.f)*P[i] *(square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching))

	dP = P_influx[i]+ system_model.s*P[i]*(1 - (P[i] / K(system_model.sigma,C[i]))) - fishing(P[i], system_model.f)/(1-system_model.m/system_model.n)*P[i] *(square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching, system_model.mgmt_strat))
	# print(P_influx)
	# print(square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching))
	dC = (system_model.i_C + system_model.r*C[i])*(1-M[i]-C[i])*(1-system_model.alpha*M[i]) - system_model.d*C[i]
	
	dM = (system_model.i_M+system_model.gamma*M[i])*(1-M[i]-C[i]) - system_model.g*M[i]*P[i]/(system_model.g*system_model.eta*M[i]+1)

	return [dP, dC, dM]
	# NEED TO CHANGE OTHER FUNCTIONS TO MATCH BLEP BLOP BLEEP METHOD
	# concatenate into 1D vector to pass to next step
	# return np.concatenate((dC, dC, dM), axis=0)


def main():

	
	yrs = 1000 #total amount of time
	t = np.linspace(0, yrs, yrs) #timestep array -- same number of timesteps as years 

	RB_yrs = 2000 #total amount of time
	RB_t = np.linspace(0, RB_yrs, RB_yrs) #timestep array -- same number of timesteps as years 


	P0, C0L, C0H, M0L, M0H, M0vL, M0vH, M0iL, M0iH = 0.1, 0.04, 0.4, 0.04, 0.4, 0.04, 0.2, 0.04, 0.2


	# create Model objects

	blackwood_mumby = Model('BM', 10, 1, mgmt_strat = 'periodic') 
	van_de_leemput = Model('vdL', 10, 1, mgmt_strat = 'periodic')
	rass_briggs = Model('RB', 10, 1, mgmt_strat = 'periodic')



	blackwood_mumby.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	van_de_leemput.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	rass_briggs.initialize_patch_model(P0, C0L, C0H, M0vL, M0vH, M0iL, M0iH)

	van_de_leemput.load_parameters()
	rass_briggs.load_parameters()
	blackwood_mumby.load_parameters()


	# rass_briggs.bistable_zone(t)
	# rass_briggs.set_mgmt_params(100, 0.2864, 3, 0)
	# rass_briggs.time_series(rass_briggs.X1, RB_t, save = False, show = True) # why does this also flatten out?


	'''

	blackwood_mumby.set_mgmt_params(25, 0.303, 1, 0)
	blackwood_mumby.time_series(blackwood_mumby.X1, t, save = False, show = True) # why does this flatten out?

	rass_briggs.set_mgmt_params(30, 0.2864, 1, 0)
	rass_briggs.time_series(rass_briggs.X1, t, save = False, show = True) # why does this also flatten out?

	van_de_leemput.set_mgmt_params(35, 0.28, 2, 0)
	van_de_leemput.time_series(van_de_leemput.X1, t, save = False, show = True) 
	'''

	'''

	van_de_leemput.coral_recovery_map(t, 0.28, van_de_leemput.X1, filename = 'vdL_heatmap_0.28_startinglow')
	blackwood_mumby.coral_recovery_map(t, 0.303, blackwood_mumby.X1, filename = 'BM_heatmap_0.303_startinglow')
	rass_briggs.coral_recovery_map(t, 0.2864, rass_briggs.X1, filename = 'RB_heatmap_0.2864_startinglow')
	van_de_leemput.coral_recovery_map(t, 0.28, van_de_leemput.X2, filename = 'vdL_heatmap_0.28_startinghigh')
	blackwood_mumby.coral_recovery_map(t, 0.303, blackwood_mumby.X2, filename = 'BM_heatmap_0.303_startinghigh')
	rass_briggs.coral_recovery_map(t, 0.2864, rass_briggs.X2, filename = 'RB_heatmap_0.2864_startinghigh')

	'''
	
	
	
	# scenario plots 
	ICs = rass_briggs.X1
	rass_briggs.scenario_plot(RB_t, 0.2864, ICs, filename = 'RB_ScenarioPlot_NoDispersal_StartingLow')
	ICs = rass_briggs.X2
	rass_briggs.scenario_plot(RB_t, 0.2864, ICs, filename = 'RB_ScenarioPlot_NoDispersal_StartingHigh')
	ICs = blackwood_mumby.X1
	blackwood_mumby.scenario_plot(t, 0.303, ICs, filename = 'BM_ScenarioPlot_NoDispersal_StartingLow')
	ICs = blackwood_mumby.X2
	blackwood_mumby.scenario_plot(t, 0.303, ICs, filename = 'BM_ScenarioPlot_NoDispersal_StartingHigh')
	ICs = van_de_leemput.X1
	van_de_leemput.scenario_plot(t, 0.28, ICs, filename = 'vdL_ScenarioPlot_NoDispersal_StartingLow')
	ICs = van_de_leemput.X2
	van_de_leemput.scenario_plot(t, 0.28, ICs, filename = 'vdL_ScenarioPlot_NoDispersal_StartingHigh')
	

	# scenario plots with poaching 
	van_de_leemput.set_mgmt_params(closure_length = 35, f = 0.303, m = 1, poaching =  0.2)
	blackwood_mumby.set_mgmt_params(closure_length = 35, f = 0.303141, m = 1, poaching =  0.2)
	rass_briggs.set_mgmt_params(closure_length = 35, f = 0.303, m = 1, poaching =  0.2)

	ICs = blackwood_mumby.X1
	blackwood_mumby.scenario_plot(t, .303, ICs, filename = 'BM_ScenarioPlot_NoDispersal_StartingLow_Poaching')
	ICs = blackwood_mumby.X2
	blackwood_mumby.scenario_plot(t, 0.303, ICs, filename = 'BM_ScenarioPlot_NoDispersal_StartingHigh_Poaching')
	ICs = van_de_leemput.X1
	van_de_leemput.scenario_plot(t, 0.28, ICs, filename = 'vdL_ScenarioPlot_NoDispersal_StartingLow_Poaching')
	ICs = van_de_leemput.X2
	van_de_leemput.scenario_plot(t, 0.28, ICs, filename = 'vdL_ScenarioPlot_NoDispersal_StartingHigh_Poaching')
	ICs = rass_briggs.X1
	rass_briggs.scenario_plot(RB_t, 0.2864, ICs, filename = 'RB_ScenarioPlot_NoDispersal_StartingLow_Poaching')
	ICs = rass_briggs.X2
	rass_briggs.scenario_plot(RB_t, 0.2864, ICs, filename = 'RB_ScenarioPlot_NoDispersal_StartingHigh_Poaching')
	

	# scenario plots for dispersal -- need to initialize new objects due to code structure 
	blackwood_mumby = Model('BM', 10, 0.95, mgmt_strat = 'periodic') 
	van_de_leemput = Model('vdL', 10, 0.95, mgmt_strat = 'periodic')
	rass_briggs = Model('RB', 10, 0.95, mgmt_strat = 'periodic')

	blackwood_mumby.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	van_de_leemput.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	rass_briggs.initialize_patch_model(P0, C0L, C0H, M0vL, M0vH, M0iL, M0iH)

	van_de_leemput.load_parameters()
	rass_briggs.load_parameters()
	blackwood_mumby.load_parameters()
	

	
	# scenario plots for 5 percent of fish moving 

	ICs = blackwood_mumby.X1
	blackwood_mumby.scenario_plot(t, 0.303, ICs, filename = 'BM_ScenarioPlot_SomeDispersal_StartingLow')
	ICs = blackwood_mumby.X2
	blackwood_mumby.scenario_plot(t, 0.303, ICs, filename = 'BM_ScenarioPlot_SomeDispersal_StartingHigh')
	ICs = van_de_leemput.X1
	van_de_leemput.scenario_plot(t, 0.28, ICs, filename = 'vdL_ScenarioPlot_SomeDispersal_StartingLow')
	ICs = van_de_leemput.X2
	van_de_leemput.scenario_plot(t, 0.28, ICs, filename = 'vdL_ScenarioPlot_SomeDispersal_StartingHigh')
	
	ICs = rass_briggs.X1
	rass_briggs.scenario_plot(RB_t, 0.2864, ICs, filename = 'RB_ScenarioPlot_SomeDispersal_StartingLow')

	ICs = rass_briggs.X2
	rass_briggs.scenario_plot(RB_t, 0.2864, ICs, filename = 'RB_ScenarioPlot_SomeDispersal_StartingHigh')


	# for 50 percent of fish moving 

	blackwood_mumby = Model('BM', 10, 0.5, mgmt_strat = 'periodic') 
	van_de_leemput = Model('vdL', 10, 0.5, mgmt_strat = 'periodic')
	rass_briggs = Model('RB', 10, 0.5, mgmt_strat = 'periodic')

	blackwood_mumby.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	van_de_leemput.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	rass_briggs.initialize_patch_model(P0, C0L, C0H, M0vL, M0vH, M0iL, M0iH)

	van_de_leemput.load_parameters()
	rass_briggs.load_parameters()
	blackwood_mumby.load_parameters()
	
 
	ICs = blackwood_mumby.X1
	blackwood_mumby.scenario_plot(t, 0.303, ICs, filename = 'BM_ScenarioPlot_FiftyPercentDispersal_StartingLow')
	ICs = blackwood_mumby.X2
	blackwood_mumby.scenario_plot(t, 0.303, ICs, filename = 'BM_ScenarioPlot_FiftyPercentDispersal_StartingHigh')
	ICs = van_de_leemput.X1
	van_de_leemput.scenario_plot(t, 0.28, ICs, filename = 'vdL_ScenarioPlot_FiftyPercentDispersal_StartingLow')
	ICs = van_de_leemput.X2
	van_de_leemput.scenario_plot(t, 0.28, ICs, filename = 'vdL_ScenarioPlot_FiftyPercentDispersal_StartingHigh')
	ICs = rass_briggs.X1
	rass_briggs.scenario_plot(RB_t, 0.2864, ICs, filename = 'RB_ScenarioPlot_FiftyPercentDispersal_StartingLow')
	ICs = rass_briggs.X2
	rass_briggs.scenario_plot(RB_t, 0.2864, ICs, filename = 'RB_ScenarioPlot_FiftyPercentDispersal_StartingHigh')




	van_de_leemput.coral_recovery_map(t, 0.28, van_de_leemput.X1, filename = 'vdL_heatmap_0.28_startinglow')
	blackwood_mumby.coral_recovery_map(t, 0.303, blackwood_mumby.X1, filename = 'BM_heatmap_0.303_startinglow')
	rass_briggs.coral_recovery_map(t, 0.2864, rass_briggs.X1, filename = 'RB_heatmap_0.2864_startinglow')
	van_de_leemput.coral_recovery_map(t, 0.28, van_de_leemput.X2, filename = 'vdL_heatmap_0.28_startinghigh')
	blackwood_mumby.coral_recovery_map(t, 0.303, blackwood_mumby.X2, filename = 'BM_heatmap_0.303_startinghigh')
	rass_briggs.coral_recovery_map(t, 0.2864, rass_briggs.X2, filename = 'RB_heatmap_0.2864_startinghigh')





	blackwood_mumby.bistable_zone(t, filename = 'BM_bistable_zone')
	rass_briggs.bistable_zone(t, filename = 'RB_bistable_zone')
	van_de_leemput.bistable_zone(t, filename = 'vdL_bistable_zone')



	'''
	blackwood_mumby.set_mgmt_params(25, 0.303, 1, 0)
	blackwood_mumby.time_series(blackwood_mumby.X2, t, save = False, show = True) # why does this flatten out?

	rass_briggs.set_mgmt_params(25, 0.303, 1, 0)
	rass_briggs.time_series(rass_briggs.X2, t, save = False, show = True) # why does this also flatten out?

	van_de_leemput.set_mgmt_params(30, 0.303, 1, 0)
	van_de_leemput.time_series(van_de_leemput.X2, t, save = False, show = True) 

	'''
	


	





	# ------------------------------------ DEMO ---------------------------------- #





	van_de_leemput = Model('vdL', n = 10, frac_nomove = 0.8, mgmt_strat = 'periodic')

	van_de_leemput.load_parameters()
	van_de_leemput.initialize_patch_model(P0, C0L, C0H, M0L, M0H)


	van_de_leemput.set_mgmt_params(closure_length = 30, f = 0.303, m = 1, poaching =  0.2)


	# van_de_leemput.coral_recovery_map(t, 0.303, IC_set = van_de_leemput.X1)


	van_de_leemput.time_series(van_de_leemput.X1, t, save = False, show = True) 



	van_de_leemput.set_mgmt_params(30, 0.303, 2, poaching = 0.1)

	van_de_leemput.time_series(van_de_leemput.X1, t, save = False, show = True) 

	van_de_leemput.set_mgmt_params(30, 0.303, 3, 0)
	van_de_leemput.time_series(van_de_leemput.X1, t, save = False, show = True) 

	

	# rass_briggs = Model('RB', 4, 1, mgmt_strat = 'periodic')


	# van_de_leemput.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	# blackwood_mumby.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	# rass_briggs.initialize_patch_model(P0, C0L, C0H, M0vL, M0vH, M0iL, M0iH)


	# van_de_leemput.load_parameters()
	# rass_briggs.load_parameters()
	# blackwood_mumby.load_parameters()

	# van_de_leemput.coral_vs_CL(t, 0.2864, IC_set = van_de_leemput.X1)


	# blackwood_mumby.set_mgmt_params(30, 0.303, 1, 0)
	# blackwood_mumby.time_series(blackwood_mumby.X1, t, save = False, show = True) # 

	# rass_briggs.set_mgmt_params(50, 0.303, 1, 0)
	# rass_briggs.time_series(rass_briggs.X1, t, save = False, show = True) # 
# 
	# van_de_leemput.set_mgmt_params(25, 0.2864, 1, 0)
# 	van_de_leemput.time_series(van_de_leemput.X1, t, save = False, show = True) 
# # 
# 	van_de_leemput.set_mgmt_params(25, 0.2864, 2, 0)
# 	van_de_leemput.time_series(van_de_leemput.X1, t, save = False, show = True) 
# # 
# 	van_de_leemput.set_mgmt_params(25, 0.2864, 3, 0)
# 	van_de_leemput.time_series(van_de_leemput.X1, t, save = False, show = True) 
# 
	# van_de_leemput.set_mgmt_params(40, 0.2864, 1, 0)
	# van_de_leemput.time_series(van_de_leemput.X2, t, save = False, show = True) 
# 
	# van_de_leemput.set_mgmt_params(40, 0.2864, 2, 0)
	# van_de_leemput.time_series(van_de_leemput.X2, t, save = False, show = True) 
# 
	# van_de_leemput.set_mgmt_params(40, 0.2864, 3, 0)
	# van_de_leemput.time_series(van_de_leemput.X2, t, save = False, show = True) 






	# blackwood_mumby.set_mgmt_params(30, 0.303, 1, 0)
	# blackwood_mumby.time_series(blackwood_mumby.X2, t, save = False, show = True) # 
# 

	# rass_briggs.set_mgmt_params(50, 0.303, 1, 0)
	# rass_briggs.time_series(rass_briggs.X2, t, save = False, show = True) # 

	# van_de_leemput.set_mgmt_params(30, 0.303, 1, 0)
	# van_de_leemput.time_series(van_de_leemput.X2, t, save = False, show = True) 


	# blackwood_mumby = Model('BM', 12,  1, mgmt_strat = 'periodic') # ISSUE WITH DISPERSAL: kP calculation or P_influx must be incorrect 
	
	# van_de_leemput = Model('vdL', 42, 1, mgmt_strat = 'periodic')

	
	# van_de_leemput.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	# blackwood_mumby.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	# rass_briggs.initialize_patch_model(P0, C0L, C0H, M0vL, M0vH, M0iL, M0iH)


	# van_de_leemput.load_parameters()
	# rass_briggs.load_parameters()
	# blackwood_mumby.load_parameters()
	# blackwood_mumby.bistable_zone(t)
	# rass_briggs.bistable_zone(t)
	# van_de_leemput.bistable_zone(t)


	# ICs = blackwood_mumby.X1
# 	blackwood_mumby.scenario_plot(t, 0.303, ICs)# 

# 	ICs = blackwood_mumby.X2# 

# 	blackwood_mumby.scenario_plot(t, 0.303, ICs)
# 	ICs = van_de_leemput.X1
# 	van_de_leemput.scenario_plot(t, 0.303, ICs)
# 	ICs = van_de_leemput.X2# 

# 	van_de_leemput.scenario_plot(t, 0.303, ICs)
# 	ICs = rass_briggs.X1# 
# 	rass_briggs.scenario_plot(t, 0.303, ICs)
# 	ICs = rass_briggs.X2# 

# 	rass_briggs.scenario_plot(t, 0.303, ICs)

	# blackwood_mumby.set_mgmt_params(30, 0.26, 0, 0)
	# blackwood_mumby.time_series(blackwood_mumby.X1, t, save = False, show = True) # 
# 

	# van_de_leemput.coral_recovery_map(t, 0.303, van_de_leemput.X1)
	# blackwood_mumby.coral_recovery_map(t, 0.303, blackwood_mumby.X1)
	# rass_briggs.coral_recovery_map(t, 0.2864, rass_briggs.X1)
	# van_de_leemput.coral_recovery_map(t, 0.303, van_de_leemput.X2)
	# blackwood_mumby.coral_recovery_map(t, 0.303, blackwood_mumby.X2)
	# rass_briggs.coral_recovery_map(t, 0.2864, rass_briggs.X2)
	




	
	# load Model parameters according to model type
	# x.load_parameters()
	# y.load_parameters()

	# set initial conditions 
	# x.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	# y.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	# y.RB_initialize_patch_model(P0, C0L, C0H, M0vH, M0vL, M0iH, M0iL)

	# z.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	# z.load_parameters() # do this inside initializer

	# print(y.get_coral_recovery_time(t))

	# y.bistable_zone(t)


	# y.set_mgmt_params(10, 0.2864, 2, 0)
	# y.time_series(y.X1, t, save = False, show = True) 

	# is first param closure length or period???
	# y.set_mgmt_params(25, 0.2864, 2, 0)
	# y.time_series(y.X1, t, save = False, show = True) 


	# y.set_mgmt_params(100, 0.2864, 2, 0)
	# y.time_series(y.X1, t, save = False, show = True) 


	# y.time_series(y.X1, t, False, True)
	# x.set_mgmt_params(20, 0.1, 1, 0) # set management parameters -- closure length, fishing effort, # of closures, poaching 
	# print(x.run_model(x.X1, t)) # IT WORKED !!!!!!

	# x = y.find_unstable_equilibrium(t)
	# print(x)

	
	# y.time_series(y.X1, t, save = False, show = True) 

	# x.set_mgmt_params(40, 0.4, 1, 0.5)
	# x.time_series(x.X1, t, False, True)
	'''
	x.coral_recovery_map(t, 0.2)
	x.coral_recovery_map(t, 0.4)
	x.coral_recovery_map(t, 0.6)
	y.set_mgmt_params(40, 0.4, 1, 0)
	y.coral_recovery_map(t, 0.15)
	y.coral_recovery_map(t, 0.20)
	y.coral_recovery_map(t, 0.2864)
	y.coral_recovery_map(t, 0.3030)
	'''
	
	
	# y.coral_recovery_map(t, 0.27)
	# y.coral_recovery_map(t, 0.3035)
	# z.coral_recovery_map(t, 0.2864)
	# z.coral_recovery_map(t, 0.3035)
	
	# print(x.run_model(parameter list))
	# x.time_series(parameter list, show = True)


if __name__ == '__main__':
	main()




