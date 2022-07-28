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

	def initialize_patch_model(self, P0, C0L, C0H, M0L, M0H):
		# do variables defined in this function only exist in the scope of this function if it is called by constructor?

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
	def RB_initialize_patch_model(self, P0, C0L, C0H, M0vH, M0vL, M0iH, M0iL):

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


		for i in range(self.n):
			ax.plot(sol[:, self.n+i],  label = 'patch % d'% (int(i) + 1))
		locs, labels = plt.xticks()
		print(locs)

		newlabels = [int(int(loc) / int(recov_time)) for loc in locs[1:]]

		print(newlabels)

		ax.set_xticks([recov_time*tick for tick in range(0, 7)])
		ax.set_xticklabels(newlabels)
		ax.set_xlabel('Time (scaled to coral recovery time)')
		
		ax.set_ylabel('Coral cover (total fraction of area)')
		ax.set_title('Coral cover in each patch over time')
		plt.legend(loc=0)
		ax.set_xlim([0, 500])
		ax.set_ylim([0, 1])
		txt="parameters" + "\nfishing when open: " + str(self.f/(1-self.m/self.n)) + "\npercent closure: " + str(self.m/self.n) +"\nclosure time: " + str(self.closure_length)
		# plt.figtext(0.7, .31, txt, wrap=True, fontsize=8)
		if save:
			plt.savefig('model_time_series1.jpg')
		if show:
			plt.show()

	# make a flag for fast or slow version -- also consider making nonspatial for faster runtime 
	# note that resolution is determined by the previously initialized patch number for spatial heatmaps
	def coral_recovery_map(self, t, fishing_intensity):


		""" MPA color bar for comparison -- really hacky right now bc it involves creating a Model object within current one """

		mako = ListedColormap(sb.color_palette('mako').as_hex())

		fig, (ax1,ax2) = plt.subplots(nrows=2, sharex=True, figsize=(10,10), gridspec_kw={'height_ratios': [1, 1.5*self.n]})


		P0, C0L, C0H, M0L, M0H, M0vH, M0vL, M0iH, M0iL = 0.1, 0.04, 0.4, 0.04, 0.4, 0.04, 0.4, 0.04, 0.4

		""" Make MPA colorbar for comparisons """
		# hacky but works 
		# need to simulate MPA for second colorbar to juxtapose with heatmap...
		# creating another model object within this one seems unnecessarily funky 
		
		z = Model(self.model_type, self.n, 1, mgmt_strat = 'MPA')
		mpa_corals = np.empty(self.n)
		# set initial conditions 
		z.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
		z.load_parameters() # do this inside initializer
		for i in range(self.n):
			z.set_mgmt_params(500, 0.25, i, 0)
			MPAsol = z.run_model(z.X2, t)
			print(MPAsol)
			extent = [0, self.n, 0, 1]
			coral_f = 0
			for j in range(self.n):
				coral_f += MPAsol[999][j+self.n]
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
		ax1.set_position([position.x0, position.y0 + position.y1, position.x1, 1.1])
		# slow version 
		IC_set = self.X1 	# default to coral-rare 
		# IC_set = self.X2
		MAX_TIME = len(t) # last year in the model run 

		coral_array =  np.zeros(shape=(self.n, self.n)) # array of long-term coral averages
		# CL_array = np.empty(int(0.75*self.n+1)) # array of closure lenghts 
		# m_array = np.empty(int(self.n / 2))  # array of number of closures 
		closure_lengths = np.empty(self.n)
		ms = np.empty(self.n)
		# iterate over  all scenarios 
		for closure_length in range(1,self.n + 1):
			period = self.n*closure_length
			for m in range(self.n):
				# set management parameters for this run 
				self.set_mgmt_params(closure_length, fishing_intensity, m, self.poaching)

				# solve ODE system 
				sol = odeint(patch_system, IC_set, t, args = (self, ))
				# average over coral cover of last two full rotations for a single patch (assumes symmetry, may fix that)
				avg = 0

				# dealing with very large period / millennium ratio
				
				for year in range(MAX_TIME - MAX_TIME % period - period, 
					MAX_TIME - MAX_TIME % period):
					avg += sol[year][self.n]
				
				avg = avg / ((period) + 1)
				print("Average coral cover for closure length: ", closure_length, " and closure num: ", m)
				print(avg)
				coral_array[closure_length-1][m] = avg
				ms[m] = m
				closure_lengths[closure_length - 1] = closure_length
				# CL_array[closure_length-1] = closure_length -- don't think this is necessary 
				# m_array[m] = m
		# plt.figure()
		
		# fig, ax = plt.subplots(figsize=(5, 5))
		# axes = plot.axes(projection='3d')
		# axes.plot_surface(X1, Y1, Z1)
		# plt.show()
		ax2.set_title('Coral cover under periodic closures', fontsize = 15)
		f = lambda y:self.n*y
		new_labels = [f(y) for y in range(1, self.n+1)]

		ax2 = sb.heatmap(coral_array, vmin = 0.0, vmax = 1.0, cmap = "jet", yticklabels = new_labels, cbar = True) #YlGnBu for original color scheme
		ax2.invert_yaxis()
		ax2.set_ylabel('Period in years', fontsize = 10) # x-axis label with fontsize 15
		ax2.set_xlabel('Number of patches closed out of {} total patches'.format(self.n), fontsize = 10) # y-axis label with fontsize 15
		# ax.yticks(ax.get_yticks(), ax.get_yticks() * 3)
		# ax.set_yticks(rotation=0)
		# plt.show()


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


		ax2.plot([self.n / 2 + 1]*100, y, 'r-')

		


		# show it all 
		plt.show()
		
		name = 'longtermcoral_fishing_' + str(fishing_intensity) + '_' + str(self.n) + '.jpg'
		# plt.savefig(name)
		plt.close()

	def bistable_zone(self, t):
		""" plot final coral cover for different values of fishing effort for two sets of initial conditions """ 
		final_coral_high = np.empty(100)
		final_coral_low = np.empty(100)


		fishing_range = np.linspace(0, 0.99, 100)

		for f in fishing_range:

			# set management parameters 
			self.set_mgmt_params(0, f, 0, self.poaching)

			# make high start solution
			high_sol = odeint(patch_system, self.X2, t, args = (self, ))

			# make low start solution 
			low_sol = odeint(patch_system, self.X1, t, args = (self, ))

			# note: this only works without periodic oscillations, which this plot assumes are not present 
			yrs = len(t)
			print(f*100)
			final_coral_high[int(f * 100)] = high_sol[yrs - 1][self.n]

			final_coral_low[int(f * 100)] = low_sol[yrs - 1][self.n]

		plt.plot(fishing_range, final_coral_low, label = 'coral starts low', color = 'blue')
		plt.plot(fishing_range, final_coral_high, label = 'coral starts high' , color = 'green')
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

	def scenario_plot(self, t, fishing_intensity, IC_set):
		P0, C0L, C0H, M0L, M0H, M0vH, M0vL, M0iH, M0iL = 0.1, 0.04, 0.4, 0.04, 0.4, 0.04, 0.4, 0.04, 0.4


		final_coral = np.empty(self.n)
		ms = np.empty(self.n)
		periods = [5, 20, 50, 100, 200]  # parametrize in terms of coral growth time? 
		# there is a cooler way to do colors than this 
		color_sequence = {5: '#1f77b4', 20: '#aec7e8', 50: '#ff7f0e', 100:'#ffbb78', 200:'#2ca02c'}

		MAX_TIME = len(t)
		

		for period in periods:
			print(period)
			for m in range(self.n):
				print(m)

				# set management parameters for this run 
				self.set_mgmt_params(period / self.n, fishing_intensity, m, self.poaching)

				# solve ODE system 
				sol = odeint(patch_system, IC_set, t, args = (self, ))

				avg = 0

				# truncate to last full period, then average over that period 
				for year in range(MAX_TIME - MAX_TIME % period - period, MAX_TIME - MAX_TIME % period):
					avg += sol[year][self.n] # only looking at one patch here
					# may amend to average over all patches for comparison with MPA 
				
				avg = avg / (period + 1)
				final_coral[m] = avg
				print(avg)
				ms[m] = m

			# plot result for this period
			
			plt.plot(ms / self.n, final_coral, label = 'period = %d' % period, color = color_sequence[period])
		plt.xlabel('Fraction closed')
		plt.ylabel('Final coral Cover')
		plt.title('Final coral state across closure scenarios')
		
		z = Model(self.model_type, self.n, 1, mgmt_strat = 'MPA')
		mpa_corals = np.empty(self.n)
		# set initial conditions 
		z.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
		z.load_parameters() # do this inside initializer
		
		# loop over frac closed 
		for i in range(self.n):
			z.set_mgmt_params(5, 0.25, i, 0)
			MPAsol = z.run_model(z.X2, t) # coral starting low 
			
			# loop over patches 
			total = 0 
			for j in range(self.n):
				coral_f = 0
				for year in range(MAX_TIME - MAX_TIME % period - period, MAX_TIME - MAX_TIME % period):
					coral_f += MPAsol[year][j+self.n]

				coral_f = coral_f / period
				print(coral_f)
				total += coral_f
			total = total / self.n
			print(total)


			# coral_f = coral_f / self.n # patch average
			
			mpa_corals[i] = total
			print(mpa_corals)

		arr = np.linspace(0, 1, self.n)
		plt.plot(arr, mpa_corals, label = 'MPA', color = 'black')
		plt.xlim([0, 0.75])
		plt.legend(loc=2)

		plt.show() 

	
	def get_coral_recovery_time(self, t):

		P0, C0L, C0H, M0L, M0H, M0vH, M0vL, M0iH, M0iL = 0.1, 0.04, 0.4, 0.04, 0.4, 0.04, 0.4, 0.04, 0.4


		z = Model(self.model_type, self.n, 1, mgmt_strat = 'MPA')
		# set initial conditions 
		z.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
		z.load_parameters() # do this inside initializer

		z.set_mgmt_params(500, 0, 0, 0)

		MPAsol = z.run_model(z.X1, t)
		# list_of_labels = ['fish']*self.n + ['coral']*self.n + ['algae']*self.n
		# plt.plot(t, MPAsol, label = list_of_labels)
		# plt.legend(loc=0)
		# plt.show()


		coral_recovery_time = -1
		for i, state in enumerate(MPAsol):
			if state[self.n] > 0.60: # what number should I put here?
				coral_recovery_time = i
				break

		return coral_recovery_time


	# modify this to take custom feedback parameters, or maybe custom anything? 
	# this doesn't work for some reason
	def load_parameters(self):
		if self.model_type == 'vdL':
			params = {
			"r": 0.3,
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
			"r": 0.3,
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
			"r": 0.3,
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
			"r": 0.3,
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
			for j in range(system_model.n):
				P_influx[i] += (system_model.kP[i][j]) * system_model.P[j]  
	
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
	dP = P_influx[i]+system_model.rH*P[i]*(1-P[i]/system_model.K) - system_model.f/(1-system_model.m/system_model.n)*P[i] *(square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching))
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
		if region >= start and region <= end:
			return poaching
		elif start + m - 1 >= n and (region >= start or region <= end):
			return poaching
		else:
			#determine whether to bake displacement into signal
			return (1-(m/n)*poaching)#/(1-(m/n))
	elif mgmt_strat == 'MPA':
		if region < m:
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

	dP = system_model.s*P[i]*(1 - (P[i] / (system_model.beta*BMK(C[i])))) - system_model.f/(1-system_model.m/system_model.n)*P[i]*square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching)
	dC = system_model.r*(1-M[i]-C[i])*C[i]-system_model.d*C[i] - system_model.a*M[i]*C[i] + system_model.i_C*(1-M[i]-C[i])
	# need to define g(P) before this model is used 
	dM = system_model.a*M[i]*C[i] - system_model.alpha*P[i]/system_model.beta*M[i] *(1/(1-C[i])) + system_model.gamma*M[i]*(1-M[i]-C[i])+system_model.i_M*(1-M[i]-C[i])

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

	dP = P_influx[i]+ system_model.s*P[i]*(1 - (P[i] / K(system_model.sigma,C[i]))) - system_model.f/(1-system_model.m/system_model.n)*P[i] *(square_signal(t, system_model.closure_length, i, system_model.m, system_model.n, system_model.poaching, system_model.mgmt_strat))
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

	P0, C0L, C0H, M0L, M0H, M0vH, M0vL, M0iH, M0iL = 0.1, 0.04, 0.4, 0.04, 0.2, 0.04, 0.2, 0.04, 0.2

	# create Model objects
	x = Model('vdL', 3, 1)

	y = Model('RB', 5,  1) # ISSUE WITH DISPERSAL: kP calculation or P_influx must be incorrect 
	
	z = Model('vdL', 20, 1)
	
	# load Model parameters according to model type
	x.load_parameters()
	y.load_parameters()

	# set initial conditions 
	x.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	# y.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	y.RB_initialize_patch_model(P0, C0L, C0H, M0vH, M0vL, M0iH, M0iL)

	z.initialize_patch_model(P0, C0L, C0H, M0L, M0H)
	z.load_parameters() # do this inside initializer

	# print(y.get_coral_recovery_time(t))

	# y.bistable_zone(t)


	y.set_mgmt_params(10, 0.25, 2, 0)
	# y.time_series(y.X1, t, save = False, show = True) 

	# is first param closure length or period???
	y.set_mgmt_params(25, 0.25, 2, 0)
	# y.time_series(y.X1, t, save = False, show = True) 


	y.set_mgmt_params(100, 0.25, 2, 0)
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
	y.coral_recovery_map(t, 0.25)
	y.coral_recovery_map(t, 0.30)
	'''
	ICs = y.X2
	y.scenario_plot(t, 0.25, ICs)


	y.coral_recovery_map(t, 0.27)
	# y.coral_recovery_map(t, 0.35)
	# z.coral_recovery_map(t, 0.25)
	# z.coral_recovery_map(t, 0.35)
	
	# print(x.run_model(parameter list))
	# x.time_series(parameter list, show = True)


if __name__ == '__main__':
	main()