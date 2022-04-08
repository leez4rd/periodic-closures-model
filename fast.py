#fast heatmap
#simplifies m over n to lowest term and does least number of patches necessary 
import numpy as np 

def fast_heatmap():
	patches = 10

	#heatmap of period vs m -- a bit messy right now 
	coral_array_HI =  np.zeros(shape=(patches,patches))
	coral_array_LO =  np.zeros(shape=(patches,patches))
	coral_array_AVG =  np.zeros(shape=(patches,patches))
	period_array = np.empty(patches)
	m_array = np.empty(patches)
	fishin = 0.35
	n = 10

	frac_nomove = 1
	for period in range(1,patches+1):
		for m in range(patches):
			GCD = np.gcd(m, n)
			m = m / GCD
			n = n / GCD
			initialize_patch_model(n, frac_nomove)
			displacer = 1/(1-m/float(n))
			final_coral_LO = odeint(patch_system, X1, t, args = (period*2, displacer*float(fishin), m, n, 0))
			avg2 = 0
			for year in range(999- (999 % (n*period*2)) - 2*(n*period*2), 999 - (999 % (n*period*2))):
				avg2 += final_coral_LO[year][n]
			avg2 = avg2 / (2*(period*n*2) + 1)
			print(avg2)
			print("________________________")
			coral_array_LO[period-1][m] = avg2
			period_array[period-1] = period
			m_array[m] = m
			show_labels = False

	plt.title('heatmap', fontsize = 20)
	sb.heatmap(coral_array_LO, vmin = 0.0, vmax = 1.0, cmap = "viridis", annot = show_labels) #YlGnBu for original color scheme
	plt.ylabel('period', fontsize = 10) # x-axis label with fontsize 15
	plt.xlabel('number of closures', fontsize = 10) # y-axis label with fontsize 15
	plt.yticks(rotation=0)
	plt.show()

