#! python2.7

#Program for calculating probability that a measured cell barcode is not compromised
import numpy as np
from binomial import poisson, poisson_distr

bd = 0.06		#bead dilution factor (mean beads/droplet)
cd = 0.05		#cell dilution factor (mean cells/droplet)
beads = poisson_distr(bd,3) #probability distr for number beads per droplet
cells = poisson_distr(cd,3) #probability distr for number cells per droplet
bead_roommates = np.array(beads[1:])/(1-beads[0]) #probability for number of 'roommates' sharing a droplet | at least one bead
cell_roommates = np.array(cells[1:])/(1-cells[0]) #probability for number of 'roommates' sharing a droplet | at least one cell


n = 10**7		#number of beads synthesized
bl = 12			#length of barcode
bp = 4**bl		#number of possible barcodes

u = (n-1)/float(bp)	#population mean for number of repeats per barcode
#nu = n/(1+u)		#number of beads containing unique barcodes
#nr = n-nu		#number of beads containing repeat barcodes
repeats = poisson_distr(u,8) #probability distr for number of repeats per barcode


p = repeats[0] #p = prob that barcode does not repeat
for k in range(1,len(repeats)):
	p+= repeats[k]*cells[0]**k
			#p = prob that given a barcode, it was unique to a cell droplet
p_valid = bead_roommates[0]*cell_roommates[0]*p	#probability that given a barcode, it came into a one-to-one relationship with a single cell

p1 = 1-p #Probability that given a barcode, it was not unique to a cell droplet
p2 = p*(1-cell_roommates[0]) #Probability that given a barcode, it was unique to a cell droplet but the droplet had cell roommates
p3 = p*cell_roommates[0]*(1-bead_roommates[0]) #Probability that given a barcode, it was unique to a droplet containing a single cell but that it had bead roommates.

p4 = p1 + p2 #Probability that given a barcode, it labeled the contents of more than one cell
p5 = 1-bead_roommates[0] #Probability that given a barcode, it had barcode roommates

print('Probability that given a cell barcode, it formed a unique and bidirectional relationship with a single cell: %f\n' %p_valid)

print('\tProbability that given a cell barcode, it labeled the contents of more than one cell: %f' %p4)
print('\t\tProbability that given a cell barcode, it was not unique to a cell droplet: %f' %p1)
print('\t\tProbability that given a cell barcode, it was unique to a cell droplet but the droplet had cell roommates: %f\n' %p2)

print('\tProbability that given a cell barcode, it had barcode roommates: %f' %p5)
print('\t\tProbability that given a cell barcode, it was unique to a cell droplet and labeled a single cell, but barcodes on its bead roommates also labeled the same cell: %f' %p3)




