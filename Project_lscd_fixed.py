

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy.stats import multinomial
import time


s = 0.004
u = 3.4e-5

def generation(cells, s_this, Print = True, stop = False, beginning = False):
    #runs one generation of the branching process
    
    #cells is a list of how many types of cells
    cells_next = cells.copy()
    for i in range(len(cells)):
        if(cells[i]> 0):
            #get the number died, born
            dj = ((1-s_this)**(i))/2
            bj = (1-dj)
            if (i == 0):
                if beginning: 
                    dj = 0
                    bj = 1
            die = np.int64(0)
            born = np.int64(0)
            mutate = np.int64(0)
            n = cells[i]
            while n > 1e9:
                change = np.random.multinomial(1e9, (dj, bj*(1-u), bj*u))
                die +=change[0]
                born = np.int64(born) + np.int64(change[1])
                mutate += change[2]
                n -= 1e9

            change = np.random.multinomial(int(cells[i]), (dj, bj*(1-u), bj*u))
            die +=change[0]
            born +=change[1]
            mutate += change[2]                
            
            #update died
            cells_next[i] -= np.int64(die)
            
            #update born
            cells_next[i] += np.int64(born)

            #make it so there's not negative
            if cells_next[i] < 0:
                cells_next[i] = 0

            #update mutate
            if mutate >0:
                try:
                    cells_next[i + 1] += np.int64(mutate)
                except:
                    if Print:
                        print("Mutate to next generation, end of list")
                    if stop:
                        return [0]
                    return cells_next
    return cells_next

ratio = 7.4e7/80
# cell_start = 7.4e7
# lscd_list = np.logspace(8.5, 11, 10)
t = 2920
lscd_list = np.logspace(6,15, 10)
n_runs = 1000

for t in [10, 80, 1000, 2900]:

    lscd_actual = []

    calculated = []
    # lscd = []
    for i in range(len(lscd_list)):
        start_time = time.time()
        lscd = lscd_list[i]
        print("lscd = ", np.log10(lscd))
    #     time = int((lscd + 2)/cell_start - 2)
        cell_start = int((lscd+2)/(2+t))
        print("time = ", t)
        print("cell_start = ", np.log10(cell_start))
        lscd_actual.append(cell_start*(2+t) - 2)


        cancer_age = []

        n_cancer_during = 0
        n_cancer = 0

        for run in range(n_runs):
            if(time.time() - start_time > 100): #stop if it takes too long
                print("Timed out." , time.time()-start_time)
                print("runs = ", run)
                break
            if(run > n_runs/10):
                if (n_cancer_during > 0.99*run):
                    break

#             if run % 999 == 0:
#                 print(run)
#                 print("n_cancer = ", n_cancer_during)

            cells = np.array([1, 0, 0, 0], dtype=np.int64)

            cells_all = []

            for gen in range(int(np.log2(cell_start))):
                cells = generation(cells,1, Print = False, beginning = True)
                cells_all.append(cells)

            cells[0] = cell_start

            for gen in range(t):
                cells = generation(cells,s, Print = False)
                if len(cells) > 1: #check if hit the end
                    cells_all.append(cells)
                else:
                    going = False #Hit end of mutation list must means it survivied. 
                    break 
                if(cells[-1]>1):
                    #coin flip
                    n = np.random.random()
                    if(n < 2*s*len(cells)):
                        n_cancer_during += 1
                        cancer_age.append(gen)
                        break
                    else:
                        cells[-1] = 0

        print("Percentage = ", n_cancer_during / run)
        calculated.append(n_cancer_during / run)
        if(n_cancer_during == run):
            break
    #     lscd.append(cell_start*(2+time) - 2)
    #     print("Expected = ", risk[i])
    print("lscd_actual = ", lscd_actual)
    print("calculated = ", calculated)
    plt.loglog(lscd_actual, calculated, '.', label = "Fixed time steps = " + str(t))
plt.xlabel("lscd")
plt.ylabel("Percent Cancer")
plt.legend()
plt.title("lscd vs Simulated Cancer incidence for a variety of fixed time steps")
plt.show()
