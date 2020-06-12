

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy.stats import multinomial
import time
import pandas as pd
from adjustText import adjust_text

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


cell_starts = [2e8, 1.35e8, 3.01e9, 6.5e7, 1.85e7]
times = [5840, 0, 88, 6, 7, 1720]
mutations = [3, 2, 8, 3, 5, 4]
name = ["Colorectal", "Glioblastoma", "Hepatocellular", "Lung (nonsmoker)", "Thyroid", "Head & Neck"]
risk = [0.048, 0.00219, 0.0071, 0.0045, 0.01026, 0.0138]
n_runs = 50000

calculated = []
lscd = []
for i in range(len(cell_starts)):
    start_time = time.time()
    #Finding age where cancer appears
    print(name[i])
    cell_start = cell_starts[i]
    t = times[i]

    cancer_age = []

    n_cancer_during = 0
    n_cancer = 0

    for run in range(n_runs):
        
        if(time.time() - start_time > 120): #stop if it takes too long
            print("Timed out." , time.time()-start_time)
            print("runs = ", run)
            break
            
        if(run > n_runs/10):
                if (n_cancer_during > 0.99*run):
                    break

        if run % 1000 == 0:
            print(run)
            print("n_cancer = ", n_cancer_during)

        cells = np.zeros(mutations[i] + 1, dtype=np.int64)
        cells[0] = 1

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

    print("Percentage = ", n_cancer_during / (run+1))
    calculated.append(n_cancer_during / (run+1))
    lscd.append(cell_start*(2+t) - 2)
#     print("Expected = ", risk[i])
    
    
    
for i in range(len(calculated)):
    if(calculated[i] == 0):
        calculated[i] = 1/(n_runs * 10000000)
        
fig, ax = plt.subplots()
for i in range(len(lscd)):
    color = next(ax._get_lines.prop_cycler)['color']
    ax.loglog(lscd[i], calculated[i],'o', color = color)
    ax.loglog(lscd[i], risk[i],'*', color = color)
# ax.set_yscale('symlog')
# ax.set_xscale('log')

# texts1 = [plt.text(lscd[i], calculated[i], name[i]) for i in range(len(lscd))]
# texts2 = [plt.text(lscd[i], risk[i], name[i]) for i in range(len(lscd))]
# adjust_text(texts1, lscd, calculated, arrowprops=dict(arrowstyle="->", color='b', lw=0.5),
#             autoalign='', only_move={'points':'y', 'text':'y'})
# adjust_text(texts2, lscd, risk, arrowprops=dict(arrowstyle="->", color='b', lw=0.5),
#             autoalign='', only_move={'points':'y', 'text':'y'})
# for i, txt in enumerate(name):
#     ax.annotate(txt, (lscd[i], calculated[i]))

for i, txt in enumerate(name):
    ax.annotate(txt, (lscd[i], risk[i]))
    
# plt.loglog(lscd, calculated, 'o')
# plt.loglog(lscd, risk, 'o')
plt.legend(["Calculated", "Risk"])
plt.title("Calculated incidence and observed risk in a number of cancers")
plt.xlabel('Total stem cell divisions')
plt.ylabel('Lifetime risk')

# x_data = lscd
# y_data = calculated

# txt_height = 0.04*(plt.ylim()[1] - plt.ylim()[0])
# txt_width = 0.02*(plt.xlim()[1] - plt.xlim()[0])
# text_positions = get_text_positions(x_data, y_data, txt_width, txt_height)
# text_plotter(x_data, y_data, text_positions, ax, txt_width, txt_height, name)

plt.show()
    # plt.semilogy(np.array(cells_all))
    # # plt.legend()
    # # plt.xlim([0,2000])
    # plt.xlabel("Number of time steps")
    # plt.ylabel("Number of cells")
    # plt.show()