

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

def get_text_positions(x_data, y_data, txt_width, txt_height):
    a = zip(y_data, x_data)
    text_positions = y_data.copy()
    for index, (y, x) in enumerate(a):
        local_text_positions = [i for i in a if i[0] > (y - txt_height)
                            and (abs(i[1] - x) < txt_width * 2) and i != (y,x)]
        if local_text_positions:
            sorted_ltp = sorted(local_text_positions)
            if abs(sorted_ltp[0][0] - y) < txt_height: #True == collision
                differ = np.diff(sorted_ltp, axis=0)
                a[index] = (sorted_ltp[-1][0] + txt_height, a[index][1])
                text_positions[index] = sorted_ltp[-1][0] + txt_height
                for k, (j, m) in enumerate(differ):
                    #j is the vertical distance between words
                    if j > txt_height * 1.5: #if True then room to fit a word in
                        a[index] = (sorted_ltp[k][0] + txt_height, a[index][1])
                        text_positions[index] = sorted_ltp[k][0] + txt_height
                        break
    return text_positions


def text_plotter(x_data, y_data, text_positions, axis,txt_width,txt_height, names):
    for x,y,t, name in zip(x_data, y_data, text_positions, names):
        axis.text(x - .03, 1.02*t, name,rotation=0, color='blue', fontsize=13)
        if y != t:
            axis.arrow(x, t+20,0,y-t, color='blue',alpha=0.2, width=txt_width*0.0,
                       head_width=.02, head_length=txt_height*0.5,
                       zorder=0,length_includes_head=True)

cell_starts = [7.4e7, 3.01e9, 1e8, 4e6, 4.18e6, 2.8e9, 1.22e9]
times = [80, 88, 2920, 1947, 5, 199, 6]
name = ["Pancreatic islet", "Hepatocellular", "Small Intestine", "Duodenum", "Osteosarcoma", "Melanoma", "Lung Adenocarcinoma"]
risk = [0.0000194, 0.0071, 0.00007, 0.0003, 0.00035, 0.0203, 0.0045]
n_runs = 50000

# calculated = []
lscd = []
for i in range(len(cell_starts)):
#     start_time = time.time()
#     #Finding age where cancer appears
#     print(name[i])
    cell_start = cell_starts[i]
    t = times[i]

#     cancer_age = []

#     n_cancer_during = 0
#     n_cancer = 0

#     for run in range(n_runs):
        
#         if(time.time() - start_time > 120): #stop if it takes too long
#             print("Timed out." , time.time()-start_time)
#             print("runs = ", run)
#             break
            
#         if(run > n_runs/10):
#                 if (n_cancer_during > 0.99*run):
#                     break

#         if run % 1000 == 0:
#             print(run)
#             print("n_cancer = ", n_cancer_during)

#         cells = np.array([1, 0, 0, 0], dtype=np.int64)

#         cells_all = []

# #         for gen in range(int(np.log2(cell_start))):
# #             cells = generation(cells,1, Print = False, beginning = True)
# #             cells_all.append(cells)

#         cells[0] = cell_start

#         for gen in range(t):
#             cells = generation(cells,s, Print = False)
#             if len(cells) > 1: #check if hit the end
#                 cells_all.append(cells)
#             else:
#                 going = False #Hit end of mutation list must means it survivied. 
#                 break 
#             if(cells[-1]>1):
#                 #coin flip
#                 n = np.random.random()
#                 if(n < 2*s*len(cells)):
#                     n_cancer_during += 1
#                     cancer_age.append(gen)
#                     break
#                 else:
#                     cells[-1] = 0

#     print("Percentage = ", n_cancer_during / (run+1))
#     calculated.append(n_cancer_during / (run+1))
    lscd.append(cell_start*(2+t) - 2)
# #     print("Expected = ", risk[i])
    
    
calculated_from_1 = [0.0012031307842826167,  0.055168485167170586, 1, 1, 0.0, 0.4338985622700201, 0.00013]
calculated_from_all = [0.00048, 0.03021951219512195, 1, 1, 0.0, 0.3238666666666667, 0.0]
for i in range(len(calculated_from_all)):
    if(calculated_from_1[i] == 0):
        calculated_from_1[i] = 1/(n_runs * 10)
        
    if(calculated_from_all[i] == 0):
        calculated_from_all[i] = 1/(n_runs * 10)
        
fig, ax = plt.subplots()
for i in range(len(lscd)):
    color = next(ax._get_lines.prop_cycler)['color']
    ax.loglog(lscd[i], calculated_from_1[i],'o', color = color)
    ax.loglog(lscd[i], calculated_from_all[i],'x', color = color)
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
plt.legend(["Calculated from 1", "Calculate from start", "Risk"])
plt.title("Calculated incidence and observed risk in a number of cancers")
plt.xlabel('Total stem cell divisions')
plt.ylabel('Lifetime risk')

# x_data = lscd
# y_data = risk

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