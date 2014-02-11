import sys
import numpy as np

filename = sys.argv[1]
f = open(filename, "r")
lines = f.readlines()
f.close()
ls_time = []
for line in lines:
    ls_line = line.split(' ')
    if len(ls_line) != 4:
        continue
    else:
        time_elips = ls_line[2]
        #print time_elips
        ls_time.append(time_elips)
print len(ls_time)

mat = np.asarray(ls_time)
x = np.array_split(mat, 5)
x = np.asmatrix(x)
print x.shape

m1 = x[0].T
m2 = x[1].T
m3 = x[2].T
m4 = x[3].T
m5 = x[4].T

res = np.vstack((m1,m2,m3,m4,m5))
print res

np.savetxt(filename+'.csv',res, delimiter=',', fmt="%s")
