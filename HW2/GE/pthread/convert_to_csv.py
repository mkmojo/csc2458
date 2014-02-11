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
 #       print time_elips
        ls_time.append(time_elips)

print len(ls_time)

mat = np.asarray(ls_time)
x = np.array_split(mat, 25)
x = np.asmatrix(x)
print x.shape
print x

m1 = x[0:5].T
m2 = x[5:10].T
m3 = x[10:15].T
m4 = x[15:20].T
m5 = x[20:25].T

res = np.vstack((m1,m2,m3,m4,m5))
print res

np.savetxt(filename+'.csv',res, delimiter=',', fmt="%s")
