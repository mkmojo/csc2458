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
x = np.array_split(mat, 30)
x = np.asmatrix(x)
print x.shape
#print x

m1 = x[0:6].T
m2 = x[6:12].T
m3 = x[12:18].T
m4 = x[18:24].T
m5 = x[24:30].T

res = np.vstack((m1,m2,m3,m4,m5))
#print res

np.savetxt(filename+'.csv',res, delimiter=',', fmt="%s")
