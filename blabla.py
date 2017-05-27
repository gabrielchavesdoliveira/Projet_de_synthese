import numpy as np
import matplotlib.pyplot as plt

x = np.array([0,1,2,3,4])
y = np.array([0,1,2,3,4])

xx,yy = np.meshgrid(x,y)
xx=xx**2
#print(xx)

M=16
a = ((np.arange(M)+0.5))
#print (a)
b= a.reshape((1, M))
#print (b)


ut = np.zeros((5,5), dtype=float)
ut[:,1:] = 1
ut[1:,:] = 1

print(ut)
