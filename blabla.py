import numpy as np
import matplotlib.pyplot as plt
import cv2


'''
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


self.stencil[0,0] = -4
self.stencil[0,1] = 1
self.stencil[0,-1] = 1
self.stencil[1,0] = 1
self.stencil[-1,0] = 1




try:
    fullscreen_flag = cv2.WINDOW_FULLSCREEN
    normal_flag = cv2.WINDOW_NORMAL
except:
    fullscreen_flag = cv2.cv.CV_WINDOW_FULLSCREEN
    normal_flag = cv2.cv.CV_WINDOW_NORMAL

cv2.namedWindow('u', cv2.WND_PROP_FULLSCREEN)
cv2.setWindowProperty("u", cv2.WND_PROP_FULLSCREEN, normal_flag)

a = np.zeros((3,3),dtype=float)
a[0,0] = 0.2
a[0,1] = 0.4
a[0,2] = 0.7
a[1,0] = 0.2
a[1,1] = 0.4
a[1,2] = 0.7
a[2,0] = 0.2
a[2,1] = 0.4
a[2,2] = 0.7

cv2.imshow('test',a)
'''



a = np.zeros((3,3),dtype = float)
b = cv2.imread(a, cv2.IMREAD_GRAYSCALE)
#essayons une modif pour voir hahahaha

#espece de connard
#pffffffffffffffffffffffffffffffffffffff
