

import numpy as np
import scipy.signal
import scipy.ndimage
import numpy as np
import time
import sys

'''
∂ₜu = Dᵤ∆u + (a-u)(u-1)u - v 
∂ₜv = Dᵥ∆v + e(bu-v)

With U= TF(U) and V = TF(V), the system reads :

∂ₜU[k1,k2] = -Dᵤ ( (2πk₁/d)^2 + (2πk₂/d)^2) U[k₁,k₂] - TF[TF⁻1(U)³] + (a+1)TF(TF⁻1(U)²)  - aU[k1,k2] - V[k1,k2]

∂ₜ V [k₁,k₂] = -Dᵥ( (2πk₁/d)² + (2πk₂/d)²)V[k₁,k₂] + e( bU[k1,k2] - V[k1,k2])

Decomposition of the equation, according to the linear and non-linear parts :

∂ₜ U[k₁, k₂] = Lᵤ U[k₁,k₂] + Nᵤ(U[k₁,k₂], V[k₁,k₂])
∂ₜ V[k₁, k₂] = Lᵥ V[k₁,k₂] + eb U[k₁,k₂])

with :   Lu =  -Dᵤ ( (2πk₁/d)² + (2πk₂/d)²) - a
         Lv =  -Dᵥ ( (2πk₁/d)² + (2πk₂/d)²) - e
         Nu(U[k1,k2], V[k1,k2]) = - TF[TF⁻1(U)³] + (a+1)TF(TF⁻1(U)²) - V[k1,k2]


'''

class SpectralModel :

   ''' Modes ? '''
   def __init__ ( self, param_name, width, height, d=1., dt= 0.1, mode = 'ETDFD'):
       self.param_name = param_name
       if (param_name = 'refractory'):
          self.a = 0.1
          self.b = 10000
          self.e = 0.0000005
          self.I = 0
       elif (self.param_name = 'propagation') :
          self.a = 0.1
          self.b = 10000
          self.e = 0.0000005
          self.I = 1

       elif (self.param_name = 'oscillation' ):
          self.a = -0.1
          self.b = 10000
          self.e = 0.0000005
          self.I = 0
         

       self.width =  width
       self.height = height
       self.h = ?
       self.d = d
       self.Du = ?
       self.Dv = ?
       self.dt = dt
       self.noise = 0.2

       self.tf_ut = np.zeros((self.height, self.width), dtype = complex)
       self.tf_vt = np.zeros((self.height, self.width), dtype = complex)

       k1
          
          
class Model :
   def __init__(self, param_name, width, height, mode, d=1., dt=0.1):
      self.param_name = param_name
       if (param_name = 'refractory'):
          self.a = 0.1
          self.b = 10000
          self.e = 0.0000005
          self.I = 0
       elif (self.param_name = 'propagation') :
          self.a = 0.1
          self.b = 10000
          self.e = 0.0000005
          self.I = 1

       elif (self.param_name = 'oscillation' ):
          self.a = -0.1
          self.b = 10000
          self.e = 0.0000005
          self.I = 0

       self.width =  width
       self.height = height
       self.h = d/self.width
       self.Du = 0.1/h**2
       self.Dv = 0.1/h**2
       self.dt = dt
       self.noise = 0.2
   
       self.ut_1 = np.zeros((self.height, self.width), dtype=float)
       self.vt_1 = np.zeros((self.height, self.width), dtype=float)
       self.ut = np.zeros((self.height, self.width), dtype=float)
       self.vt = np.zeros((self.height, self.width), dtype=float)

       self.mode = mode
       if (self.mode == 0):
          self.stencil[0,0] = -4
	  self.stencil[0,1] = 1
	  self.stencil[0,-1] = 1
	  self.stencil[1,0] = 1
	  self.stencil[-1,0] = 1
	  self.fft_mask = np.fft.rfft2(self.stencil)
       elif(self.mode == 1):
          self.stencil = np.array([[0, 1., 0], [1., -4., 1.], [0, 1., 0]], dtype=float)

   def init(self):
      dN = min(self.width, self.height)/4
      self.ut_1[:,:] = 1
      self.ut_1[:,:] = 1
      self.ut_1[(self.height/2 - dN/2): (self.height/2+dN/2+1), (self.width/2 - dN/2) : (self.width/2+dN/2+1)] = 0.5
      self.ut_1 += self.noise * (2 * np.random.random((self.height, self.width)) - 1)
      self.ut_1[self.ut_1 <= 0] = 0
                
      self.vt_1[:,:] = 0
      self.vt_1[(self.height/2 - dN/2): (self.height/2+dN/2+1), (self.width/2 - dN/2) : (self.width/2+dN/2+1)] = 0.25
      self.vt_1 += self.noise * (2 * np.random.random((self.height, self.width)) - 1)
      self.vt_1[self.vt_1 <= 0] = 0

      self.vt[:,:] = self.vt_1[:,:]
      self.ut[:,:] = self.ut_1[:,:]
      

   def laplacian(self, x):
      if(self.mode == 0):
	 return np.fft.irfft2(np.fft.rfft2(x)*self.fft_mask)
      elif(self.mode == 1):
	 return scipy.ndimage.convolve(x, self.stencil, mode='wrap')
      elif(self.mode == 2):
	 return scipy.ndimage.laplace(x, mode='wrap')

   def get_ut(self):
      return self.ut

   def erase_reactant(self, center, radius):
      pass

   def step(self):

      lu = self.laplacian(ut_1)
      lv = self.laplacian(vt_1)
      self.ut[:,:] = self.ut_1 + self.dt*(self.Du * lu + self.ut_1*(self.a - self.ut_1)*(self.ut_1 - 1) - self.vt_1 + self.I)

      self.vt[:,:] = self.vt_1 + self.dt*(self.Dv * lv + self.e * (self.b*self.ut_1 - vt_1) + self.I)

      self.ut_1, self.vt_1 = self.ut, self.vt



if (__name__ == '__main__'):
   if(len(sys.argv) <= 1):
      print("Usage : %s mode "% sys.argv[0])
      print("With mode : ")
      print("   0 : spatial model with FFT convolution in python, forward euler") # 100 fps
      print("   1 : spatial model with ndimage.convolve in python, forward euler") # 165 fps
      print("   2 : spatial model with ndimage.laplace in python, forward euler") # 150 fps
      print("   3 : spatial model with fast laplacian in C++, forward euler LOL") # 400 fps
      print("   4 : spectral model in python using ETDRK4")
      sys.exit(-1)

   mode = int(sys.argv[1])
    
   height = 100
   width = 100
   pattern = 'oscillation'      
   d = 1.
   dt = 1.

   if(mode <= 2):
      model = Model(pattern, width, height, mode=mode)
   elif mode == 3:
      model = None #libgrayscott.GrayScott(pattern, width, height, d, dt)
   elif mode == 4:
      #model = SpectralModel(pattern, height=height, width=width)
      print("   LOL : on n'en est pas là haha")
        
   model.init()
    
   epoch = 0
   t0 = time.time()
    
   while True:
      model.step()
      epoch += 1
      if(epoch % 500 == 0):
         t1 = time.time()
         print("FPS : %f f/s" % (500 / (t1 - t0)))
         t0 = t1
