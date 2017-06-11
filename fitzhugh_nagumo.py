

import numpy as np
import scipy.signal
import scipy.ndimage
import numpy as np
import time
import sys



class SpectralModel :

   
   def __init__ ( self, param_name, width, height, d=1., dt= 0.1, mode = 'ETDFD'):
       self.param_name = param_name
       self.width =  width
       self.height = height
       self.h = d/self.width
       self.d = d
       self.Du = 0.00001/self.h**2
       self.Dv = 0.000001/self.h**2
       self.dt = dt
       self.noise = 0.2
       if (param_name == 'refractory'):
          self.a = 0.1
          self.b = 3000
          self.e = 0.000003
          self.I = 0
       elif (self.param_name == 'propagation') :
          self.a = 0.1
          self.b = 1000
          self.e = 0.000003
          dN = min(self.height, self.width)/2
       
          self.I = np.zeros((self.height, self.width), dtype=float)
          self.I[:,:] = 0
          #self.I[32,64]=1
          self.I[(self.height/2 - dN/2): (self.height/2+dN/2+1), (self.width/2 - dN/2) : (self.width/2+dN/2+1)] = 1
          

       elif (self.param_name == 'oscillation' ):
          self.a = -0.02
          self.b = 1000
          self.e = 1e-5
          self.I = 0
         

       

       self.tf_ut = np.zeros((self.height, self.width), dtype = complex)
       self.tf_vt = np.zeros((self.height, self.width), dtype = complex)

       self.mode = mode      
       if(not self.mode in ['ETDFD','ETDRK4']):
          print("The numerical scheme you mentioned is not implemented")
          raise Exception("Unknown numerical scheme, must be ETDFD or ETDRK4")

       # Precompute various ETDRK4 scalar quantities
       k1, k2 = np.meshgrid(np.arange(self.width).astype(float), np.arange(self.height).astype(float))
       k1[:,self.width/2+1:] -= self.width
       k2[self.height/2+1:,:] -= self.height

       k1 *= 2.0 * np.pi / self.width
       k2 *= 2.0 * np.pi / self.height
       self.Lu = -(self.Du*(k1**2 + k2**2)) - self.a
       self.Lv = -(self.Dv*(k1**2 + k2**2)) - self.e
       self.E2u = np.exp(self.dt * self.Lu/2.)
       self.Eu = self.E2u ** 2
   
       self.E2v = np.exp(self.dt * self.Lv/2.)
       self.Ev = self.E2v ** 2

       M = 16 # Nb of points for complex means
       r = (np.exp(1j * np.pi * (np.arange(M)+0.5)/M)).reshape((1, M))
       # Generate the points along the unit circle contour over which to compute the mean
       LRu = (self.dt * self.Lu).reshape((self.width*self.height, 1)) + r
       LRv = (self.dt * self.Lv).reshape((self.width*self.height, 1)) + r
       
       # The matrix for integrating the constant F term in the equation of u
       self.F2u = -np.real(np.mean(self.dt * (1. - np.exp(LRu/2.))/LRu, axis=1).reshape((self.height, self.width)))
       self.F2u[1:,:] = 0
       self.F2u[:,1:] = 0
       self.Fu = -np.real(np.mean(self.dt * (1. - np.exp(LRu))/LRu, axis=1).reshape((self.height, self.width)))
       self.Fu[1:,:] = 0
       self.Fu[:,1:] = 0
       if(mode == 'ETDFD'):
           self.FNu = -np.real(np.mean(self.dt * (1. - np.exp(LRu))/LRu, axis=1).reshape((self.height, self.width)))
           self.FNv = -np.real(np.mean(self.dt * (1. - np.exp(LRv))/LRv, axis=1).reshape((self.height, self.width)))
       elif(mode == 'ETDRK4'): 
           LRu_2 = LRu**2.
           LRu_3 = LRu**3.
           self.Qu = np.real(np.mean(self.dt * (np.exp(LRu/2.) - 1.) / LRu, axis=1).reshape((self.height, self.width)))
           self.f1u = np.real(np.mean(self.dt * (-4. - LRu + np.exp(LRu) * (4. - 3 * LRu + LRu_2)) / LRu_3, axis=1).reshape((self.height, self.width)))
           self.f2u = np.real(np.mean(self.dt * 2. * (2. + LRu + np.exp(LRu) * (-2. + LRu)) / LRu_3, axis=1).reshape((self.height, self.width)))
           self.f3u = np.real(np.mean(self.dt * (-4. - 3 * LRu - LRu_2 + np.exp(LRu) * (4. - LRu)) / LRu_3, axis=1).reshape((self.height, self.width)))

           LRv_2 = LRv**2.
           LRv_3 = LRv**3.
           self.Qv = np.real(np.mean(self.dt * (np.exp(LRv/2.) - 1.) / LRv, axis=1).reshape((self.height, self.width)))
           self.f1v = np.real(np.mean(self.dt * (-4. - LRv + np.exp(LRv) * (4. - 3 * LRv + LRv_2)) / LRv_3, axis=1).reshape((self.height, self.width)))
           self.f2v = np.real(np.mean(self.dt * 2. * (2. + LRv + np.exp(LRv) * (-2. + LRv)) / LRv_3, axis=1).reshape((self.height, self.width)))
           self.f3v = np.real(np.mean(self.dt * (-4. - 3 * LRv - LRv_2 + np.exp(LRv) * (4. - LRv)) / LRv_3, axis=1).reshape((self.height, self.width)))

   def init(self):
      dN = min(self.height, self.width)/4


      #Carre blanc centre
      #ut = np.zeros((self.height, self.width), dtype=float)
      #ut[:,:] = 0
      #ut[(self.height/2 - dN/2): (self.height/2+dN/2+1), (self.width/2 - dN/2) : (self.width/2+dN/2+1)] = 1
      #ut += self.noise * (2 * np.random.random((self.height, self.width)) - 1)
      #ut[ut <= 0] = 0

      #Bande  gauche
      ut = np.zeros((self.height, self.width), dtype=float)
      ut[:,:] = 0
      ut[1:128,1] = 1

      #Point haut gauche
      #ut = np.zeros((self.height, self.width), dtype=float)
      #ut[:,:] = 0
      #ut[20:22,20:22] = 1

      #Contour
      #ut = np.zeros((self.height, self.width), dtype=float)
      #ut[:,:] = 0
      #ut[0,:] = 1
      # ut[127,:] = 1
      # ut[:,0] = 1
      # ut[:,255] = 1
      
      vt = np.zeros((self.height, self.width), dtype=float)
      vt[:,:] = 0
      #vt[(self.height/2 - dN/2): (self.height/2+dN/2+1), (self.width/2 - dN/2) : (self.width/2+dN/2+1)] = 0.25
      #vt += self.noise * (2 * np.random.random((self.height, self.width)) - 1)
      #vt[vt <= 0] = 0
        
      self.tf_ut = np.fft.fft2(ut)
      self.tf_vt = np.fft.fft2(vt)
       
   def get_ut(self):
       return np.real(np.fft.ifft2(self.tf_ut))
    
   def get_vt(self):
       return np.real(np.fft.ifft2(self.tf_vt))

   def compute_Nu_Nv(self,tf_u,tf_v):
      Nu = -np.fft.fft2(np.fft.ifft2(tf_u).real**3) + (self.a +1)*np.fft.fft2(np.fft.ifft2(tf_u).real**2) - tf_v
      Nv = self.e*self.b*tf_u
      return Nu,Nv
   # Erase the reactant in a box
   def erase_reactant(self, center, radius):
       vt =np.real(np.fft.ifft2(self.tf_vt))
       vt[(center[0]-radius):(center[0]+radius), (center[1]-radius):(center[1]+radius)] = 0
       self.tf_vt = np.fft.fft2(vt)      
   # Mask the reactant,
   # mask.shape = self.height, self.width
   # mask.dtype = float
   # mask_ij in [0, 1]
   def mask_reactant(self, mask):
       #vt =np.real(np.fft.ifft2(self.tf_vt))
       #vt = vt * mask
       #self.tf_vt = np.fft.fft2(vt)


       vt = np.real(np.fft.ifft2(self.tf_vt))
       vt[mask >= 0.5] = 1.0
       self.tf_vt = np.fft.fft2(vt)
        
       ut = np.real(np.fft.ifft2(self.tf_ut))
       ut[mask >= 0.5] = 0.0
       self.tf_ut = np.fft.fft2(ut)


   def step(self):
       if(self.mode == 'ETDFD'):
           Nu, Nv = self.compute_Nu_Nv(self.tf_ut, self.tf_vt)
           self.tf_ut = self.Eu * self.tf_ut + self.Fu * self.I * self.width * self.height + self.FNu * Nu #+ np.fft.fft2(0.025*np.random.random((self.height, self.width))) 
           self.tf_vt = self.Ev * self.tf_vt + self.FNv * Nv  
       elif(self.mode == 'ETDRK4'):
           Nu, Nv = self.compute_Nu_Nv(self.tf_ut, self.tf_vt)
           au = self.E2u * self.tf_ut + self.F2u * self.I *self.width*self.height+ self.Qu * Nu
           av = self.E2v * self.tf_vt + self.Qv * Nv
           Nau, Nav = self.compute_Nu_Nv(au, av)
           bu = self.E2u * self.tf_ut + self.F2u * self.I * self.width * self.height + self.Qu * Nau
           bv = self.E2v * self.tf_vt + self.Qv * Nav
           Nbu, Nbv = self.compute_Nu_Nv(bu, bv)
           cu = self.E2u * au + self.F2u * self.I * self.width * self.height + self.Qu * (2. * Nbu - Nu)
           cv = self.E2v * av + self.Qv * (2. * Nbv - Nv)
           Ncu, Ncv = self.compute_Nu_Nv(cu, cv)
           
           self.tf_ut = self.Eu * self.tf_ut + self.Fu * self.I * self.width * self.height + self.f1u * Nu + self.f2u * (Nau + Nbu) + self.f3u * Ncu #+ np.fft.fft2(np.random.random((self.height, self.width)))
           self.tf_vt = self.Ev * self.tf_vt + self.f1v * Nv + self.f2v * (Nav + Nbv) + self.f3v * Ncv 


        
class Model :
   def __init__(self, param_name, width, height, mode, d=1., dt=0.1):
      self.param_name = param_name
      self.width =  width
      self.height = height
      self.h = d/self.width
      self.Du = 0.00001/self.h**2
      self.Dv = 0.000001/self.h**2
      self.dt = dt
      self.noise = 0.2
      if (param_name == 'refractory'):
          self.a = 0.1
          self.b = 3000
          self.e = 0.000003
          self.I = 0
      elif (self.param_name == 'propagation') :
         self.a = 0.1
         self.b = 1000
         self.e = 0.000003
         dN = min(self.height, self.width)/2
       
         self.I = np.zeros((self.height, self.width), dtype=float)
         self.I[:,:] = 0
         #self.I[32,64]=1
         self.I[(self.height/2 - dN/2): (self.height/2+dN/2+1), (self.width/2 - dN/2) : (self.width/2+dN/2+1)] = 1
          

      elif (self.param_name == 'oscillation' ):
         self.a = -0.02
         self.b = 1000
         self.e = 0.00001
         self.I = 0

      
      
      self.ut_1 = np.zeros((self.height, self.width), dtype=float)
      self.vt_1 = np.zeros((self.height, self.width), dtype=float)
      self.ut = np.zeros((self.height, self.width), dtype=float)
      self.vt = np.zeros((self.height, self.width), dtype=float)

      self.mode = mode
      if (self.mode == 0):
         self.stencil = np.zeros((self.height, self.width))
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
      
      #Carre blanc centre
      #ut = np.zeros((self.height, self.width), dtype=float)
      #ut[:,:] = 0
      #ut[(self.height/2 - dN/2): (self.height/2+dN/2+1), (self.width/2 - dN/2) : (self.width/2+dN/2+1)] = 1
      #ut += self.noise * (2 * np.random.random((self.height, self.width)) - 1)
      #ut[ut < 0] = 0

      #Bande  gauche
      ut = np.zeros((self.height, self.width), dtype=float)
      ut[:,:] = 0
      ut[1:128,1] = 1

      #Point haut gauche
      #ut = np.zeros((self.height, self.width), dtype=float)
      #ut[:,:] = 0
      #ut[20,20] = 1

      #Contour
      #ut = np.zeros((self.height, self.width), dtype=float)
      #ut[:,:] = 0
      #ut[0,:] = 1
      #ut[127,:] = 1
      #ut[:,0] = 1
      #ut[:,255] = 1
      
      vt = np.zeros((self.height, self.width), dtype=float)
      vt[:,:] = 0
      #vt[(self.height/2 - dN/2): (self.height/2+dN/2+1), (self.width/2 - dN/2) : (self.width/2+dN/2+1)] = 0.25
      #vt += self.noise * (2 * np.random.random((self.height, self.width)) - 1)
      #vt[vt <= 0] = 0

      
      #self.vt_1[:,:] = 0
      #self.vt_1[(self.height/2 - dN/2): (self.height/2+dN/2+1), (self.width/2 - dN/2) : (self.width/2+dN/2+1)] = 0
      #self.vt_1 += self.noise * (2 * np.random.random((self.height, self.width)) - 1)
      #self.vt_1[self.vt_1 <= 0] = 0

      self.vt_1[:,:] = vt[:,:]
      self.ut_1[:,:] = ut[:,:]
      

   def laplacian(self, x):
      if(self.mode == 0):
	 return np.fft.irfft2(np.fft.rfft2(x)*self.fft_mask)
      elif(self.mode == 1):
	 return scipy.ndimage.convolve(x, self.stencil, mode='wrap')
      elif(self.mode == 2):
	 return scipy.ndimage.laplace(x, mode='wrap')

   def get_ut(self):
      return self.ut
   def get_vt(self):
      return self.vt

   def erase_reactant(self, center, radius):
      pass

   def step(self):

      lu = self.laplacian(self.ut_1)
      lv = self.laplacian(self.vt_1)
      self.ut[:,:] = self.ut_1 + self.dt*(self.Du * lu + self.ut_1*(self.a - self.ut_1)*(self.ut_1 - 1) - self.vt_1 + self.I)
      #print(self.ut)
      self.vt[:,:] = self.vt_1 + self.dt*(self.Dv * lv + self.e * (self.b*self.ut_1 - self.vt_1) + self.I)

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
   pattern = 'refractory'      
   d = 1.
   dt = 1.

   if(mode <= 2):
      model = Model(pattern, width, height, mode=mode)
   elif mode == 3:
      model = None #libgrayscott.GrayScott(pattern, width, height, d, dt)
   elif mode == 4:
      #model = SpectralModel(pattern, height=height, width=width)
      print("on n en est pas la haha")
        
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
