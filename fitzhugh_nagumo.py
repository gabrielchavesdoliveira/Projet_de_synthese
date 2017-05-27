

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
       if (param_name = 'nerves'):
          self.a =
          self.b =
          self.e =
       elif (self.param_name = 'oscillator') :
          self.a =
          self.b =
          self.e =

       else :
          self.a =
          self.b =
          self.e =

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
       if (param_name = 'nerves'):
          self.a =
          self.b =
          self.e =
       elif (self.param_name = 'oscillator') :
          self.a =
          self.b =
          self.e =

       else :
          self.a =
          self.b =
          self.e =

       self.width =  width
       self.height = height
       self.h = d/self.width
       self.Du = ?
       self.Dv = ?
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
      
