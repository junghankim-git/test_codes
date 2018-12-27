import os
import sys
import numpy as np
import math
import matplotlib.pyplot as plt


class beta_dist:
   def __init__(self,nsize,alpha,beta):
      self.nsize   = nsize
      self.alpha   = alpha
      self.beta    = beta
      self.PDF     = [0.0 for i in range(nsize)]
      self.CDF     = [0.0 for i in range(nsize)]
      self.x       = [float(i)/float(nsize-1) for i in range(nsize)]
      self.get_pdf()


   # Probability density function
   def get_pdf(self):
      alpha = self.alpha
      beta  = self.beta
      const = self.beta_func(alpha, beta)

      for i in range(self.nsize):
         self.PDF[i] = 1.0/const*self.x[i]**(alpha-1)*(1.0-self.x[i])**(beta-1)

      return self.x, self.PDF


   # Cumulative density function
   def get_cdf(self):
      alpha = self.alpha
      beta  = self.beta
      const = self.beta_func(alpha, beta)

      for i in range(self.nsize):
         self.PDF[i] = 1.0/const*self.x[i]**(alpha-1)*(1.0-self.x[i])**(beta-1)

      return self.x, self.PDF


   def draw_pdf(self):
      pltfig = plt.figure(figsize=(8,6))
      pltfig.suptitle('Beta Distribution',fontsize=18)
      axis = pltfig.add_subplot(1,1,1)
      axis.clear()
      axis.set_title('Beta Distribution')
      axis.set_xlabel('x')
      axis.set_ylabel('PDF')
      axis.set_xlim(0.0,1.0)
      axis.set_ylim(0.0,2.5)
      axis.grid(True)
      axis.plot(self.x,self.PDF,'r-')
      plt.show()


   # private
   def beta_func(self, alpha_in, beta_in):
      return float(math.factorial(alpha_in-1))*float(math.factorial(beta_in-1))/float(math.factorial(alpha_in+beta_in-1))

   # private: must be developed
   def incomplete_beta_func(self, x_in, alpha_in, beta_in):
      print math.factorial(alpha_in-1), math.factorial(beta_in-1), math.factorial(alpha_in+beta_in-1)
      return float(math.factorial(alpha_in-1))*float(math.factorial(beta_in-1))/float(math.factorial(alpha_in+beta_in-1))
