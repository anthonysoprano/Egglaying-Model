#!/usr/bin/python

#Simulate egg-laying in hermaphrodites

import math

class Hermaphrodite : 

    # do/dt = k_o-k_c*O
    # egg-laying rage = min(k_f*O,k_s*S)

    maturation_time = 10 # Initial time to mat
    bin_size = 4 # bin size for egg array in hours
    delta_t = .1

    def __init__(self, k_o, k_c, k_f, k_s): 
        self.sperm = 300
        self.oocytes = 0
        self.eggs = []
        self.total_t = 0
        self.k_o = k_o  #oocyte generation rate
        self.k_c = k_c #related to carrying capacity of the animals
        self.k_f = k_f #fertilization rate
        self.k_s = k_s #MSP feedback

    def change_parameters(self, k_o, k_c, k_f, k_s):
        self.sperm = 300
        self.oocytes = 0
        self.eggs = []
        self.total_t = 0
        self.k_o = k_o  #oocyte generation rate
        self.k_c = k_c #related to carrying capacity of the animals
        self.k_f = k_f #fertilization rate
        self.k_s = k_s #MSP feedback

    def increment_time(self) : # must be run by simulate. Holds logic of model

        self.oocytes += self.delta_t*(self.k_o-self.k_c*self.oocytes)
        new_eggs = self.delta_t*min(self.k_f*self.oocytes,self.k_s*self.sperm)
        self.sperm -= new_eggs
        self.oocytes -= new_eggs
        index = int(self.total_t/self.bin_size)
       

        if(len(self.eggs)<=index) : 
            self.eggs.append(new_eggs)
        else :
            self.eggs[index] += new_eggs
        self.total_t += self.delta_t
        print (self.total_t)

    def simulate(self) : #simulates a life of egg-laying

        while(self.sperm>1) :
            self.increment_time()
        x = [i*self.bin_size+self.bin_size/2 for i in range(0,len(self.eggs))]
        for i in range(0,len(self.eggs)) :
            print (str(self.eggs[i]/self.bin_size))

x = Hermaphrodite(12, .02,.0266,.0625)
x.simulate()
#x.change_parameters(10, .02,.0266,.0625)
#x.simulate()
#x.change_parameters(8, .02,.0266,.0625)
#x.simulate()

    
