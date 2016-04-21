'''
Created on 20 Apr 2016

@author: mathii
'''
from __future__ import division, print_function
import numpy as np

class hmm2(object):
    '''
    HMM-solving object
    '''

    def __init__(self, pair, obs, pos):
        '''
        Constructor
        '''
        self.pair = pair
        self.obs = obs
        self.pos = pos
        self.n_obs = len(obs)
        self.p = np.mean(self.obs) #Estimate the sharing prob by genome-wide Problem
        self._em = self.emission_matrix(self.p)
        self.n_states = self._em.shape[0]
        self._vit = None
        self._tb = None
        
    def emission_matrix(self, p):
        """
        indexed by IBD_state,shared
        """
        matrix = np.zeros((2,2), dtype=np.float)
        matrix[0,0] = 1-p
        matrix[0,1] = p
        matrix[1,0] = 0.75*(1-p)
        matrix[1,1] = 0.25*(1+3*p)
        
        return matrix
    
    def transition_matrix(self, r):
        """
        r is the transition probability per base
        """
        return 1-np.exp(-r*np.diff(self.pos))
    
    def traceback(self):
        if not hasattr(self, "_tb"):
            self.Viterbi()
            
        trace = np.zeros((self.n_obs), dtype=np.int)
        best = np.argmax(self._vit[:,self.n_obs-1])
        trace[self.n_obs-1]=best
        for i in xrange(1, self.n_obs):
            best=self._tb[trace[self.n_obs-i],self.n_obs-i]
            trace[self.n_obs-1-i]=best
        
        return trace    
    
    def get_chunks(self):
        """
        return a list of the chunks in each state
        """
        trace=self.traceback()
        
        chunks =  [[] for _ in range(self.n_states)]
        
        trace_enum=enumerate(trace)
        i,last_st=trace_enum.next()
        last_st_start=self.pos[i]
        for i,st in trace_enum:
            if st != last_st:
                chunks[last_st].append(self.pos[i]-last_st_start)
                last_st = st
                last_st_start=self.pos[i]
        
        chunks[last_st].append(self.pos[i]-last_st_start)
        return chunks
                        
    def Viterbi(self):
        """
        Run the Viterbi algorithm
        """
        trans = self.transition_matrix(1e-8)
        vit = np.zeros((2, len(self.obs)), dtype=np.float) #Viterbi matrix 
        tb = np.zeros((2,len(self.obs)), dtype=np.int) #Traceback matrix
        
        vit[:,0] = self._em[:,self.obs[0]] 

        for i in xrange(1,len(self.obs)):
            for jj in range(self.n_states):
                col = trans[i-1]+np.zeros((self.n_states))
                col[jj] = 1-col[jj]
                col = col*vit[:,i-1]*self._em[:,self.obs[i]] 
                best = np.argmax(col)
                vit[jj,i] = col[best]
                tb[jj,i] = best
                
            vit[:,i] = vit[:,i]/sum(vit[:,i])
                
            self._vit = vit
            self._tb = tb
        