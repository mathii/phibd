'''
Created on 20 Apr 2016

@author: mathii
'''
from __future__ import division, print_function
import numpy as np

################################################################################

class hmm2(object):
    '''
    HMM-solving object
    '''

    def __init__(self, pair, obs, pos):
        '''
        Constructor
        '''
<<<<<<< HEAD
        self.pair = pair
        self.obs = obs
        self.pos = pos
        self.n_obs = len(obs)
        self.p = np.mean(self.obs) #Estimate the sharing prob by genome-wide Problem
        self._em = self.emission_matrix(self.p)
        self.n_states = self._em.shape[0]
        self._vit = None
        self._tb = None
=======
        self.pair=pair
        self.obs=obs
        self.pos=pos
        self.n_obs=len(obs)
        self.p=np.mean(self.obs) #Estimate the sharing prob by genome-wide Problem
        self.n_states=self.emission_matrix().shape[0]
>>>>>>> 312d45dbe73029a5a3790915a50b08c08ae578d8
        
                
    def emission_matrix(self):
        """
        indexed by IBD_state,shared
        """
<<<<<<< HEAD
        matrix = np.zeros((2,2), dtype=np.float)
        matrix[0,0] = 1-p
        matrix[0,1] = p
        matrix[1,0] = 0.75*(1-p)
        matrix[1,1] = 0.25*(1+3*p)
=======
        matrix=np.zeros((3,2), dtype=np.float)
        matrix[0,0]=1-self.p
        matrix[0,1]=self.p
        matrix[1,0]=0.75*(1-self.p)
        matrix[1,1]=0.25*(1+3*self.p)
        matrix[2,0]=0.5*(1-self.p)
        matrix[2,1]=0.5*(1+self.p)
>>>>>>> 312d45dbe73029a5a3790915a50b08c08ae578d8
        
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
<<<<<<< HEAD
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
=======
        trans=self.transition_matrix(1e-8)
        vit=np.zeros((self.n_states, len(self.obs)), dtype=np.float) #Viterbi matrix 
        tb=np.zeros((self.n_states,len(self.obs)), dtype=np.short) #Traceback matrix
        em=self.emission_matrix()
        
        vit[:,0]=em[:,self.obs[0]] 

        for i in xrange(1,len(self.obs)):
            for jj in range(self.n_states):
                tp=trans[i-1]
                #col is column of Viterbi matrix
                #Only transition from adjacent states. 
                col=np.zeros((self.n_states))
                col[jj]=1-col[jj]
                if jj>0:
                    col[jj-1]=tp
                if jj<self.n_states-1:
                    col[jj+1]=tp
                col=col*vit[:,i-1]*em[:,self.obs[i]] 
                best=np.argmax(col)
                vit[jj,i]=col[best]
                tb[jj,i]=best
>>>>>>> 312d45dbe73029a5a3790915a50b08c08ae578d8
                
            vit[:,i] = vit[:,i]/sum(vit[:,i])
                
<<<<<<< HEAD
            self._vit = vit
            self._tb = tb
        
=======
            self._vit=vit
            self._tb=tb

#END Class      
################################################################################
  
class multi_hmm(object):
    """  
    Contains multiple hmms
    """
    
    def __init__(self, hmms, tolerance=0.001, max_iters=15):
        self.hmms=hmms
        #Set the p on each of the sub hmms to be equal to the mean
        self.base_p=np.mean([hmm.p for hmm in self.hmms])
        self.p=self.base_p
        for hmm in self.hmms:
            hmm.p=self.p
        self.max_iters=max_iters
        self.tol=tolerance

    def get_chunks(self):
        """
        get the chunk proportions
        """
        return [hmm.get_chunks() for hmm in self.hmms]
        
    def get_proportions(self):
        """
        get the proportion of the genome in each IBD state
        """
        chunks= self.get_chunks()
        lengths=np.array([[sum(x) for x in y] for y in chunks])
        total=np.sum(lengths)
        proportions=np.sum(lengths, axis=0)/total
        return proportions

    def update_p(self):
        """
        Update p when a proportion x of the genome is IBD 
        Also update the sub hmms
        """
        prop=self.get_proportions()
        x=prop[1]
        y=prop[2]
#        self.p=(self.base_p-0.25*x)/(1-0.25*x)
        self.p=(self.base_p-0.25*x-0.5*y)/(1-0.25*x-0.5*y)
        for hmm in self.hmms:
            hmm.p=self.p

       
    def train(self):
        """
        Viterbi training
        """
        old_p=-1
        iter=0
        while abs(old_p-self.p)>self.tol and iter<self.max_iters:
            old_p=self.p
            [hmm.Viterbi() for hmm in self.hmms]
            prop=self.get_proportions()
            #print("%1.4f\t%1.4f\t%1.4f\t%1.4f"%(self.p, prop[1], prop[2], prop[1]+prop[2]))
            self.update_p()
            iter+=1
            
>>>>>>> 312d45dbe73029a5a3790915a50b08c08ae578d8
