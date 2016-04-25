#!/usr/local/bin/python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import *
import sys
import os
import pdb
import yaml

'''
Code for running multiple 2D IsingModels for Statistical Mechanics

Inputs: Parameter yaml file

Ouputs: Graphs/txt files with various thermodynamic properties

'''

class IsingModelEnsemble:
    def __init__(self, param_path):
        self.param_path =param_path

        with open(param_path) as f:
            self.params = yaml.load(f)
        self.T_list = self.params["T_list"]
        # self.H_list = self.params["H_list"]
        # self.L_list = self.params["L_list"]
        self.graph_dict = {}

    def RunSystem(self, model):
        # if self.params["param1"] == "T": model.Equilibriate()

        alg = self.params["sys_params"]["sim_type"]
        if alg == "Heat Bath":
            model.RunHeatBath()
        elif alg == "Metropolis":
            model.RunMetropolis()
        else:
            print "Did you check your spelling of the algorithms?"
        return

    def GraphvsParam(self):
        print "Not implemented yet"

    #Graphs changes in a thermodynamic variable(s) (graph_type) with respect to time(step)
    #Multiple values of a parameter (param2) can be given to get multiple graphs
    def Graph(self):
        var = self.params["param2"]
        orig = float(self.params["sys_params"][var])
        Torig = float(self.params["sys_params"]["T"])

        #Determine whether this is a multivalue run or not
        if var+"_list" not in self.params.keys():
            var_list = [orig]
            # colors = [(0,0,1)]
            colors = ['b']
        else:
            var_list = self.params[var+"_list"]
            colors = mpl.cm.rainbow(np.linspace(0,1,len(var_list)))

        #Number of graph types
        n_graph = len(self.params["graph_type"])

        #Make number of subfigs required
        self.fig, axarr = plt.subplots(n_graph, sharex=True)

        #list of different Ising models
        IM_list = []

        #Run Ising Models with different values of param2
        for v in var_list:
            self.params["sys_params"][var] = v 
            #Make a 2D array if you are running it with temperature
            if self.params["param1"] == "T":
                IM_Tlist = []
                for T in self.T_list:
                    self.params["sys_params"]["T"] = T
                    IM_Tlist.append(IsingModel(self.params))
                    self.RunSystem(IM_Tlist[-1])
                IM_list.append(IM_Tlist)
            else:
                IM_list.append(IsingModel(self.params))
                self.RunSystem(IM_list[-1])
        
        #loop through subfigs for each graph type wanted
        for ax, gt in zip(axarr, self.params["graph_type"]):
            #Make plots on each subfig 
            if self.params['param1'] == 'step':
                self.GraphVsTime( ax, IM_list, colors, gt, var_list, var)
            elif self.params['param1'] == 'T':
                self.GraphVsTemp( ax, IM_list, colors, gt, var_list, var)
            else:
                print "param1 must be either T or step"
        
        self.params["sys_params"][var] = orig
        self.params["sys_params"]["T"] = Torig
        return
        
    def GraphVsTime(self, ax, IM_list, colors, gt, var_list, var):
        for m, color, v in zip(IM_list, colors, var_list):
            y_array = self.GetValArray(m, gt)
            ax.plot(m.step_array, y_array, color=color, label = '{0} = {1}'.format(var, v)) 

        ax.set_title("{} vs Time".format(gt))
        ax.set_xlabel("steps")
        ax.set_ylabel("{} / {}_max ".format(gt,gt))
        ax.set_ylim([-1.1,1.1])
        ax.legend(loc='center left', bbox_to_anchor=(1,.5))
        return

    def GraphVsTemp(self, ax, IM_list, colors, gt, var_list, var):
        for IM_Tlist, color, v in zip(IM_list, colors, var_list):
            y_array=[]
            for m in IM_Tlist:
                y_array.append(self.GetValAvg(m, gt))  

            ax.plot(self.T_list, y_array, color=color, marker='o', label = '{0} = {1}'.format(var, v)) 
        ax.set_title("{} vs Temperature".format(gt))
        ax.set_xlabel("T")
        ax.set_ylabel("{} / {}_max ".format(gt,gt))
        ax.legend(loc='center left', bbox_to_anchor=(1,.5))
        return


    def ShowPlots(self):
        #Make title
        title = ""
        for key, param in self.params["sys_params"].items():
            if key != self.params["param1"] and key != self.params["param2"]:
                title += "{} = {}, ".format(key, param)
            
        self.fig.suptitle(title)
        self.fig.tight_layout()
        self.fig.subplots_adjust(right=.75, top=.87)
        plt.show()
        return

    def SavePlots(self):
        self.fig.savefig(os.path.splitext("Graphs/"+self.param_path)[0]+".png")


    def GetValArray(self, model, Val):
        if Val == 'E':
            return model.E_array
        elif Val == 'M':
            return model.M_array
        elif Val == '|M|':
            return np.absolute(model.M_array)
        elif Val == 'C':
            return model.C_array
        elif Val == 'X':
            return model.S_array

    def GetValAvg(self, model, Val):
        if Val == 'E':
            return model.GetAvgEnergy()
        elif Val == 'M':
            return model.GetAvgMag()
        elif Val == '|M|':
            return model.GetAvgAbsMag()
        elif Val == 'C':
            return model.GetAvgC()
        elif Val == 'X':
            return model.GetAvgX()

class IsingModel:
    def __init__(self, param_dict=None, param_path=None):

        if param_dict != None:
            self.params = param_dict

        elif param_path != None:
            with open(param_path) as f:
                self.params = yaml.load(f)

        else:
            print "Need either a list or a dictionary of parameters"

        self.seed = self.params["seed"]
        np.random.seed(self.seed)
        
        #Get variables used in algorithm
        self.beta = 1/float(self.params["sys_params"]["T"])
        self.H = float(self.params["sys_params"]["H"])
        self.J = float(self.params["sys_params"]["J"])
        self.L = int(self.params["sys_params"]["L"])
        self.V = float(self.L*self.L)
        self.step_array = range(self.params["sys_params"]["n_steps"]) #Array of all steps

        #Create empty lists for data output
        self.M_array = [] #Magnetization at every step
        self.E_array = [] #Energy at every step
        self.C_array = [] #Specific Heat at every step
        self.X_array = [] #Susceptibility at every step
        
        #Initialize system with random orientations of spin
        self.system = self.GenerateRandomMatrix(self.L)

    def Equilibriate(self):
        for i_step in range(int(self.V*100*self.beta)):
           i, j = self.ChooseRandomParticle() 
           self.ThermalizeParticleMet(i, j)
        return

    def RunHeatBath(self):
        for i_step in self.step_array:
           i, j = self.ChooseRandomParticle() 
           self.ThermalizeParticleHB(i, j)
           self.M_array.append(self.GetMagnetization())
           self.E_array.append(self.GetTotalEnergy())

        print self.GetAvgEnergy()
        print self.GetAvgMag()
        return


    def RunMetropolis(self):
        for i_step in self.step_array:
           i, j = self.ChooseRandomParticle() 
           self.ThermalizeParticleMet(i, j)
           self.M_array.append(self.GetMagnetization())
           self.E_array.append(self.GetTotalEnergy())

        print self.GetAvgEnergy()
        print self.GetAvgMag()
        return

    def RunCluster(self):
        print "Not implemented yet"


    def ThermalizeParticleHB(self, i, j):
        E_up = self.GetParticleEnergy(i,j,1)
        E_down = -E_up
        # E_down = self.GetParticleEnergy(i,j,-1)

        P = self.ProbUp(E_up, E_down)
        r = np.random.ranf()

        if (P > r):
            self.system[i,j] = 1
        else:
            self.system[i,j] = -1

        return

    def ThermalizeParticleMet(self, i, j):
        p = self.ProbFlip(i, j)

        if ( p == 1 or p > np.random.ranf()):
            self.system[i,j] *= -1

        return

    def GenerateRandomMatrix(self, N):
        M = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                M[i][j] = np.random.choice([-1,1])
        return M

    def ChooseRandomParticle(self):
        i = np.random.randint(0, self.L)
        j = np.random.randint(0, self.L)
        return i, j

    #Probability calculating algorithms
    def ProbUp(self, E_up, E_down):
        P = np.exp(-self.beta*E_up)/(np.exp(-self.beta*E_up) + np.exp(-self.beta*E_down))
        return P

    def ProbFlip(self, i, j):
        E_c = self.GetParticleEnergy(i, j)
        E_f = -E_c
        if (E_f-E_c) < 0:
            return 1
        else:
            return np.exp(-(E_f-E_c)*self.beta)
        

    # We are going to assume periodic boundary conditions
    def GetParticleEnergy(self, i, j, s=0):
        M = self.system
        if s == 0:
            s = M[i,j] 
        
        i_u = i-1
        i_d = i+1
        j_l = j-1
        j_r = j+1

        if i == 0 :
            i_u = self.L-1
        elif i == self.L-1 :
            i_d = 0
        if j == 0:
            j_l = self.L-1
        elif j == self.L-1 :
            j_r = 0

        E = -s*self.J*(M[i_u,j] + M[i_d,j] + M[i,j_l] + M[i,j_r]) - self.H*s

        return E

    def PrintSystem(self):
        print self.system

    def GetMagnetization(self):
        self.M = np.sum(self.system)
        return self.M/self.V

    def GetTotalEnergy(self):
        E = 0.0
        for i in range(self.L):
            for j in range(self.L):
                E += self.GetParticleEnergy(i, j)

        return E/(4*self.V)
    
    def GetAvgEnergy(self):
        return np.mean(self.E_array)
    
    def GetAvgAbsMag(self):
        return np.mean(np.absolute(self.M_array))

    def GetAvgMag(self):
        return np.mean(self.M_array)

    def GetAvgC(self):
        return np.mean(self.C_array)

    def GetAvgX(self):
        return np.mean(self.X_array)


if __name__ == "__main__":
    ens = IsingModelEnsemble(sys.argv[1])
    ens.Graph()
    ens.ShowPlots()
    ens.SavePlots()







