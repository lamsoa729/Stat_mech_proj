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

def GetStats(arr, bins):
    mean_array = []
    std_array = []
    var_array = []
    n_steps = int(floor(len(arr)/bins))

    for i in range(1,bins):
        mean_array.append(np.mean(arr[i*n_steps:((i+1)*n_steps-1)]))
        var_array.append(np.var(arr[i*n_steps:((i+1)*n_steps-1)]))

    mean = np.mean(mean_array)
    mean_var = np.var(mean_array)
    var = np.mean(var_array)
    var_var = np.var(var_array)

    return mean, mean_var, var, var_var


class IsingModelEnsemble:
    def __init__(self, param_path):
        self.param_path =param_path

        with open(param_path) as f:
            self.params = yaml.load(f)

    def RunSystem(self, model):
        # if self.params["param1"] != 'step' : model.Equilibriate()

        alg = self.params["sys_params"]["sim_type"]
        if alg == "Heat Bath":
            model.RunHeatBath()
        elif alg == "Metropolis":
            model.RunMetropolis()
        else:
            print "Did you check your spelling of the algorithms?"
        return


    #Graphs changes in a thermodynamic variable(s) (graph_type) with respect to time(step)
    #Multiple values of a parameter (param2) can be given to get multiple graphs
    def Graph(self):
        P1 = self.params["param1"]
        P2 = self.params["param2"]
        if P1 != "step":
            P1_list = self.params[P1+"_list"]
            if len(P1_list) == 2:
                temp = np.linspace(P1_list[0], P1_list[1], 50)
                P1_list = temp
        P2orig = float(self.params["sys_params"][P2])


        #Determine whether this is a multivalue run or not
        if P2+"_list" not in self.params.keys():
            P2_list = [P2orig]
            # colors = [(0,0,1)]
            colors = ['b']
        else:
            P2_list = self.params[P2+"_list"]
            colors = mpl.cm.rainbow(np.linspace(0,1,len(P2_list)))

        #Number of graph types
        n_graph = len(self.params["graph_type"])

        #Make number of subfigs required
        self.fig, axarr = plt.subplots(n_graph, sharex=True)

        #list of different Ising models
        IM_list = []

        #Run Ising Models with different values of param2
        for val2 in P2_list:
            self.params["sys_params"][P2] = val2 
            if self.params["param1"] == "step":
                print 'check'
                IM_list.append(IsingModel(self.params))
                self.RunSystem(IM_list[-1])
            #Make a 2D array if you are running it with anything but step
            else:
                IM_Plist = []
                for val1 in P1_list:
                    self.params["sys_params"][P1] = val1
                    IM_Plist.append(IsingModel(self.params))
                    self.RunSystem(IM_Plist[-1])
                IM_list.append(IM_Plist)
        
        #loop through subfigs for each graph type wanted
        for ax, gt in zip(axarr, self.params["graph_type"]):
            #Make plots on each subfig 
            if self.params['param1'] == 'step':
                self.GraphVsTime( ax, IM_list, colors, gt, P2_list, P2)
            else:
                self.GraphVsParam( ax, IM_list, colors, gt, P1_list, P1, P2_list, P2)
        
        # self.params["sys_params"][P1] = P1orig
        # self.params["sys_params"][P2] = P2orig
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

    def GraphVsParam(self, ax, IM_list, colors, gt, P1_list, P1, P2_list, P2):
        for IM_Plist, color, val2 in zip(IM_list, colors, P2_list):
            y_array=[]
            for m in IM_Plist:
                y_array.append(self.GetValStat(m, gt))  

            a = np.asarray(y_array)
            ax.errorbar(P1_list, a[:,0], yerr=np.sqrt(a[:,1]), color=color, marker='o', label = '{0} = {1}'.format(P2, val2)) 
        ax.set_title("{} vs {}".format(gt, P1))
        ax.set_xlabel(P1)
        if gt == 'C' or gt == 'X':
            ax.set_ylabel("{} / V ".format(gt))
        else:
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
        # elif Val == 'C':
            # return model.C_array
        # elif Val == 'X':
            # return model.S_array

    def GetValStat(self, model, Val):
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

        self.seed = int(self.params['sys_params']["seed"])
        np.random.seed(self.seed)
        
        #Get variables used in algorithm
        self.beta = 1/float(self.params["sys_params"]["T"])
        self.H = float(self.params["sys_params"]["H"])
        self.J = float(self.params["sys_params"]["J"])
        self.L = int(self.params["sys_params"]["L"])
        self.V = float(self.L*self.L)
        self.step_array = range(self.params["sys_params"]["n_steps"]) #Array of all steps
        self.bins = self.params["sys_params"]["n_bins"]

        #Create empty lists for data output
        self.M_array = [] #Magnetization at every step
        self.E_array = [] #Energy at every step
        # self.C_array = [] #Specific Heat at every step
        # self.X_array = [] #Susceptibility at every step
        
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

           # self.M_array.append(self.GetMagnetization())
           # self.E_array.append(self.GetTotalEnergy())

        print self.GetAvgEnergy()
        print self.GetAvgMag()
        print self.GetAvgC()
        print self.GetAvgX()
        return


    def RunMetropolis(self):
        for i_step in self.step_array:
           i, j = self.ChooseRandomParticle() 
           self.ThermalizeParticleMet(i, j)

           # self.M_array.append(self.GetMagnetization())
           # self.E_array.append(self.GetTotalEnergy())

        print self.GetAvgEnergy()
        print self.GetAvgMag()
        print self.GetAvgC()
        print self.GetAvgX()
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
            self.SetTotalEnergy()
            self.SetMagnetization()
        else:
            self.system[i,j] = -1
            self.SetTotalEnergy()
            self.SetMagnetization()
        return

    def ThermalizeParticleMet(self, i, j):
        p, E = self.ProbFlip(i, j)

        if ( p == 1 or p > np.random.ranf()):
            self.system[i,j] *= -1
            self.SetTotalEnergy(E)
            self.SetMagnetization(self.system[i,j])
        else:
            self.SetTotalEnergy(0)
            self.SetMagnetization(0)

        return

    def GenerateRandomMatrix(self, N):
        M = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                M[i][j] = np.random.choice([-1,1])
        return M

    def ChooseRandomParticle(self):
        i = np.random.randint(self.L)
        j = np.random.randint(self.L)
        return i, j

    #Probability calculating algorithms
    def ProbUp(self, E_up, E_down):
        P = np.exp(-self.beta*E_up)/(np.exp(-self.beta*E_up) + np.exp(-self.beta*E_down))
        return P

    def ProbFlip(self, i, j):
        E_c = self.GetParticleEnergy(i, j)
        E_f = -E_c
        if (E_f-E_c) < 0:
            return 1, E_f
        else:
            return np.exp(-(E_f-E_c)*self.beta), E_f

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

    #Can make this more efficient
    def GetMagnetization(self):
        return self.M_array[-1]

    #Can make this more efficient
    def GetTotalEnergy(self):
        return self.E_array[-1]
        # E = 0.0
        # for i in range(self.L):
            # for j in range(self.L):
                # E += self.GetParticleEnergy(i, j)
     
    def SetTotalEnergy(self, E=None):
        if self.E_array == [] or E==None:
            E = 0.0
            for i in range(self.L):
                for j in range(self.L):
                    E += self.GetParticleEnergy(i, j)
            for i in range(self.L):
                for j in range(self.L):
                    E += self.GetParticleEnergy(i, j)
            self.E_array.append(E/(4*self.V))
        else:
            self.E_array.append(self.E_array[-1]+(E/(self.V)))

    def SetMagnetization(self, s=None):
        if s == None or self.M_array == []:
            self.M_array.append(np.sum(self.system)/self.V)
        else:
            self.M_array.append((self.M_array[-1])+(2*s/self.V))
    
    def GetAvgEnergy(self):
        mean, mean_var, var, var_var = GetStats(self.E_array, self.bins)
        return [mean, mean_var]
    
    def GetAvgAbsMag(self):
        mean, mean_var, var, var_var = GetStats(np.absolute(self.M_array), self.bins)
        return [mean, mean_var]

    def GetAvgMag(self):
        mean, mean_var, var, var_var = GetStats(self.M_array, self.bins)
        return [mean, mean_var]

    def GetAvgC(self):
        mean, mean_var, var, var_var = GetStats(self.E_array, self.bins)
        return [self.V*var, self.V*var_var]

    def GetAvgX(self):
        mean, mean_var, var, var_var = GetStats(self.M_array, self.bins)
        return [self.V*self.beta*var, self.V*self.beta*mean_var]


if __name__ == "__main__":
    ens = IsingModelEnsemble(sys.argv[1])
    ens.Graph()
    ens.ShowPlots()
    ens.SavePlots()







