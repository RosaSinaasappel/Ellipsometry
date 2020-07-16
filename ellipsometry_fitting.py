
"""
@author: Rosa Sinaasappel 2020
rosa.sinaasappel@gmail.com
"""
"""
#Run this to initialize a file with material parameters in the right format. You should work in the resulting txt file after this
datas = np.array([[1.0,0.0,1000],[3.88,0.019,10000],[1.45,0.0,2]])
materials = pd.DataFrame(data = datas, index = ['Air', 'Si', 'SiO2'], columns = ['n','k','d'])
materials.to_csv(r'C:/Users/Rosa/Documents/ellipsometer/materials.txt',sep = '\t')
"""


from matplotlib import pyplot  as plt
import numpy as np
import math
import cmath
import pandas as pd
import lmfit as lm
from lmfit import Model

#-------------------------------------------------------------------------------------------------
"""user defined input"""

#path to measured data
locdata = 'D:/small lab project/metingen/0,5condens.phi'
#path to matrials file
locmaterials = 'C:/Users/Rosa/OneDrive/Documents/studie/master-citrus-sinensis/ellipsometer/data_analysis/materials.txt'
wavelength = 632.8 #in nanometer
#--------------------------------------------------------------------------------------------------


""" Initialize the problem"""
#Import measurement and the matrials file
def importdata(locdata, locmaterials):
    data = pd.read_csv(locdata ,sep = '\s+',skiprows = 2, names = ['phi', 'Delta', 'Psi'])
    materials = pd.read_csv(locmaterials ,sep = '\t',index_col = 0)
    return data, materials

#this function will inquire the user for information about the measured system
#through prompts you can input the material of the layers, what quantities to fit
#and what the bounds are for the fit.
#the output is used is the fitting function.
def problem_parameters(materials):
    layers = []
    fitn = []
    fitd = []
    fitk = []
    nrange =[]
    drange = []
    krange = []
    number_of_layers = int(input('Number of layers? '))
    #count from baselayer up. A bit confusing, sorry
    for i in range (0,number_of_layers+1):
        l = '%i'%i
        if i == 0:
            l = 'base'
        if i == number_of_layers:
            l = 'ambiant'
        layer = input('material layer %s '%l)#ask for name of layer
        n  = input('fit real index of refraction? [y/n] ')
        if n == 'y':
            fitn.append(True)
            nr = input('error on initial guess? ')
            nr = float(nr)
            nrange.append(nr)
        if n == 'n':
            fitn.append(False)
            nrange.append(0)
        k  = input('fit complex index of refraction? [y/n] ')
        if k == 'y':
            fitk.append(True)
            kr = input('error on initial guess? ')
            kr = float(kr)
            krange.append(kr)
        if k == 'n':
            fitk.append(False)
            krange.append(0)

        d = input('fit layer thickness? [y/n] ')
        if d == 'y':
            fitd.append(True)
            dr = input('error on initial guess? (in nm) ')
            dr = float(dr)
            drange.append(dr)
        if d == 'n':
            fitd.append(False)
            drange.append(0)

        layers.append(layer)
    fit_parameters = materials.loc[layers]
    fit_parameters.insert(1,'Error on init. guess n',nrange)
    fit_parameters.insert(2,'Fit n?',fitn)
    fit_parameters.insert(4,'Error on init. guess k',krange)
    fit_parameters.insert(5,'Fit k?',fitk)
    fit_parameters.insert(7,'Error on init. guess d',drange)
    fit_parameters.insert(8,'Fit d?', fitd)
    return fit_parameters, number_of_layers
#------------------------------------------------------------------------------
"""possible models"""

#
def Rsp(r01,r12,beta): #calculate total s and p refelaction
    x = r01+r12*cmath.exp(-1j*2*beta)
    y = 1+r01*r12*cmath.exp(-1j*2*beta)
    return (x/y)
def r_ab(na,nb,phia,phib): #calculate reflection coefficient between layer a and b
    x = na*cmath.cos(phia)-nb*cmath.cos(phib)
    y = na*cmath.cos(phia)+nb*cmath.cos(phib)
    return (x/y)
def snell(n0,ni,phi0): #snell's law
    phii = cmath.asin((n0*cmath.sin(phi0))/ni)
    return phii

#model in the case of two layers, taken from "ellipsometry and polarized light"
#by Azzam and Bashara  
def model1layer(phi, wavelength, n1, k1, n2, k2):
    n_0 = n1-1j*k1 #air
    n_1 = n2-1j*k2 #substrate
    rholist = []
    Psilist = []
    Deltalist = []
    for phi0 in phi:
        phi0 = math.radians(phi0)
        phi1 = snell(n_0,n_1,phi0)
        Rp = cmath.tan(phi0-phi1)/cmath.tan(phi0+phi1)
        Rs = (-cmath.sin(phi0-phi1))/cmath.sin(phi0+phi1)
        rho = (Rp/Rs)
        rholist.append(rho)
        #rholistim.append(rho.imag)
        Psilist.append(math.degrees(math.atan(abs(rho))))
        Deltalist.append(math.degrees(cmath.phase(rho)))
    return rholist
    
#model in the case of two layers, taken from "ellipsometry and polarized light"
#by Azzam and Bashara    
def model2layer(phi, wavelength, n1, k1, n2, k2, n3, k3, d ):
    #substrate = n3-->n_0, film = n2-->n_1, ambiant = n1 --> n_2
    n_0 = n1-1j*k1 #air
    n_1 = n2-1j*k2 #film
    n_2 = n3-1j*k3 #substrate
    rholist =[]
    #rholistim = []
    Psilist = []
    Deltalist = []
    
    for phi0 in phi:
        phi0 = math.radians(phi0)
        #x = cmath.sqrt((n_1**2)-(n_0**2)*(math.sin(phi0)**2))
        phi1 = snell(n_0,n_1,phi0) 
        phi2 = snell(n_0,n_2,phi0)
        beta = 2*math.pi*(d/wavelength)*n_1*cmath.cos(phi1)
        r01s = r_ab(n_0,n_1,phi0, phi1)
        r12s = r_ab(n_1,n_2,phi1, phi2)
        r01p = r_ab(n_1,n_0,phi0, phi1)
        r12p = r_ab(n_2,n_1,phi1, phi2)
        Rs = Rsp(r01s,r12s,beta)
        Rp = Rsp(r01p,r12p,beta)
        rho = (Rp/Rs)
        rholist.append(rho)
        #rholistim.append(rho.imag)
        Psilist.append(math.degrees(math.atan(abs(rho))))
        Deltalist.append(math.degrees(cmath.phase(rho)))
    return rholist
    #return Psilist, Deltalist

#model in the case of three layers, taken from "ellipsometry and polarized light"
#by Azzam and Bashara  
def model3layer(phi, wavelength, n1, k1, n2, k2, n3, k3, n4, k4, d1, d2 ):
    #substrate = n3-->n_0, film = n2-->n_1, ambiant = n1 --> n_2
    n_0 = n1-1j*k1 #air
    n_1 = n2-1j*k2 #film 1
    n_2 = n3-1j*k3 #film 1    
    n_3 = n4-1j*k4 #substrate
    rholist =[]
    #rholistim = []
    Psilist = []
    Deltalist = []
    
    for phi0 in phi:
        phi0 = math.radians(phi0)
        #x = cmath.sqrt((n_1**2)-(n_0**2)*(math.sin(phi0)**2))
        phi1 = snell(n_0,n_1,phi0) 
        phi2 = snell(n_0,n_2,phi0)
        phi3 = snell(n_0,n_3,phi0)
        beta1 = 2*math.pi*(d1/wavelength)*n_1*cmath.cos(phi1)
        beta2 = 2*math.pi*(d2/wavelength)*n_2*cmath.cos(phi2)
        r01s = r_ab(n_0,n_1,phi0, phi1)
        r12s = r_ab(n_1,n_2,phi1, phi2)
        r23s = r_ab(n_2,n_3,phi2, phi3)
        r123s = Rsp(r12s,r23s,beta2)
        r01p = r_ab(n_1,n_0,phi0, phi1)
        r12p = r_ab(n_2,n_1,phi1, phi2)
        r23p = r_ab(n_3,n_2,phi2, phi3)
        r123p = Rsp(r12p,r23p,beta2)
        
        Rs = Rsp(r01s,r123s,beta1)
        Rp = Rsp(r01p,r123p,beta1)
        rho = (Rp/Rs)
        rholist.append(rho)
        #rholistim.append(rho.imag)
        Psilist.append(math.degrees(math.atan(abs(rho))))
        Deltalist.append(math.degrees(cmath.phase(rho)))
    return rholist
               

def model2layerBrewster(phi, wavelength, n1, k1, n2, k2, n3, k3, d ):
    """find rho from these variables. n0,n1,n2 are complex numbers, d is real
    and phi is a list"""
    n_1 = n1-1j*k1
    n_2 = n2-1j*k2
    n_3 = n3-1j*k3
    x = (((n_1**2)-(n_2**2))*(n_3**2))
    y = (cmath.sqrt((n_1**2)+(n_2**2))*((n_3**2)-(n_2**2))*((n_3**2)-(n_1**2)))
    rho = (y/x)*(math.pi/wavelength)*d
    rholist = []
    for i in range(0,len(phi)):
        rholist.append(rho)
    
    #return np.append(np.real(rholist), np.imag(rholist))
    return rholist
        
#------------------------------------------------------------------------------
"""find values to fit to"""
def measuredRho(data):
    """find rho from measurements. data should be a pandas dataframe"""
    rho = []
    for i in range (0,len(data.phi)):
        #convert delta and psi to radians
        psir = math.radians(data.Psi[i])
        delta = data.Delta[i]
        if delta > 300:
            delta = delta-360
        deltar = math.radians(data.Delta[i])
        #find rho
        rhoi = cmath.tan(psir)*cmath.exp(deltar*1j)
        rho.append(rhoi)
    """ To Do: find error on rho"""
    data.insert(3,'rho', rho)
    return data


#------------------------------------------------------------------------------
def fitfunction(layermodel,fit_parameters,parameterlist, data):
    lmodel = Model(layermodel)
    
    #set parameters including initial values, errors, and which parameters to fit
    lmodel.set_param_hint('wavelength', value=wavelength,  vary = False)

    for i in range(0,len(parameterlist[:,0])):
        s = parameterlist[i,0]
        layer = parameterlist[i,1]
        if 'n' in s:
            Value = fit_parameters.at[layer,'n']
            Vary = fit_parameters.at[layer,'Fit n?']
            Min = Value-fit_parameters.at[layer, 'Error on init. guess n']
            Max = Value+fit_parameters.at[layer, 'Error on init. guess n']
        if 'k' in s:
            Value = fit_parameters.at[layer,'k']
            Vary = fit_parameters.at[layer,'Fit k?']
            Min = Value-fit_parameters.at[layer, 'Error on init. guess k']
            Max = Value+fit_parameters.at[layer, 'Error on init. guess k']
        if 'd' in s: 
            Value = fit_parameters.at[layer,'d']
            Vary = fit_parameters.at[layer,'Fit d?']
            Min = Value-fit_parameters.at[layer, 'Error on init. guess d']
            Max = Value+fit_parameters.at[layer, 'Error on init. guess d']
        if Min == Max:
            Min = -np.inf
            Max = np.inf
        lmodel.set_param_hint('%s'%s, value=Value, vary=Vary, min=Min, max=Max)
        
    parameters = lmodel.make_params()
    parameters.pretty_print()
    
    rho = data['rho'].to_numpy()
    #rhoim = rho.imag
    #rhore = rho.real
    phi = data['phi'].to_numpy()
    #resultre = lmodel.fit(rhore, parameters, phi = phi)
    #resultim = lmodel.fit(rhoim, parameters, phi = phi)
    result = lmodel.fit(rho, parameters, phi = phi)#, method = 'nelder')
    #return resultim, resultre, rhoim, rhore, phi
    return result, phi, rho
    
    
        


#------------------------------------------------------------------------------
""" aquire all data"""

#import data and assamble dataframe with the details of the problem
data, materials = importdata(locdata, locmaterials)
fit_parameters, number_of_layers = problem_parameters(materials)

#add measured rho to data
data = measuredRho(data)



#------------------------------------------------------------------------------
""" pick the right model and fit"""
if number_of_layers == 1 :
    
    #link layer names to their parameters
    layers = fit_parameters.index
    parameterlist = np.array([['n1',layers[1]], ['k1',layers[1]], 
                              ['n2',layers[0]], ['k2',layers[0]]])
    
    #fit
    fit, phi, rho = fitfunction(model1layer, fit_parameters, parameterlist, data)


if number_of_layers == 2 :
    
    #link layer names to their parameters
    layers = fit_parameters.index
    parameterlist = np.array([['n1',layers[2]], ['k1',layers[2]],
                              ['n2',layers[1]], ['k2',layers[1]],
                              ['n3',layers[0]], ['k3',layers[0]],
                              ['d',layers[1]]])
    
    #fit
    fit, phi, rho = fitfunction(model2layer, fit_parameters, parameterlist, data)
    
  
    
if number_of_layers == 3 :
        #link layer names to their parameters
    layers = fit_parameters.index
    parameterlist = np.array([['n1',layers[3]], ['k1',layers[3]],
                              ['n2',layers[2]], ['k2',layers[2]],
                              ['n3',layers[1]], ['k3',layers[1]],
                              ['n4',layers[0]], ['k4',layers[0]],
                              ['d1',layers[2]], ['d2',layers[1]]])
    
    #fit
    #fitre, fitim, rhoim, rhore, phi = fitfunction(model2layer, fit_parameters, parameterlist, data)
    fit, phi, rho = fitfunction(model3layer, fit_parameters, parameterlist, data)


#------------------------------------------------------------------------------
"""plot results"""

#plot rho, real and imaginary part
plt.figure()
plt.title('Real fit')
plt.plot(phi, rho, 'bo')
plt.plot(phi, fit.init_fit, 'k--', label='initial fit')
plt.plot(phi, fit.best_fit, 'r-', label='best fit')
plt.legend(loc='best')
plt.xlabel('Angle of incidence (phi)')
plt.ylabel('Re[rho]')

plt.figure()
plt.title('imaginary fit')
plt.plot(phi, np.imag(rho), 'bo')
plt.plot(phi, np.imag(fit.init_fit), 'k--', label='initial fit')
plt.plot(phi, np.imag(fit.best_fit), 'r-', label='best fit')
plt.legend(loc='best')
plt.xlabel('Angle of incidence (phi)')
plt.ylabel('Im[rho]')

print(fit.fit_report())



#plot Psi and Delta
Psifit = []
Deltafit = []
rhofit = fit.best_fit
for r in rhofit:
    Psifit.append(math.degrees(math.atan(abs(r))))
    Deltafit.append(math.degrees(cmath.phase(r)))
Psifitinit = []
Deltafitinit = []
rhofitinit = fit.init_fit
for r in rhofitinit:
    Psifitinit.append(math.degrees(math.atan(abs(r))))
    Deltafitinit.append(math.degrees(cmath.phase(r)))

plt.figure()
plt.title('Delta')
plt.xlabel('Angle of incidence (phi)')
plt.ylabel('Delta')
plt.plot(phi,data.Delta, 'bo', label = 'data')
plt.plot(phi,Deltafit, 'r-', label = 'fit')
plt.plot(phi,Deltafitinit, 'k--', label = 'initial fit')
plt.legend()

plt.figure()
plt.title('Psi')
plt.xlabel('Angle of incidence (phi)')
plt.ylabel('Psi')
plt.plot(phi,data.Psi, 'bo', label = 'data')
plt.plot(phi,Psifit, 'r-', label = 'fit')
plt.plot(phi,Psifitinit, 'k--', label = 'initial fit')
plt.legend()


#-----------------------------------------------------------------------------
#diagnostic sripts
"""
n1 = 1.00
k1 = 0
n2 = 1.457
k2 = 0.0
n3 = 3.8827
k3 = 0.019626
d = 6

Psifit, Deltafit = model2layer(data.phi, 633, n1, k1, n2, k2, n3, k3, d)
plt.figure()
plt.plot(data.phi, Deltafit)
plt.plot(data.phi, data.Delta)

plt.figure()
plt.plot(data.phi, Psifit)
plt.plot(data.phi, data.Psi)
"""
"""
n1 = 1.00
k1 = 0
n2 = 1.457
k2 = 0.0
n3 = 3.8827
k3 = 0.2
d = 6.5

plt.figure()
plt.title('n3')
for n3 in np.arange(3.5,4.5,0.1):
    rhos = np.imag(model2layer(phi, wavelength, n1, k1, n2, k2, n3, k3, d ))
    plt.plot(phi,rhos, label = n3)
    plt.legend()
plt.show()
"""
"""
datatest = data[{'phi','Delta','Psi'}]
plt.figure()
plt.title('Delta')
for i in np.arange(-10,10,1):
    datatest.Delta = data.Delta + i
    datatest = measuredRho(datatest)
    plt.plot(datatest.phi,datatest.rho, label = i)
    datatest = datatest.drop(columns = 'rho')
plt.legend()
plt.show()
"""  
"""
rslist = []
rplist = [] 
for a in np.arange(30,50,0.1):
    n_air = 1
    n_water = 1.331
    phi_water = snell(n_air,n_water, a)
    rs01 = r_ab(n_air, n_water, a, phi_water)
    rslist.append(rs01)
    rp01 = r_ab(n_air, n_water, phi_water, a)
    rplist.append(rp01)

plt.figure()
plt.plot(np.arange(30,50,0.1), rslist)
plt.plot(np.arange(30,50,0.1), rplist)
"""
"""

locdata2 = 'D:/small lab project/metingen/glass.phi'
data2, materials = importdata(locdata2, locmaterials)
data2 = measuredRho(data2)


plt.figure()
plt.title('Delta')
plt.xlabel('phi')
plt.ylabel('Delta')
plt.plot(data2.phi,data2.Delta, 'ro', label = 'glass')
plt.plot(data.phi,data.Delta, 'bo', label = 'glass+ propanol')
plt.legend()

plt.figure()
plt.title('Psi')
plt.xlabel('phi')
plt.ylabel('Psi')
plt.plot(data2.phi,data2.Psi, 'ro', label = 'glass')
plt.plot(data.phi,data.Psi, 'bo', label = 'glass+propanol')
plt.legend()


plt.figure()
plt.title('Real fit')
plt.plot(data2.phi,data2.rho, 'ro', label = 'glass')
plt.plot(data.phi,data.rho, 'bo', label = 'glass+propanol')
plt.legend(loc='best')
plt.xlabel('phi')
plt.ylabel('Re[rho]')

plt.figure()
plt.title('imaginary fit')
plt.plot(data2.phi,np.imag(data2.rho), 'ro', label = 'glass')
plt.plot(data.phi,np.imag(data.rho), 'bo', label = 'glass+propanol')
plt.xlabel('phi')
plt.ylabel('Im[rho]')
plt.legend()
"""