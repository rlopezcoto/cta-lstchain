"""Module for plotting results from reconstruction.

Usage:
"import plot_dl2"
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
from scipy.stats import norm
from matplotlib import gridspec
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_curve

def plot_features(data,truehadroness=False):
    """Plot the distribution of different features that characterize
    events, such as hillas parameters or MC data.

    Parameters:
    -----------
    data: pandas DataFrame
    
    truehadroness:
    True: True gammas and proton events are plotted (they are separated using true hadroness). 
    False: Gammas and protons are separated using reconstructed hadroness (Hadrorec)
    """
    hadro = "Hadrorec"
    if truehadroness:
        hadro = "hadroness"

    #Energy distribution
    plt.subplot(331)
    plt.hist(data[data[hadro]<1]['mc_energy'],
             histtype=u'step',bins=100,
             label="Gammas")
    plt.hist(data[data[hadro]>0]['mc_energy'],
             histtype=u'step',bins=100,
             label="Protons")
    plt.ylabel(r'# of events',fontsize=15)
    plt.xlabel(r"$log_{10}E$(GeV)")
    plt.legend()

    #Disp distribution
    plt.subplot(332)
    plt.hist(data[data[hadro]<1]['disp'],
             histtype=u'step',bins=100,
             label="Gammas")
    plt.hist(data[data[hadro]>0]['disp'],
             histtype=u'step',bins=100,
             label="Protons")
    plt.ylabel(r'# of events',fontsize=15)
    plt.xlabel(r"Disp (m)")

    #Intensity distribution
    plt.subplot(333)
    plt.hist(data[data[hadro]<1]['intensity'],
             histtype=u'step',bins=100,
             label="Gammas")
    plt.hist(data[data[hadro]>0]['intensity'],
             histtype=u'step',bins=100,
             label="Protons")
    plt.ylabel(r'# of events',fontsize=15)
    plt.xlabel(r"$log_{10}Intensity$")

    dataforwl = data[data['intensity']>np.log10(200)]
    #Width distribution
    plt.subplot(334)
    plt.hist(dataforwl[dataforwl[hadro]<1]['width'],
             histtype=u'step',bins=100,
             label="Gammas")
    plt.hist(dataforwl[dataforwl[hadro]>0]['width'],
             histtype=u'step',bins=100,
             label="Protons")
    plt.ylabel(r'# of events',fontsize=15)
    plt.xlabel(r"Width (º)")

    #Length distribution
    plt.subplot(335)
    plt.hist(dataforwl[dataforwl[hadro]<1]['length'],
             histtype=u'step',bins=100,
             label="Gammas")
    plt.hist(dataforwl[dataforwl[hadro]>0]['length'],
             histtype=u'step',bins=100,
             label="Protons")
    plt.ylabel(r'# of events',fontsize=15)
    plt.xlabel(r"Length (º)")

    #r distribution
    plt.subplot(336)
    plt.hist(data[data[hadro]<1]['r'],
             histtype=u'step',bins=100,
             label="Gammas")
    plt.hist(data[data[hadro]>0]['r'],
             histtype=u'step',bins=100,
             label="Protons")
    plt.ylabel(r'# of events',fontsize=15)
    plt.xlabel(r"r (m)")

    #psi distribution

    plt.subplot(337)
    plt.hist(data[data[hadro]<1]['psi'],
             histtype=u'step',bins=100,
             label="Gammas")
    plt.hist(data[data[hadro]>0]['psi'],
             histtype=u'step',bins=100,
             label="Protons")
    plt.ylabel(r'# of events',fontsize=15)
    plt.xlabel(r"psi angle(rad)")
    
    #psi distribution

    plt.subplot(338)
    plt.hist(data[data[hadro]<1]['phi'],
             histtype=u'step',bins=100,
             label="Gammas")
    plt.hist(data[data[hadro]>0]['phi'],
             histtype=u'step',bins=100,
             label="Protons")
    plt.ylabel(r'# of events',fontsize=15)
    plt.xlabel(r"phi angle(m)")
    
    #Time gradient

    plt.subplot(339)
    plt.hist(data[data[hadro]<1]['time_gradient'],
             histtype=u'step',bins=100,
             label="Gammas")
    plt.hist(data[data[hadro]>0]['time_gradient'],
             histtype=u'step',bins=100,
             label="Protons")
    plt.ylabel(r'# of events',fontsize=15)
    plt.xlabel(r"Time gradient")


def plot_E(data,truehadroness=False):
    
    """Plot the performance of reconstructed Energy. 

    Parameters:
    -----------
    data: pandas DataFrame
    
    truehadroness:
    True: True gammas and proton events are plotted (they are separated using true hadroness). 
    False: Gammas and protons are separated using reconstructed hadroness (Hadrorec)
    
    """
    hadro = "Hadrorec"
    if truehadroness:
        hadro = "hadroness"
    
    gammas = data[data[hadro]==0] 
    
    plt.subplot(221)
    difE = ((gammas['mc_energy']-gammas['Erec'])*np.log(10))
    section = difE[abs(difE) < 1.5]
    mu,sigma = norm.fit(section)
    print(mu,sigma)
    n, bins, patches = plt.hist(difE,100,
                                density=1,alpha=0.75)
    y = norm.pdf( bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=2)
    plt.xlabel('$(log_{10}(E_{gammas})-log_{10}(E_{rec}))*log_{N}(10)$',
               fontsize=10)
    plt.figtext(0.15,0.7,'Mean: '+str(round(mu,4)),
                fontsize=10)
    plt.figtext(0.15,0.65,'Std: '+str(round(sigma,4)),
                fontsize=10)
    
    plt.subplot(222)
    hE = plt.hist2d(gammas['mc_energy'],gammas['Erec'],
                    bins=100)
    plt.colorbar(hE[3])
    plt.xlabel('$log_{10}E_{gammas}$',
               fontsize=15)
    plt.ylabel('$log_{10}E_{rec}$',
               fontsize=15)
    plt.plot(gammas['mc_energy'],gammas['mc_energy'],
             "-",color='red')

    #Plot a profile
    subplot = plt.subplot(223)
    means_result = scipy.stats.binned_statistic(
        gammas['mc_energy'],[difE,difE**2],
        bins=50,range=(1,6),statistic='mean')
    means, means2 = means_result.statistic
    standard_deviations = np.sqrt(means2 - means**2)
    bin_edges = means_result.bin_edges
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
    
    #fig = plt.figure()
    gs = gridspec.GridSpecFromSubplotSpec(2, 1,
                                          height_ratios=[2, 1],
                                          subplot_spec=subplot)
    ax0 = plt.subplot(gs[0])
    plot0 = ax0.errorbar(x=bin_centers, y=means, yerr=standard_deviations, 
                         linestyle='none', marker='.')
    plt.ylabel('$(log_{10}(E_{true})-log_{10}(E_{rec}))*log_{N}(10)$',
               fontsize=10)
    
    ax1 = plt.subplot(gs[1], sharex = ax0)
    plot1 = ax1.plot(bin_centers,standard_deviations,
                     marker='+',linestyle='None')
    plt.setp(ax0.get_xticklabels(), 
             visible=False)
    yticks = ax1.yaxis.get_major_ticks()
    yticks[-1].label1.set_visible(False)
    plt.ylabel("Std",fontsize=10)
    plt.xlabel('$log_{10}E_{true}$',
               fontsize=10)
    plt.subplots_adjust(hspace=.0)

    
def plot_Disp(data,truehadroness=False):
    
    """Plot the performance of reconstructed position

    Parameters:
    -----------
    data: pandas DataFrame
    
    truehadroness: boolean
    True: True gammas and proton events are plotted (they are separated
    using true hadroness). 
    False: Gammas and protons are separated using reconstructed
    hadroness (Hadrorec)
    """
    hadro = "Hadrorec"
    if truehadroness:
        hadro = "hadroness"

    gammas = data[data[hadro]==0] 

    plt.subplot(221)
    difD = ((gammas['disp']-gammas['Disprec'])/gammas['disp'])
    section = difD[abs(difD) < 0.5]
    mu,sigma = norm.fit(section)
    print(mu,sigma)
    n,bins,patches = plt.hist(difD,100,density=1,
                              alpha=0.75,range=[-2,1.5])
    y = norm.pdf( bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=2)
    plt.xlabel('$\\frac{Disp_{gammas}-Disp_{rec}}{Disp_{gammas}}$',
               fontsize=15)
    plt.figtext(0.15,0.7,'Mean: '+str(round(mu,4)),
                fontsize=12)
    plt.figtext(0.15,0.65,'Std: '+str(round(sigma,4)),
                fontsize=12)
                
    plt.subplot(222)
    hD = plt.hist2d(gammas['disp'],gammas['Disprec'],
                    bins=100,range=([0,1.1],[0,1.1]))
    plt.colorbar(hD[3])
    plt.xlabel('$Disp_{gammas}$',
               fontsize=15)
    plt.ylabel('$Disp_{rec}$',
               fontsize=15)
    plt.plot(gammas['disp'],gammas['disp'],
             "-",color='red')   
 
    plt.subplot(223)
    theta2 = (gammas['src_x']-gammas['src_xrec'])**2
    +(gammas['src_y']-gammas['src_y'])**2
    plt.hist(theta2,bins=100,
             range=[0,0.1],histtype=u'step')
    plt.xlabel(r'$\theta^{2}(º)$',
               fontsize=15)
    plt.ylabel(r'# of events',
               fontsize=15)
    

def plot_pos(data,truehadroness=False):
    
    """Plot the performance of reconstructed position

    Parameters:
    data: pandas DataFrame
    
    truehadroness: boolean
    True: True gammas and proton events are plotted (they are separated
    using true hadroness). 
    False: Gammas and protons are separated using reconstructed
    hadroness (Hadrorec)
    """
    hadro = "Hadrorec"
    if truehadroness:
        hadro = "hadroness"

    #True position

    trueX = data[data[hadro]==0]['src_x']
    trueY = data[data[hadro]==0]['src_y']
    trueXprot = data[data[hadro]==1]['src_x']
    trueYprot = data[data[hadro]==1]['src_y']

    #Reconstructed position

    recX = data[data[hadro]==0]['src_xrec']
    recY = data[data[hadro]==0]['src_yrec']
    recXprot = data[data[hadro]==1]['src_xrec']
    recYprot = data[data[hadro]==1]['src_yrec']

    plt.subplot(221)
    plt.hist2d(trueXprot,trueYprot,
               bins=100,label="Protons")
    plt.colorbar()
    plt.title("True position Protons")
    plt.xlabel("x(m)")
    plt.ylabel("y (m)")
    plt.subplot(222)
    plt.hist2d(trueX,trueY,
               bins=100,label="Gammas",
               range=np.array([(-1, 1), (-1, 1)]))
    plt.colorbar()
    plt.title("True position Gammas")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.subplot(223)
    plt.hist2d(recXprot,recYprot,
               bins=100,label="Protons",
               range=np.array([(-1, 1), (-1, 1)]))
    plt.colorbar()
    plt.title("Reconstructed position Protons")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.subplot(224)
    plt.hist2d(recX,recY,
               bins=100,label="Gammas",
               range=np.array([(-1, 1), (-1, 1)]))
    plt.colorbar()
    plt.title("Reconstructed position Gammas ")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")


def plot_importances(clf,features):

    importances = clf.feature_importances_
    std = np.std([tree.feature_importances_ for tree in clf.estimators_],
                 axis=0)
    indices = np.argsort(importances)[::-1]

    print("Feature importances (gini index)")
    for f in range(len(features)):
        print("%d. %s (%f)" % (f + 1, 
                               features[indices[f]], 
                               importances[indices[f]]))

    ordered_features=[]
    for index in indices:
        ordered_features=ordered_features+[features[index]]
    
    plt.title("Feature importances for G/H separation",
              fontsize=15)
    plt.bar(range(len(features)), 
            importances[indices],
            color="r", yerr=std[indices], align="center")
    plt.xticks(range(len(features)), 
               ordered_features)
    plt.xlim([-1, 
              len(features)])

def plot_ROC(clf,data,features, Energy_cut):
    # Plot ROC curve:
    check = clf.predict_proba(data[features])[0:,1]
    accuracy = accuracy_score(data['hadroness'], 
                              data['Hadrorec'])
    print(accuracy)
    
    fpr_rf, tpr_rf, _ = roc_curve(data['hadroness'],
                                  check)
    
    plt.plot(fpr_rf, tpr_rf, 
             label='Energy Cut: '+'%.3f'%(pow(10,Energy_cut)/1000)+' TeV')
    plt.xlabel('False positive rate',
               fontsize=15)
    plt.ylabel('True positive rate',
               fontsize=15)
    plt.legend(loc='best')
