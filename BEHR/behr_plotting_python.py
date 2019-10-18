"""
This code is the original work of Simon Wright used in their 2017 Wesleyan Undergratuate Thesis
"Exploring the High Energy Sky: A Survey of X-ray Sources in Galaxies Within 15 Megaparsecs of the Milky Way".
This thesis was submitted to the faculty of Wesleyan University in partial fulfillment of the requirements for the Degree of Bachelor of Arts
with Departmental Honors in Astronomy.

This version of the code is now modifed by Anthony Santini (Wesleyan 2019 Master's alumnus) for use in current research.


"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import readsav

def Plot(Error_Bar_B=True):
    data = readsav('BEHR_results_allsources.sav')

    # print data.keys()

    hardmean = data['hr2mean']
    softmean = data['hr1mean']
    harderror = [data['hr2lower'],data['hr2upper']]
    softerror = [data['hr1lower'],data['hr1upper']]

    fig = plt.figure(figsize=(10,8))
    if(Error_Bar_B):
        plt.suptitle('Bayesian Estimated Hardness Ratios with Errors',fontsize=16)
    else:
        plt.suptitle('Bayesian Estimated Hardness Ratios',fontsize=16)
    ax1 = plt.subplot2grid((4,4),(0,0),rowspan=3)
    ax2 = plt.subplot2grid((4,4),(0,1),colspan=3, rowspan=3, sharey=ax1)
    ax4 = plt.subplot2grid((4,4),(3,1),colspan=3, sharex=ax2)
    ax4.yaxis.tick_right()
    ax4.set_xlabel(r'$\mathrm{H}-\mathrm{M}$ $/$ $\mathrm{H}+\mathrm{M}$')
    ax1.set_ylabel(r'$\mathrm{M}-\mathrm{S}$ $/$ $\mathrm{M}+\mathrm{S}$')

    fig.subplots_adjust(hspace=0,wspace=0,top=0.89)

    ax2.set_xlim(-1,1)
    ax2.set_ylim(-1,1)
    if(Error_Bar_B):
        ax2.errorbar(hardmean,softmean,xerr=harderror,yerr=softerror,alpha=.4,fmt='+',markersize=0.0001)#,fmt='o')
    ax2.plot(hardmean,softmean,'+g')
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax1.xaxis.tick_top()

    ax4.hist(hardmean,bins=np.linspace(-1,1,40),color='g')
    ax1.hist(softmean,bins=np.linspace(-1,1,40),orientation=u'horizontal',color='g')

    # plt.show()
    if(Error_Bar_B):
        plt.savefig('behr_histplot_with_EB.pdf')
        plt.savefig('behr_histplot_with_EB.png',dpi=1500)
    else:
        plt.savefig('behr_histplot.pdf')
        plt.savefig('behr_histplot.png',dpi=1500)
    plt.clf()

Plot()
Plot(Error_Bar_B=False)
