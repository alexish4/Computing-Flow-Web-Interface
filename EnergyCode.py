import numpy as np
import scipy as sp
import pandas as pd
import sklearn as skl
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
import ipywidgets as widgets
from ipywidgets import interact, interact_manual
import gc
import os
import tqdm
import collections
import re
import correlation_data_utilities
import networkx as nx
from itertools import islice
import copy
import MDAnalysis as mda
import nglview as nv
from IPython.display import display, HTML

energyDataDir='energy_topology'
systems=['igps']
variants=['apo','holo']

summaryDataSuffix='decomp_summary.csv'
seriesDataSuffix='energy_series.csv'

nSkip_Summary=6
headerRows_Summary=[0,1]

subDirStr='{system}.{variant}'

dataPrefixStr='{system}.prot_lig_{variant}.mmpbsa.pairwise_decomp'

#The header for the energy summary is rather messy and needs to
#be reformatted to work properly.
def parse_summary_header(headerColumns):
    summaryHead=[]

    head1=''
    for headEntry in headerColumns: #summaryData.columns:
        if not ('Unnamed' in headEntry[0]):
            head1=headEntry[0].replace(' ','_').replace('-','_').replace('.','')
        if 'Unnamed' in headEntry[1]:
            head2=''
        else:
            head2=headEntry[1].replace(' ','_').replace('-','_').replace('.','')
        if head2=='':
            summaryHead.append(head1)
        else:
            summaryHead.append('.'.join([head1,head2]))
    summaryHead=np.array(summaryHead)
    return summaryHead

#lets test it and see if it worked.
#this will also give us a list of columns so we can figure out which ones we need
parse_summary_header(
pd.read_csv(
    'energy_topology/igps.apo/igps.prot_lig_apo.mmpbsa.pairwise_decomp.decomp_summary.csv',
    skiprows=nSkip_Summary,
    header=headerRows_Summary).columns
)

outputFileName='joint_interaction_energy_summary.TOTAL.csv'
outputFilePath='/'.join([energyDataDir,outputFileName])

system=systems[0]
variant=variants[0]
summaryTables=[]
for system in systems:
    for variant in variants:
        print('--- --- ---')
        print('- Loading',end=', ')
        strDict={
            "system":system,
            "variant":variant,}
        inputPrefix=dataPrefixStr.format(**strDict)
        summaryFileName='.'.join([inputPrefix,summaryDataSuffix])
        print('fileName:',summaryFileName)
        subDir=subDirStr.format(**strDict)
        fileDir='/'.join([energyDataDir,subDir])
        summaryFilePath='/'.join([fileDir,summaryFileName])
        print('filePath:',summaryFilePath)
        summaryData=pd.read_csv(summaryFilePath,skiprows=nSkip_Summary,header=headerRows_Summary)
        print('parsing column names',end=', ')
        summaryData.columns=parse_summary_header(summaryData.columns)
        print('extracting interaction terms',end=', ')
        summaryData['System']=system
        summaryData['Variant']=variant
        summaryData['E_Interact.Mean']=summaryData['TOTAL.Avg'] #(summaryData['van_der_Waals.Avg']+\
                                        # summaryData['Electrostatic.Avg'])
        summaryData['E_Interact.Std_Err']=summaryData['TOTAL.Std_Err_of_Mean'] #(summaryData['van_der_Waals.Std_Err_of_Mean']+\
                                          # summaryData['Electrostatic.Std_Err_of_Mean']).abs()
        summaryData=summaryData[np.concatenate([
            ['System','Variant','Resid_1','Resid_2'],
            [colName for colName in summaryData.columns if 'E_Interact' in colName]
        ])]
        print('formatting residue data columns',end=' -\n')
        summaryData['ResNum_1']=summaryData['Resid_1'].map(
            lambda x: re.sub('[A-Za-z ]','',x))
        summaryData['ResNum_2']=summaryData['Resid_2'].map(
            lambda x: re.sub('[A-Za-z ]','',x))
        summaryData['ResName_1']=summaryData['Resid_1'].map(
            lambda x: re.sub('[0-9 ]','',x))
        summaryData['ResName_2']=summaryData['Resid_2'].map(
            lambda x: re.sub('[0-9 ]','',x))
        summaryData=summaryData.drop(columns=['Resid_1','Resid_2'])
        summaryData=summaryData[np.concatenate([
            ['System','Variant'],
            [colName for colName in summaryData.columns if 'Res' in colName],
            [colName for colName in summaryData.columns if 'Interact' in colName]
        ])]
        summaryTables.append(summaryData.copy())
print('--- --- ---')
print('Done')
summaryData=pd.concat(summaryTables)
del(summaryTables)
summaryData.to_csv(outputFilePath,index=False)
summaryData.head()

energyDataDir='energy_topology'
summaryFileName='joint_interaction_energy_summary.TOTAL.csv'
summaryFilePath='/'.join([energyDataDir,summaryFileName])
summaryData=pd.read_csv(summaryFilePath)
summaryData.head()

#Use a cutoff of .5 kb*T energy (note, energies from mmpbsa.py are in kcal/mol)
#This is a somewhat arbitrary choice
T_sim=300.15 #room temperature in Kelvins
kb=.0019872041 #Boltzmann's constant in kcal/(mol*K)

thermalEnergy=kb*T_sim
cutoffEnergy=.5*thermalEnergy

#now we want to group by residue pair and filter out any pairs for which
#the absolute value of the interaction mean is not at least 1 standard error
#above the cutoff in any of the systems.
#Conversely if this holds for either system we will keep the pair for ALL systems
#This will allow us to see how much the energy shifted in cases where one system
#has an interaction and the other normally would not.

tqdm.tqdm.pandas() #progress bar
#compute a filter value column to speed filtering
print("Generating Grouping Column")
summaryData['ResNumPair']=summaryData[
        ['ResNum_1','ResNum_2']
    ].progress_apply(lambda x: '_'.join(map(str,x)),axis=1)
print("Computing Filter Values")
sig_pairs=summaryData['ResNumPair'][
        (summaryData['E_Interact.Mean'].abs() - \
         summaryData['E_Interact.Std_Err']) > cutoffEnergy
    ].unique()
print('Significant Pairs:',sig_pairs)
sigData=summaryData[summaryData['ResNumPair'].isin(sig_pairs)].drop(columns=['ResNumPair'])
sigData.head()

#Save list of significant residue pairs
outSigPairName='significant_pairs.TOTAL.txt'
outSigPairPath='/'.join([energyDataDir,outSigPairName])
np.savetxt(fname=outSigPairPath,X=sig_pairs,fmt='%s')
outSigDataName='significant_interactions_summary.TOTAL.csv'
outSigDataPath='/'.join([energyDataDir,outSigDataName])
sigData.to_csv(outSigDataPath,index=False)

#Save filtered energy data in matrix format for each system
nRes=454
outTopMatNameStr='{system}.{variant}.energy_topology_matrix.TOTAL.txt'
outValMatNameStr='{system}.{variant}.energy_value_matrix.TOTAL.txt'
headerComment='# cutoff energy = %.4e kcal/mol'%cutoffEnergy
for groupName,groupData in sigData.groupby(['System','Variant']):
    outMatName=outTopMatNameStr.format(
        system=groupName[0],variant=groupName[1])
    print('saving {matName} to disk'.format(matName=outMatName))
    outMatPath='/'.join([energyDataDir,outMatName])
    outMat=np.savetxt(
        fname=outMatPath,
        X=sp.sparse.coo_matrix(
            (np.array(groupData['E_Interact.Mean'].abs()>0,dtype=float),
             (np.array(groupData['ResNum_1'],dtype=int)-1,
              np.array(groupData['ResNum_2'],dtype=int)-1)),
            shape=(nRes,nRes)).todense())
    
    outMatName=outValMatNameStr.format(
        system=groupName[0],variant=groupName[1])
    print('saving {matName} to disk'.format(matName=outMatName))
    outMatPath='/'.join([energyDataDir,outMatName])
    outMat=np.savetxt(
        fname=outMatPath,
        X=sp.sparse.coo_matrix(
            (np.array(groupData['E_Interact.Mean'].abs()>0,dtype=float),
             (np.array(groupData['ResNum_1'],dtype=int)-1,
              np.array(groupData['ResNum_2'],dtype=int)-1)),
            shape=(nRes,nRes)).todense())
print('done')

matNameSearchStr='energy_topology_matrix.TOTAL'

matDir='energy_topology'

matFiles=[fileName for fileName in \
          os.listdir(matDir) if \
          matNameSearchStr in fileName]

topologyMatrixDict={}

print('loading topology matrices:',end=' ')
for matFileName in matFiles:
    print(matFileName,end=', ')
    entryKey='.'.join(list(np.array(matFileName.split('.'))[:2]))
    matFilePath='/'.join([matDir,matFileName])
    topologyMatrixDict[entryKey]=np.loadtxt(matFilePath)
print('')

matNameSearchStr='energy_value_matrix.TOTAL'
matFiles=[fileName for fileName in \
          os.listdir(matDir) if \
          matNameSearchStr in fileName]
valueMatrixDict={}
print('loading value matrices:',end=' ')
for matFileName in matFiles:
    print(matFileName,end=', ')
    entryKey='.'.join(list(np.array(matFileName.split('.'))[:2]))
    matFilePath='/'.join([matDir,matFileName])
    valueMatrixDict[entryKey]=np.loadtxt(matFilePath)
print('')

sigPairName='significant_pairs.TOTAL.txt'
sigPairPath='/'.join([energyDataDir,sigPairName])
sigDataName='significant_interactions_summary.TOTAL.csv'
sigDataPath='/'.join([energyDataDir,sigDataName])


sig_pairs=np.loadtxt(sigPairPath,dtype=str)
sigData=pd.read_csv(sigDataPath)

print('sig. pairs:',sig_pairs)
print('Sig. Data:')

for matKey in valueMatrixDict:
    print(matKey,":",
          np.nonzero(valueMatrixDict[matKey])[0].shape,'entries')
    
nMats=len(topologyMatrixDict)
fig,axs=plt.subplots(nMats,nMats)

fig.set_figwidth(12)
fig.set_figheight(11)

for iMat,matName in enumerate(topologyMatrixDict.keys()):
    print('Plotting {matName}'.format(matName=matName))
    ax=axs[iMat,0]
    sns.heatmap(topologyMatrixDict[matName],ax=ax)
    ax.set_title(matName)
    
    ax=axs[iMat,1]
    sns.heatmap(valueMatrixDict[matName],ax=ax)
    ax.set_title(matName)
    
plt.tight_layout()

strucDataDir='energy_topology'
systems=['igps']
variants=['apo','holo']

strucFileNameStr="{system}.prot_lig_{variant}.pdb"

subDirStr='{system}.{variant}'

strucDict={}
for system in systems:
    for variant in variants:
        subDir=subDirStr.format(system=system,variant=variant)
        strucFileName=strucFileNameStr.format(system=system,variant=variant)
        strucFilePath='/'.join([strucDataDir,subDir,strucFileName])
        strucKey=subDir
        print("loading",strucKey, " ", strucFilePath)
        
        strucDict[strucKey] = mda.Universe(strucFilePath)

print(strucDict)

system='igps'
variant='apo'
nRes=454
Tsim=310.15
kb=0.0019872041
magCut=kb*Tsim*1.0
seqDelta=0
#filter out desired interaction energy edges based
#for the given system based upon energy and
#sequence delta cutoffs
matData=sigData[
    (sigData['System']==system) & \
    (sigData['Variant']==variant) & \
    ((sigData['ResNum_1']-sigData['ResNum_2'])<seqDelta) &
    ((sigData['E_Interact.Mean'].abs() - sigData['E_Interact.Std_Err']) > magCut)]
tempMat=sp.sparse.coo_matrix(
    (matData['E_Interact.Mean'],
     (matData['ResNum_1']-1,matData['ResNum_2']-1)),
    shape=(nRes,nRes))
plotMat=np.abs(tempMat.todense())

print("Nedges=",len(np.nonzero(tempMat)[0]))

nzInds=np.nonzero(plotMat)

#Create a color map for edges based on log(abs(E_Interact.Mean))
tempCmap=matplotlib.cm.get_cmap('Spectral',2048)
tempCmat=np.array((tempMat.todense()))
vMin=np.min(tempCmat)
vCenter=0
vMax=np.max(tempCmat)
cNorm=matplotlib.colors.TwoSlopeNorm(vmin=vMin,vcenter=vCenter,vmax=vMax)
#tempCmat[nzInds]=np.log(tempCmat[nzInds])
edgeColors=correlation_data_utilities.getCorrNetEdgeColors(tempCmat,maskInds=nzInds,cmap=tempCmap,
                                cNorm=cNorm)

#Compute widths for edges based on lob(abs(E_Interact.Mean))
#rMin and rMax set minimum and maximum edgewidths
#edgwidth will then interpolate linearly between those two
#bounds.
eMin=.125
eMax=0.625
radiiMat=correlation_data_utilities.getCorrNetEdgeRadii(plotMat,maskInds=nzInds,eMin=eMin,eMax=eMax)

#Plot distribution histograms for edge widths and edge color rgb values
sns.distplot(radiiMat[np.nonzero(radiiMat)].flatten(),
             kde=False)
plt.title('Edge Width distribution')
plt.show()

clist=['red','green','blue']
for ii in np.arange(3):
    sns.distplot(edgeColors[:,:,ii][np.nonzero(edgeColors[:,:,ii])].flatten(),
                 kde=False,color=clist[ii])
plt.title('Edge rgb color value distribution')
plt.show()

sortArr=np.argsort(tempCmat[np.nonzero(tempCmat)])
sortedInds=np.array(list(zip(nzInds[0],nzInds[1])))[sortArr]
sortedNzInds=(sortedInds[:,0],sortedInds[:,1])

clist=['red','green','blue']
for ii in np.arange(3):
    plt.plot(edgeColors[:,:,ii][sortedNzInds],c=clist[ii])
plt.title('Value Sorted Edge rgb color values')
plt.show()

plt.plot(tempCmat[sortedNzInds])
plt.title('Sorted Interaction_Energy.Mean')
plt.show()

CSS = """
.output {
    display: flex;
    align-items: center;
    text-align: center;
}
"""

HTML('<style>{}</style>'.format(CSS))

#Draw a color bar for the edge coloring,
#(nglview does not seem to support this yet)
correlation_data_utilities.drawEdgeColorBar(tempCmat,maskInds=nzInds,
                 cmap=tempCmap,label='E_Interaction.Mean (kcal/mol)',
                 cNorm=cNorm,
                 orientation='horizontal')

#Render the filtered interaction network using nglview with
#the edge width and colormaps generated above
struc=strucDict[list(strucDict.keys())[0]]
view=nv.show_mdanalysis(struc,overwrite=True)
view.clear_representations()
view.add_representation('cartoon',alpha=.5)
correlation_data_utilities.drawProtCorrMat(protStruc=struc,corrMat=plotMat,ngViewOb=view,
                    frame=0,colorsArray=edgeColors,radiiMat=radiiMat,
                    undirected=True)
view