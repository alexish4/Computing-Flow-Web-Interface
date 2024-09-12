import numpy as np
import scipy as sp
import pandas as pd
import sklearn as skl
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
# import ipywidgets as widgets
# from ipywidgets import interact, interact_manual
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
import tempfile
from io import StringIO

def visualizeBetweenness():
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

    #In order to apply current flow betweenness methods, we need one or more source
    #and target residues.
    #For IGPS, two well known interaction residues are LEU50 (Allosteric Ligand biniding pocket)
    #and GLU 180 (active site residue)
    #Our structure here has 2 chains, the first is the allosteric region and the
    #second contains the active site.
    #The allosteric region is 253 residues long... and numbering in python starts at 0
    #so our two residues of interest are 50-1 = 49
    #and 253+180-1=432
    #So residues 49 and 432 will serve as the source and target residues repsectively
    print("Keys in strucDict:", strucDict.keys())

    struc=strucDict[list(strucDict.keys())[0]]
    print("Selected structure:", struc)

    print (struc.select_atoms("resid 50").residues[0])
    print (struc.select_atoms("resid 433").residues[0])
    sourceSet=[49]
    targetSet=[432]

    def calculatePathLength(pathGraph,path,weight='weight'):
        return(np.sum([pathGraph.edges()[(edge[0],edge[1])][weight] \
                    for edge in zip(path[:-1],path[1:])]))

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
        shape=(454,454))
    wMat=tempMat.todense()
    nzInds=np.nonzero(tempMat)
    wMat[nzInds]=1./(np.abs(wMat[nzInds]))
    wMat=np.array(wMat)
    wMatGraph=nx.from_numpy_array(wMat)

    pathList,alphas=correlation_data_utilities.converge_subopt_paths_betweenness(
        wMatGraph,source=49,target=432,maxPaths=10000,verbose=True,giveAlphas=True)
    print(pathList)
    pathLengths=[calculatePathLength(wMatGraph,path) for path in pathList]
    print(pathLengths) #Got path lengths

    plt.plot(alphas)
    plt.loglog()
    plt.show()
    plt.plot(pathLengths)
    plt.show()

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
        shape=(454,454))
    datMat=np.abs(tempMat.todense())
    nzInds=np.nonzero(datMat)

    btwMat=np.array(correlation_data_utilities.getBtwMat(
        datMat,sources=sourceSet,targets=targetSet,verbose=False,verboseLevel=0,
        useProgressBar=True,pbarFun=tqdm.tqdm_notebook))
    btwMat

    #assigning weights and using -ln?
    btwWeightMat=copy.deepcopy(btwMat+btwMat.T)
    btwWeightMat=np.array(btwWeightMat)
    nzInds=np.nonzero(btwWeightMat)
    btwWeightMat[nzInds]=-np.log(btwWeightMat[nzInds])
    btwWeightGraph=nx.from_numpy_array(btwWeightMat)

    btwWeights=[btwWeightGraph.edges()[edge]['weight'] for edge in btwWeightGraph.edges()]
    sns.distplot(btwWeights)
    plt.show()
    print(np.min(btwWeights),np.mean(btwWeights),np.max(btwWeights))

    #Use 0 indexing when working with correlation_data_utilities and 1 indexing when working with MDAnalysis
    btwPathList,btwAlphas=correlation_data_utilities.converge_subopt_paths_betweenness(
        btwWeightGraph,source=49,target=432,maxPaths=10000,verbose=True,giveAlphas=True)
    print(btwPathList)
    btwPathLengths=[calculatePathLength(btwWeightGraph,path) for path in btwPathList]
    print(btwPathLengths)

    plt.figure(figsize=(6.4,4.8))
    plt.plot(alphas,label='weight=1/abs(E)')
    plt.plot(btwAlphas,label='weight=CurrentFlow(1/abs(E))')
    plt.axvline(100,color='teal',label='top 100 paths')
    plt.axhline(alphas[99],color='teal')
    plt.axhline(btwAlphas[99],color='teal')
    plt.loglog()
    plt.legend(loc="lower left")
    plt.xlabel('Number of Paths')
    plt.ylabel('Suboptimal Path Betweennes Convergence')
    plt.title("Converge of Suboptimal Path Usage")
    plt.grid()
    plt.show()

    plt.figure(figsize=(6.4,4.8))
    plt.plot(np.arange(len(pathLengths))+1,
            pathLengths/np.min(pathLengths),
            label='weight=1/abs(E)')
    plt.plot(np.arange(len(btwPathLengths))+1,
            btwPathLengths/np.min(btwPathLengths),
            label='weight=CurrentFlow(1/abs(E))')
    plt.axvline(100,color='teal',label='top 100 paths')
    plt.axhline(pathLengths[99]/pathLengths[0],color='teal')
    plt.axhline(btwPathLengths[99]/btwPathLengths[0],color='teal')
    plt.loglog()
    plt.legend(loc="upper left")
    plt.xlabel('Number of Paths')
    plt.ylabel('Relative Path Length')
    plt.title("Scaling of Relative Path Length")
    plt.grid()
    plt.show()

    path=btwPathList[0]
    pathNz=(path[1:],path[:-1])

    # Function to get the chain ID from the residue object
    def get_chain_id(residue):
        return residue.segid

    # Function to get the residue number from the residue object
    def get_residue_number(residue):
        return residue.resid

    #Shortest path!!!

    resNumMap=lambda x: x-253*(x>=253)+1
    resChainMap=lambda x: ['A','B'][int(x>=253)]
    print("Shortest path length (LEU_50 - GLU180):",btwPathLengths[0])
    print('Path (edge list):',[tuple(np.flip(edge,axis=0)) for edge in zip(*pathNz)])
    print('Path (Node Sequence Residue Names):',
        [struc.residues[iRes].resname + "_" + \
        get_chain_id(struc.residues[iRes]) + "_" + \
        str(get_residue_number(struc.residues[iRes])) \
        for iRes in np.flip(np.array(path),axis=0)])

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
                    cmap=tempCmap,label='Betweenness',
                    cNorm=cNorm,
                    orientation='horizontal')

    #Render the filtered interaction network using nglview with
    #the edge width and colormaps generated above
    struc=strucDict[list(strucDict.keys())[0]]
    view=nv.show_mdanalysis(struc)
    view.clear_representations()
    for res in sourceSet:
        view.add_representation('spacefill',selection=str(res+1)+' and .CA')
    for res in targetSet:
        view.add_representation('spacefill',selection=str(res+1)+' and .CA')
    view.add_representation('cartoon',selection='backbone',alpha=.5)
    correlation_data_utilities.drawProtCorrMat(protStruc=struc,corrMat=plotMat,ngViewOb=view,
                        frame=0,colorsArray=edgeColors,radiiMat=radiiMat,
                        undirected=True)
    buffer = StringIO()
    nv.write_html(buffer, [view]) 

    return buffer.getvalue()

if __name__ == '__main__':
    visualizeBetweenness()