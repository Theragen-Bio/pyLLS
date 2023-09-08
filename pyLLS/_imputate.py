'''
This script contains functions to estimate missing values by pyLLS.
pyLLS is composed of three steps.
In the first step, it calculates pair-wise similarity between missing genes (targets) and other genes (probes).
The similarity can be obtained by following metrices;
a) Pearson Correlation Coefficients (PCCs, default)
b) L1-norm (Manhattan distance)
c) L2-norm (Euclidean distance)
Next, it finds optimal k for each target by Kneedle algorithm.
Finally, missing values are estimated for each target using individual linear regression models modeled with k features.
'''
import pandas as _pd
import numpy as _np
from tqdm import tqdm as _tqdm
import warnings
import os
import sys

from sklearn.metrics import pairwise_distances as _pairwise_distances
from sklearn.metrics import r2_score as _r2_score 
from kneed import KneeLocator as _KneeLocator


#----------------------------------------------------
# Similarity Calculation
# _pairwise_correlation : Pearson correlation coefficient
# _L1norm : L1-norm (Manhattan distance)
# _L2norm : L2-nrom (Euclidean distance)
#----------------------------------------------------
def _pairwise_correlation(A, B,n_jobs=-1):
    '''
    Calculate Pearsons Correlation Coefficient (PCC) between two genes.
    A : missing gene matrix
    B : reference gene matrix
    n_jobs = The number of jobs to use for the computation. 'all' means using all processors.
    n_jobs is ignored if neither A nor B has np.nan. (i.e single-core process)
    return MG x Probes matrix containing 1-pearson correlation
    '''
    n_jobs = n_jobs if type(n_jobs)==int else -1
    # Check NANs
    if _np.isnan(A).sum()+_np.isnan(B).sum()==0 and n_jobs==-1:
        A, B = A.T, B.T # Transpose matrixes
        am = A - _np.mean(A, axis=0, keepdims=True)
        bm = B - _np.mean(B, axis=0, keepdims=True)
        return 1 - abs(am.T @ bm /  (_np.sqrt(_np.sum(am**2, axis=0,keepdims=True)).T * _np.sqrt(_np.sum(bm**2, axis=0, keepdims=True))))
    else:
        # Pairwise_distances function of scipy returns 1-Pearson correlation distance
        # therefore, we need to recalculate 1-PCC for reverse direction.
        cor_dist=_pairwise_distances(X=A,Y=B,metric='correlation',n_jobs=n_jobs,force_all_finite=False)
        return _np.minimum(cor_dist,2-cor_dist)

def _L1norm(A, B, n_jobs='all'):
    '''
    Calculate Manhattan distance (L1-norm) between two genes 
    A : missing gene matrix
    B : reference gene matrix
    n_jobs = The number of jobs to use for the computation. 'all' means using all processors.
    return MG x Probes matrix containing L1-distance
    '''
    n_jobs = n_jobs if type(n_jobs)==int else -1
    return _pairwise_distances(X=A,Y=B,metric='l1',n_jobs=n_jobs,force_all_finite=False)

def _L2norm(A, B,n_jobs='all'):
    '''
    Calculate Euclidean distance (L2-norm) between two genes 
    n_jobs = The number of jobs to use for the computation. 'all' means using all processors.
    return MG x Probes matrix containing L2-distance
    '''
    n_jobs = n_jobs if type(n_jobs)==int else -1
    return _pairwise_distances(X=A,Y=B,metric='l2',n_jobs=n_jobs,force_all_finite=False)

#----------------------------------------------------
# Select neighbors
# _identify_neighbors : Return a table containing k probes for each target feature.
#----------------------------------------------------
def _identify_neighbors(ref=None,missing_genes=None,maxK=100,metric='correlation',n_jobs='all'):
    '''
    This function sorts the specified number of k probes correlative with each missing gene.
    # parameter
    ref : reference pd.DataFrame. Index is gene-name and Column is sample.
    missing_genes : list of missing genes (ex. missing_genes=['GeneA','GeneB','GeneC'])
    maxK : The number of probes to be evaluated.
    metric : Metric to calculate distance between missing gene and probes.
    Available options are ['correlation','L1','L2']. 'correlation' is default.
    # return
    pd.DataFrame consists of two columns ('Missing Gene' and 'Candidate Probes').
    Each candidate probe for an individual target is listed according to similarity (from highest to lowest).
    '''
    # Split matrix 
    probes = list(set(ref.index)-set(missing_genes))
    probe_mat, mg_mat = ref.loc[probes,:], ref.loc[missing_genes,:] #probe_mat : probe, #mg_mat : missing gene
    # Calculate distance
    metric=metric.lower()
    if metric in ['correlation','cor']:
        print('Calculating 1-Pearson Correlation distance')
        distMat = _pairwise_correlation(A=mg_mat.values, B=probe_mat.values,n_jobs=n_jobs)
    elif metric in ['l1','manhattan']:
        print('Calculating L1-norm similarity')
        distMat = _L1norm(A=mg_mat.values, B=probe_mat.values,n_jobs=n_jobs)
    elif metric in ['l2','euclidean']:
        print('Calculating L2-norm similarity')
        distMat = _L2norm(A=mg_mat.values, B=probe_mat.values,n_jobs=n_jobs)
    # Assigning candidate probes to the target genes individually.
    result=_pd.DataFrame([])
    result['Missing Gene']=missing_genes
    result['Candidate Probes']=''
    print(f'Selecting k = {maxK} probes for each missing genes')
    for row in range(result.shape[0]):
        idx=_np.argsort(distMat[row,:])[:maxK]
        result.loc[row,'Candidate Probes']=','.join([probes[i] for i in idx])
    print('Selection is finished')
    return result

#----------------------------------------------------
# Imputation
# imputating missing values of missing entries.
# a : reference table
# b : target table with missing values.
# probe : neighbor genes list
# gene : missing gene
#----------------------------------------------------
def _impute(a,b,probe,gene):
    Apart = a.loc[:,probe].values
    Bpart = b.loc[:,probe].values
    X = _np.linalg.pinv(Apart) @ a.loc[:,gene]
    return Bpart @ X

#----------------------------------------------------
# Select best probes
# _select_best_probes : Return a table containing k probes for each target feature after probe selection by Kneedle algorithm.
#----------------------------------------------------
def _select_best_probes(ref=None,mgcp=None,r2_cut=None,addK=1):
    '''
    This function selects k probes using the Kneedle algorithm.
    The increase of k may cause overfitting in linear regression.
    Therefore, we implemented kneedle algorithm to locate K where accuracy is not significantly increased.
    # parameters
    ref = reference table
    mgcp = missing gene candidate probe table obtained by _identify_neighbors()
    r2_cut = R-squared cutoff (0 < r2_cut < 1). Default is None.
             if 0.5 is set, then the minimum number of neighbors to achieve 0.5 will be selected.
             if there is no combination above the cutoff, Kneedle's K will be used.
    addK = intenger that added to Kneedle's K to prevent underfit. This will use K+addK probes to estimate missing values of a gene.
    # return
    pd.DataFrame which is similar to output from _identify_neighbors()
    '''
    # Start
    mgcp_fin=mgcp.copy()
    print('Finding optimal K')
    ref_T=ref.T
    for idx in _tqdm(range(mgcp.shape[0])): #idx=0
        mg = mgcp.loc[idx,'Missing Gene']
        y_answer = ref_T.loc[:,mg].values
        probe=mgcp.loc[idx,'Candidate Probes'].split(',')
        scores=[]
        for k in range(1,len(probe)): # calculate R-squared score with linear regression models with k=1 to k=maxK probes.
            x=ref_T.loc[:,probe[:k]]
            y_pred = _impute(a=ref_T,b=x,probe=probe[:k],gene=mg)
            scores.append(_r2_score(y_answer,y_pred))
        if scores[0]!=1:
            Kneedle_K = _KneeLocator(range(1,len(probe)),scores, S=1.0, curve="concave", direction="increasing").knee
            try:
                bestK = min([Kneedle_K+addK,len(scores)]) # Kneedle_K+addK or maxK
            except:
                bestK=len(probe)-1 # If there is no kneedle, then assign maxK
        else:
            bestK=1 # If R^2==1 at k=1, then return bestK=1
        if r2_cut is not None and scores[bestK-1]<r2_cut and any(_np.array(scores)>=r2_cut):
            bestK=_np.where(_np.array(scores)>=r2_cut)[0][0]+1
        mgcp_fin.loc[idx,'Candidate Probes']=','.join(probe[:bestK])
        mgcp_fin.loc[idx,'R-square']=scores[bestK-1]
    mgcp_fin=mgcp_fin.rename(columns={'Candidate Probes':'Final Probes'})
    return mgcp_fin

#----------------------------------------------------
# Imputation
# impute_missing_gene : This function is not hidden function.
#----------------------------------------------------
def impute_missing_gene(
    ref=None,target=None,metric='correlation',maxK=100,r2_cut=None,useKneedle=True,addK=1,verbose=True,n_jobs='all',
    return_probes=False
):
    if ref is None:
        print(
            '''
            This function estimates missing values of the specified target probes.
            # parameters
            ref (pd.DataFrame): reference data. gene x sample (n x p) DataFrame.
            target (pd.DataFrame) : Target table containing missing values. gene x sample (i x k) DataFrame.
            metric (str) : ['correlation'(default),'L1','L2']
                           Similarity metric to prioritize the probes for each target.
            maxK : maximum number of probe to be used in missing value estimation.
            useKneedle : It determines whether Kneedle algorithm should be used (True) or not (False).
                         If useKneedle==False, then maxK probes will be used to estimate missing values.
            verbose : If True, progress is reported. Otherwise, no progress is reported.
            n_jobs : Use all threads ('all') or speicified number of threads (int)
            addK (int) = Intenger that will be added to Kneedle's K to prevent underfitting.
                   This function will use K+addK probes to estimate missing values of a gene. (default = 1)
            return_probes = if true, 'target-table and mgcp' will be returned else 'target' will be returned.
            # Return
            * target : table with estimated values of missing genes that are not present in original target table.
            matrix shape will be (n x k).
            * mgcp : missing gene correlative probes. If useKneedle == True, mgcp will have R2-square column.
            # tutorial
            import pandas as pd
            import numpy as np
            import random
            tmp=pd.DataFrame(np.array(random.sample(range(1000),1000)).reshape(100,10))
            tmp.index=['g'+str(i) for i in tmp.index]
            tmp.columns=['s'+str(i) for i in tmp.columns]
            tmp2=tmp.iloc[:90,:5]
            tmp3=pyLLS.impute_missing_gene(ref=tmp,target=tmp2)
            ''')
        return
    if verbose==False:
        sys.stdout=open('/dev/null','w')
    # Idenfitying missing genes
    missing_genes=list(set(ref.index)-set(target.index))
    # Finding correlative neighbors for each target
    mgcp=_identify_neighbors(ref=ref,missing_genes=missing_genes,maxK=maxK,metric=metric,n_jobs=n_jobs)
    # Finding best K neighbors for each target if useKneedle is True.
    if useKneedle: # use kneedle algorithm
        mgcp=_select_best_probes(ref=ref,mgcp=mgcp,r2_cut=r2_cut,addK=addK)
    else: # use maxK probes.
        mgcp=mgcp.rename(columns={'Candidate Probes':'Final Probes'})
    # Estimating missing value in target-table.
    pred=_pd.DataFrame([])
    warnings.simplefilter('ignore')
    target_T=target.T
    ref_T=ref.T
    for idx,gene in enumerate(mgcp['Missing Gene']): #idx,gene=0,mgcp.loc[0,'Missing Gene']
        try:
            probes = mgcp.loc[idx,'Final Probes'].split(',')
            pred[gene] = _impute(a=ref_T,b=target_T,probe=probes,gene=gene)
        except:
            if verbose==False:
                sys.stdout=sys.__stdout__
            print('Error occured. Please check outputs and whether there are duplicated genes in the target table.')
            return (gene, probes), mgcp
    pred.index=target.columns
    target_fin=_pd.concat([target,pred.T],axis=0)
    print('Imputation finished. New target table shape is',target_fin.shape)
    if verbose==False:
        sys.stdout=sys.__stdout__
    warnings.resetwarnings()
    if return_probes:
        return target_fin,mgcp
    else:
        return target_fin