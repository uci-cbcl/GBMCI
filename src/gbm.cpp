//------------------------------------------------------------------------------
//
//  GBM by Greg Ridgeway  Copyright (C) 2003
//  File:       gbm.cpp
//
//------------------------------------------------------------------------------

#include<vector>
#include "gbm.h"

#include <iostream>
using namespace std;

unsigned long gbm_setup
(
    double *adY,
    double *adOffset,
    double *adX,
    int *aiXOrder,
    double *adWeight,
    double *adMisc,
    int cRows,
    int cCols,
    int *acVarClasses,
    int *alMonotoneVar,
    const char *pszFamily,
    int cTrees,
    int cDepth,
    int cMinObsInNode,
    double dShrinkage,
    double dBagFraction,
    int cTrain,
    CDataset *pData,
    PCDistribution &pDist
)
{
    unsigned long hr = 0;

    hr = pData->SetData(adX,aiXOrder,adY,adOffset,adWeight,adMisc,
                        cRows,cCols,acVarClasses,alMonotoneVar);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }

    // set the distribution
    if(strncmp(pszFamily,"bernoulli",2) == 0)
    {
        pDist = new CBernoulli();
        if(pDist==NULL)
        {
            hr = GBM_OUTOFMEMORY;
            goto Error;
        }
    }
    else if(strncmp(pszFamily,"gaussian",2) == 0)
    {
        pDist = new CGaussian();
        if(pDist==NULL)
        {
            hr = GBM_OUTOFMEMORY;
            goto Error;
        }
    }
    else if(strncmp(pszFamily,"poisson",2) == 0)
    {
        pDist = new CPoisson();
        if(pDist==NULL)
        {
            hr = GBM_OUTOFMEMORY;
            goto Error;
        }
    }
    else if(strncmp(pszFamily,"adaboost",2) == 0)
    {
        pDist = new CAdaBoost();
        if(pDist==NULL)
        {
            hr = GBM_OUTOFMEMORY;
            goto Error;
        }
    }
    else if(strncmp(pszFamily,"coxph",2) == 0)
    {
        pDist = new CCoxPH();
        if(pDist==NULL)
        {
            hr = GBM_OUTOFMEMORY;
            goto Error;
        }
        //cout<<"OK coxph initialization!"<<endl;
    }
	else if(strncmp(pszFamily,"sci",2) == 0)
	{
        pDist = new CSCI();
        if(pDist==NULL)
        {
            hr = GBM_OUTOFMEMORY;
            goto Error;
        }
        //cout<<"OK sci initialization!"<<endl;
	}
    else if(strncmp(pszFamily,"laplace",2) == 0)
    {
        pDist = new CLaplace();
        if(pDist==NULL)
        {
            hr = GBM_OUTOFMEMORY;
            goto Error;
        }
    }
    else if(strncmp(pszFamily,"quantile",2) == 0)
    {
        pDist = new CQuantile(adMisc[0]);
        if(pDist==NULL)
        {
            hr = GBM_OUTOFMEMORY;
            goto Error;
        }
    }
    if(pDist==NULL)
    {
        hr = GBM_INVALIDARG;
        goto Error;
    }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}



GBMRESULT gbm_transfer_to_R
(
    CGBM *pGBM,
    VEC_VEC_CATEGORIES &vecSplitCodes,
    int *aiSplitVar,
    double *adSplitPoint,
    int *aiLeftNode,
    int *aiRightNode,
    int *aiMissingNode,
    double *adErrorReduction,
    double *adWeight,
    double *adPred,
    int cCatSplitsOld
)
{
    GBMRESULT hr = GBM_OK;

    hr = pGBM->TransferTreeToRList(aiSplitVar,
                                   adSplitPoint,
                                   aiLeftNode,
                                   aiRightNode,
                                   aiMissingNode,
                                   adErrorReduction,
                                   adWeight,
                                   adPred,
                                   vecSplitCodes,
                                   cCatSplitsOld);
    if(GBM_FAILED(hr)) goto Error;

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


GBMRESULT gbm_transfer_catsplits_to_R
(
    int iCatSplit,
    VEC_VEC_CATEGORIES &vecSplitCodes,
    int *aiSplitCodes
)
{
    unsigned long i=0;

    for(i=0; i<vecSplitCodes[iCatSplit].size(); i++)
    {
        aiSplitCodes[i] = vecSplitCodes[iCatSplit][i];
    }

    return GBM_OK;
}


int size_of_vector
(
    VEC_VEC_CATEGORIES &vec,
    int i
)
{
    return vec[i].size();
}



