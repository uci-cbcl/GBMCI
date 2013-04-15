//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "gbm_engine.h"


CGBM::CGBM()
{
    adFadj = NULL;
    adZ = NULL;
    afInBag = NULL;
    aiNodeAssign = NULL;
    aNodeSearch = NULL;
    
    cDepth = 0;
    cMinObsInNode = 0;
    dBagFraction = 0.0;
    dLambda = 0.0;
    fInitialized = false;
    cTotalInBag = 0;
    cTrain = 0;
    cValid = 0;

    pData = NULL;
    pDist = NULL;
    pNodeFactory = NULL;
    ptreeTemp = NULL;
}


CGBM::~CGBM()
{
    if(adFadj != NULL)
    {
        delete [] adFadj;
        adFadj = NULL;
    }
    if(adZ != NULL)
    {
        delete [] adZ;
        adZ = NULL;
    }
    if(afInBag != NULL)
    {
        delete [] afInBag;
        afInBag = NULL;
    }
    if(aiNodeAssign != NULL)
    {
        delete [] aiNodeAssign;
        aiNodeAssign = NULL;
    }
    if(aNodeSearch != NULL)
    {
        delete [] aNodeSearch;
        aNodeSearch = NULL;
    }
    if(ptreeTemp != NULL)
    {
        delete ptreeTemp;
        ptreeTemp = NULL;
    }
    // must delete the node factory last!!! at least after deleting trees
    if(pNodeFactory != NULL)
    {
        delete pNodeFactory;
        pNodeFactory = NULL;
    }
}


GBMRESULT CGBM::Initialize
(
    CDataset *pData,
    CDistribution *pDist,
    double dLambda,
    unsigned long cTrain,
    double dBagFraction,
    unsigned long cDepth,
    unsigned long cMinObsInNode
)
{
    GBMRESULT hr = GBM_OK;
    unsigned long i=0;

    if(pData == NULL)
    {
        hr = GBM_INVALIDARG;
        goto Error;
    }
    if(pDist == NULL)
    {
        hr = GBM_INVALIDARG;
        goto Error;
    }

    this->pData = pData;
    this->pDist = pDist;
    this->dLambda = dLambda;
    this->cTrain = cTrain;
    this->dBagFraction = dBagFraction;
    this->cDepth = cDepth;
    this->cMinObsInNode = cMinObsInNode;

    // allocate the tree structure
    ptreeTemp = new CCARTTree;
    if(ptreeTemp == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }

    cValid = pData->cRows - cTrain;
    cTotalInBag = (unsigned long)(dBagFraction*cTrain);

    adZ = new double[cTrain];
    if(adZ == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    adFadj = new double[pData->cRows];
    if(adFadj == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }

    pNodeFactory = new CNodeFactory();
    if(pNodeFactory == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    hr = pNodeFactory->Initialize(cDepth);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }
    ptreeTemp->Initialize(pNodeFactory);

    // array for flagging those observations in the bag
    afInBag = new bool[cTrain];
    if(afInBag==NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    // aiNodeAssign tracks to which node each training obs belongs
    aiNodeAssign = new ULONG[cTrain];
    if(aiNodeAssign==NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    // NodeSearch objects help decide which nodes to split
    aNodeSearch = new CNodeSearch[2*cDepth+1];
    if(aNodeSearch==NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    for(i=0; i<2*cDepth+1; i++)
    {
        aNodeSearch[i].Initialize(cMinObsInNode);
    }
    vecpTermNodes.resize(2*cDepth+1,NULL);

    fInitialized = true;

Cleanup:
    return hr;
Error:
    goto Cleanup;
}




GBMRESULT CGBM::Predict
(
    unsigned long iVar,
    unsigned long cTrees,
    double *adF,
    double *adX,
    unsigned long cLength
)
{
    GBMRESULT hr = GBM_OK;


    return hr;
}


GBMRESULT CGBM::Predict
(
    double *adX,
    unsigned long cRow,
    unsigned long cCol,
    unsigned long cTrees,
    double *adF
)
{
    GBMRESULT hr = GBM_OK;


    return hr;
}



GBMRESULT CGBM::GetVarRelativeInfluence
(
    double *adRelInf,
    unsigned long cTrees
)
{
    GBMRESULT hr = GBM_OK;
    int iVar=0;

    for(iVar=0; iVar<pData->cCols; iVar++)
    {
        adRelInf[iVar] = 0.0;
    }

    return hr;
}


GBMRESULT CGBM::PrintTree()
{
    GBMRESULT hr = GBM_OK;

    hr = ptreeTemp->Print();
    if(GBM_FAILED(hr)) goto Error;

Cleanup:
    return hr;
Error:
    goto Cleanup;
}




GBMRESULT CGBM::iterate
(
    double *adF,
    double &dTrainError,
    double &dValidError,
    double &dOOBagImprove,
    int &cNodes
)
{
    GBMRESULT hr = GBM_OK;
    unsigned long i = 0;
    unsigned long cBagged = 0;

    if(!fInitialized)
    {
        hr = GBM_FAIL;
        goto Error;
    }

    dTrainError = 0.0;
    dValidError = 0.0;
    dOOBagImprove = 0.0;

    vecpTermNodes.assign(2*cDepth+1,NULL);

    // randomly assign observations to the Bag
    cBagged = 0;
    for(i=0; i<cTrain; i++)
    {
        if(unif_rand()*(cTrain-i) < cTotalInBag-cBagged)
        {
            afInBag[i] = true;
            cBagged++;
        }
        else
        {
            afInBag[i] = false;
        }
    }

    #ifdef NOISY_DEBUG
    Rprintf("Compute working response\n");
    #endif

    // CYF
    /*Rprintf("# training data: %d \n", cTrain);
    Rprintf("one iteration \n");
    for(i=0; i < cTrain; i++)
    {
    	Rprintf("%d  ", afInBag[i]); //adF[i]
    }
    Rprintf("\n");*/

    hr = pDist->ComputeWorkingResponse(pData->adY, 
                                       pData->adMisc,
                                       pData->adOffset,
                                       adF, 
                                       adZ,
                                       pData->adWeight,
                                       afInBag,
                                       cTrain);

    // CYF
    /*for(i=0; i < cTrain; i++)
    {
    	Rprintf("%f\n", adZ[i]); //adF[i]
    }*/
    //Rprintf("\n");



    if(GBM_FAILED(hr))
    {
        goto Error;
    }

    #ifdef NOISY_DEBUG
    Rprintf("Reset tree\n");
    #endif
    hr = ptreeTemp->Reset();
    #ifdef NOISY_DEBUG
    Rprintf("grow tree\n");
    #endif
    hr = ptreeTemp->grow(adZ,pData,pData->adWeight,adFadj,
                         cTrain,cTotalInBag,dLambda,cDepth,
                         cMinObsInNode,
                         afInBag,
                         aiNodeAssign,aNodeSearch,vecpTermNodes);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }

    #ifdef NOISY_DEBUG
    Rprintf("get node count\n");
    #endif
    hr = ptreeTemp->GetNodeCount(cNodes);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }


    // CYF: ??!! only node 0 1 3 5 has data
    //Rprintf("# training data: %d \n", cTrain);
    /*Rprintf("# total nodes: %d \n", cNodes);
    for(i=0; i < cTrain; i++)
    {
    	Rprintf("%d, ", aiNodeAssign[i]); //adF[i]
    }
    Rprintf("\n");*/

    // Now I have adF, adZ, and vecpTermNodes (new node assignments)
    // Fit the best constant within each terminal node
    #ifdef NOISY_DEBUG
    Rprintf("fit best constant\n");
    #endif
    hr = pDist->FitBestConstant(pData->adY,
                                pData->adMisc,
                                pData->adOffset,
                                pData->adWeight,
                                adF,
                                adZ,
                                aiNodeAssign,
                                cTrain,
                                vecpTermNodes,
                                (2*cNodes+1)/3, // number of terminal nodes
                                cMinObsInNode,
                                afInBag,
                                adFadj);
    //Rprintf("FitBestConstant finished \n");

    if(GBM_FAILED(hr))
    {
        goto Error;
    }

    // CYF
    // Rprintf("2*cDepth+1, (2*cNodes+1)/3: %d %d \n", 2*cDepth+1, (2*cNodes+1)/3);

    // update training predictions
    // fill in missing nodes where N < cMinObsInNode
    hr = ptreeTemp->Adjust(aiNodeAssign,adFadj,cTrain,
                           vecpTermNodes,cMinObsInNode);
/*
    for (i=0; i<cTrain; i++)
    {
    	Rprintf("%f, ", adFadj[i]);
    }
    Rprintf("\n");
*/
    //Rprintf("Tree adjust finished \n");


    /*Rprintf("# total nodes: %d \n", cNodes);
    for(i=0; i < cTrain; i++)
    {
    	Rprintf("%d, ", aiNodeAssign[i]); //  adF[i]
    }
    Rprintf("\n");*/


    if(GBM_FAILED(hr))
    {
        goto Error;
    }
    ptreeTemp->SetShrinkage(dLambda);

    dOOBagImprove = pDist->BagImprovement(pData->adY,
                                          pData->adMisc,
                                          pData->adOffset,
                                          pData->adWeight,
                                          adF,
                                          adFadj,
                                          afInBag,
                                          dLambda,
                                          cTrain);

    //Rprintf("BagImprovement finished \n");

    // update the training predictions
    for(i=0; i < cTrain; i++)
    {
        adF[i] += dLambda*adFadj[i];
    }
    dTrainError = pDist->Deviance(pData->adY,
                                  pData->adMisc,
                                  pData->adOffset,
                                  pData->adWeight,
                                  adF,
                                  cTrain);

    // update the validation predictions
    hr = ptreeTemp->PredictValid(pData,cValid,adFadj);
    for(i=cTrain; i < cTrain+cValid; i++)
    {
        adF[i] += adFadj[i];
    }
    if(pData->fHasOffset)
    {
        dValidError =
            pDist->Deviance(&(pData->adY[cTrain]),
                            &(pData->adMisc[cTrain]),
                            &(pData->adOffset[cTrain]),
                            &(pData->adWeight[cTrain]),
                            &(adF[cTrain]),
                            cValid);
    	//dValidError = 0; // temp
    }
    else
    {
        dValidError = pDist->Deviance(&(pData->adY[cTrain]),
                                      &(pData->adMisc[cTrain]),
                                      NULL,
                                      &(pData->adWeight[cTrain]),
                                      &(adF[cTrain]),
                                      cValid);
    	//dValidError = 0;
    }
    //Rprintf("%f  ", dValidError);

    //Rprintf("Deviance computing finished \n");

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


GBMRESULT CGBM::TransferTreeToRList
(
    int *aiSplitVar,
    double *adSplitPoint,
    int *aiLeftNode,
    int *aiRightNode,
    int *aiMissingNode,
    double *adErrorReduction,
    double *adWeight,
    double *adPred,    
    VEC_VEC_CATEGORIES &vecSplitCodes,
    int cCatSplitsOld
)
{
    GBMRESULT hr = GBM_OK;

    hr = ptreeTemp->TransferTreeToRList(pData,
                                        aiSplitVar,
                                        adSplitPoint,
                                        aiLeftNode,
                                        aiRightNode,
                                        aiMissingNode,
                                        adErrorReduction,
                                        adWeight,
                                        adPred,
                                        vecSplitCodes,
                                        cCatSplitsOld,
                                        dLambda);

    return hr;
}


