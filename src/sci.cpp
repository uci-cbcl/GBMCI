#include "sci.h"

#include <iostream>
#include <cmath>
using namespace std;


CSCI::CSCI()
{
	alpha = 1.0;  // 1.0 (default), 2.0, 5.0, 10.0
}

CSCI::~CSCI()
{
}


GBMRESULT CSCI::ComputeWorkingResponse
(
    double *adT,
    double *adDelta,
    double *adOffset,
    double *adF, 
    double *adZ, 
    double *adWeight,
    bool *afInBag,
    unsigned long nTrain
)
{
    unsigned long i = 0;
    unsigned long j = 0;

    // compute total number of pairs
    unsigned long nPair = 0;
    unsigned long nRisk = 0;
    for (int i=0; i<nTrain; i++)
    {
    	if (!afInBag[i]) continue;

    	if (adDelta[i])
    		nPair += nRisk;
    	nRisk++;
    }
    // nPair = 1; // for debugging show 0.25 step! OK


    // compute gradient
    double tmp;
    for (i=0; i<nTrain; i++)
    {
    	adZ[i] = 0.0; // initialization

    	if (afInBag[i])
    	{
        	for (j=i+1; j<nTrain; j++)	// {(j,i)|t_j<t_i}
        	{
        		if (afInBag[j] && adDelta[j])
        		{
        			//tmp = exp(alpha*(adF[j]-adF[i]));
        			//adZ[i] += alpha*tmp/((1+tmp)*(1+tmp)*nPair);
        			tmp = 1/(1+exp(alpha*(adF[j]-adF[i])));
        			adZ[i] += tmp*(1-tmp);
        		}
        	}

        	if (adDelta[i])
        	{
            	for (j=0; j<i; j++)			// {(i,j)|t_i<t_j}
            	{
            		if (afInBag[j])
            		{
            			//tmp = exp(alpha*(adF[i]-adF[j]));
            			//adZ[i] -= alpha*tmp/((1+tmp)*(1+tmp)*nPair);
            			tmp = 1/(1+exp(alpha*(adF[i]-adF[j])));
            			adZ[i] -= tmp*(1-tmp);
            		}
            	}
        	}
    	}
    	adZ[i] *= alpha/nPair;
    	//Rprintf("%.10f ", adZ[i]);
    }
    //Rprintf("\n");

    return GBM_OK;
}



GBMRESULT CSCI::InitF
(
    double *adY,
    double *adMisc,
    double *adOffset, 
    double *adWeight,
    double &dInitF, 
    unsigned long cLength
)
{
    dInitF = 0.0;

    return GBM_OK;
}


double CSCI::Deviance
(
    double *adT,
    double *adDelta,
    double *adOffset, 
    double *adWeight,
    double *adF,
    unsigned long cLength
)
{
    unsigned long i=0;
    unsigned long j=0;

    // compute total number of pairs
    unsigned long nPair = 0;
    unsigned long nRisk = 0;
    for (i=0; i<cLength; i++)
    {
       	if (adDelta[i])
    		nPair += nRisk;
    	nRisk++;
    }

    // compute smoothed concordance index
    unsigned long nP_check = 0;
    double sci = 0.0;;
    for (i=0; i<cLength; i++)
    {
    	for (j=i+1; j<cLength; j++)	// {(j,i)|t_j<=t_i}
    	{
    		if (adT[j]>adT[i])
    			Rprintf("The survival time is not in decreasing order!!!\n");
    		if (!adDelta[j]) continue;
    		sci += 1/(1+exp(alpha*(adF[j]-adF[i])));
    		nP_check++;
    	}

    	// CYF: big bug fixed
    	/*if (!adDelta[i]) continue;

    	for (j=0; j<i; j++)			// {(i,j)|t_i<t_j}
    	{
    		sci += 1/(1+exp(alpha*(adF[i]-adF[j])));
    		nP_check++;
    	}*/
    }


    // compute concordance index
    double count = 0;
    double score = 0;
    for (i=0; i<cLength; i++)
    {
    	for (j=i+1; j<cLength; j++)	// actually, {(j,i)|t_j<=t_i}}
    	{
    		if ((adDelta[i]==1.0 && adDelta[j]==1.0) ||
    				(adDelta[i]==1.0 && adT[i]<=adT[j]) ||
    				(adDelta[j]==1.0 && adT[j]<=adT[i]))
    		{
    			count +=2; 		// This pair is comparable.
    			if (adT[i] < adT[j]) {
    				if (adF[i] < adF[j]) {
    					score = score + 2;
    				}
    				if (adF[i] == adF[j]) {
    					score = score + 1;
    				}
    			}
    			if (adT[i] == adT[j]) {
    				if (adF[i] == adF[j]) {
    					score = score + 2;
    				}
    			}
    			if (adT[i] > adT[j]) {
    				if (adF[i] > adF[j]) {
    					score = score + 2;
    				}
    				if (adF[i] == adF[j]) {
    					score = score + 1;
    				}
    			}
    		}

    	}
    }
    score /= count;

    //Rprintf("sci: %f, ci: %f , nPair: %d, count: %f\n", sci/nPair, score, nPair, count);
    //Rprintf("%f    %f\n", sci/nPair, score);
    //Rprintf("sci: %d  %f  %d %d\n", cLength, sci/nPair, nPair, nP_check); //sci/nPair

    //cout<<fixed
    return -sci/nPair;	// bug fixing: minimize energy
}



double CSCI::BagImprovement
(
    double *adT,
    double *adDelta,
    double *adOffset,
    double *adWeight,
    double *adF,
    double *adFadj,
    bool *afInBag,
    double dStepSize,
    unsigned long nTrain
)
{
    unsigned long i=0;
    unsigned long j=0;

    // compute total number of pairs
    unsigned long nPair = 0;
    unsigned long nRisk = 0;
    for (i=0; i<nTrain; i++)
    {
    	if (afInBag[i])	// bug fixing: only consider OOB
    		continue;

       	if (adDelta[i])
    		nPair += nRisk;
    	nRisk++;
    }

    // compute update of smoothed concordance index
    double dSci = 0.0;;
    for (i=0; i<nTrain; i++)
    {
    	if (afInBag[i])	// bug fixing: only consider OOB
    		continue;

    	for (j=i+1; j<nTrain; j++)	// {(j,i)|t_j<t_i}
    	{
        	if (afInBag[j])	// bug fixing: only consider OOB
        		continue;
    		if (!adDelta[j]) continue;
    		dSci += 1/(1+exp(alpha*(adF[j]+dStepSize*adFadj[j]-adF[i]-dStepSize*adFadj[i])))
    				- 1/(1+exp(alpha*(adF[j]-adF[i])));
    	}

    	if (!adDelta[i]) continue;

    	for (j=0; j<i; j++)			// {(i,j)|t_i<t_j}
    	{
        	if (afInBag[j])	// bug fixing: only consider OOB
        		continue;
    		dSci += 1/(1+exp(alpha*(adF[i]+dStepSize*adFadj[i]-adF[j]-dStepSize*adFadj[j])))
    				- 1/(1+exp(alpha*(adF[i]-adF[j])));
    	}
    }

    return dSci/nPair;
}


// 1D fitting /rho with line-search, all leaves will be updated

GBMRESULT CSCI::FitBestConstant
(
    double *adT,
    double *adDelta,
    double *adOffset,
    double *adW,
    double *adF,
    double *adZ,
    unsigned long *aiNodeAssign,
    unsigned long nTrain,
    VEC_P_NODETERMINAL vecpTermNodes,
    unsigned long cTermNodes,
    unsigned long cMinObsInNode,
    bool *afInBag,
    double *adFadj
)
{
	//cout<<"in SCI constant fitting"<<endl;
    GBMRESULT hr = GBM_OK;

    unsigned long nPair = 0;
    unsigned long nRisk = 0;
    for (int i=0; i<nTrain; i++)
    {
    	if (afInBag[i])
    	{
        	if (adDelta[i])
        		nPair += nRisk;
        	nRisk++;
    	}
    }


	double rho, sci;
	double rho_optimal, sci_optimal=-1;
	double f_i, f_j, w_ji;
	double rho_0=0.1, rho_max=100, k=3;	// line-search parameter 100, 1000, 10000
	int n = ceil(log(rho_max/rho_0)/log(k));
	for (int ii=n; ii>=0; ii--)
	{
		rho = rho_max/pow(k,ii);
	    sci = 0;
	    f_i=f_j=w_ji=0.0;
	    for (int i=0; i<nTrain; i++)
	    {
	    	if (afInBag[i])
	    	{
	        	for (int j=i+1; j<nTrain; j++)	// {(j,i)|t_j<=t_i}
	        	{
	        		if (afInBag[j] && adDelta[j])
	        		{
	            		f_i = vecpTermNodes[aiNodeAssign[i]]->dPrediction;
	            		f_j = vecpTermNodes[aiNodeAssign[j]]->dPrediction;
	            		w_ji = 1/(1+exp(alpha*(adF[j]-adF[i]+rho*f_j-rho*f_i)));
	            		sci += w_ji;;
	        		}
	        	}
	    	} // if (afInBag[i])
	    } // for (int i=0; i<nTrain; i++)
	    sci = sci/nPair;
	    //Rprintf("%f  ",sci); // the obj is monotonously increasing, if rho_max is not too big
	    if (sci>sci_optimal)
	    {
	    	rho_optimal = rho;
	    	sci_optimal = sci;
	    }
	}
	//Rprintf("%f    %f    ", rho_optimal, sci_optimal);

	for(int iNode=0; iNode<cTermNodes; iNode++)
	{
		if(vecpTermNodes[iNode]!=NULL)
		{
			//Rprintf("%f  ", vecpTermNodes[iNode]->dPrediction);
			vecpTermNodes[iNode]->dPrediction *= rho_optimal;
		}
	}
	//Rprintf("\n");


    return hr;
}



// 1D fitting /rho with line-search, k-1 leaves will be updated -----> paused, not as expected before
/*
GBMRESULT CSCI::FitBestConstant
(
    double *adT,
    double *adDelta,
    double *adOffset,
    double *adW,
    double *adF,
    double *adZ,
    unsigned long *aiNodeAssign,
    unsigned long nTrain,
    VEC_P_NODETERMINAL vecpTermNodes,
    unsigned long cTermNodes,
    unsigned long cMinObsInNode,
    bool *afInBag,
    double *adFadj
)
{
	//cout<<"in SCI constant fitting"<<endl;
    GBMRESULT hr = GBM_OK;

    unsigned long nPair = 0;
    unsigned long nRisk = 0;
    for (int i=0; i<nTrain; i++)
    {
    	if (afInBag[i])
    	{
        	if (adDelta[i])
        		nPair += nRisk;
        	nRisk++;
    	}
    }


    unsigned long K = 0;
    veciK2Node.resize(cTermNodes);
    for(int iNode=0; iNode<cTermNodes; iNode++)
    {
        if(vecpTermNodes[iNode]!=NULL && vecpTermNodes[iNode]->cN >= cMinObsInNode)
        {
            veciK2Node[K] = i;
            K++;
        }
    }


	double rho, sci;
	double rho_optimal, sci_optimal=-1;
	double f_i, f_j, w_ji;
	double rho_0=0.1, rho_max=100, k=3;	// line-search parameter 100, 1000, 10000
	int n = ceil(log(rho_max/rho_0)/log(k));
	for (int ii=n; ii>=0; ii--)
	{
		rho = rho_max/pow(k,ii);
	    sci = 0;
	    f_i=f_j=w_ji=0.0;
	    for (int i=0; i<nTrain; i++)
	    {
	    	if (afInBag[i])
	    	{
	        	for (int j=i+1; j<nTrain; j++)	// {(j,i)|t_j<=t_i}
	        	{
	        		if (afInBag[j] && adDelta[j])
	        		{
	            		f_i = vecpTermNodes[aiNodeAssign[i]]->dPrediction;
	            		f_j = vecpTermNodes[aiNodeAssign[j]]->dPrediction;
	            		w_ji = 1/(1+exp(alpha*(adF[j]-adF[i]+rho*f_j-rho*f_i)));
	            		sci += w_ji;
	        		}
	        	}
	    	} // if (afInBag[i])
	    } // for (int i=0; i<nTrain; i++)
	    sci = sci/nPair;
	    // Rprintf("%f  ",sci); // It's funny the obj is monotonously increasing
	    if (sci>sci_optimal)
	    {
	    	rho_optimal = rho;
	    	sci_optimal = sci;
	    }
	}
	//Rprintf("%f    %f    ", rho_optimal, sci_optimal);

    for(int k=0; k<K-1; k++)
    {
    	if (k<K-2)
    		vecpTermNodes[veciK2Node[k]]->dPrediction *= rho_optimal;
    	else
    		vecpTermNodes[veciK2Node[k]]->dPrediction = 0.0;	// keep the last entry to 0 !
    }


    return hr;
}
*/

GBMRESULT CSCI::FitBestConstant_Cuo
(
    double *adT,
    double *adDelta,
    double *adOffset,
    double *adW,
    double *adF,
    double *adZ,
    unsigned long *aiNodeAssign,
    unsigned long nTrain,
    VEC_P_NODETERMINAL vecpTermNodes,
    unsigned long cTermNodes,
    unsigned long cMinObsInNode,
    bool *afInBag,
    double *adFadj
)
{
	//cout<<"in SCI constant fitting"<<endl;
    GBMRESULT hr = GBM_OK;

    //return hr; // no stupid Newton's currently

    unsigned long nPair = 0;
    unsigned long nRisk = 0;
    for (int i=0; i<nTrain; i++)
    {
    	if (!afInBag[i]) continue;

    	if (adDelta[i])
    		nPair += nRisk;
    	nRisk++;
    }


    unsigned long K = 0;
    veciK2Node.resize(cTermNodes);
    veciNode2K.resize(cTermNodes);
    for(int i=0; i<cTermNodes; i++)
    {
        veciNode2K[i] = 0;
        if(vecpTermNodes[i]->cN >= cMinObsInNode)
        {
            veciK2Node[K] = i;
            veciNode2K[i] = K;
            K++;
        }
    }


    double w_ij = 0.0;
    matH.setactualsize(K-1);
    vecdG.resize(K-1);

    // gradient
    for (int k=0; k<K-1; k++)
    {
    	vecdG[k] = 0.0;

    	for (int j=0; j<nTrain; j++)
    	{
        	if (!afInBag[j] || (vecpTermNodes[aiNodeAssign[j]]->cN < cMinObsInNode)) continue;

        	for (int i=j+1; i<nTrain; i++)	// {(i,j)|t_i<t_j}
        	{
        		if (!afInBag[i] || (vecpTermNodes[aiNodeAssign[i]]->cN < cMinObsInNode)
        				|| !adDelta[i]) continue;
        		w_ij = 1/(1+exp(alpha*(adF[i]-adF[j])));
        		vecdG[k] += (-(aiNodeAssign[i]==veciK2Node[k]) + (aiNodeAssign[j]==veciK2Node[k]))*w_ij*(1-w_ij);
        	}
    	}
    	vecdG[k] *= alpha/nPair;
    }

    // Hessian
    double h_mn = 0.0;
    //unsigned long I_mn = 0;
    for (int m=0; m<K-1; m++)
    {
    	for (int n=0; n<=m; n++)
    	{
    		h_mn = 0.0;
    		//I_mn = 0;
        	for (int j=0; j<nTrain; j++)
        	{
            	if (!afInBag[j] || (vecpTermNodes[aiNodeAssign[j]]->cN < cMinObsInNode)) continue;

            	for (int i=j+1; i<nTrain; i++)	// {(i,j)|t_i<t_j}
            	{
            		if (!afInBag[i] || (vecpTermNodes[aiNodeAssign[i]]->cN < cMinObsInNode)
            				|| !adDelta[i]) continue;
            		w_ij = 1/(1+exp(alpha*(adF[i]-adF[j])));
            		h_mn += (-(aiNodeAssign[i]==veciK2Node[m])+(aiNodeAssign[j]==veciK2Node[m]))
            				* (-(aiNodeAssign[i]==veciK2Node[n])+(aiNodeAssign[j]==veciK2Node[n]))
            				* (1-2*w_ij)*w_ij*(1-w_ij);
            		//I_mn += (-(aiNodeAssign[i]==veciK2Node[m])+(aiNodeAssign[j]==veciK2Node[m]))
                    	//			* (-(aiNodeAssign[i]==veciK2Node[n])+(aiNodeAssign[j]==veciK2Node[n]));
            	}
        	}
        	//cout<<I_mn<<" ";
        	h_mn *= alpha*alpha/nPair;
        	matH.setvalue(m,n,h_mn);
        	matH.setvalue(n,m,h_mn);
    	}
    	//cout<<endl;
    }
    //cout<<"================="<<endl;



    double dTemp = 0.0;
    bool fTemp = false;

    // CYF
    /*
    for (int k=0; k<cTermNodes; k++)
    	cout<<veciK2Node[k]<<" ";
    cout<<endl<<"==============================="<<endl;
    */
    /*Rprintf("\n gradient:\n");
    for (int k=0; k<K-1; k++)
    	Rprintf("%f \n", vecdG[k]);
    Rprintf("----Hessian:----\n");
    for(int k=0; k<K-1; k++)
    {
        for(int m=0; m<K-1; m++)
        {
            matH.getvalue(k,m,dTemp,fTemp);
            Rprintf("%f ",dTemp);
        }
        Rprintf("\n");
    }*/


    // one step to get leaf predictions
    matH.invert();

    /*
    Rprintf("----Inverse Hessian:----\n");
    for(int k=0; k<K-1; k++)
    {
        for(int m=0; m<K-1; m++)
        {
            matH.getvalue(k,m,dTemp,fTemp);
            Rprintf("%f ",dTemp);
        }
        Rprintf("\n");
    }
    Rprintf("=============================\n");
    */


    for(int k=0; k<cTermNodes; k++)
    {
        vecpTermNodes[k]->dPrediction = 0.0;	// keep the last entry to 0 !
    }

    // Switch 1: try Newton's method first
    for(int m=0; m<K-1; m++)
    {
        for(int k=0; k<K-1; k++)
        {
            matH.getvalue(m,k,dTemp,fTemp);
            if(!R_FINITE(dTemp)) // occurs if matH was not invertible
            {
            	cerr<<"Hessian not invertable!"<<endl;
                vecpTermNodes[veciK2Node[m]]->dPrediction = 0.0;
                break;
            }
            else
            {
                vecpTermNodes[veciK2Node[m]]->dPrediction -= dTemp*vecdG[k];	// CYF bugging
            }
        }
    }

    // check if S.C.I. is improving with Newton's method
    for(int iObs=0; iObs<nTrain; iObs++)
    {
        adFadj[iObs] = vecpTermNodes[aiNodeAssign[iObs]]->dPrediction;
    }

    double dSci = 0.0;;
    for (int i=0; i<nTrain; i++)
    {
    	if (afInBag[i])
    	{
        	for (int j=i+1; j<nTrain; j++)	// {(j,i)|t_j<t_i}
        	{
            	if (afInBag[j] && adDelta[j])
            	{
            		dSci += 1/(1+exp(alpha*(adF[j]+adFadj[j]-adF[i]-adFadj[i])))
            				- 1/(1+exp(alpha*(adF[j]-adF[i])));
            	}
        	}
    	}
    }

    if (dSci>0)
    	return hr;

    // Switch 2: try line-search & gradient descend, if Newton's method doesn't work
    for(int m=0; m<K-1; m++)
    {
        vecpTermNodes[veciK2Node[m]]->dPrediction = vecdG[m];	// simple gradient ascend for debugging
    }

    return hr;
}


GBMRESULT CSCI::Compute_SCInew(double *adT,
		double *adDelta,
		double *adOffset,
		double *adW,
		double *adF,
		double *adZ,
		unsigned long *aiNodeAssign,
		unsigned long nTrain,
		VEC_P_NODETERMINAL vecpTermNodes,
		unsigned long cTermNodes,
		unsigned long cMinObsInNode,
		bool *afInBag,
		double *adFadj,
		unsigned long nPair,
		double rho, double &sci, double &dsci, double &ddsci)
{
    // compute S.C.I.new
    double f_i=0.0, f_j=0.0, w_ji=0.0;
    sci=dsci=ddsci=0.0;
    for (int i=0; i<nTrain; i++)
    {
    	if (afInBag[i])
    	{
        	for (int j=i+1; j<nTrain; j++)	// {(j,i)|t_j<=t_i}
        	{
        		if (adT[j]>adT[i])
        			Rprintf("The survival time is not in decreasing order!!!\n");
        		if (afInBag[j] && adDelta[j])
        		{
            		f_i = vecpTermNodes[aiNodeAssign[i]]->dPrediction;
            		f_j = vecpTermNodes[aiNodeAssign[j]]->dPrediction;
            		w_ji = 1/(1+exp(alpha*(adF[j]-adF[i]+rho*f_j-rho*f_i)));
            		sci += w_ji;;
            		dsci += (f_j-f_i)*w_ji*(1-w_ji);
            		ddsci += (f_j-f_i)*(f_j-f_i)*(1-2*w_ji)*w_ji*(1-w_ji);
        		}
        	}
    	} // if (afInBag[i])
    } // for (int i=0; i<nTrain; i++)

    sci = sci/nPair;
    dsci = -alpha*dsci/nPair;
    ddsci = alpha*alpha*ddsci/nPair;

    return GBM_OK;
}
