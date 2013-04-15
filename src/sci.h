//------------------------------------------------------------------------------
//  An expansion model to GBM by Yifei Chen  Copyright (C) 20012
//
//  File:       sci.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   Smoothed Concordance Index Object
//        	  
//  Owner:      yifeic@uci.edu
//
//  History:    11/25/2012   yifeic created
//
//------------------------------------------------------------------------------



#ifndef SCI_H
#define SCI_H

#include "distribution.h"
#include "matrix.h"

class CSCI : public CDistribution
{

public:

    CSCI();

    virtual ~CSCI();

    GBMRESULT ComputeWorkingResponse(double *adT,
                                   double *adDelta,
                                   double *adOffset,
                                   double *adF, 
                                   double *adZ, 
                                   double *adWeight,
                                   bool *afInBag,
                                   unsigned long nTrain);

    GBMRESULT InitF(double *adT,
                  double *adDelta,
                  double *adOffset,
                  double *adWeight,
                  double &dInitF, 
                  unsigned long cLength);

    GBMRESULT FitBestConstant_Cuo(double *adT,
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
                            double *adFadj);

    double Deviance(double *adT,
                    double *adDelta,
                    double *adOffset,
                    double *adWeight,
                    double *adF,
                    unsigned long cLength);

    double BagImprovement(double *adT,
                          double *adDelta,
                          double *adOffset,
                          double *adWeight,
                          double *adF,
                          double *adFadj,
                          bool *afInBag,
                          double dStepSize,
                          unsigned long nTrain);

    GBMRESULT FitBestConstant(double *adT,
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
                            double *adFadj);


private:
    double alpha;
    vector<double> vecdG;
    matrix<double> matH;

    vector<unsigned long> veciK2Node;
    vector<unsigned long> veciNode2K;

    GBMRESULT Compute_SCInew(double *adT,
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
    		double rho, double &sci, double &dsci, double &ddsci);
};

#endif // SCI_H

