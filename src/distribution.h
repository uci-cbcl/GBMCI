//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       distribution.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   distribution object
//        	  
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include "node_terminal.h"

class CDistribution
{

public:

    CDistribution();
    virtual ~CDistribution();

    virtual GBMRESULT ComputeWorkingResponse(double *adY,
                                           double *adMisc,
                                           double *adOffset,
                                           double *adF, 
                                           double *adZ,
                                           double *adWeight,
                                           bool *afInBag,
                                           unsigned long nTrain) = 0;

    virtual GBMRESULT InitF(double *adY,
                          double *adMisc,
                          double *adOffset,
                          double *adWeight,
                          double &dInitF, 
                          unsigned long cLength) = 0;

    virtual double Deviance(double *adY,
                            double *adMisc,
                            double *adOffset,
                            double *adWeight,
                            double *adF,
                            unsigned long cLength) = 0;

    virtual GBMRESULT FitBestConstant(double *adY,
                                    double *adMisc,
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
                                    double *adFadj) = 0;

    virtual double BagImprovement(double *adY,
                                  double *adMisc,
                                  double *adOffset,
                                  double *adWeight,
                                  double *adF,
                                  double *adFadj,
                                  bool *afInBag,
                                  double dStepSize,
                                  unsigned long nTrain) = 0;
};
typedef CDistribution *PCDistribution;
#endif // DISTRIBUTION_H



