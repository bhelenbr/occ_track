//
//  timeDeriv.h
//  occ_track
//
//  Created by Brian Helenbrook on 10/27/13.
//  Copyright (c) 2013 Brian Helenbrook. All rights reserved.
//

#ifndef __occ_track__
#define __occ_track__

#include <vector>
#include <Geom_BSplineSurface.hxx>
#include <Geom_BSplineCurve.hxx>
#include <iostream>

/* Routine to make spline of track spine */
void makeSpine(const std::vector<double>& arcLength, const std::vector<gp_Pnt>& points, Handle(Geom_BSplineCurve)& C);
void makeTrackSpline2D(std::string const sourceFolder, std::string const destinationFolder);
void makeTrackLoft(std::string const sourceFolder, std::string const destinationFolder);
void makeStraightSpline2D(std::string const destinationFolder);
void outputSpline1D(const std::string filename, int nPts, Handle(const Geom_BSplineCurve) const C);
void outputSpline2D(const std::string filename, int nPtsU, int nPtsV, Handle(Geom_BSplineSurface) const C);
void makeStraightProfile(std::string const filename,double WL, double WR, double HL, double HR, double RL, double RR);
void convertProfiles2BRep(std::string const sourceFolder);
void step2BRep(std::string const filename);
void timeIntegrate(std::string destinationFolder, Handle(Geom_BSplineSurface const) const trackSurface);
void timeDeriv(std::vector<double> y, std::vector<double> dydt, Handle(Geom_BSplineSurface) trackSurface);


#endif /* defined(__occ_track__timeDeriv__) */
