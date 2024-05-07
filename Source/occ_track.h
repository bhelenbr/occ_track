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

/* Routine to make spline of track spine & output*/
void makeTrackSpine(std::string const sourceFolder, std::string const destinationFolder, Handle(Geom_BSplineCurve)& C);
void makeSpine(const std::vector<double>& arcLength, const std::vector<gp_Pnt>& points, Handle(Geom_BSplineCurve)& C);
void outputSpline1D(const std::string filename, int nPts, Handle(const Geom_BSplineCurve) const C);

/* Routines to make track surfaces of different types and output surface */
struct profileData {
    std::string filename;
    double sx;
};
void loadProfileData(std::string const filename, std::vector<profileData>& profiles);
void makeTrackLoft(std::string const sourceFolder, std::string const destinationFolder);
void makeTrackSpline2D(std::string const sourceFolder, std::string const destinationFolder);
void makeStraightSpline2D(std::string const destinationFolder);
void outputSpline2D(const std::string filename, int nPtsU, int nPtsV, Handle(Geom_BSplineSurface) const C);

/* Routines to make cross-sections profiles and manipulate*/
void makeProfiles(std::string const sourceFile, std::string const destinationFolder);
void makeStraightProfile(std::string const filename, std::string const Dir, double WL, double WR, double HL, double HR, double RL, double RR);
void makeStraightTransitionProfile(std::string const filename, std::string const Dir, double A, double B, double R1, double KH, double HZ, double R2, double BW, double Rr, double BH, double CE);
void makeCurveTransitionProfile(std::string const filename, std::string const Dir, double A, double B, double HZ, double R1, double KH, double BW, double R2, double BH, double B1, double Rr, double CE);
void makeCurveProfile(std::string const filename, std::string const Dir, double A, double B, double HZ, double WZ, double KH, double BW, double R2, double BH, double B1, double Rr, double CE);
void convertProfiles2BRep(std::string const sourceFolder);
void step2BRep(std::string const filename);


/* Routines to calculate trajectory down track */
void timeIntegrate(std::string destinationFolder, Handle(Geom_BSplineSurface const) const trackSurface);
void timeDeriv(std::vector<double> y, std::vector<double> dydt, Handle(Geom_BSplineSurface) trackSurface);



#endif /* defined(__occ_track__timeDeriv__) */
