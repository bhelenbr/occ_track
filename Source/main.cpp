//
//  main.cpp
//  occ_track
//
//  Created by Brian Helenbrook on 4/25/24.
//

#include "occ_track.h"
#include <unistd.h>

int main(int argc, char **argv) {
    
    enum run {straightProfile,straightTransitionProfile,curveProfile,curveTransitionProfile,
        straightTrackSpline2D,straightTrackEdges,LakePlacid,Sochi,ParkCity,ConvertProfiles,makeProfilesParkCity};
    
    run myRun;
    myRun = straightTrackSpline2D;
    myRun = Sochi;
    myRun = ConvertProfiles;
    myRun = straightTrackEdges;
    myRun = LakePlacid;
    myRun = ParkCity;
    myRun = straightProfile;
    myRun = straightTransitionProfile;
    myRun = curveProfile;
    myRun = curveTransitionProfile;
    myRun = straightTrackEdges;


    switch (myRun) {
        case(straightProfile): {
            makeStraightProfile("Straight/Results/straight.brep", 750., 750., 150., 150., 25., 25.);
            break;
        }
        case(straightTransitionProfile): {
            makeStraightTransitionProfile("Straight/Results/straightTransition.brep",11., 1008., 150., 672., 1798., 125., 600.);
            break;
        }
        case(curveProfile): {
            makeCurveProfile("Straight/Results/curve.brep",1350., 1182., 518., 13.8, 2601., 753., 1300., 600., 518., 125.);
            break;
        }
        case(curveTransitionProfile): {
            makeCurveTransitionProfile("Straight/Results/curveTransition.brep",472., 915., 73., 481., 1279., 598., 4129., 600., 473., 125.);
            break;
        }
        case(straightTrackSpline2D): {
            makeStraightSpline2D("Straight/Results");
            break;
        }
        case(straightTrackEdges): {
            makeTrackLoft("Straight","Straight/Results");
            break;
        }
        case(LakePlacid): {
            makeTrackSpline2D("TrackData_LKP", "TrackData_LKP/Results");
            break;
        }
        case(Sochi): {
            makeTrackSpline2D("TrackData_Sochi", "TrackData_Sochi/Results");
            break;
        }
        case(ParkCity): {
            makeTrackLoft("TrackData_ParkCity", "TrackData_ParkCity/Results");
            break;
        }
        case(ConvertProfiles): {
            convertProfiles2BRep("TrackData_LKP");
            break;
        }
        case(makeProfilesParkCity): {
            makeProfiles("TrackData_ParkCity/XSECTMOD.csv","TrackData_ParkCity/Results");
        }
        default:
            std::cout << "These aren't done yet\n";
    }
    
    return 0;
}
