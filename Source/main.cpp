//
//  main.cpp
//  occ_track
//
//  Created by Brian Helenbrook on 4/25/24.
//

#include "occ_track.h"
#include <unistd.h>

int main(int argc, char **argv) {
    
    enum run {straightTrackSpline2D,straightTrackEdges,LakePlacid,Sochi,ParkCity,ConvertProfiles,makeProfilesParkCity};
    
    run myRun;
    myRun = straightTrackSpline2D;
    myRun = straightTrackEdges;
    myRun = straightTrackEdges;
    myRun = ConvertProfiles;
    myRun = Sochi;
    myRun = makeProfilesParkCity;
    myRun = ParkCity;
    myRun = LakePlacid;


    switch (myRun) {
        case(straightTrackSpline2D): {
            makeStraightSpline2D("Straight/Results");
            break;
        }
        case(straightTrackEdges): {
            makeTrackLoft("Straight","Straight/Results",0.001);
            break;
        }
        case(LakePlacid): {
            makeTrackSpline2D("TrackData_LKP", "TrackData_LKP/Results",0.01);
            //makeTrackLoft("TrackData_LKP", "TrackData_LKP/Results",0.01);
            const double offset = 0.25/(4.26-2.85)*60.0/100.0;  // Taken from graphical measurements from curve 1.  Includes ice thickness
            offsetSurface("TrackData_LKP/Results/theTrack.brep","TrackData_LKP/Results/withIce.brep", -offset);
            break;
        }
        case(Sochi): {
            makeTrackSpline2D("TrackData_Sochi", "TrackData_Sochi/Results",0.001);
            break;
        }
        case(ParkCity): {
            makeTrackLoft("TrackData_ParkCity", "TrackData_ParkCity/Results",0.001);
            break;
        }
        case(ConvertProfiles): {
            convertProfiles2BRep("TrackData_Sochi");
//            convertProfiles2BRep("TrackData_LKP");
            break;
        }
        case(makeProfilesParkCity): {
            makeProfiles("TrackData_ParkCity/profiledata.csv","TrackData_ParkCity");
            break;
        }
        default:
            std::cout << "These aren't done yet\n";
    }
    
    return 0;
}
