//
//  main.cpp
//  occ_track
//
//  Created by Brian Helenbrook on 4/25/24.
//

#include "occ_track.h"
#include <unistd.h>

int main(int argc, char **argv) {
    
    enum run {straightProfile,straightTrack,LakePlacid,Sochi,ParkCity,ConvertProfiles};
    
    run myRun;
    // myRun = straightProfile;
    // myRun = straightTrack;
    // myRun = LakePlacid;
    // myRun = Sochi;
    // myRun = ParkCity;
    myRun = ConvertProfiles;
    
    switch (myRun) {
        case(straightProfile): {
            makeStraightProfile("Straight/Results/straight.brep",750., 750., 150., 150., 25., 25.);
            break;
        }
        case(straightTrack): {
            makeStraightSpline2D("Straight/Results");
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
            //makeTrackSpline2D("TrackData_ParkCity", "TrackData_ParkCity/Results");
            makeTrackLoft("TrackData_ParkCity", "TrackData_ParkCity/Results");
            break;
        }
        case(ConvertProfiles): {
            convertProfiles2BRep("TrackData_LKP");
            break;
        }
        default:
            std::cout << "These aren't done yet\n";
    }
    
    return 0;
}
