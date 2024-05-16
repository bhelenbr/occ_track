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
    myRun = LakePlacid;
    myRun = ConvertProfiles;
    myRun = Sochi;
    myRun = makeProfilesParkCity;
    myRun = ParkCity;


    switch (myRun) {
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
