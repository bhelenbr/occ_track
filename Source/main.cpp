//
//  main.cpp
//  occ_track
//
//  Created by Brian Helenbrook on 4/25/24.
//

#include "occ_track.h"
#include <unistd.h>

int main(int argc, char **argv) {
    
    enum location {StraightTrack,LakePlacid,Sochi,ParkCity};
    enum operation {loft,BSpline,convertProfiles,createProfiles};
    
    location myLocation;
    myLocation = StraightTrack;
    myLocation = LakePlacid;
    myLocation = Sochi;
    myLocation = ParkCity;

    operation myOp;
    myOp = createProfiles;
    myOp = convertProfiles;
    myOp = loft;
    myOp = BSpline;

    double offset;  // Offset from defined surface to top of ice surface
    double scale; // scaling factor for profiles
    std::string Folder;

    switch(myLocation) {
        case(StraightTrack): {
            offset = -0.03;
            scale = 0.001;
            Folder = "Straight";
            break;
        }
        case(LakePlacid): {
            offset = -0.25/(4.26-2.85)*60.0/100.0;  // Taken from graphical measurements from curve 1.  Includes ice thickness
            scale = 0.01;
            Folder = "TrackData_LKP";
            break;
        }
        case(Sochi): {
            offset = -0.03;  // Taken from graphical measurements from curve 1.  Includes ice thickness
            scale = 0.001;
            Folder = "TrackData_Sochi";
            break;
        }
        case(ParkCity): {
            offset = -1.5*2.54/100;  // From Jon Owen
            scale = 0.001;
            Folder = "TrackData_ParkCity";
            break;
        }
    }
    
    switch (myOp) {
        case(loft): {
            makeTrackLoft(Folder, Folder +"/Results",scale);
            offsetSurface(Folder +"/Results/theTrack.brep",Folder +"/Results/withIce.brep", offset);
            break;
        }
        case(BSpline): {
            makeTrackBSpline(Folder, Folder +"/Results",scale);
            offsetSurface(Folder +"/Results/theTrack.brep",Folder +"/Results/withIce.brep", offset);
            break;
        }
        case(convertProfiles): {
            convertProfiles2BRep(Folder);
            break;
        }
        case(createProfiles): {
            makeProfiles(Folder +"/profiledata.csv",Folder);
        }
        default:
            std::cout << "Don't know how to do that" << std::endl;
    }
    
    return 0;
}
