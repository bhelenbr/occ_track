//
//  main.cpp
//  occ_track
//
//  Created by Brian Helenbrook on 4/25/24.
//

#include "occ_track.h"
#include "ExportSTL.hpp"
#include <unistd.h>
#include <filesystem>
namespace fs = std::filesystem;

int main(int argc, char **argv) {
    
    enum location {StraightTrack,LakePlacid,Sochi,ParkCity,ParkCityNoEntries,ParkCityLugeStart};
    enum operation {loft,BSpline,convertProfiles,createProfiles};
    
    location myLocation;
    myLocation = StraightTrack;
    myLocation = LakePlacid;
    myLocation = Sochi;
    myLocation = ParkCity;
    myLocation = ParkCityLugeStart;
    myLocation = ParkCityNoEntries;


    operation myOp;
    myOp = loft;
    myOp = BSpline;

    double offset;  // Offset from defined surface to top of ice surface
    double scale; // scaling factor for profiles
    std::string Folder;
    
    /* to update include file list cd /opt/homebrew/include/opencascade and ls */
    /* and library names cd /opt/homebrew/lib and ls libTK* | grep -v 7.9 */
    
    /* The makeTrack<Loft,BSpline> routines calls makeTrackSpine to load the file <source folder>/track.txt and use it to make the spine of the track using the routine makeSpine.  makeSpineTrack then outputs points on the spine curve to the file "<destination folder>/track_spine_pts.txt" makeTrackXXX then opens the profile brep curves which are at the location specified in the file <source folder>/profiles.txt.  It outputs two files:  <destinationFolder>/theTrack.brep and <destinationFolder>/track_surface_pts.txt.  The routine offset surface loads the track file and then offsets a surface and outputs a brep file "withice.brep" and a surface definition withice.txt".  It finally calls timeIntegrate to calculate a gravity path, but it doesn't let V change do that is kind of funky.  This ouputs a file called <destination folder>/path.txt */
    
    /* The makeProfile routine reads the input file (typically <Folder>/provile.csv and for every line in that file creates a profile called Profiles/XX.brep where XX is the arc location with an ending describing the type of profile and whether it is left or right. It then outputs a file <Destination Folder>/profiles.txt that contains a list of locations and filenames  */
    
    switch(myLocation) {
        case(StraightTrack): {
            offset = -0.03;
            scale = 0.001;
            Folder = "Straight";
            makeProfiles(Folder +"/Inputs/profiledata.csv",Folder);
            break;
        }
        case(LakePlacid): {
            offset = -0.25/(4.26-2.85)*60.0/100.0;  // Taken from graphical measurements from curve 1.  Includes ice thickness
            scale = 0.01;
            Folder = "TrackData_LKP";
            convertProfiles2BRep(Folder);
            fs::copy_file(Folder + "/Inputs/profiles.txt", Folder + "/Results/profiles.txt", fs::copy_options::overwrite_existing);
            break;
        }
        case(Sochi): {
            offset = -0.03;  // Taken from graphical measurements from curve 1.  Includes ice thickness
            scale = 0.001;
            Folder = "TrackData_Sochi";
            convertProfiles2BRep(Folder);
            fs::copy_file(Folder + "/Inputs/profiles.txt", Folder + "/Results/profiles.txt", fs::copy_options::overwrite_existing);
            break;
        }
        case(ParkCity): {
            /* The file profiledata_luge.csv is the men's luge start I think not used currently */
            offset = -1.5*2.54/100;  // From Jon Owen
            scale = 0.001;
            Folder = "TrackData_ParkCity";
            makeProfiles(Folder +"/Inputs/profiledata.csv",Folder);
            break;
        }
        case(ParkCityNoEntries): {
            /* The file profiledata_luge.csv is the men's luge start I think not used currently */
            offset = -1.5*2.54/100;  // From Jon Owen
            scale = 0.001;
            Folder = "TrackData_ParkCity_noentries";
            makeProfiles(Folder +"/Inputs/profiledata.csv",Folder);
            break;
        }
        case(ParkCityLugeStart): {
            /* The file profiledata_luge.csv is the men's luge start I think not used currently */
            offset = -1.5*2.54/100;  // From Jon Owen
            scale = 0.001;
            Folder = "TrackData_ParkCity_lugestart";
            makeProfiles(Folder +"/Inputs/profiledata.csv",Folder);
            break;
        }
    }
    
    switch (myOp) {
        case(loft): {
            makeTrackLoft(Folder, Folder +"/Results",scale);
            offsetSurface(Folder +"/Results/theTrack.brep",Folder +"/Results/withIce", offset);
            break;
        }
        case(BSpline): {
            makeTrackBSpline(Folder, Folder +"/Results",scale);
            offsetSurface(Folder +"/Results/theTrack.brep",Folder +"/Results/withIce", offset);
            break;
        }
        default:
            std::cout << "Don't know how to do that" << std::endl;
    }
    
    std::string brepTrack   = Folder + "/Results/theTrack.brep";
    std::string brepWithIce = Folder + "/Results/withIce.brep";

    // if withIce is written into a subfolder, find it
    if (!fs::exists(brepWithIce) && fs::exists(Folder + "/Results/withIce"))
    {
        for (auto& p : fs::directory_iterator(Folder + "/Results/withIce"))
            if (p.path().extension() == ".brep") { brepWithIce = p.path().string(); break; }
    }

    ExportBREPFileToSTL(brepTrack,   Folder + "/Results/theTrack.stl", 0.05, 0.35, false);
    ExportBREPFileToSTL(brepWithIce, Folder + "/Results/withIce.stl",  0.05, 0.35, false);
    return 0;
}
