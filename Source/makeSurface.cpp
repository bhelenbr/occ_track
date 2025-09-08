//
//  makeSurface.cpp
//  occ_track
//
//  Created by Brian Helenbrook on 5/4/24.
//

#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Ax1.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <GeomTools.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <GC_MakeSegment.hxx>
#include <GC_MakeCircle.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GC_MakeEllipse.hxx>
#include <GC_MakeArcOfEllipse.hxx>

#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_CompCurve.hxx>
//#include <BRepAdaptor_Curve2D.hxx>
#include <BRepOffset.hxx>
// #include <BRepPrimAPI_MakeBox.hxx>
// #include <BRepAlgoAPI_Cut.hxx>
// #include <BRepAlgoAPI_Fuse.hxx>
//#include <BRepPrimAPI_MakePrism.hxx>


#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Builder.hxx>
//#include <TopoDS_CompSolid.hxx>
#include <TopExp_Explorer.hxx>

#include <TColgp_HArray1OfPnt.hxx>
#include <TColgp_HArray2OfPnt.hxx>
#include <TColStd_HArray1OfReal.hxx>
#include <TColStd_HArray1OfInteger.hxx>

#include <STEPControl_Reader.hxx>
// #include <BOPTools_DSFiller.hxx>

#include <array>
#include <fstream>
#include <sstream>
// #include <assert>

#include "occ_track.h"
//#define TEST_PROFILES

void loadProfileData(std::string const filename, std::vector<profileData>& profiles) {
    profileData myProfile;
    std::ifstream profileFile;
    profileFile.open(filename);
    while(profileFile >> myProfile.filename >> myProfile.sx >> myProfile.BRK) {
        profiles.push_back(myProfile);
    }
    profileFile.close();
}

void transformProfile(Handle(Geom_BSplineCurve) C, double sx, double scale, std::string filename, TopoDS_Wire& aWire) {
    BRep_Builder aBuilder;
    TopoDS_Shape aShape;
    std::ifstream fs;
    fs.open(filename);
    if ( (fs.rdstate() & std::ifstream::failbit ) != 0 ) {
        std::cout << "Trouble reading " << filename << std::endl;
        exit(1);
    }
    BRepTools::Read(aShape,fs,aBuilder);
    
    /* Just exploring this surface a bit */
    TopExp_Explorer aWireExplorer(aShape, TopAbs_WIRE);
    if (aWireExplorer.More()) {
        /* If profile is a wire use that */
        aWire = TopoDS::Wire(aWireExplorer.Current());

        /* Check if orientation is correct */
        BRepAdaptor_CompCurve myComposite(aWire);
        Standard_Real begin = myComposite.FirstParameter();
        Standard_Real end = myComposite.LastParameter();
        gp_Pnt beginPnt = myComposite.Value(begin);
        gp_Pnt endPnt = myComposite.Value(end);
        if (beginPnt.X() > endPnt.X()) {
            /* Need to reverse */
            TopExp_Explorer edgeExplorer(aWire, TopAbs_EDGE);
            struct curveData {
                Handle(Geom_Curve) aCurve;
                double begin;
                double end;
            } cData;
            std::vector<curveData> storeCurves;
            for(;edgeExplorer.More();edgeExplorer.Next()){
                TopoDS_Edge anEdge = TopoDS::Edge(edgeExplorer.Current());
                cData.aCurve = BRep_Tool::Curve(anEdge, cData.begin, cData.end);
                storeCurves.push_back(cData);
            }
            BRepBuilderAPI_MakeWire mkWire;
            for (auto riter = storeCurves.rbegin(); riter != storeCurves.rend(); ++riter)
            {
                Handle(Geom_Curve) aCurve = riter->aCurve;
                begin = riter->begin;
                end = riter->end;
                begin = aCurve->ReversedParameter(begin);
                end = aCurve->ReversedParameter(end);
                TopoDS_Edge anEdge = BRepBuilderAPI_MakeEdge(aCurve->Reversed(), end, begin);
                mkWire.Add(anEdge);
            }
            aWire = mkWire.Wire();
        }
        
    }
    else {
        /* Convert profile to wire and use */
        TopExp_Explorer anEdgeExplorer(aShape, TopAbs_EDGE);
        TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer.Current());
        
        /* This is a fix for the reversed direction curves of Lake Placid */
        Standard_Real begin, end;
        Handle(Geom_Curve) aCurve = BRep_Tool::Curve(anEdge, begin, end);
        gp_Pnt beginPnt = aCurve->Value(begin);
        gp_Pnt endPnt = aCurve->Value(end);
        if (beginPnt.X() > endPnt.X()) {
            begin = aCurve->ReversedParameter(begin);
            end = aCurve->ReversedParameter(end);
            anEdge = BRepBuilderAPI_MakeEdge(aCurve->Reversed(), end, begin);
        }
        aWire = BRepBuilderAPI_MakeWire(anEdge);
    }
 
    TopExp_Explorer aVertexExplorer(aShape, TopAbs_VERTEX, TopAbs_WIRE);
    TopoDS_Vertex aVertex = TopoDS::Vertex(aVertexExplorer.Current());
    gp_Pnt CE = BRep_Tool::Pnt(aVertex);
    gp_Pnt origin(0.,0.,0.);
    
    /* Find transformation from profile section to position on track */
    gp_Pnt loc;
    gp_Vec dir;
    C->D1(sx, loc, dir);
#ifdef TEST_PROFILES
    loc = gp_Pnt(0.0,sx,0.0);
#endif
    gp_Trsf trackScale;
    trackScale.SetScale(origin,scale);  // Convert to meters
    
    CE.Scale(origin,scale);
    gp_Trsf profileTrans;
    profileTrans.SetTranslation(CE, loc);
    
    gp_Trsf profileRot;
    gp_Dir xAxis(1.0,0.0,0.0);
    gp_Ax1 rotAx(CE, xAxis);
    profileRot.SetRotation(rotAx, M_PI/2);
    
    gp_Trsf trackRot;
    Standard_Real angle = atan2(dir.Y(),dir.X())-M_PI/2;
    gp_Dir zAxis(0.0,0.0,1.0);
    gp_Ax1 rotAz(CE, zAxis);
    trackRot.SetRotation(rotAz, angle);
    
#ifdef TEST_PROFILES
    gp_Trsf Total = profileTrans*profileRot*trackScale;
#else
    gp_Trsf Total = profileTrans*trackRot*profileRot*trackScale;
#endif
    BRepBuilderAPI_Transform myTransform(aWire,Total);
    aWire = TopoDS::Wire(myTransform.ModifiedShape(aWire));
}

void makeTrackLoft(std::string const sourceFolder, std::string const destinationFolder, double scale) {
    Handle(Geom_BSplineCurve) C;
    makeTrackSpine(sourceFolder, destinationFolder, C);
    
    /* Load profiles */
    std::vector<profileData> profiles;
    loadProfileData(sourceFolder +"/profiles.txt", profiles);
    const int nProfiles = int(profiles.size());
    
    int BRKcount = 0;
    int row = 0;
    while (row < nProfiles) {
        BRepOffsetAPI_ThruSections trackLoft(false,true);
        do {
            std::string filename = sourceFolder +"/" +profiles[row].filename;
            double sx = profiles[row].sx;
            TopoDS_Wire newWire;
            transformProfile(C, sx, scale, filename, newWire);
#ifdef TEST_PROFILES
            BRepTools::Write(newWire,(destinationFolder +"/" +profiles[row].filename).c_str());
#endif
            trackLoft.AddWire(newWire);
            ++row;
        } while(profiles[row-1].BRK == 0 && row < nProfiles);
        
        TopoDS_Shape track = trackLoft.Shape();
        std::ostringstream nstr;
        nstr << "/theTrack" << BRKcount << ".brep";
        /* Write brep file */
        BRepTools::Write(track,(destinationFolder +nstr.str()).c_str());
        trackLoft.Init();
        ++BRKcount;
    }
}

void makeTrackBSpline(std::string const sourceFolder, std::string const destinationFolder, double scale) {
    Handle(Geom_BSplineCurve) C;
    makeTrackSpine(sourceFolder, destinationFolder, C);
    
    /* Load profiles */
    std::vector<profileData> profiles;
    loadProfileData(sourceFolder +"/Results/profiles.txt", profiles);
    const int nProfiles = int(profiles.size());
    const int nCrossPts = 100;
    
    TColgp_Array2OfPnt myPoints(1,nProfiles+1,1,nCrossPts);
    TColgp_Array1OfPnt myProfilePoints(1,nCrossPts);
    TColStd_Array1OfReal UKnots(1,nProfiles);
    TColStd_Array1OfInteger UMults(1,nProfiles);
    TColStd_Array1OfReal VKnots(1,nCrossPts);
    TColStd_Array1OfInteger VMults(1,nCrossPts);
    
    /* Set up cross point spline */
    for (int col=1;col<=nCrossPts;++col) {
        VKnots.SetValue(col,static_cast<double>(col-1)/(nCrossPts-1));
        VMults.SetValue(col,1);
    }
    VMults.SetValue(1,3);
    VMults.SetValue(nCrossPts,3);
    
    for (int row = 0; row < nProfiles; ++row) {
        std::string filename = sourceFolder +"/" +profiles[row].filename;
        double sx = profiles[row].sx;
        TopoDS_Wire aWire;
        transformProfile(C, sx, scale, filename, aWire);
#ifdef TEST_PROFILES
        BRepTools::Write(aWire,(destinationFolder +"/" +profiles[row].filename).c_str());
#endif
        
        UKnots.SetValue(row+1,sx);
        UMults.SetValue(row+1,1);
        wireToPoints(aWire, nCrossPts, myProfilePoints);
        for(int col=0;col < nCrossPts;++col) {
            gp_Pnt myPnt = myProfilePoints(col+1);
#ifdef UNITY
            gp_Pnt unitySucks(myPnt);
            // Going to switch to Unity's stupid coordinate system;
            unitySucks.SetY(myPnt.Z());
            unitySucks.SetZ(-myPnt.Y());
            myPnt = unitySucks;
#endif
            myPoints.SetValue(row+2,col+1,myPnt);
        }
    }
    
    /* Add additional control point */
    for(int col=0;col < nCrossPts;++col) {
        gp_Pnt P1 = myPoints(2,col+1);
        myPoints.SetValue(1,col+1,P1);
        gp_Pnt P2 = myPoints(3,col+1);
        P1.SetX(0.5*(P1.X() +P2.X()));
        P1.SetY(0.5*(P1.Y() +P2.Y()));
        P1.SetZ(0.5*(P1.Z() +P2.Z()));
        myPoints.SetValue(2,col+1,P1);
    }
    UMults.SetValue(1,4);
    UMults.SetValue(nProfiles,3);
    
    Handle(Geom_BSplineSurface) trackSurface = new Geom_BSplineSurface(myPoints, UKnots, VKnots, UMults, VMults, 3, 3);
    BRepBuilderAPI_MakeFace FaceMaker(trackSurface, 1.0e-7);
    if (!FaceMaker.IsDone()) {
        std::cout << "FaceMaker Error " << FaceMaker.Error() << std::endl;
    }
    TopoDS_Face F = FaceMaker.Face();
    
    /* Write brep file */
    BRepTools::Write(F,(destinationFolder +"/theTrack.brep").c_str());
    outputBSplineSurface(destinationFolder +"/track_surface_pts.txt", nProfiles, nCrossPts, trackSurface);
    timeIntegrate(destinationFolder,trackSurface);

}


void offsetSurface(std::string input, std::string output, double offset) {
    std::ifstream fs;
    fs.open(input);
    if ( (fs.rdstate() & std::ifstream::failbit ) != 0 ) {
        std::cout << "Trouble reading " << input << std::endl;
        exit(1);
    }
    BRep_Builder aBuilder;
    TopoDS_Shape aShape;
    BRepTools::Read(aShape,fs,aBuilder);
    TopExp_Explorer aFaceExplorer(aShape, TopAbs_FACE);
    
    if (!aFaceExplorer.More()) {
        /* If profile is a wire use that */
        std::cout << "couldn't find face to offset " << input << std::endl;
        exit(1);
    }
    TopoDS_Face trackFace = TopoDS::Face(aFaceExplorer.Current());
    Handle(Geom_Surface) trackSurface = BRep_Tool::Surface(trackFace);
    //You probably want to create a BRepAdaptor_Surface

    // Offset surface to account for cement & ice thickness
    BRepOffset_Status status;
    
    Handle(Geom_Surface) offsetSurface = BRepOffset::Surface(trackSurface,offset,status);
    if (status != BRepOffset_Good) {
        std::cout << "Offset did not work: " << status << :: std::endl;
    }
    
    /* Make a topological face from the spline surface */
    BRepBuilderAPI_MakeFace FaceMaker(offsetSurface, 1.0e-7);
    if (!FaceMaker.IsDone()) {
        std::cout << "FaceMaker Error " << FaceMaker.Error() << std::endl;
    }
    
    TopoDS_Face F = FaceMaker.Face();
    
    /* Write brep file */
    std::string brepname = output +".brep";
    BRepTools::Write(F,brepname.c_str());
  
	/* Write surface file */
	std::ofstream mySurfaceFile;
	mySurfaceFile.open(output +".txt");
	GeomTools::Write(trackSurface,mySurfaceFile);
	mySurfaceFile.close();
}

void outputBSplineSurface(const std::string filename, int nPtsU, int nPtsV, Handle(Geom_BSplineSurface) const C) {
    Standard_Real U1,U2,V1,V2;
    C->Bounds (U1, U2, V1, V2);
    
    Standard_Real du = (U2-U1)/nPtsU;
    Standard_Real dv = (V2-V1)/nPtsV;
    
    /* Test track curve */
    std::ofstream track_test;
    track_test.open(filename);
    track_test << nPtsU << ' ' << nPtsV << std::endl;
    for (int i=0; i<nPtsU; ++i) {
        Standard_Real U = U1+du*i;
        for (int j=0; j<nPtsV; ++j) {
            Standard_Real V = V1+dv*j;
            gp_Pnt P;
            C->D0 (U, V, P);
#ifdef UNITY
            track_test << P.X() << ' ' << P.Z() << ' ' << -P.Y() << std::endl;
#else
            track_test << P.X() << ' ' << P.Y() << ' ' << P.Z() << std::endl;
#endif
            
        }
    }
    track_test.close();
}

