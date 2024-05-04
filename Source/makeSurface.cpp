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

void loadProfileData(std::string const filename, std::vector<profileData>& profiles) {
    profileData myProfile;
    std::ifstream profileFile;
    profileFile.open(filename);
    while(profileFile >> myProfile.filename >> myProfile.sx) {
        profiles.push_back(myProfile);
    }
    profileFile.close();
}

void makeTrackLoft(std::string const sourceFolder, std::string const destinationFolder) {
    Handle(Geom_BSplineCurve) C;
    makeTrackSpine(sourceFolder, destinationFolder, C);
    
    /* Load profiles */
    std::vector<profileData> profiles;
    loadProfileData(sourceFolder +"/profiles.txt", profiles);
    const int nProfiles = int(profiles.size());
    
    BRepOffsetAPI_ThruSections trackLoft;
    for(int row = 0; row < nProfiles; ++ row) {
        BRep_Builder aBuilder;
        TopoDS_Shape aShape;
        Standard_Boolean result = BRepTools::Read(aShape,(sourceFolder +"/" +profiles[row].filename).c_str(),aBuilder);
        if (!result) {
            std::cout << row << std::endl;
            std::cout << "Trouble reading " << sourceFolder +"/" +profiles[row].filename << std::endl;
            exit(1);
        }
        /* Just exploring this surface a bit */
        TopExp_Explorer aWireExplorer(aShape, TopAbs_WIRE);
        TopoDS_Wire aWire = TopoDS::Wire(aWireExplorer.Current());
        
        TopExp_Explorer aVertexExplorer(aShape, TopAbs_VERTEX);
        TopoDS_Vertex aVertex = TopoDS::Vertex(aVertexExplorer.Current());
        gp_Pnt origin = BRep_Tool::Pnt(aVertex);
        
        /* Find transformation from profile section to position on track */
        double sx = profiles[row].sx;
        gp_Pnt loc;
        gp_Vec dir;
        C->D1(sx, loc, dir);
        
        gp_Trsf trackScale;
        trackScale.SetScale(origin, 0.001);  // Convert to meters
        
        gp_Trsf profileTrans;
        profileTrans.SetTranslation(origin, loc);
        
        gp_Trsf profileRot;
        gp_Dir xAxis(1.0,0.0,0.0);
        gp_Ax1 rotAx(origin, xAxis);
        profileRot.SetRotation(rotAx, M_PI/2);
        
        gp_Trsf trackRot;
        Standard_Real angle = atan2(dir.Y(),dir.X())-M_PI/2;
        gp_Dir zAxis(0.0,0.0,1.0);
        gp_Ax1 rotAz(origin, zAxis);
        trackRot.SetRotation(rotAz, angle);

        gp_Trsf Total = profileTrans*trackRot*profileRot*trackScale;
        BRepBuilderAPI_Transform myTransform(aWire,Total);
        TopoDS_Wire newWire = TopoDS::Wire(myTransform.ModifiedShape(aWire));

        trackLoft.AddWire(newWire);
    }
    TopoDS_Shape track = trackLoft.Shape();
    /* Write brep file */
    BRepTools::Write(track,(destinationFolder +"/theTrack.brep").c_str());
}

void makeTrackSpline2D(std::string const sourceFolder, std::string const destinationFolder) {
    Handle(Geom_BSplineCurve) C;
    makeTrackSpine(sourceFolder, destinationFolder, C);
    
    /* Load profiles */
    std::vector<profileData> profiles;
    loadProfileData(sourceFolder +"/profiles.txt", profiles);
    const int nProfiles = int(profiles.size());
    const int nCrossPts = 100;
    
    TColgp_Array2OfPnt myPoints(1,nProfiles+1,1,nCrossPts);
    TColStd_Array1OfReal UKnots(1,nProfiles);
    TColStd_Array1OfInteger UMults(1,nProfiles);
    TColStd_Array1OfReal VKnots(1,nCrossPts);
    TColStd_Array1OfInteger VMults(1,nCrossPts);
    
    /* Now load in profiles and use to make 2D array of points for B-Spline Surface */
    for(int row=0;row<nProfiles;++row) {
        BRep_Builder aBuilder;
        TopoDS_Shape aShape;
        Standard_Boolean result = BRepTools::Read(aShape,(sourceFolder +"/" +profiles[row].filename).c_str(),aBuilder);
        if (!result) {
            std::cout << row << std::endl;
            std::cout << "Trouble reading " << sourceFolder +"/" +profiles[row].filename << std::endl;
            exit(1);
        }
        /* Just exploring this surface a bit */
        TopExp_Explorer anEdgeExplorer(aShape, TopAbs_EDGE);
        TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer.Current());
                
        Standard_Real begin, end;
        Handle(Geom_Curve) aCurve = BRep_Tool::Curve(anEdge, begin, end);
        gp_Pnt beginPnt = aCurve->Value(begin);
        gp_Pnt endPnt = aCurve->Value(end);
        
        if (beginPnt.X() < endPnt.X()) {
            aCurve->Reverse();
        }
        
        TopExp_Explorer aVertexExplorer(aShape, TopAbs_VERTEX);
        TopoDS_Vertex aVertex = TopoDS::Vertex(aVertexExplorer.Current());
        gp_Pnt PV = BRep_Tool::Pnt(aVertex);
        
        /* Find transformation from profile section to position on track */
        double sx = profiles[row].sx;
        gp_Pnt loc;
        gp_Vec dir;
        C->D1(sx, loc, dir);
        
        gp_Trsf trackScale;
        trackScale.SetScale(PV, 0.01);  // Convert to meters
        
        gp_Trsf profileTrans;
        profileTrans.SetTranslation(PV, loc);
        
        gp_Trsf profileRot;
        gp_Dir xAxis(1.0,0.0,0.0);
        gp_Ax1 rotAx(PV, xAxis);
        profileRot.SetRotation(rotAx, M_PI/2);
        
        gp_Trsf trackRot;
        Standard_Real angle = atan2(dir.Y(),dir.X())-M_PI/2;
        gp_Dir zAxis(0.0,0.0,1.0);
        gp_Ax1 rotAz(PV, zAxis);
        trackRot.SetRotation(rotAz, angle);
        
        
        gp_Trsf Total = profileTrans*trackRot*profileRot*trackScale;
        
        //            std::ofstream fout;
        //            std::ostringstream nstr;
        //            nstr.str("");
        //            nstr << destinationFolder +"/profile"  << row << ".pts";
        //            fout.open(nstr.str());
        
        UKnots.SetValue(row+1,sx);
        UMults.SetValue(row+1,1);
        for(int col=0;col < nCrossPts;++col) {
            Standard_Real u = begin +col*(end-begin)/(nCrossPts-1);
            gp_Pnt PV = aCurve->Value(u);
            gp_Pnt P2 = PV.Transformed(Total);
            
            // Going to switch to Unity's stupid coordinate system;
            PV = P2;
            P2.SetY(PV.Z());
            P2.SetZ(-PV.Y());
            
            //                fout << P2.X() << ' ' << P2.Y() << ' ' << P2.Z() << std::endl;
            myPoints.SetValue(row+2,col+1,P2);
        }
        //            fout.close();
        
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
    
    
    for (int col=1;col<=nCrossPts;++col) {
        VKnots.SetValue(col,static_cast<double>(col-1)/(nCrossPts-1));
        VMults.SetValue(col,1);
    }
    VMults.SetValue(1,3);
    VMults.SetValue(nCrossPts,3);
    
    Handle(Geom_BSplineSurface) trackSurface = new Geom_BSplineSurface(myPoints, UKnots, VKnots, UMults, VMults, 3, 3);
    outputSpline2D(destinationFolder +"/track_surface_pts.txt", nProfiles, nCrossPts, trackSurface);

    
    // Offset surface to account for cement & ice thickness
    BRepOffset_Status Status;
    const Standard_Real Offset = 0.25/(4.26-2.85)*60.0/100.0;  // Taken from graphical measurements from curve 1.  Includes ice thickness
    Handle(Geom_Surface) offsetSurface = BRepOffset::Surface(trackSurface,Offset,Status);
    if (Status != BRepOffset_Good) {
        std::cout << "Offset did not work: " << Status << :: std::endl;
    }
    
    /* Make a topological face from the spline surface */
    BRepBuilderAPI_MakeFace FaceMaker(offsetSurface, 1.0e-7);
    if (!FaceMaker.IsDone()) {
        std::cout << "FaceMaker Error " << FaceMaker.Error() << std::endl;
    }
    
    TopoDS_Face F = FaceMaker.Face();
    
    
    /* Just exploring this surface a bit */
    TopExp_Explorer anEdgeExplorer(F, TopAbs_EDGE);
    
    while(anEdgeExplorer.More()) {
        TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer.Current());
        std::cout << "I found an edge" << std::endl;
        anEdgeExplorer.Next();
    }
    
    /* Write brep file */
    BRepTools::Write(F,(destinationFolder +"/theTrack.brep").c_str());
    
    /* Write surface file */
    std::ofstream mySurfaceFile;
    mySurfaceFile.open(destinationFolder +"/surface.txt");
    GeomTools::Write(trackSurface,mySurfaceFile);
    mySurfaceFile.close();
    
    timeIntegrate(destinationFolder,trackSurface);
}

void makeStraightSpline2D(std::string const destinationFolder) {
    
    /* Make cross sectional profile of straight track test */
    /* I'm guessing 2 meters wide and 0.5 m high */
    
    int Npts = 300;
    double ds = 3.0/Npts;
    double radius = 0.1;
    double width = 2;
    double height = 0.5;
    
    int NptsBottom = (width-2.*radius)/ds;
    int NptsArc = M_PI/2.*radius/ds;
    int NptsSide = (height-radius)/ds;
    
    std::cout << NptsBottom << ' ' << NptsArc << ' ' << NptsSide << std::endl;
    std::cout << height << ' ' << radius << ' ' << ds << std::endl;
    
    std::vector<gp_Pnt2d> profile_pnts(NptsBottom+2*NptsArc +2*NptsSide);
    
    
    /* First side */
    int pntcnt = 0;
    for (int pnt=0; pnt < NptsSide; ++pnt) {
        profile_pnts[pntcnt++].SetCoord(-width*0.5,height- (height-radius)/NptsSide*pnt);
    }
    /* Arc */
    for (int pnt=0; pnt < NptsArc; ++pnt) {
        profile_pnts[pntcnt++].SetCoord(-width*0.5+radius -radius*cos(M_PI*pnt/(2.*NptsArc)),radius - radius*sin(M_PI*pnt/(2.*NptsArc)));
    }
    /* Bottom */
    for (int pnt=0; pnt < NptsBottom; ++pnt) {
        profile_pnts[pntcnt++].SetCoord(-width*0.5+radius +pnt*(width-2*radius)/NptsBottom,0.0);
    }
    /* Arc */
    for (int pnt=0; pnt < NptsArc; ++pnt) {
        profile_pnts[pntcnt++].SetCoord(width*0.5-radius +radius*sin(M_PI*pnt/(2.*NptsArc)),radius - radius*cos(M_PI*pnt/(2.*NptsArc)));
    }
    /* Final side */
    for (int pnt=0; pnt < NptsSide; ++pnt) {
        profile_pnts[pntcnt++].SetCoord(width*0.5,radius +(height-radius)/NptsSide*pnt);
    }
    
    
    
    /* Make a B-spline Surface */
    double L = 500.0;
    double grade = -0.17;
    double dL = 1.0;
    
    int nrow = L/dL;
    int ncol = NptsBottom+2*(NptsArc +NptsSide);
    std::cout << ncol << ' ' << pntcnt << std::endl;
    TColgp_Array2OfPnt myPoints(1,nrow,1,ncol);
    
    gp_Pnt P1;
    for (int row=1;row<=nrow;++row) {
        double y = dL*row;
        double z = y*grade;
        for(int col=1;col<=ncol;++col) {
            P1 = gp_Pnt(profile_pnts[col-1].X(),y,-z-profile_pnts[col-1].Y());
            myPoints.SetValue(row,col,P1);
        }
    }
    
    
    TColStd_Array1OfReal UKnots(1,nrow);
    TColStd_Array1OfInteger UMults(1,nrow);
    for (int row=1;row<=nrow;++row) {
        UKnots.SetValue(row,static_cast<double>(row-1)/(nrow-1));
        UMults.SetValue(row,1);
    }
    UMults.SetValue(1,3);
    UMults.SetValue(nrow,3);
    
    TColStd_Array1OfReal VKnots(1,ncol);
    TColStd_Array1OfInteger VMults(1,ncol);
    for (int col=1;col<=ncol;++col) {
        VKnots.SetValue(col,static_cast<double>(col-1)/(ncol-1));
        VMults.SetValue(col,1);
    }
    VMults.SetValue(1,3);
    VMults.SetValue(ncol,3);
    
    Handle(Geom_BSplineSurface) C = new Geom_BSplineSurface(myPoints, UKnots, VKnots, UMults, VMults, 3, 3);
    
    
    // Lets check surface to make sure it makes sense
    C->D0 (0.5, 0.5, P1);
    std::cout << " Mid Point of Surface " << P1.X() << ' ' << P1.Y() << ' ' << P1.Z() << std::endl;
    
    Handle(Geom_Curve) U1 =    C->UIso(0.0);
    Handle(Geom_Curve) U2 =    C->UIso(1.0);
    U2->Reverse();
    Handle(Geom_Curve) V1 =    C->VIso(0.0);
    V1->Reverse();
    Handle(Geom_Curve) V2 =    C->VIso(1.0);
    
    U1->D0 (0.0, P1);
    std::cout << "First Point of U1 " << P1.X() << ' ' << P1.Y() << ' ' << P1.Z() << std::endl;
    U1->D0 (1.0, P1);
    std::cout << "Last Point of U1 " << P1.X() << ' ' << P1.Y() << ' ' << P1.Z() << std::endl;
    
    V2->D0 (0.0, P1);
    std::cout << "First Point of V2 " << P1.X() << ' ' << P1.Y() << ' ' << P1.Z() << std::endl;
    V2->D0 (1.0, P1);
    std::cout << "Last Point of V2 " << P1.X() << ' ' << P1.Y() << ' ' << P1.Z() << std::endl;
    
    U2->D0 (0.0, P1);
    std::cout << "First Point of U2 " << P1.X() << ' ' << P1.Y() << ' ' << P1.Z() << std::endl;
    U2->D0 (1.0, P1);
    std::cout << "Last Point of U2 " << P1.X() << ' ' << P1.Y() << ' ' << P1.Z() << std::endl;
    
    V1->D0 (0.0, P1);
    std::cout << "First Point of V1 " << P1.X() << ' ' << P1.Y() << ' ' << P1.Z() << std::endl;
    V1->D0 (1.0, P1);
    std::cout << "Last Point of V1 " << P1.X() << ' ' << P1.Y() << ' ' << P1.Z() << std::endl;
    
    
    /* This was just poking around trying to figure out how things worked */
    //    TopoDS_Edge anEdge1OnSurf1 = BRepBuilderAPI_MakeEdge(U1, C);
    //    TopoDS_Edge aEdgeU1 = BRepBuilderAPI_MakeEdge(U1);
    //    TopoDS_Edge aEdgeU2 = BRepBuilderAPI_MakeEdge(U2);
    //    TopoDS_Edge aEdgeV1 = BRepBuilderAPI_MakeEdge(V1);
    //    TopoDS_Edge aEdgeV2 = BRepBuilderAPI_MakeEdge(V2);
    
    //    BRepBuilderAPI_MakeWire MakeWire(aEdgeU1,aEdgeV2);
    //    MakeWire.Add(aEdgeU2);
    //    MakeWire.Add(aEdgeV1);
    //    TopoDS_Wire myWireProfile = MakeWire.Wire();
    
    //TopoDS_Face myFaceProfile = BRepBuilderAPI_MakeFace(myWireProfile);
    //gp_Vec aPrismVec(0, 0, 1.0);
    //TopoDS_Shape myBody = BRepPrimAPI_MakePrism(myFaceProfile, aPrismVec);
    //BRepTools::Write(myBody,"myBody.brep");
    
    
    /* Make a topological face from the spline surface */
    BRepBuilderAPI_MakeFace FaceMaker(C, 0.0, 1, 0.0, 1, 1.0e-7);
    if (!FaceMaker.IsDone()) {
        std::cout << "FaceMaker Error " << FaceMaker.Error() << std::endl;
    }
    
    TopoDS_Face F = FaceMaker.Face();

    /* Write brep file */
    BRepTools::Write(F,(destinationFolder +"/spline.brep").c_str());
    
    return;
}

void outputSpline2D(const std::string filename, int nPtsU, int nPtsV, Handle(Geom_BSplineSurface) const C) {
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
            track_test << P.X() << ' ' << P.Y() << ' ' << P.Z() << std::endl;
        }
    }
    track_test.close();
}

