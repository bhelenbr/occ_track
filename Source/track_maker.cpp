#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <TColgp_HArray2OfPnt.hxx>
#include <TColStd_HArray1OfReal.hxx>
#include <TColStd_HArray1OfInteger.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp_Explorer.hxx>
#include <STEPControl_Reader.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Curve2D.hxx>
#include <Geom_Curve.hxx>
#include <GeomTools.hxx>
#include <gp_Ax1.hxx>
#include <BRepOffset.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <Geom_Surface.hxx>
#include <GC_MakeSegment.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <array>
#include <fstream>
#include <sstream>

#include "occ_track.h"

// #include <assert.h>
// #include <BRepPrimAPI_MakeBox.hxx>
// #include <BRepAlgoAPI_Cut.hxx>
// #include <BRepAlgoAPI_Fuse.hxx>
// #include <BOPTools_DSFiller.hxx>
// #include <TopoDS_Builder.hxx>
// #include <TopoDS_CompSolid.hxx>

// To compile
// g++ spline.cpp -I${HOME}/Packages/include -I${HOME}/Packages/include/oce -L${HOME}/Packages/lib -lTKBRep -lTKPrim -lTKernel -lTKMath -lTKTopAlgo -lTKBO -lTKGeomAlgo -lTKG3D

// To find what library has a routine
// nm -o libTK* | grep GeomTools | grep Read | grep ' T '

// Tutorial is in opencascade directory under doc / tutorial

// Doxygen is in that folder as well

// Use spotlight to find header files or search in XCode File

void makeTrackSpline2D(std::string const sourceFolder, std::string const destinationFolder) {
    const int nCrossPts = 100;
    
    /* Read file to see how many points there are */
    std::ifstream trackFile;
    trackFile.open(sourceFolder +"/track.txt");
    int nPoints = 0;
    double sx,x,y,z;
    while (trackFile >> sx >> x >> y >> z)
        ++nPoints;
    trackFile.close();
    
    /* load arcLengths & points */
    std::vector<double> arcLength(nPoints);
    std::vector<gp_Pnt> points(nPoints);
    trackFile.open(sourceFolder +"/track.txt");
    
    for (int n=0; n < nPoints; ++n) {
        trackFile >> sx >> x >> y >> z;
        arcLength[n] = sx;
        points[n].SetCoord(x,y,z);
    }
    trackFile.close();
    
    /* Make spine of track */
    Handle(Geom_BSplineCurve) C;
    makeSpine(arcLength,points,C);
    outputSpline1D(destinationFolder +"/track_test_pts.txt", nPoints, C);
    
    /* Load names of profile files & arcLength of profile */
    struct profileData {
        std::string filename;
        double sx;
    };
    
    /* Count number of profiles */
    std::ifstream  profileFile;
    profileFile.open(sourceFolder +"/profiles.txt");
    int nProfiles = 0;
    std::string FileName;
    while (profileFile >> FileName >> sx) {
        ++nProfiles;
    }
    profileFile.close();
    
    /* Load profiles */
    std::vector<profileData> profiles(nProfiles);
    profileFile.open(sourceFolder +"/profiles.txt");
    for (int n=0; n < nProfiles; ++n)	{
        profileFile >> profiles[n].filename >> profiles[n].sx;
    }

    
    TColgp_Array2OfPnt myPoints(1,nProfiles+1,1,nCrossPts);
    TColStd_Array1OfReal UKnots(1,nProfiles);
    TColStd_Array1OfInteger UMults(1,nProfiles);
    TColStd_Array1OfReal VKnots(1,nCrossPts);
    TColStd_Array1OfInteger VMults(1,nCrossPts);
    
    /* Now load in profiles and use to make 2D array of points for B-Spline Surface */
    for(int row=0;row<nProfiles;++row) {
        /* Load in the step file */
        STEPControl_Reader reader;
        IFSelect_ReturnStatus stat = reader.ReadFile((sourceFolder +"/" +profiles[row].filename).c_str());
        if (!stat) {
            std::cout << "There was a problem " << stat << std::endl;
            std::cout << sourceFolder +"/" +profiles[row].filename << std::endl;
            exit(1);
        }
        // ItemsByEntity - CountByItem - ListByItem -
        IFSelect_PrintCount mode = IFSelect_ItemsByEntity;
        reader.PrintCheckLoad(false,mode);
        // Handle_StepData_StepModel MyStepModel = reader.StepModel ();
        
        
        Standard_Integer num = reader.TransferRoots();
        if (num != 1) {
            std::cout << "Uh-Oh more than 1 shape in step file\n";
        }
        for (int rank = 1; rank <= num; ++rank) {
            TopoDS_Shape shape = reader.Shape(rank);
            
            /* Just exploring this surface a bit */
            TopExp_Explorer anEdgeExplorer(shape, TopAbs_EDGE);
            
            TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer.Current());
            
            Standard_Real begin, end;
            Handle(Geom_Curve) aCurve = BRep_Tool::Curve(anEdge, begin, end);
            gp_Pnt beginPnt = aCurve->Value(begin);
            gp_Pnt endPnt = aCurve->Value(end);
            
            if (beginPnt.X() < endPnt.X()) {
                aCurve->Reverse();
            }
            
            /* Find transformation from profile section to position on track */
            sx = profiles[row].sx;
            gp_Pnt loc;
            gp_Vec dir;
            C->D1(sx, loc, dir);
            
            TopExp_Explorer aVertexExplorer(shape, TopAbs_VERTEX);
            TopoDS_Vertex aVertex = TopoDS::Vertex(aVertexExplorer.Current());
            gp_Pnt PV = BRep_Tool::Pnt(aVertex);
            
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

void makeTrackLoft(std::string const sourceFolder, std::string const destinationFolder) {
    
    BRep_Builder aBuilder;
    TopoDS_Wire aWire;
    Standard_Boolean result = BRepTools::Read(aWire,"Straight/Results/straight.brep",aBuilder);
    std::cout << result << ' ' << !aWire.IsNull() << std::endl;
    /* Just exploring this surface a bit */
    TopExp_Explorer anEdgeExplorer(aWire, TopAbs_EDGE);
    
    while(anEdgeExplorer.More()) {
        TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer.Current());
        std::cout << "I found an edge" << std::endl;
        anEdgeExplorer.Next();
    }
    
    TopExp_Explorer aVertexExplorer(aWire, TopAbs_VERTEX);
    while(aVertexExplorer.More()) {
        TopoDS_Vertex aVertex = TopoDS::Vertex(aVertexExplorer.Current());
        std::cout << "I found an vertex" << std::endl;
        anEdgeExplorer.Next();
    }

//    std::assert(result);
//    std::assert(!aShape.IsNull());
    
    
    
    
    // BRepOffsetAPI_ThruSections

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
    
    Handle(Geom_Curve) U1 =	C->UIso(0.0);
    Handle(Geom_Curve) U2 =	C->UIso(1.0);
    U2->Reverse();
    Handle(Geom_Curve) V1 =	C->VIso(0.0);
    V1->Reverse();
    Handle(Geom_Curve) V2 =	C->VIso(1.0);
    
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
    //	TopoDS_Edge anEdge1OnSurf1 = BRepBuilderAPI_MakeEdge(U1, C);
    //	TopoDS_Edge aEdgeU1 = BRepBuilderAPI_MakeEdge(U1);
    //	TopoDS_Edge aEdgeU2 = BRepBuilderAPI_MakeEdge(U2);
    //	TopoDS_Edge aEdgeV1 = BRepBuilderAPI_MakeEdge(V1);
    //	TopoDS_Edge aEdgeV2 = BRepBuilderAPI_MakeEdge(V2);
    
    //	BRepBuilderAPI_MakeWire MakeWire(aEdgeU1,aEdgeV2);
    //	MakeWire.Add(aEdgeU2);
    //	MakeWire.Add(aEdgeV1);
    //	TopoDS_Wire myWireProfile = MakeWire.Wire();
    
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

void outputSpline1D(const std::string filename, int nPts, Handle(const Geom_BSplineCurve) const C) {
    /* Test track curve */
    Standard_Real sBgn = C->FirstParameter();
    Standard_Real sEnd = C->LastParameter();
    Standard_Real ds = (sEnd-sBgn)/(nPts-1);
    
    std::ofstream track_test;
    track_test.open(filename);
    for (int i=0;i<nPts;++i) {
        double sx = sBgn +ds*i;
        gp_Pnt P;
        C->D0(sx, P);
        track_test << P.X() << ' ' << P.Z() << ' ' << -P.Y() << std::endl;
    }
    track_test.close();
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

void makeStraightProfile(std::string const filename,double WL, double WR, double HL, double HR, double RL, double RR) {
    double drl = RL*(1.-1./sqrt(2.));
    double drr = RR*(1.-1./sqrt(2.));
    
    gp_Pnt aPnt1(-WL, HL, 0);
    gp_Pnt aPnt2(-WL, RL, 0);
    gp_Pnt aPnt3(-WL+drl,drl,0);
    gp_Pnt aPnt4(-WL+RL, 0.0, 0);
    gp_Pnt aPnt5(WR-RR, 0.0, 0);
    gp_Pnt aPnt6(WR-drr,drr,0);
    gp_Pnt aPnt7(WR, RR, 0);
    gp_Pnt aPnt8(WR, HR, 0);
    
    Handle(Geom_TrimmedCurve) aSegment1    = GC_MakeSegment(aPnt1, aPnt2);
    Handle(Geom_TrimmedCurve) aArcOfCircle1 = GC_MakeArcOfCircle(aPnt2,aPnt3,aPnt4);
    Handle(Geom_TrimmedCurve) aSegment2    = GC_MakeSegment(aPnt4, aPnt5);
    Handle(Geom_TrimmedCurve) aArcOfCircle2 = GC_MakeArcOfCircle(aPnt5,aPnt6,aPnt7);
    Handle(Geom_TrimmedCurve) aSegment3    = GC_MakeSegment(aPnt7, aPnt8);
    
    TopoDS_Edge anEdge1 = BRepBuilderAPI_MakeEdge(aSegment1);
    TopoDS_Edge anEdge2 = BRepBuilderAPI_MakeEdge(aArcOfCircle1);
    TopoDS_Edge anEdge3 = BRepBuilderAPI_MakeEdge(aSegment2);
    TopoDS_Edge anEdge4 = BRepBuilderAPI_MakeEdge(aArcOfCircle2);
    TopoDS_Edge anEdge5 = BRepBuilderAPI_MakeEdge(aSegment3);
    
    TopoDS_Wire aWire1 = BRepBuilderAPI_MakeWire(anEdge1, anEdge2, anEdge3);
    TopoDS_Wire aWire2 = BRepBuilderAPI_MakeWire(anEdge4, anEdge5);
    BRepBuilderAPI_MakeWire mkWire;
    mkWire.Add(aWire1);
    mkWire.Add(aWire2);
    TopoDS_Wire myWireProfile = mkWire.Wire();
    
    std::ofstream of;
    of.open(filename);
    BRepTools::Write(myWireProfile,of);
    
    gp_Pnt Origin (0., 0., 0.);
    BRep_Builder aBuilder;
    TopoDS_Vertex V = BRepBuilderAPI_MakeVertex(Origin);
    BRepTools::Write(V,of);
    of.close();
}

void makeSpine(const std::vector<double>& arcLength, const std::vector<gp_Pnt>& points, Handle(Geom_BSplineCurve)& C) {
    const int nPoints = static_cast<int>(arcLength.size());
    
    /* Make spine of track */
    /* Array of points */
    Handle(TColgp_HArray1OfPnt) myPointsHandle = new TColgp_HArray1OfPnt(1,nPoints);
    /* array for arcLengths */
    Handle(TColStd_HArray1OfReal) mySxHandle = new TColStd_HArray1OfReal(1,nPoints);
    
    for (int n=0;n<nPoints;++n) {
        myPointsHandle->SetValue(n+1,points[n]);
        mySxHandle->SetValue(n+1,arcLength[n]);
    }
    GeomAPI_Interpolate Interp(myPointsHandle,mySxHandle,false,0.0005);
    Interp.Perform();
    if (!Interp.IsDone()) {
        std::cout << "Uh-Oh Spine of track not made\n";
        exit(1);
    }
    C = Interp.Curve();
}

void step2BRep(std::string const filename) {
    
    std::string outfile = filename.substr(0,filename.rfind(".")) +".brep";
    
    /* Load in the step file */
    STEPControl_Reader reader;
    IFSelect_ReturnStatus stat = reader.ReadFile(filename.c_str());
    if (!stat) {
        std::cout << "There was a problem " << stat << std::endl;
        std::cout << filename << std::endl;
        exit(1);
    }
    // ItemsByEntity - CountByItem - ListByItem -
    IFSelect_PrintCount mode = IFSelect_ItemsByEntity;
    reader.PrintCheckLoad(false,mode);
    // Handle_StepData_StepModel MyStepModel = reader.StepModel ();
    
    
    Standard_Integer num = reader.TransferRoots();
    if (num != 1) {
        std::cout << "Uh-Oh more than 1 shape in step file\n";
    }
    for (int rank = 1; rank <= num; ++rank) {
        TopoDS_Shape shape = reader.Shape(rank);
        BRepTools::Write(shape,outfile.c_str());
    }
}

void convertProfiles2BRep(std::string const sourceFolder) {
    /* Load names of profile files & arcLength of profile */
    struct profileData {
        std::string filename;
        double sx;
    };
    
    /* Count number of profiles */
    std::ifstream  profileFile;
    profileFile.open(sourceFolder +"/profiles.txt");
    int nProfiles = 0;
    std::string FileName;
    double sx;
    while (profileFile >> FileName >> sx) {
        ++nProfiles;
    }
    profileFile.close();
    
    /* Load profiles */
    std::vector<profileData> profiles(nProfiles);
    profileFile.open(sourceFolder +"/profiles.txt");
    for (int n=0; n < nProfiles; ++n)    {
        profileFile >> profiles[n].filename >> profiles[n].sx;
    }
    
    for(int row=0;row<nProfiles;++row) {
        step2BRep(sourceFolder +"/" +profiles[row].filename);
    }
}


void experiments() {
    //    gp_Vec aPrismVec(0, 0, 1.0);
    //    TopoDS_Shape myBody = BRepPrimAPI_MakePrism(F, aPrismVec);
    //    BRepTools::Write(myBody,"myBody.brep");
    
    
    //     GeomAPI_PointsToBSplineSurface mySurface;
    //     mySurface.Init(myPoints);
    //     if (!mySurface.IsDone()) {
    //           std::cout << "Uh-Oh\n";
    //           exit(1);
    //      }
    //     Handle(Geom_BSplineSurface) C = mySurface.Surface();
    //      BRep_Builder FaceMaker;
    //      TopoDS_Face F;
    //      FaceMaker.MakeFace(F,C,0.1);
    //     BRepTools::Write(F,"spline.brep");
    
    //     gp_Pnt P1(0.0,2.,2.);
    //     gp_Pnt P2;
    //
    // Without Handle
    //     TColgp_HArray1OfPnt points(0,10);
    //     points.Init(P1);
    //     P2 = points.Value(0);
    //     std::cout << P2.Y() << std::endl;
    //
    //
    //     // With Handle
    //     Handle(TColgp_HArray1OfPnt) myPointsHandle = new TColgp_HArray1OfPnt(1,3);
    //     Handle(TColStd_HArray1OfReal) myRealsHandle = new TColStd_HArray1OfReal(1,3);
    //
    //
    //
    //     myPointsHandle->Init(P1) ;
    //     P1.SetX(1.0);
    //     myPointsHandle->SetValue(1,P1);
    //     P1.SetY(0.1);
    //     myPointsHandle->SetValue(2,P1);
    //     for (int n=1;n<=3;++n) {
    //         myRealsHandle->SetValue(n,n);
    //         P2 = myPointsHandle->Value(n);
    //         std::cout << P2.X() << ' ' << P2.Y() << ' ' << P2.Z() << std::endl;
    //     }
    //
    //      GeomAPI_Interpolate Interp(myPointsHandle,myRealsHandle,false,0.0005);
    //      Interp.Perform();
    //      if (!Interp.IsDone()) {
    //          std::cout << "Uh-Oh\n";
    //          exit(1);
    //      }
    //      Handle(Geom_BSplineCurve) C = Interp.Curve();
    
    
    
    
    //
    //     BRepPrimAPI_MakeSphere my_sphere(50.);
    //     my_sphere.Build();
    //     ColladaAPI_Writer myColladaWriter;
    //     myColladaWriter.Write(my_sphere.Shape(), "sphere.dae");
    
    //     BRepPrimAPI_MakeBox mainBox(10.,10.,10.);
    //     mainBox.Build();
    //     assert(mainBox.IsDone());
    //     gp_Pnt P1(10.,2.,2.);
    //     BRepPrimAPI_MakeBox rightBox(P1,6.,6.,6.);
    //      rightBox.Build();
    //      assert(rightBox.IsDone());
    //      //TopoDS_Shape shp_result = BRepAlgoAPI_Fuse(mainBox,rightBox);
    //
    //     BRepTools::Write(mainBox.Solid(),"mainBoxSolid.brep");
    //     //BRepTools::Write(mainBox.Shape(),"mainBoxShape.brep");
    //
    // //     BRep_Builder aBuilder;
    // //  TopoDS_Shape aShape;
    // //     Standard_Boolean result = BRepTools::Read(aShape,"mainBoxSolid.brep",aBuilder);
    // //     std::cout << result << ' ' << !aShape.IsNull() << std::endl;
    // //     assert(result);
    // //     assert(!aShape.IsNull());
    // //
    //     TopoDS_Solid box1 = (TopoDS_Solid) BRepPrimAPI_MakeBox(10.,10.,10.);
    //     assert(!box1.IsNull());
    //     // the second shape: a smaller cube, a corner at the origin
    //     //TopoDS_Solid box2 = (TopoDS_Solid) BRepPrimAPI_MakeBox(5.,5.,5.);
    //     TopoDS_Solid box2 = (TopoDS_Solid) BRepPrimAPI_MakeBox(P1,6.,6.,6.);
    //     assert(!box2.IsNull());
    //     // boolean cut
    //     BOPTools_DSFiller notsure;
    //     TopoDS_Shape shp_result = BRepAlgoAPI_Fuse(box1,box2,notsure);
    //
    //     TopoDS_Builder C;
    //     TopoDS_CompSolid D;
    //     C.MakeCompSolid(D);
    //     C.Add(D,box1);
    //     C.Add(D,box2);
    //
    //
    //
    //
    //     TopTools_ListOfShape aRL;
    //     aRL.Append(box1);
    //     aRL.Append(box2);
    //
    //     TopoDS_Shape composite;
    //
    //
    //     assert(!shp_result.IsNull());
    // //     TopoDS_Shape shp_result2 = BRepAlgoAPI_Cut(box2,box1);
    // //     assert(!shp_result2.IsNull());
    //     BRepTools::Write(D,"shp_result.brep");
    // //     BRepTools::Write(shp_result2,"shp_result2.brep");
    
    return;
}
