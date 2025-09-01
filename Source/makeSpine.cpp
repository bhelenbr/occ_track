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

//#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_CompCurve.hxx>
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

// To compile
// g++ spline.cpp -I${HOME}/Packages/include -I${HOME}/Packages/include/oce -L${HOME}/Packages/lib -lTKBRep -lTKPrim -lTKernel -lTKMath -lTKTopAlgo -lTKBO -lTKGeomAlgo -lTKG3D

// To find what library has a routine
// nm -o libTK* | grep GeomTools | grep Read | grep ' T '

// Tutorial is in opencascade directory under doc / tutorial

// Doxygen is in that folder as well

// Use spotlight to find header files or search in XCode File

/* Read file to see how many points there are */
void makeTrackSpine(std::string const sourceFolder, std::string const destinationFolder, Handle(Geom_BSplineCurve)& C) {
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
    makeSpine(arcLength,points,C);
    outputBSplineCurve(destinationFolder, nPoints, C);
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

void outputBSplineCurve(const std::string destinationFolder, int nPts, Handle(const Geom_BSplineCurve) const C) {
    /* Test track curve */
    Standard_Real sBgn = C->FirstParameter();
    Standard_Real sEnd = C->LastParameter();
    Standard_Real ds = (sEnd-sBgn)/(nPts-1);
    
    std::ofstream track_test;
    track_test.open(destinationFolder +"/track_spine_pts.txt");
    for (int i=0;i<nPts;++i) {
        double sx = sBgn +ds*i;
        gp_Pnt P;
        gp_Vec T, d2T;
        C->D2(sx, P, T, d2T);
#ifdef UNITY
        track_test << sx << ' ' << P.X() << ' ' << P.Z() << ' ' << -P.Y() << ' ' << T.X() << ' ' << T.Z() << ' ' << -T.Y() << ' ' << d2T.X() << ' ' << d2T.Z() << ' ' <<  -d2T.Y() <<<< std::endl;
#else
        track_test << sx << ' ' << P.X() << ' ' << P.Y() << ' ' << P.Z() << ' ' << T.X() << ' ' << T.Y() << ' ' <<  T.Z() << ' ' << d2T.X() << ' ' << d2T.Y() << ' ' <<  d2T.Z() << std::endl;
#endif
    }
    track_test.close();
    
    /* Write spline file */
    std::ofstream mySplineFile;
    mySplineFile.open(destinationFolder +"/spine_spline.txt");
    GeomTools::Write(C,mySplineFile);
    mySplineFile.close();
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
