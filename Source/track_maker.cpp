//#define COLLADA

#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <TColgp_HArray2OfPnt.hxx>
#include <TColStd_HArray1OfReal.hxx>
#include <TColStd_HArray1OfInteger.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_BSplineCurve.hxx>
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
#ifdef COLLADA
#include <ColladaAPI_Writer.hxx>
#endif
#include <STEPControl_Reader.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Curve2D.hxx>
#include <Geom_Curve.hxx>
#include <GeomTools.hxx>
#include <gp_Ax1.hxx>
#include <BRepOffset.hxx>
#include <Geom_Surface.hxx>

#include <blitz/array.h>
#include <fstream>
#include <sstream>

#include "timeDeriv.h"

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

#define MAKE_TRACK

void timeDeriv(blitz::array<double,1> y, Handle(Geom_BSplineSurface) trackSurface);

int main(int argc, char **argv) {
	
	
#ifdef MAKE_TRACK
	const int nCrossPts = 100;
	
	/* Read file to see how many points there are */
	std::ifstream  arclengthFile;
	arclengthFile.open("arclength.txt");
	int nPoints = 0;
	double sx;
	while (arclengthFile >> sx)
		++nPoints;
	arclengthFile.close();
	
	/* load arclengths */
	arclengthFile.open("arclength.txt");
	blitz::Array<double,1> arclength(nPoints);
	for (int n=0; n < nPoints; ++n) {
		arclengthFile >> arclength(n);
	}
	arclengthFile.close();
	
	std::ifstream  pointsFile;
	pointsFile.open("track.txt");
	blitz::Array<double,2> points(nPoints,3);
	for (int n=0; n < nPoints; ++n) {
		pointsFile >> points(n,0) >> points(n,1) >> points(n,2);
	}
	pointsFile.close();

	std::ifstream  profileFile;
	struct profileData {
		std::string filename;
		double sx;
	};
	
	profileFile.open("profiles.txt");
	int nProfiles = 0;
	std::string FileName;
	while (profileFile >> FileName >> sx) {
		++nProfiles;
	}
	profileFile.close();
	
	blitz::Array<profileData,1> profiles(nProfiles);
	profileFile.open("profiles.txt");
	for (int n=0; n < nProfiles; ++n)	{
		profileFile >> profiles(n).filename >> profiles(n).sx;
	}

	/* Make spine of track */
	// With Handle
	Handle(TColgp_HArray1OfPnt) myPointsHandle = new TColgp_HArray1OfPnt(1,nPoints);
	Handle(TColStd_HArray1OfReal) mySxHandle = new TColStd_HArray1OfReal(1,nPoints);
	

	for (int n=0;n<nPoints;++n) {
		gp_Pnt P1(points(n,0),points(n,1),points(n,2));
		myPointsHandle->SetValue(n+1,P1);
		mySxHandle->SetValue(n+1,arclength(n));
	}
	GeomAPI_Interpolate Interp(myPointsHandle,mySxHandle,false,0.0005);
	Interp.Perform();
	if (!Interp.IsDone()) {
		std::cout << "Uh-Oh Spine of track not made\n";
		exit(1);
	}
	Handle(Geom_BSplineCurve) C = Interp.Curve();
	
	/* Test track curve */
	std::ofstream track_test;
	track_test.open("Results/track_test_pts");
	for (sx = arclength(0); sx < arclength(nPoints-1);sx += 0.2) {
		gp_Pnt P2;
		gp_Vec dir;
		C->D1(sx, P2, dir);
		track_test << P2.X() << ' ' << P2.Z() << ' ' << -P2.Y() << std::endl;
	}
	track_test.close();
	
	TColgp_Array2OfPnt myPoints(1,nProfiles+1,1,nCrossPts);
	TColStd_Array1OfReal UKnots(1,nProfiles);
	TColStd_Array1OfInteger UMults(1,nProfiles);
	TColStd_Array1OfReal VKnots(1,nCrossPts);
	TColStd_Array1OfInteger VMults(1,nCrossPts);
	
	/* Now load in profiles and use to make 2D array of points for B-Spline Surface */
	for(int row=0;row<nProfiles;++row) {
		/* Load in the step file */
		STEPControl_Reader reader;
		IFSelect_ReturnStatus stat = reader.ReadFile(profiles(row).filename.c_str());
		if (!stat) {
			std::cout << "There was a problem " << stat << std::endl;
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
			sx = profiles(row).sx;
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

			
//			ofstream fout;
//			std::ostringstream nstr;
//			nstr.str("");
//			nstr << "Results/profile"  << row << ".pts";
//			fout.open(nstr.str());

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
				
				// fout << P2.X() << ' ' << P2.Y() << ' ' << P2.Z() << std::endl;
				myPoints.SetValue(row+2,col+1,P2);
			}
//			fout.close();
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
	BRepTools::Write(F,"Results/theTrack.brep");
	
	/* Write dae file */
#ifdef COLLADA
	ColladaAPI_Writer myColladaWriter;
	myColladaWriter.SetCoefficient(1e-4);
	myColladaWriter.Write(F, "Results/theTrack.dae");
#endif
	
	/* Write surface file */
	std::ofstream mySurfaceFile;
	mySurfaceFile.open("Results/surface.txt");
	GeomTools::Write(trackSurface,mySurfaceFile);
	mySurfaceFile.close();

	timeIntegrate(trackSurface);
	
#endif
	
#ifdef STRAIGHT
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
	
	blitz::Array<double,2> profile_pnts(NptsBottom+2*NptsArc +2*NptsSide,2);

	
	/* First side */
	int pntcnt = 0;
	for (int pnt=0; pnt < NptsSide; ++pnt) {
		profile_pnts(pntcnt,0) = -width*0.5 ;
		profile_pnts(pntcnt++,1) = height- (height-radius)/NptsSide*pnt;
	}
	/* Arc */
	for (int pnt=0; pnt < NptsArc; ++pnt) {
		profile_pnts(pntcnt,0) = -width*0.5+radius -radius*cos(M_PI*pnt/(2.*NptsArc));
		profile_pnts(pntcnt++,1) = radius - radius*sin(M_PI*pnt/(2.*NptsArc));
	}
	/* Bottom */
	for (int pnt=0; pnt < NptsBottom; ++pnt) {
		profile_pnts(pntcnt,0) = -width*0.5+radius +pnt*(width-2*radius)/NptsBottom;
		profile_pnts(pntcnt++,1) = 0.0;
	}
	/* Arc */
	for (int pnt=0; pnt < NptsArc; ++pnt) {
		profile_pnts(pntcnt,0) = width*0.5-radius +radius*sin(M_PI*pnt/(2.*NptsArc));
		profile_pnts(pntcnt++,1) = radius - radius*cos(M_PI*pnt/(2.*NptsArc));
	}
	/* Final side */
	for (int pnt=0; pnt < NptsSide; ++pnt) {
		profile_pnts(pntcnt,0) = width*0.5 ;
		profile_pnts(pntcnt++,1) = radius +(height-radius)/NptsSide*pnt;
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
			P1 = gp_Pnt(profile_pnts(col-1,0),y,-z-profile_pnts(col-1,1));
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
	
	
	/* Just exploring this surface a bit */
	TopExp_Explorer anEdgeExplorer(F, TopAbs_EDGE);
	
	while(anEdgeExplorer.More()) {
		TopoDS_Edge anEdge = TopoDS::Edge(anEdgeExplorer.Current());
		std::cout << "I found an edge" << std::endl;
		anEdgeExplorer.Next();
	}
	
	/* Write brep file */
	BRepTools::Write(F,"spline.brep");

	/* Write dae file */
#ifdef COLLADA
	ColladaAPI_Writer myColladaWriter;
	myColladaWriter.SetCoefficient(1e-5);
	myColladaWriter.Write(F, "spline.dae");
#endif
#endif

	
//	gp_Vec aPrismVec(0, 0, 1.0);
//	TopoDS_Shape myBody = BRepPrimAPI_MakePrism(F, aPrismVec);
//	BRepTools::Write(myBody,"myBody.brep");
	

// 	GeomAPI_PointsToBSplineSurface mySurface;
// 	mySurface.Init(myPoints);
// 	if (!mySurface.IsDone()) {
//   		std::cout << "Uh-Oh\n";
//   		exit(1);
//  	}
// 	Handle(Geom_BSplineSurface) C = mySurface.Surface(); 
//  	BRep_Builder FaceMaker;
//  	TopoDS_Face F;
//  	FaceMaker.MakeFace(F,C,0.1);
// 	BRepTools::Write(F,"spline.brep"); 






// 	gp_Pnt P1(0.0,2.,2.);
// 	gp_Pnt P2;
// 
	// Without Handle
// 	TColgp_HArray1OfPnt points(0,10);
// 	points.Init(P1);
// 	P2 = points.Value(0);
// 	std::cout << P2.Y() << std::endl;
// 
// 	
// 	// With Handle
// 	Handle(TColgp_HArray1OfPnt) myPointsHandle = new TColgp_HArray1OfPnt(1,3);
// 	Handle(TColStd_HArray1OfReal) myRealsHandle = new TColStd_HArray1OfReal(1,3);
// 	
// 	
// 	
// 	myPointsHandle->Init(P1) ;
// 	P1.SetX(1.0);
// 	myPointsHandle->SetValue(1,P1);
// 	P1.SetY(0.1);
// 	myPointsHandle->SetValue(2,P1);
// 	for (int n=1;n<=3;++n) {
// 		myRealsHandle->SetValue(n,n);
// 		P2 = myPointsHandle->Value(n);
// 		std::cout << P2.X() << ' ' << P2.Y() << ' ' << P2.Z() << std::endl;
// 	}
// 	
//  	GeomAPI_Interpolate Interp(myPointsHandle,myRealsHandle,false,0.0005);
//  	Interp.Perform();
//  	if (!Interp.IsDone()) {
//  		std::cout << "Uh-Oh\n";
//  		exit(1);
//  	}
//  	Handle(Geom_BSplineCurve) C = Interp.Curve();
 	
 	
 	
 	
//   
// 	BRepPrimAPI_MakeSphere my_sphere(50.);
// 	my_sphere.Build();
// 	ColladaAPI_Writer myColladaWriter;
// 	myColladaWriter.Write(my_sphere.Shape(), "sphere.dae");
  
//     BRepPrimAPI_MakeBox mainBox(10.,10.,10.);
//     mainBox.Build();
//     assert(mainBox.IsDone());
// 	gp_Pnt P1(10.,2.,2.);
// 	BRepPrimAPI_MakeBox rightBox(P1,6.,6.,6.);
//  	rightBox.Build();
//  	assert(rightBox.IsDone());    
//  	//TopoDS_Shape shp_result = BRepAlgoAPI_Fuse(mainBox,rightBox);
//  	
// 	BRepTools::Write(mainBox.Solid(),"mainBoxSolid.brep"); 
// 	//BRepTools::Write(mainBox.Shape(),"mainBoxShape.brep"); 
// 		
// // 	BRep_Builder aBuilder;
// //  TopoDS_Shape aShape;
// // 	Standard_Boolean result = BRepTools::Read(aShape,"mainBoxSolid.brep",aBuilder);
// // 	std::cout << result << ' ' << !aShape.IsNull() << std::endl;
// // 	assert(result);
// //     assert(!aShape.IsNull());
// // 	    
// 	TopoDS_Solid box1 = (TopoDS_Solid) BRepPrimAPI_MakeBox(10.,10.,10.);
// 	assert(!box1.IsNull());
// 	// the second shape: a smaller cube, a corner at the origin
// 	//TopoDS_Solid box2 = (TopoDS_Solid) BRepPrimAPI_MakeBox(5.,5.,5.);
// 	TopoDS_Solid box2 = (TopoDS_Solid) BRepPrimAPI_MakeBox(P1,6.,6.,6.);
// 	assert(!box2.IsNull());
// 	// boolean cut
// 	BOPTools_DSFiller notsure;	
// 	TopoDS_Shape shp_result = BRepAlgoAPI_Fuse(box1,box2,notsure);
// 	
// 	TopoDS_Builder C;
// 	TopoDS_CompSolid D;
// 	C.MakeCompSolid(D);
// 	C.Add(D,box1);
// 	C.Add(D,box2);
// 	
// 	
// 	
// 	
// 	TopTools_ListOfShape aRL;
// 	aRL.Append(box1);
// 	aRL.Append(box2);
// 	
// 	TopoDS_Shape composite;
// 	
// 	
// 	assert(!shp_result.IsNull());
// //     TopoDS_Shape shp_result2 = BRepAlgoAPI_Cut(box2,box1);
// //     assert(!shp_result2.IsNull());
//     BRepTools::Write(D,"shp_result.brep"); 
// //     BRepTools::Write(shp_result2,"shp_result2.brep");

	return 0;
}
