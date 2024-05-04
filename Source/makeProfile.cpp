//
//  makeProfile.cpp
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

void makeProfiles(std::string const sourceFile, std::string const destinationFolder) {
    
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
    
    
    /* Playing around with this */
//    BRepAdaptor_CompCurve myComposite(myWireProfile);
//    Standard_Real U1 = myComposite.FirstParameter();
//    Standard_Real U2 = myComposite.LastParameter();
//    gp_Pnt aPnt;
//
//    const int nPts = 100;
//    double du = (U2-U1)/(nPts-1);
//    for (int i = 0; i < 100; ++i) {
//        Standard_Real U = U1 +du*i;
//        myComposite.D0(U,aPnt);
//        std::cout << aPnt.X() << ' ' << aPnt.Y() << ' ' << aPnt.Z() << std::endl;
//    }
    
    std::ofstream of;
    of.open(filename);
    BRepTools::Write(myWireProfile,of);
    
    gp_Pnt Origin (0., 0., 0.);
    TopoDS_Vertex V = BRepBuilderAPI_MakeVertex(Origin);
    BRepTools::Write(V,of);
    of.close();

//    If you need the real Geom_Curve then you have to approximate CompCurve into bspline using available approximation tool
//    (e.g. AdvApprox_ApproxAFunction, see GeomLib::BuildCurve3d).
        
        
}

void makeStraightTransitionProfile(std::string const filename,double A, double B, double R1, double KH, double BW, double Rr, double BH) {
    const gp_Dir xAxis(1.0,0.0,0.0);
    const gp_Dir yAxis(0.0,1.0,0.0);
    const gp_Dir zAxis(0.0,0.0,1.0);
    
    gp_Pnt aPnt, aPnt1, aPnt2;
    gp_Ax2 anAx2;
    gp_Elips anEllipse;
    
    aPnt = gp_Pnt(0.0,A,0.0);
    if (B > A) {
        anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
        anEllipse = gp_Elips(anAx2,B,A);
    }
    else {
        anAx2 = gp_Ax2(aPnt, zAxis, yAxis);
        anEllipse = gp_Elips(anAx2,A,B);
    }
    
    double BT=0.0,AT=0.0,MY;
    double dx=0.0;

    /* Find tangency point between ellipse and circle */
    BT = 0.9*B;
    for (int iter = 0; iter < 100; ++iter) {
        // x^2/B^2 +(y-A)^2/A^2 = 1
        AT = A*(1.-sqrt(1.-pow(BT/B,2)));
        // 2x/B^2 +2*(y-A)*dy/dx/A^2 = 0
        double slope = -A*A*(2*BT/(B*B))/(2*(AT-A));
        double inv = -1/slope;
        // R1^2 = dy^2+dx^2
        // R1^2/(1+dy/dx^2) = dx^2;
        dx = sqrt(R1*R1/(1+inv*inv));
        double Biter = BT -dx +R1;
        double BT_old = BT;
        BT = (B-Biter)+BT;
        if (abs(BT_old-BT) < 1.0e-14) {
            break;
        }
        std::cout << BT << std::endl;
    }
    BT = -BT;
    MY = AT+sqrt(R1*R1-dx*dx);
    std::cout << BT << ' ' << AT << ' ' << MY << std::endl;
       
    aPnt1 = gp_Pnt(-B,KH,0.0);
    aPnt2 = gp_Pnt(-B,MY,0.0);
    Handle(Geom_TrimmedCurve) aSegment1    = GC_MakeSegment(aPnt1, aPnt2);
    TopoDS_Edge lEdge = BRepBuilderAPI_MakeEdge(aSegment1);
    
    aPnt = gp_Pnt(-B+R1,MY,0.0);
    anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
    gp_Circ aCirc(anAx2,R1);
    aPnt1 = aPnt2;
    aPnt2 = gp_Pnt(BT,AT,0.0);
    Handle(Geom_TrimmedCurve) aArcOfCircle1 = GC_MakeArcOfCircle(aCirc,aPnt1,aPnt2,Standard_True);
    TopoDS_Edge r1Edge = BRepBuilderAPI_MakeEdge(aArcOfCircle1);
              
    aPnt1 = aPnt2;
    gp_Pnt origin(0.0,0.0,0.0);
    GC_MakeArcOfEllipse eMaker(anEllipse,aPnt1,origin,true);
    Handle(Geom_TrimmedCurve) eHandle = eMaker.Value();
    TopoDS_Edge eEdge = BRepBuilderAPI_MakeEdge(eHandle);
    
    aPnt2 = gp_Pnt(BW-Rr,0.0,0.0);
    Handle(Geom_TrimmedCurve) aSegment2    = GC_MakeSegment(origin, aPnt2);
    TopoDS_Edge bEdge = BRepBuilderAPI_MakeEdge(aSegment2);
    
    aPnt = gp_Pnt(BW-Rr,Rr,0.0);
    anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
    aCirc = gp_Circ(anAx2,Rr);
    Handle(Geom_TrimmedCurve) aArcOfCircle2 = GC_MakeArcOfCircle(aCirc,3.*M_PI/2.,2.*M_PI,true);
    TopoDS_Edge r2Edge = BRepBuilderAPI_MakeEdge(aArcOfCircle2);
    
    aPnt1 = gp_Pnt(BW,Rr,0.0);
    aPnt2 = gp_Pnt(BW,BH,0.0);
    Handle(Geom_TrimmedCurve) aSegment3    = GC_MakeSegment(aPnt1, aPnt2);
    TopoDS_Edge rEdge = BRepBuilderAPI_MakeEdge(aSegment3);
    
    
    BRepBuilderAPI_MakeWire mkWire;
    mkWire.Add(lEdge);
    mkWire.Add(r1Edge);
    mkWire.Add(eEdge);
    mkWire.Add(bEdge);
    mkWire.Add(r2Edge);
    mkWire.Add(rEdge);
    TopoDS_Wire myWireProfile = mkWire.Wire();
    
    std::ofstream of;
    of.open(filename);
    BRepTools::Write(myWireProfile,of);
    
    TopoDS_Vertex V = BRepBuilderAPI_MakeVertex(origin);
    BRepTools::Write(V,of);
    of.close();
}

void makeCurveProfile(std::string const filename,double A, double B, double HZ, double WZ, double KH, double BW, double R2, double BH, double B1, double Rr) {
    const gp_Dir xAxis(1.0,0.0,0.0);
    const gp_Dir yAxis(0.0,1.0,0.0);
    const gp_Dir zAxis(0.0,0.0,1.0);
    
    gp_Pnt aPnt, aPnt1, aPnt2;
    gp_Ax2 anAx2;
    gp_Elips anEllipse;
    
    aPnt = gp_Pnt(0.0,A,0.0);
    if (B > A) {
        anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
        anEllipse = gp_Elips(anAx2,B,A);
    }
    else {
        anAx2 = gp_Ax2(aPnt, zAxis, yAxis);
        anEllipse = gp_Elips(anAx2,A,B);
    }
    
    double BT=0.0,AT=0.0;
    double dx,dy;
    double slope_target = WZ/HZ;
    /* Find tangency point between ellipse and circle */
    AT = 0.95*A;
    for (int iter = 0; iter < 100; ++iter) {
        // x^2/B^2 +(y-A)^2/A^2 = 1
        BT = B*sqrt(1-pow((AT-A)/A,2));
        std::cout << BT/B << std::endl;
        // 2x dx/dy/B^2 +2*(y-A)/A^2 = 0
        double slope = -2*(AT-A)/(A*A*2*BT)*B*B;
        double AT_old = AT;
        AT = (slope-slope_target)*(A*A*2*BT)/(B*B*2) +AT;
        if (abs(AT_old-AT) < 1.0e-14) {
            break;
        }
        std::cout << AT << std::endl;
    }
    std::cout << BT << ' ' << AT << std::endl;
       
    aPnt = gp_Pnt(-BT-WZ+R2,AT+HZ,0.0);
    anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
    gp_Circ aCirc(anAx2,R2);
    dy = KH-AT-HZ;
    dx = sqrt(R2*R2-dy*dy);
    aPnt1 =  gp_Pnt(-BT-WZ+R2-dx,KH,0.0);
    aPnt2 =  gp_Pnt(-BT-WZ,AT+HZ,0.0);
    Handle(Geom_TrimmedCurve) aArcOfCircleR2 = GC_MakeArcOfCircle(aCirc,aPnt1,aPnt2,Standard_True);
    TopoDS_Edge r2Edge = BRepBuilderAPI_MakeEdge(aArcOfCircleR2);
  
    aPnt1 = aPnt2;
    aPnt2 = gp_Pnt(-BT,AT,0.0);
    Handle(Geom_TrimmedCurve) aSegment1    = GC_MakeSegment(aPnt1, aPnt2);
    TopoDS_Edge lEdge = BRepBuilderAPI_MakeEdge(aSegment1);
                                                                 
    GC_MakeArcOfEllipse eMaker(anEllipse,aPnt2,M_PI,true);
    Handle(Geom_TrimmedCurve) eHandle = eMaker.Value();
    TopoDS_Edge eEdge = BRepBuilderAPI_MakeEdge(eHandle);
    
    gp_Pnt origin(0.0,0.0,0.0);
    aPnt2 = gp_Pnt(B1,0.0,0.0);
    Handle(Geom_TrimmedCurve) aSegment2    = GC_MakeSegment(origin, aPnt2);
    TopoDS_Edge bEdge = BRepBuilderAPI_MakeEdge(aSegment2);
    
    aPnt = gp_Pnt(B1,Rr,0.0);
    anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
    aCirc = gp_Circ(anAx2,Rr);
    
    /* Find tangent point */
    dx = sqrt(pow(BW-B1,2)+pow(BH-Rr,2));
    dy = sqrt(dx*dx-Rr*Rr);
    double beta = atan((BW-B1)/(BH-Rr));
    double alpha = atan(dy/Rr);
    double theta = M_PI -alpha -beta;
    aPnt1 = aPnt2;
    aPnt2 = gp_Pnt(B1+Rr*sin(theta),Rr-Rr*cos(theta),0.0);

    Handle(Geom_TrimmedCurve) aArcOfCircleRr = GC_MakeArcOfCircle(aCirc,aPnt1,aPnt2,true);
    TopoDS_Edge rrEdge = BRepBuilderAPI_MakeEdge(aArcOfCircleRr);
    
    aPnt1 = aPnt2;
    aPnt2 = gp_Pnt(BW,BH,0.0);
    Handle(Geom_TrimmedCurve) aSegment3    = GC_MakeSegment(aPnt1, aPnt2);
    TopoDS_Edge rEdge = BRepBuilderAPI_MakeEdge(aSegment3);
    
    
    BRepBuilderAPI_MakeWire mkWire;
    mkWire.Add(r2Edge);
    mkWire.Add(lEdge);
    mkWire.Add(eEdge);
    mkWire.Add(bEdge);
    mkWire.Add(rrEdge);
    mkWire.Add(rEdge);
    TopoDS_Wire myWireProfile = mkWire.Wire();
    
    std::ofstream of;
    of.open(filename);
    BRepTools::Write(myWireProfile,of);
    
    TopoDS_Vertex V = BRepBuilderAPI_MakeVertex(origin);
    BRepTools::Write(V,of);
    of.close();
}

void makeCurveTransitionProfile(std::string const filename,double A, double B, double HZ, double R1, double KH, double BW, double R2, double BH, double B1, double Rr) {
    const gp_Dir xAxis(1.0,0.0,0.0);
    const gp_Dir yAxis(0.0,1.0,0.0);
    const gp_Dir zAxis(0.0,0.0,1.0);
    
    gp_Pnt aPnt, aPnt1, aPnt2;
    gp_Ax2 anAx2;
    gp_Elips anEllipse;
    
    aPnt = gp_Pnt(0.0,A,0.0);
    if (B > A) {
        anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
        anEllipse = gp_Elips(anAx2,B,A);
    }
    else {
        anAx2 = gp_Ax2(aPnt, zAxis, yAxis);
        anEllipse = gp_Elips(anAx2,A,B);
    }
    
    double BT=0.0,AT=0.0,MY=0.0;
    double dx=0.0,dy=0.0;

    /* Find tangency point between ellipse and circle */
    BT = 0.9*B;
    for (int iter = 0; iter < 100; ++iter) {
        // x^2/B^2 +(y-A)^2/A^2 = 1
        AT = A*(1.-sqrt(1.-pow(BT/B,2)));
        // 2x/B^2 +2*(y-A)*dy/dx/A^2 = 0
        double slope = -A*A*(2*BT/(B*B))/(2*(AT-A));
        double inv = -1/slope;
        // R1^2 = dy^2+dx^2
        // R1^2/(1+dy/dx^2) = dx^2;
        dx = sqrt(R1*R1/(1+inv*inv));
        double Biter = BT -dx +R1;
        double BT_old = BT;
        BT = (B-Biter)+BT;
        if (abs(BT_old-BT) < 1.0e-14) {
            break;
        }
        std::cout << BT << std::endl;
    }
    BT = -BT;
    MY = AT+sqrt(R1*R1-dx*dx);
    std::cout << BT << ' ' << AT << ' ' << MY << std::endl;
       
    aPnt = gp_Pnt(-B+R2,MY+HZ,0.0);
    anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
    gp_Circ aCirc(anAx2,R2);
    dy = KH-MY-HZ;
    dx = sqrt(R2*R2-dy*dy);
    aPnt1 =  gp_Pnt(-B+R2-dx,KH,0.0);
    aPnt2 =  gp_Pnt(-B,MY+HZ,0.0);
    Handle(Geom_TrimmedCurve) aArcOfCircleR2 = GC_MakeArcOfCircle(aCirc,aPnt1,aPnt2,Standard_True);
    TopoDS_Edge r2Edge = BRepBuilderAPI_MakeEdge(aArcOfCircleR2);
  
    aPnt1 = aPnt2;
    aPnt2 = gp_Pnt(-B,MY,0.0);
    Handle(Geom_TrimmedCurve) aSegment1    = GC_MakeSegment(aPnt1, aPnt2);
    TopoDS_Edge lEdge = BRepBuilderAPI_MakeEdge(aSegment1);
    
    aPnt = gp_Pnt(-B+R1,MY,0.0);
    anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
    aCirc = gp_Circ(anAx2,R1);
    aPnt1 = aPnt2;
    aPnt2 = gp_Pnt(BT,AT,0.0);
    Handle(Geom_TrimmedCurve) aArcOfCircleR1 = GC_MakeArcOfCircle(aCirc,aPnt1,aPnt2,Standard_True);
    TopoDS_Edge r1Edge = BRepBuilderAPI_MakeEdge(aArcOfCircleR1);
                                                                 
    GC_MakeArcOfEllipse eMaker(anEllipse,aPnt2,3./2.*M_PI,true);
    Handle(Geom_TrimmedCurve) eHandle = eMaker.Value();
    TopoDS_Edge eEdge = BRepBuilderAPI_MakeEdge(eHandle);
    
    gp_Pnt origin(0.0,0.0,0.0);
    aPnt2 = gp_Pnt(B1,0.0,0.0);
    Handle(Geom_TrimmedCurve) aSegment2    = GC_MakeSegment(origin, aPnt2);
    TopoDS_Edge bEdge = BRepBuilderAPI_MakeEdge(aSegment2);
    
    aPnt = gp_Pnt(B1,Rr,0.0);
    anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
    aCirc = gp_Circ(anAx2,Rr);
    
    /* Find tangent point */
    dx = sqrt(pow(BW-B1,2)+pow(BH-Rr,2));
    dy = sqrt(dx*dx-Rr*Rr);
    double beta = atan((BW-B1)/(BH-Rr));
    double alpha = atan(dy/Rr);
    double theta = M_PI -alpha -beta;
    aPnt1 = aPnt2;
    aPnt2 = gp_Pnt(B1+Rr*sin(theta),Rr-cos(theta),0.0);

    Handle(Geom_TrimmedCurve) aArcOfCircleRr = GC_MakeArcOfCircle(aCirc,aPnt1,aPnt2,true);
    TopoDS_Edge rrEdge = BRepBuilderAPI_MakeEdge(aArcOfCircleRr);
    aPnt1 = aPnt2;
    aPnt2 = gp_Pnt(BW,BH,0.0);
    Handle(Geom_TrimmedCurve) aSegment3    = GC_MakeSegment(aPnt1, aPnt2);
    TopoDS_Edge rEdge = BRepBuilderAPI_MakeEdge(aSegment3);
    
    
    BRepBuilderAPI_MakeWire mkWire;
    mkWire.Add(r2Edge);
    mkWire.Add(lEdge);
    mkWire.Add(r1Edge);
    mkWire.Add(eEdge);
    mkWire.Add(bEdge);
    mkWire.Add(rrEdge);
    mkWire.Add(rEdge);
    TopoDS_Wire myWireProfile = mkWire.Wire();
    
    std::ofstream of;
    of.open(filename);
    BRepTools::Write(myWireProfile,of);
    
    TopoDS_Vertex V = BRepBuilderAPI_MakeVertex(origin);
    BRepTools::Write(V,of);
    of.close();
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
