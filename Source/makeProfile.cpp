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
#include <GeomAPI_PointsToBSpline.hxx>
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
#include <TColgp_Array1OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>

#include <STEPControl_Reader.hxx>
// #include <BOPTools_DSFiller.hxx>

#include <array>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
// #include <assert>
#include <filesystem>
namespace fs = std::filesystem;

#include "occ_track.h"



#define VERBOSE

void makeProfiles(std::string const sourceFile, std::string const destinationFolder) {
    std::ifstream input;
    std::string line;
    std::string word;
    
    std::filesystem::create_directory(destinationFolder +"/Results");
    std::filesystem::create_directory(destinationFolder +"/Results/Profiles");
    
    input.open(sourceFile);
    
    std::getline(input, line);
    /* Read headers */
    std::map<std::string,int> hcols;
    std::stringstream hstream(line);
    int cols = 0;
    while (hstream >> word) {
        hcols[word] = cols++;
    }
    
    std::ofstream of;
    of.open(destinationFolder +"/Results/dims.txt");

    /* Input line and make profile*/
    std::vector<profileData> profiles;
    profileData myProfile;
    std::string Type,Dir,num;
    char com;
    int BRK;
    
    const double RR = 125.0; // Right Radius
    double Sx,Az,A,B,AT,BT,H1,MY,MN,R1,INCL,HZ,WZ,H2,R2,KH,UEB,B1,AK1,BW1,AK3,HW2,BHS,BW,BGES,CE,DFA,R1x,R1y,R2x,R2y,Offset,WR2;
    int rows = 0;
    int prevProfileType = 0, profileType;
#ifdef VERBOSE
           std::cout << "Sx Type Dir A B HZ WZ KH BW R2 BH B1 R1 RR CE" << std::endl;
#endif
    
    while (std::getline(input, line)) {
        std::istringstream dstream(line);
        std::getline(dstream,Type,',');
        std::getline(dstream,Dir,',');
        dstream >> Sx >> com >> Az >> com >> A >> com >> B >> com >> AT >> com >> BT >> com >> H1 >> com >> MY >> com >> MN >> com >> R1 >> com >> INCL >> com;
        dstream >> HZ >> com >> WZ >> com >> H2 >> com >> R2 >> com >> KH >> com >> UEB >> com >> B1 >> com >> AK1 >> com >> BW1 >> com >> AK3 >> com >> HW2 >> com >> BHS >> com;
        dstream >> BW >> com >> BGES >> com >> CE >> com >> DFA >> com >> R1x >> com >> R1y >> com >> R2x >> com >> R2y >> com >> Offset >> com >> WR2 >> com >> BRK;
        double BH = AK3+HW2;
        
        int dir = 0;
        if (Dir == "L") {
            dir = -1;
            of << Sx << ' ' << -BW-CE  << ' ' << B-CE << ' ' << BH << ' ' << KH << ' ' << BRK << std::endl;
        }
        else {
            of << Sx << ' ' << -B+CE << ' ' << BW+CE << ' ' << KH << ' ' << BH << ' ' << BRK << std::endl;
            if (Dir == "R") {
                dir = 1;
            }
        }

        /* Decide what kind of profile it is */
        std::ostringstream nstr;
        nstr << Sx << '_' << Type << '_' << Dir;
        myProfile.filename = "Results/Profiles/" +nstr.str() +".brep";
        myProfile.sx = Sx;
        myProfile.BRK = BRK;
        std::string filename = destinationFolder +"/" +myProfile.filename;

        if (abs(AT) < 1.0e-6) {
            /* If AT,BT = 0.0 then it is a straight section */
#ifdef VERBOSE
            std::cout << Sx << " 0 0 " << -1 << ' ' << B << " -1 -1 " << KH << ' ' << BW << " -1 " << BH << " -1 " << R1 << ' ' << RR << ' ' << CE << std::endl;
#endif
            makeStraightProfile(filename, Dir, B, BW, KH, BH, R1, RR, CE);
            profileType = 0;
        }
        else if (R1 < 1.0e-6) {
            /* If R1 = 0.0 then it is a curve */
#ifdef VERBOSE
           std::cout << Sx << " 3 " << dir << ' ' << A << ' ' << B << ' ' << HZ << ' ' << WZ << ' ' << KH << ' ' << BW << ' ' << R2 << ' ' << BH << ' ' << B1 << ' ' << -1 << ' ' << RR << ' ' << CE << std::endl;
#endif
            makeCurveProfile(filename, Dir, A, B, HZ, WZ, KH, BW, R2, BH, B1, RR, CE);
            profileType = 1;
        }
        else if (BW1 < 1.0e-6) {
            /* then if BW1 = 0 it is a straight transiton */
#ifdef VERBOSE
            std::cout << Sx << " 1 " << dir << ' ' << A << ' ' << B << ' ' << HZ << ' ' << -1 << ' ' << KH << ' ' << BW << ' ' << -1 << ' ' << BH << ' ' << -1 << ' ' << R1 << ' ' << RR << ' ' << CE << std::endl;
#endif
            makeStraightTransitionProfile(filename, Dir, A, B, R1, KH, HZ, R2, BW, RR, BH, CE);
            profileType = 2;
        }
        else {
            /* else it is a curve transition */
#ifdef VERBOSE
            std::cout << Sx << " 2 " << dir << ' ' << A << ' ' << B << ' ' << HZ << ' ' << -1 << ' ' << KH << ' ' << BW << ' ' << R2 << ' ' << BH << ' ' << B1 << ' ' << R1 << ' ' << RR << ' ' << CE << std::endl;
#endif
            makeCurveTransitionProfile(filename, Dir, A, B, HZ, R1, KH, BW, R2, BH, B1, RR, CE);
            profileType = 3;
        }
    
//        if (profileType != prevProfileType) {
//            myProfile.BRK = 1;
//        }
        profiles.push_back(myProfile);
        prevProfileType = profileType;
        ++rows;
    }
    of.close();
    
    of.open(destinationFolder +"/Results/profiles.txt");
    /* Output profile file*/
    for (int i=0; i < rows; ++i) {
        of << profiles[i].filename << ' ' << profiles[i].sx << ' ' << profiles[i].BRK << std::endl;
    }
    of.close();
}


void makeStraightProfile(std::string const filename, std::string const Dir, double B, double BW, double KH, double BH, double R1, double RR, double CE) {
    const gp_Dir xAxis(1.0,0.0,0.0);
    const gp_Dir yAxis(0.0,1.0,0.0);
    const gp_Dir zAxis(0.0,0.0,1.0);
    
    double drl = R1*(1.-1./sqrt(2.));
    double drr = RR*(1.-1./sqrt(2.));
    
    gp_Pnt aPnt1(-B, KH, 0);
    gp_Pnt aPnt2(-B, R1, 0);
    gp_Pnt aPnt3(-B+drl,drl,0);
    gp_Pnt aPnt4(-B+R1, 0.0, 0);
    gp_Pnt aPnt5(BW-RR, 0.0, 0);
    gp_Pnt aPnt6(BW-drr,drr,0);
    gp_Pnt aPnt7(BW, RR, 0);
    gp_Pnt aPnt8(BW, BH, 0);
    
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
    
    gp_Pnt origin(0.0,0.0,0.0);
    gp_Pnt aPnt;
    TopoDS_Wire finalWire;
    if (Dir == "L") {
        /* Mirror Edge */
        gp_Trsf Mirror;
        gp_Ax2 anAx2(origin,xAxis,zAxis);
        Mirror.SetMirror(anAx2);
        BRepBuilderAPI_Transform myTransform(myWireProfile,Mirror);
        finalWire = TopoDS::Wire(myTransform.Shape());
        aPnt = gp_Pnt(CE,0.0,0.0);
    }
    else {
        finalWire = myWireProfile;
        aPnt = gp_Pnt(-CE,0.0,0.0);
    }
    
    TopoDS_Vertex V = BRepBuilderAPI_MakeVertex(aPnt);
    BRep_Builder builder;
    TopoDS_Compound Comp;
    builder.MakeCompound(Comp);
    builder.Add(Comp,finalWire);
    builder.Add(Comp,V);
    BRepTools::Write(Comp,filename.c_str());
}

void makeStraightTransitionProfile(std::string const filename, std::string const Dir, double A, double B, double R1, double KH, double HZ, double R2, double BW, double Rr, double BH, double CE) {
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
    
    double BT=0.0,AT=0.0,MY,BT_old=0.0;
    double dx=0.0,dy=0.0;

    /* Find tangency point between ellipse and circle */
    BT = 0.9*B;
    bool converged = false;
    for (int iter = 0; iter < 100000; ++iter) {
        // x^2/B^2 +(y-A)^2/A^2 = 1
        AT = A*(1.-sqrt(1.-pow(BT/B,2)));
        // 2x/B^2 +2*(y-A)*dy/dx/A^2 = 0
        double slope = -A*A*(2*BT/(B*B))/(2*(AT-A));
        double inv = -1/slope;
        // R1^2 = dy^2+dx^2
        // R1^2/(1+dy/dx^2) = dx^2;
        dx = sqrt(R1*R1/(1+inv*inv));
        double Biter = BT -dx +R1;
        BT_old = BT;
        BT = (B-Biter)+BT;
        if (abs(BT_old-BT) < 1.0e-12) {
            converged = true;
            break;
        }
    }
    if (!converged) {
        std::cout << "difficulty converging ST" << std::endl;
        std::cout << BT -BT_old << std::endl;
    }
    BT = -BT;
    MY = AT+sqrt(R1*R1-dx*dx);
    
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
    TopoDS_Edge lEdge;
    if (HZ > 1.0e-6) {
        Handle(Geom_TrimmedCurve) aSegment1    = GC_MakeSegment(aPnt1, aPnt2);
        lEdge = BRepBuilderAPI_MakeEdge(aSegment1);
    }
    
    aPnt = gp_Pnt(-B+R1,MY,0.0);
    anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
    aCirc = gp_Circ(anAx2,R1);
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
    TopoDS_Edge rrEdge = BRepBuilderAPI_MakeEdge(aArcOfCircle2);
    
    aPnt1 = gp_Pnt(BW,Rr,0.0);
    aPnt2 = gp_Pnt(BW,BH,0.0);
    Handle(Geom_TrimmedCurve) aSegment3    = GC_MakeSegment(aPnt1, aPnt2);
    TopoDS_Edge rEdge = BRepBuilderAPI_MakeEdge(aSegment3);
    
    
    BRepBuilderAPI_MakeWire mkWire;
    mkWire.Add(r2Edge);
    if (HZ > 1.0e-6) {
        mkWire.Add(lEdge);
    }
    mkWire.Add(r1Edge);
    mkWire.Add(eEdge);
    mkWire.Add(bEdge);
    mkWire.Add(rrEdge);
    mkWire.Add(rEdge);
    TopoDS_Wire myWireProfile = mkWire.Wire();

    TopoDS_Wire finalWire;
    if (Dir == "L") {
        /* Mirror Edge */
        gp_Trsf Mirror;
        anAx2 = gp_Ax2(origin,xAxis,zAxis);
        Mirror.SetMirror(anAx2);
        BRepBuilderAPI_Transform myTransform(myWireProfile,Mirror);
        finalWire = TopoDS::Wire(myTransform.Shape());
        aPnt = gp_Pnt(CE,0.0,0.0);
    }
    else {
        finalWire = myWireProfile;
        aPnt = gp_Pnt(-CE,0.0,0.0);
    }

    TopoDS_Vertex V = BRepBuilderAPI_MakeVertex(aPnt);
    BRep_Builder builder;
    TopoDS_Compound Comp;
    builder.MakeCompound(Comp);
    builder.Add(Comp,finalWire);
    builder.Add(Comp,V);
    BRepTools::Write(Comp,filename.c_str());
}

void makeCurveProfile(std::string const filename, std::string const Dir, double A, double B, double HZ, double WZ, double KH, double BW, double R2, double BH, double B1, double RR, double CE) {
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
    
    double BT=0.0,AT=0.0,AT_old=0.0;
    double dx,dy;
    if (HZ > 1.0e-6) {
        double slope_target = WZ/HZ;
        /* Find tangency point between ellipse and circle */
        AT = 0.95*A;
        bool converged = false;
        for (int iter = 0; iter < 100; ++iter) {
            // x^2/B^2 +(y-A)^2/A^2 = 1
            BT = B*sqrt(1-pow((AT-A)/A,2));
            // 2x dx/dy/B^2 +2*(y-A)/A^2 = 0
            double slope = -2*(AT-A)/(A*A*2*BT)*B*B;
            AT_old = AT;
            AT = (slope-slope_target)*(A*A*2*BT)/(B*B*2) +AT;
            if (abs(AT_old-AT) < 1.0e-12) {
                converged = true;
                break;
            }
        }
        if (!converged) {
            std::cout << "difficulty converging CP" << std::endl;
            std::cout << AT -AT_old << std::endl;
        }
    }
    else {
        AT = A;
        BT = B;
    }
    aPnt = gp_Pnt(-BT-WZ+R2,AT+HZ,0.0);
    anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
    gp_Circ aCirc(anAx2,R2);
    dy = KH-AT-HZ;
    dx = sqrt(R2*R2-dy*dy);
    aPnt1 =  gp_Pnt(-BT-WZ+R2-dx,KH,0.0);
    aPnt2 =  gp_Pnt(-BT-WZ,AT+HZ,0.0);
    Handle(Geom_TrimmedCurve) aArcOfCircleR2 = GC_MakeArcOfCircle(aCirc,aPnt1,aPnt2,Standard_True);
    TopoDS_Edge r2Edge = BRepBuilderAPI_MakeEdge(aArcOfCircleR2);
  
    TopoDS_Edge lEdge;
    if (HZ > 1.0e-6) {
        aPnt1 = aPnt2;
        aPnt2 = gp_Pnt(-BT,AT,0.0);
        Handle(Geom_TrimmedCurve) aSegment1    = GC_MakeSegment(aPnt1, aPnt2);
        lEdge = BRepBuilderAPI_MakeEdge(aSegment1);
    }
                                                                 
    GC_MakeArcOfEllipse eMaker(anEllipse,aPnt2,M_PI,true);
    Handle(Geom_TrimmedCurve) eHandle = eMaker.Value();
    TopoDS_Edge eEdge = BRepBuilderAPI_MakeEdge(eHandle);
    
    gp_Pnt origin(0.0,0.0,0.0);
    aPnt2 = gp_Pnt(B1,0.0,0.0);
    Handle(Geom_TrimmedCurve) aSegment2    = GC_MakeSegment(origin, aPnt2);
    TopoDS_Edge bEdge = BRepBuilderAPI_MakeEdge(aSegment2);
    
    aPnt = gp_Pnt(B1,RR,0.0);
    anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
    aCirc = gp_Circ(anAx2,RR);
    
    /* Find tangent point */
    dx = sqrt(pow(BW-B1,2)+pow(BH-RR,2));
    dy = sqrt(dx*dx-RR*RR);
    double beta = atan((BW-B1)/(BH-RR));
    double alpha = atan(dy/RR);
    double theta = M_PI -alpha -beta;
    aPnt1 = aPnt2;
    aPnt2 = gp_Pnt(B1+RR*sin(theta),RR-RR*cos(theta),0.0);

    Handle(Geom_TrimmedCurve) aArcOfCircleRr = GC_MakeArcOfCircle(aCirc,aPnt1,aPnt2,true);
    TopoDS_Edge rrEdge = BRepBuilderAPI_MakeEdge(aArcOfCircleRr);
    
    aPnt1 = aPnt2;
    aPnt2 = gp_Pnt(BW,BH,0.0);
    Handle(Geom_TrimmedCurve) aSegment3    = GC_MakeSegment(aPnt1, aPnt2);
    TopoDS_Edge rEdge = BRepBuilderAPI_MakeEdge(aSegment3);
    
    
    BRepBuilderAPI_MakeWire mkWire;
    mkWire.Add(r2Edge);
    if (HZ > 1.0e-6) {
        mkWire.Add(lEdge);
    }
    mkWire.Add(eEdge);
    mkWire.Add(bEdge);
    mkWire.Add(rrEdge);
    mkWire.Add(rEdge);
    TopoDS_Wire myWireProfile = mkWire.Wire();
    
    TopoDS_Wire finalWire;
    if (Dir == "L") {
        /* Mirror Edge */
        gp_Trsf Mirror;
        anAx2 = gp_Ax2(origin,xAxis,zAxis);
        Mirror.SetMirror(anAx2);
        BRepBuilderAPI_Transform myTransform(myWireProfile,Mirror);
        finalWire = TopoDS::Wire(myTransform.Shape());
        aPnt = gp_Pnt(CE,0.0,0.0);
    }
    else {
        finalWire = myWireProfile;
        aPnt = gp_Pnt(-CE,0.0,0.0);

    }

    TopoDS_Vertex V = BRepBuilderAPI_MakeVertex(aPnt);
    BRep_Builder builder;
    TopoDS_Compound Comp;
    builder.MakeCompound(Comp);
    builder.Add(Comp,finalWire);
    builder.Add(Comp,V);
    BRepTools::Write(Comp,filename.c_str());
}

void makeCurveTransitionProfile(std::string const filename, std::string const Dir, double A, double B, double HZ, double R1, double KH, double BW, double R2, double BH, double B1, double Rr, double CE) {
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
    
    double BT=0.0,AT=0.0,MY=0.0,BT_old=0.0;
    double dx=0.0,dy=0.0;

    /* Find tangency point between ellipse and circle */
    BT = 0.9*B;
    bool converged = false;
    for (int iter = 0; iter < 50000; ++iter) {
        // x^2/B^2 +(y-A)^2/A^2 = 1
        AT = A*(1.-sqrt(1.-pow(BT/B,2)));
        // 2x/B^2 +2*(y-A)*dy/dx/A^2 = 0
        double slope = -A*A*(2*BT/(B*B))/(2*(AT-A));
        double inv = -1/slope;
        // R1^2 = dy^2+dx^2
        // R1^2/(1+dy/dx^2) = dx^2;
        dx = sqrt(R1*R1/(1+inv*inv));
        double Biter = BT -dx +R1;
        BT_old = BT;
        BT = (B-Biter)+BT;
        if (abs(BT_old-BT) < 1.0e-12) {
            converged = true;
            break;
        }
    }
    if (!converged) {
        std::cout << "difficulty converging CT" << std::endl;
        std::cout << BT -BT_old << std::endl;
    }
    BT = -BT;
    MY = AT+sqrt(R1*R1-dx*dx);
       
    aPnt = gp_Pnt(-B+R2,MY+HZ,0.0);
    anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
    gp_Circ aCirc(anAx2,R2);
    dy = KH-MY-HZ;
    dx = sqrt(R2*R2-dy*dy);
    aPnt1 =  gp_Pnt(-B+R2-dx,KH,0.0);
    aPnt2 =  gp_Pnt(-B,MY+HZ,0.0);
    Handle(Geom_TrimmedCurve) aArcOfCircleR2 = GC_MakeArcOfCircle(aCirc,aPnt1,aPnt2,Standard_True);
    TopoDS_Edge r2Edge = BRepBuilderAPI_MakeEdge(aArcOfCircleR2);
  
    
    TopoDS_Edge lEdge;
    if (HZ > 1.0e-6) {
        aPnt1 = aPnt2;
        aPnt2 = gp_Pnt(-B,MY,0.0);
        Handle(Geom_TrimmedCurve) aSegment1    = GC_MakeSegment(aPnt1, aPnt2);
        lEdge = BRepBuilderAPI_MakeEdge(aSegment1);
    }
    
    aPnt = gp_Pnt(-B+R1,MY,0.0);
    anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
    aCirc = gp_Circ(anAx2,R1);
    aPnt1 = aPnt2;
    aPnt2 = gp_Pnt(BT,AT,0.0);
    Handle(Geom_TrimmedCurve) aArcOfCircleR1 = GC_MakeArcOfCircle(aCirc,aPnt1,aPnt2,Standard_True);
    TopoDS_Edge r1Edge = BRepBuilderAPI_MakeEdge(aArcOfCircleR1);
               
    gp_Pnt origin(0.0,0.0,0.0);
    GC_MakeArcOfEllipse eMaker(anEllipse,aPnt2,origin,true);
    Handle(Geom_TrimmedCurve) eHandle = eMaker.Value();
    TopoDS_Edge eEdge = BRepBuilderAPI_MakeEdge(eHandle);
    
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
    double thetaf = -M_PI/2+theta;
    aPnt1 = aPnt2;
    aPnt2 = gp_Pnt(B1+Rr*cos(thetaf),Rr+Rr*sin(thetaf),0.0);
    Handle(Geom_TrimmedCurve) aArcOfCircleRr = GC_MakeArcOfCircle(aCirc,aPnt1,-M_PI/2+theta,true);
    TopoDS_Edge rrEdge = BRepBuilderAPI_MakeEdge(aArcOfCircleRr);
    
    aPnt1 = aPnt2;
    aPnt2 = gp_Pnt(BW,BH,0.0);
    Handle(Geom_TrimmedCurve) aSegment3    = GC_MakeSegment(aPnt1, aPnt2);
    TopoDS_Edge rEdge = BRepBuilderAPI_MakeEdge(aSegment3);
        
    BRepBuilderAPI_MakeWire mkWire;
    mkWire.Add(r2Edge);
    if (HZ > 1.0e-6)
        mkWire.Add(lEdge);
    mkWire.Add(r1Edge);
    mkWire.Add(eEdge);
    mkWire.Add(bEdge);
    mkWire.Add(rrEdge);
    mkWire.Add(rEdge);
    TopoDS_Wire myWireProfile = mkWire.Wire();
    
    TopoDS_Wire finalWire;
    if (Dir == "L") {
        /* Mirror Edge */
        gp_Trsf Mirror;
        anAx2 = gp_Ax2(origin,xAxis,zAxis);
        Mirror.SetMirror(anAx2);
        BRepBuilderAPI_Transform myTransform(myWireProfile,Mirror);
        finalWire = TopoDS::Wire(myTransform.Shape());
        aPnt = gp_Pnt(CE,0.0,0.0);
    }
    else {
        finalWire = myWireProfile;
        aPnt = gp_Pnt(-CE,0.0,0.0);
    }
    TopoDS_Vertex V = BRepBuilderAPI_MakeVertex(aPnt);
    BRep_Builder builder;
    TopoDS_Compound Comp;
    builder.MakeCompound(Comp);
    builder.Add(Comp,finalWire);
    builder.Add(Comp,V);
    BRepTools::Write(Comp,filename.c_str());
}

void step2BRep(std::string const inputfile, std::string const outputfile) {

    /* Load in the step file */
    STEPControl_Reader reader;
    IFSelect_ReturnStatus stat = reader.ReadFile(inputfile.c_str());
    if (stat != IFSelect_RetDone) {
        std::cout << "There was a problem " << stat << std::endl;
        std::cout << inputfile << std::endl;
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
    else {
        TopoDS_Shape shape = reader.Shape(1);
        TopExp_Explorer aVertexExplorer(shape, TopAbs_VERTEX, TopAbs_EDGE);
        if (!aVertexExplorer.More()) {
            /* This is a hack to add a centerline point for now */
            gp_Pnt origin(0.0,0.0,0.0);
            TopoDS_Vertex V = BRepBuilderAPI_MakeVertex(origin);
            
            BRep_Builder builder;
            TopoDS_Compound Comp;
            builder.MakeCompound(Comp);
            builder.Add(Comp,shape);
            builder.Add(Comp,V);
            BRepTools::Write(Comp,outputfile.c_str());
        }
        else {
            BRepTools::Write(shape,outputfile.c_str());
        }
    }
}

void findAndReplaceAll(std::string& str, const std::string& from, const std::string& to) {
    if (from.empty())
        return; // avoid infinite loop

    size_t startPos = 0;
    while ((startPos = str.find(from, startPos)) != std::string::npos) {
        str.replace(startPos, from.length(), to);
        startPos += to.length(); // advance to avoid replacing inside the new text
    }
}

void convertProfiles2BRep(std::string const sourceFolder) {
    /* Load profiles */
    std::vector<profileData> profiles;
    loadProfileData(sourceFolder +"/Inputs/profiles.txt", profiles);
    const int nProfiles = int(profiles.size());
    
    std::filesystem::create_directory(sourceFolder +"/Results");
    std::filesystem::create_directory(sourceFolder +"/Results/Profiles");
    
    for(int row=0;row<nProfiles;++row) {
        
        std::string outputfile = sourceFolder +"/" +profiles[row].filename;
        fs::path po(outputfile);
        std::string dir = po.parent_path().string();  // "Results/Profiles/C01"
        std::filesystem::create_directory(dir);
        
        std::string inputfile = outputfile;
        findAndReplaceAll(inputfile, "Results","Inputs");
        findAndReplaceAll(inputfile, ".brep",".stp");

        
        step2BRep(inputfile,outputfile);
        
    }
}

void wireToPoints(TopoDS_Wire aWire, int nPts, TColgp_Array1OfPnt& profilePoints) {
    /* Playing around with this */
    
    TopoDS_Edge anEdge;
    TopExp_Explorer edgeExplorer(aWire, TopAbs_EDGE);
    
    double length[10];
    int edgeCounter = 0;
    for(;edgeExplorer.More();edgeExplorer.Next()){
        anEdge = TopoDS::Edge(edgeExplorer.Current());
        Standard_Real begin, end;
        Handle(Geom_Curve) aCurve = BRep_Tool::Curve(anEdge, begin, end);
        gp_Pnt pnt;
        gp_Vec vec;
        length[edgeCounter] = 0.0;
        double du = (begin-end)/(nPts-1.);
        for (int i=0;i<nPts;++i) {
            Standard_Real U1 = i*du +begin;
            Standard_Real U2 = U1 +du;
            aCurve->D1(U1,pnt,vec);
            length[edgeCounter] += 0.5*vec.Magnitude()*du;
            aCurve->D1(U2,pnt,vec);
            length[edgeCounter] += 0.5*vec.Magnitude()*du;
        }
        ++edgeCounter;
    }
    
    double totalLength = 0.0;
    for(int i=0;i<edgeCounter;++i) {
        totalLength += length[i];
    }
    
    int ptsArray[10];
    int remPts = nPts-1;
    for(int i=0;i<edgeCounter;++i) {
        ptsArray[i] = round(length[i]/totalLength*remPts);
        totalLength -= length[i];
        remPts -= ptsArray[i];
    }
    
    edgeCounter = 0;
    edgeExplorer.ReInit();
    Standard_Real begin, end = 1.0;
    gp_Pnt pnt;
    Handle(Geom_Curve) aCurve;
    int pntCount = 1;
    for(;edgeExplorer.More();edgeExplorer.Next()){
        anEdge = TopoDS::Edge(edgeExplorer.Current());
        aCurve = BRep_Tool::Curve(anEdge, begin, end);
        double du = (end-begin)/ptsArray[edgeCounter];
        for (int i=0;i<ptsArray[edgeCounter];++i) {
            Standard_Real U1 = i*du +begin;
            aCurve->D0(U1,pnt);
            profilePoints.SetValue(pntCount++,pnt);
//            std::cout << U1 << ' ' << pnt.X() << ' ' << pnt.Y() << ' ' << pnt.Z() << std::endl;
        }
        ++edgeCounter;
    }
    /* Add Last Point */
    aCurve->D0(end,pnt);
    profilePoints.SetValue(pntCount++,pnt);
//    std::cout << end << ' ' << pnt.X() << ' ' << pnt.Y() << ' ' << pnt.Z() << std::endl;

//    BRepAdaptor_CompCurve myComposite(aWire);
//    Standard_Real U1 = myComposite.FirstParameter();
//    Standard_Real U2 = myComposite.LastParameter();
//    gp_Pnt aPnt;
//
//    double du = (U2-U1)/(nPts-1);
//    for (int i = 0; i < 100; ++i) {
//        Standard_Real U = U1 +du*i;
//        gp_Vec aVec;
//        myComposite.D1(U,aPnt,aVec);
//        std::cout << "comp " << U << ' ' << aPnt.X() << ' ' << aPnt.Y() << ' ' << aPnt.Z() << std::endl;
//    }
/*    You can use the class BRepAdaptor_CompCurve if the interface of Adaptor3d_Curve is sufficient for you. If you need the real Geom_Curve then you have to approximate CompCurve into bspline using available approximation tool (e.g. AdvApprox_ApproxAFunction, see GeomLib::BuildCurve3d). */

}



Handle(Geom_BSplineCurve) makeCubicSplineWithSlopes(const gp_Pnt& p0, const gp_Vec& v0, const gp_Pnt& p1, const gp_Vec& v1) {
    // We’ll approximate a cubic Hermite segment using a 4-point BSpline
    // Equivalent to a cubic with clamped ends and specified derivatives.

    // Compute chord length between endpoints
    gp_Vec dxy(p0, p1);
    double L = dxy.Magnitude();

    // Defensive check
    if (L < 1e-12) {
        std::cerr << "Error: endpoints are coincident.\n";
        return nullptr;
    }

    // Choose handle length proportional to chord length
    // For a cubic Bézier with clamped tangents, the typical choice is L/3
    double handle = L / 3.0;
    
    // Create intermediate control points based on the Hermite formulation:
    // P0 = p0
    // P1 = p0 + v0 / 3
    // P2 = p1 - v1 / 3
    // P3 = p1
    gp_Pnt P0 = p0;
    gp_Pnt P1 = gp_Pnt(p0.XYZ() + handle * v0.XYZ()/v0.Magnitude());
    gp_Pnt P2 = gp_Pnt(p1.XYZ() - handle * v1.XYZ()/v1.Magnitude());
    gp_Pnt P3 = p1;

    // Define control points
    TColgp_Array1OfPnt poles(1, 4);
    poles(1) = P0;
    poles(2) = P1;
    poles(3) = P2;
    poles(4) = P3;

    // Uniform knot vector for a cubic Bézier (two knots with multiplicity 4)
    TColStd_Array1OfReal knots(1, 2);
    knots(1) = 0.0;
    knots(2) = 1.0;

    TColStd_Array1OfInteger mults(1, 2);
    mults(1) = 4;  // clamped start
    mults(2) = 4;  // clamped end

    const Standard_Integer degree = 3;

    Handle(Geom_BSplineCurve) spline =
        new Geom_BSplineCurve(poles, knots, mults, degree);

    return spline;
}

Handle(Geom_BSplineCurve) makeCubicSpline(const gp_Pnt& p0, const gp_Pnt& p1, const gp_Pnt& p2, const gp_Pnt& p3) {
    // We’ll approximate a cubic Hermite segment using a 4-point BSpline
    // Equivalent to a cubic with clamped ends and specified derivatives.

    // Define control points
    TColgp_Array1OfPnt poles(1, 4);
    poles(1) = p0;
    poles(2) = p1;
    poles(3) = p2;
    poles(4) = p3;

    // Uniform knot vector for a cubic Bézier (two knots with multiplicity 4)
    TColStd_Array1OfReal knots(1, 2);
    knots(1) = 0.0;
    knots(2) = 1.0;

    TColStd_Array1OfInteger mults(1, 2);
    mults(1) = 4;  // clamped start
    mults(2) = 4;  // clamped end

    const Standard_Integer degree = 3;

    Handle(Geom_BSplineCurve) spline =
        new Geom_BSplineCurve(poles, knots, mults, degree);

    return spline;
}

void makeProfilesCortina(std::string const sourceFile, std::string const destinationFolder) {
    
#ifdef TEST_CORTINA
    /* Curve 1 Parameters */
    double Rr = 212.5;
    double BW = 842+Rr;
    double BH = 496+Rr;
    /* % 0,0 is left of arc
     % 1274.1937,-1433.959 is near bottom
     a = 1274.1937;
     b = -1433.959;
     % xc == sqrt((xc-a)^2 +b^2))
     % xc^2 = xc^2-2*xc*a +a^2+b^2;
     xc = (a^2+b^2)/(2*a)

     xy = [0,0; 0,-709.3066; 538.365,-1315.175; 1274.1937,-1433.959; xc,0.0];
     scatter(xy(:,1),xy(:,2))
     axis equal */
    double R1 = 1444;
    double MY = 1274.1937;
    double CE = 15535-15129;
    double UEB = -88.0;
    double DXDY = -6.0/(836.0-739.0);
    double KH = MY+836.0;
    makeCortinaProfile("CortinaTest.brep", "L", DXDY, UEB, MY, R1, KH, BW, BH, Rr, CE);
    return 0;
#endif // TEST_CORTINA
    
    std::ifstream input;
    std::string line;
    std::string word;
    
    std::filesystem::create_directory(destinationFolder +"/Results");
    std::filesystem::create_directory(destinationFolder +"/Results/Profiles");
    
    input.open(sourceFile);
    
    std::getline(input, line);
    /* Read headers */
    std::map<std::string,int> hcols;
    std::stringstream hstream(line);
    int cols = 0;
    while (hstream >> word) {
        hcols[word] = cols++;
    }
    /* Type    Dir    Sx    X1    X2    Y2    X3    Y3    BT    UEB    KH    B1    BH    Rr    CE */
    
    /* Input line and make profile*/
    std::vector<profileData> profiles;
    profileData myProfile;
    std::string Type,Dir,num;
    char com;
    int prevProfileType = 0, profileType;
    
    double Sx, X1, X2, Y2, X3, Y3, BT, UEB, KH, B1, BH, Rr, CE;
    int rows = 0;
    while (std::getline(input, line)) {
        std::istringstream dstream(line);
        std::getline(dstream,Type,',');
        std::getline(dstream,Dir,',');
        dstream >> Sx >> com >> X1 >> com >> X2 >> com >> Y2 >> com >> X3 >> com >> Y3 >> com >> BT >> com >> UEB >> com >> KH >> com >> B1 >> com >> BH >> com >> Rr >> com >> CE;
        
        /* Decide what kind of profile it is */
        std::ostringstream nstr;
        nstr << Sx << '_' << Type << '_' << Dir;
        myProfile.filename = "Results/Profiles/" +nstr.str() +".brep";
        myProfile.sx = Sx;
        myProfile.BRK = 0;
        std::string filename = destinationFolder +"/" +myProfile.filename;
        
        if (Dir == "T") {
            double BW = B1+Rr;
            double B = BT+Rr; /* Going to assume radii are the same rather than using spline definition */
            makeStraightProfile(filename, Dir, B, BW, KH, BH, Rr, Rr, CE);
            profileType = 0;
        }
        else {
            makeCortinaProfile(filename, Dir, X1, X2, Y2, X3, Y3, BT, UEB, KH, B1, BH, Rr, CE);
            profileType = 1;
        }
//        if (profileType != prevProfileType) {
//            myProfile.BRK = 1;
//        }
        profiles.push_back(myProfile);
        prevProfileType = profileType;
        ++rows;
    }
    std::ofstream of;
    of.open(destinationFolder +"/Results/profiles.txt");
    /* Output profile file*/
    for (int i=0; i < rows; ++i) {
        of << profiles[i].filename << ' ' << profiles[i].sx << ' ' << profiles[i].BRK << std::endl;
    }
    of.close();
}

void makeCortinaProfile(std::string const filename, std::string const Dir, double X1, double X2, double Y2, double X3, double Y3, double BT, double UEB, double KH, double B1, double BH, double Rr, double CE) {
    const gp_Dir xAxis(1.0,0.0,0.0);
    const gp_Dir yAxis(0.0,1.0,0.0);
    const gp_Dir zAxis(0.0,0.0,1.0);
    const gp_Pnt origin(0.0,0.0,0.0);

    gp_Pnt aPnt, aPnt1, aPnt2;
    gp_Ax2 anAx2;
    
    aPnt1 = gp_Pnt(-X3-UEB-BT,KH,0.0);
    aPnt2 = gp_Pnt(-X3-BT,Y3,0.0);
    double theta = atan((X3-X2)/(Y3-Y2));
    double UEH = KH-Y3;
    double R2 = abs(0.5*(UEB*UEB+UEH*UEH)/(UEB*cos(theta)+UEH*sin(theta)));
    aPnt = gp_Pnt(R2*cos(theta) +aPnt2.X(),R2*sin(theta)+aPnt2.Y(),0.0);
    anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
    gp_Circ aCirc(anAx2,R2);
    Handle(Geom_TrimmedCurve) aArcOfCircleR1 = GC_MakeArcOfCircle(aCirc,aPnt1,aPnt2,Standard_True);
    TopoDS_Edge r1Edge = BRepBuilderAPI_MakeEdge(aArcOfCircleR1);
    
    aPnt1 = gp_Pnt(-BT,0.0,0.0);
    Handle(Geom_BSplineCurve) spline = makeCubicSpline(aPnt2,gp_Pnt(-X2-BT,Y2,0.0),gp_Pnt(-X1-BT,0.0,0.0),aPnt1);
    TopoDS_Edge splineEdge = BRepBuilderAPI_MakeEdge(spline);

    aPnt2 = gp_Pnt(B1,0.0,0.0);
    Handle(Geom_TrimmedCurve) aSegment2    = GC_MakeSegment(aPnt1, aPnt2);
    TopoDS_Edge bEdge = BRepBuilderAPI_MakeEdge(aSegment2);
    
    aPnt = gp_Pnt(B1,Rr,0.0);
    anAx2 = gp_Ax2(aPnt, zAxis, xAxis);
    aCirc = gp_Circ(anAx2,Rr);
    Handle(Geom_TrimmedCurve) aArcOfCircle2 = GC_MakeArcOfCircle(aCirc,3.*M_PI/2.,2.*M_PI,true);
    TopoDS_Edge rrEdge = BRepBuilderAPI_MakeEdge(aArcOfCircle2);
    
    double BW = B1+Rr;
    aPnt1 = gp_Pnt(BW,Rr,0.0);
    aPnt2 = gp_Pnt(BW,BH,0.0);
    Handle(Geom_TrimmedCurve) aSegment3    = GC_MakeSegment(aPnt1, aPnt2);
    TopoDS_Edge rEdge = BRepBuilderAPI_MakeEdge(aSegment3);
        
    BRepBuilderAPI_MakeWire mkWire;
    mkWire.Add(r1Edge);
    mkWire.Add(splineEdge);
    mkWire.Add(bEdge);
    mkWire.Add(rrEdge);
    mkWire.Add(rEdge);
    TopoDS_Wire myWireProfile = mkWire.Wire();
    
    TopoDS_Wire finalWire;
    if (Dir == "L") {
        /* Mirror Edge */
        gp_Trsf Mirror;
        anAx2 = gp_Ax2(origin,xAxis,zAxis);
        Mirror.SetMirror(anAx2);
        BRepBuilderAPI_Transform myTransform(myWireProfile,Mirror);
        finalWire = TopoDS::Wire(myTransform.Shape());
        aPnt = gp_Pnt(CE,0.0,0.0);
    }
    else {
        finalWire = myWireProfile;
        aPnt = gp_Pnt(-CE,0.0,0.0);
    }
    TopoDS_Vertex V = BRepBuilderAPI_MakeVertex(aPnt);
    BRep_Builder builder;
    TopoDS_Compound Comp;
    builder.MakeCompound(Comp);
    builder.Add(Comp,finalWire);
    builder.Add(Comp,V);
    BRepTools::Write(Comp,filename.c_str());
}
