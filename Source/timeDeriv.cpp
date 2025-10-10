//
//  timeDeriv.cpp
//  occ_track
//
//  Created by Brian Helenbrook on 10/27/13.
//  Copyright (c) 2013 Brian Helenbrook. All rights reserved.
//

#include "occ_track.h"

#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

void timeDeriv(std::array<double,4>& y, std::array<double,4>& dydt, Handle(Geom_Surface) trackSurface) {
    
    //  Calculate time derivative at a point.
    //  y(0) = u coordinate
    //  y(1) = v coordinate
    //  y(2) = du/dt coordinate velocity
    //  y(3) = dv/dt coordinate velocity
    
    const double rho_CdA_2m = 1.2*(50/(0.5*1.2*35*35))/(2.*100.0);
    const double g = 9.81;
    const double mu = 0.015;
    
    Standard_Real U, V, dUdt, dVdt;
    
    U = y[0];
    V = y[1];
    dUdt = y[2];
    dVdt = y[3];
    
    gp_Pnt loc;
    gp_Vec coU, coV, D2U, D2V, D2UV;
    
    // Calculate tangents at point & curvatures
    trackSurface->D2(U,V,loc,coU,coV,D2U,D2V,D2UV);
    
    // Calculate normal to track
    gp_Vec normal = coU^coV;
    normal /= normal.Magnitude();
    
    // Now calculate contra-variant vectors;
    gp_Vec conU, conV;
    
    conU = normal^coU;
    conU /= conU*coV;
    conV = coV^normal;
    conV /= conV*coU;
    
    
    // Calculate normal force
    double trackForce;
    gp_Vec accelCurv = (D2U*dUdt*dUdt +2*D2UV*dUdt*dVdt +D2V*dVdt*dVdt);
#ifdef UNITY
    trackForce = accelCurv*normal  +g*normal.Y();
#else
    trackForce = accelCurv*normal  +g*normal.Z();
#endif
    
    // Calculate tangent forces
    // these act in a direction opposite to the velocity
    
    gp_Vec physVel = coU*dUdt +coV*dVdt;
    Standard_Real speed = physVel.Magnitude();
    
    gp_Vec drag = -rho_CdA_2m*physVel*speed;
    gp_Vec friction = -trackForce*mu*physVel/speed;
#ifdef UNITY
    gp_Vec tangentForces = drag +friction -g*gp_Vec(0,1,0);
#else
    gp_Vec tangentForces = drag +friction -g*gp_Vec(0,0,1);
#endif
    
    double d2Udt2 = (tangentForces -accelCurv)*conV;
    double d2Vdt2 = (tangentForces -accelCurv)*conU;
    
    dydt[0] = dUdt;
    dydt[1] = 0*dVdt;
    dydt[2] = d2Udt2;
    dydt[3] = 0*d2Vdt2;
    return;
}
	

void timeIntegrate(std::string filename, Handle(Geom_Surface const) const trackSurface) {
	
	
	Standard_Real U1,U2,V1,V2;  // Get bounds of track
	trackSurface->Bounds(U1,U2,V1,V2);
	
	std::array<double,4> y, y0, k1, k2, k3, k4;

	y[0] = U1;	// Beginning of track
	y[1] = 0.5*(V1+V2); // in the middle
	y[2] = 27*1000/3600;  // If U is arc-length then this will be 27 km/h otherwise have to adjust
	y[3] = 0.0; // No side ways velocity to start
	
//	std::cout << U1 << ' ' << U2 << ' ' << V1 << ' ' << V2 << std::endl;
	
	double dt = 1.0/30;  // Should move about 1/10 meter per time step
	double time = 0.0;
	
	std::ofstream path_file;
	path_file.open(filename);
	
	while (y[0] <= U2 && y[2] > 0.0) { // && y(1) >= V1 && y(1) <= V2) {
		y0 = y;
		timeDeriv(y,k1,trackSurface);
        for (int n=0;n<4;++n)
            y[n] = y0[n]+0.5*dt*k1[n];
        
		timeDeriv(y,k2,trackSurface);
        for (int n=0;n<4;++n)
            y[n] = y0[n]+0.5*dt*k2[n];
        
		timeDeriv(y,k3,trackSurface);
        for (int n=0;n<4;++n)
            y[n] = y0[n]+dt*k3[n];
        
		timeDeriv(y,k4,trackSurface);
        
        for (int n=0;n<4;++n)
            y[n] = y0[n] +dt/6.*(k1[n] +2.*k2[n] +2.*k3[n] +k4[n]);
		time += dt;
		
		gp_Pnt loc;
		gp_Vec coU, coV;
		trackSurface->D1(y[0],y[1],loc,coU,coV);
		
		gp_Vec normal = coU^coV;
		normal /= normal.Magnitude();
		
		gp_Vec physVel = coU*y[2] +coV*y[3];
		
#ifdef UNITY
        path_file << time << ' ' << loc.X() << ' ' << loc.Z() << ' ' << -loc.Y() << ' ' << normal.X() << ' ' << normal.Z() << ' ' << -normal.Y() << ' ' << physVel.X() << ' ' << physVel.Z() << ' ' << -physVel.Y() << std::endl;
#else
        path_file << time << ' ' << loc.X() << ' ' << loc.Y() << ' ' << loc.Z() << ' ' << normal.X() << ' ' << normal.Y() << ' ' << normal.Z() << ' ' << physVel.X() << ' ' << physVel.Y() << ' ' << physVel.Z() << std::endl;
#endif
	}
	path_file.close();
}
