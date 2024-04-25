//
//  timeDeriv.h
//  occ_track
//
//  Created by Brian Helenbrook on 10/27/13.
//  Copyright (c) 2013 Brian Helenbrook. All rights reserved.
//

#ifndef __occ_track__timeDeriv__
#define __occ_track__timeDeriv__

#include <blitz/array.h>
#include <Geom_BSplineSurface.hxx>
#include <iostream>

void timeIntegrate(Handle(Geom_BSplineSurface) trackSurface);
void timeDeriv(blitz::array<double,1> y, blitz::Array<double,1> dydt, Handle(Geom_BSplineSurface) trackSurface);



#endif /* defined(__occ_track__timeDeriv__) */
