//
//  ExortSTL.hpp
//  occ_track
//
//  Created by Poorna Raavi on 9/10/25.
//

// ExportSTL.hpp
#pragma once
#include <TopoDS_Shape.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <StlAPI.hxx>
#include <filesystem>
#include <string>

// 1) Shape -> STL
inline bool ExportBREPtoSTL(const TopoDS_Shape& shape,
                            const std::string& outPath,
                            double linearDefl = 0.05,
                            double angularDefl = 0.35,
                            bool ascii = false)
{
    BRepMesh_IncrementalMesh mesher(shape, linearDefl, /*isRelative*/false,
                                    angularDefl, /*parallel*/true);
    std::filesystem::create_directories(std::filesystem::path(outPath).parent_path());
    return StlAPI::Write(shape, outPath.c_str(), ascii);
}

// 2) BREP file -> STL file (convenience)
inline bool ExportBREPFileToSTL(const std::string& brepPath,
                                const std::string& stlPath,
                                double linearDefl = 0.05,
                                double angularDefl = 0.35,
                                bool ascii = false)
{
    TopoDS_Shape s;
    BRep_Builder b;
    if (!BRepTools::Read(s, brepPath.c_str(), b)) return false;
    return ExportBREPtoSTL(s, stlPath, linearDefl, angularDefl, ascii);
}
