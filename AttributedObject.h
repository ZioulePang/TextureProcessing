///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  AttributedObject.h
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//	Variant on TexturedObject that stores explicit RGB
//	values for each vertex
//  
///////////////////////////////////////////////////

// include guard for AttributedObject
#ifndef _ATTRIBUTED_OBJECT_H
#define _ATTRIBUTED_OBJECT_H

// include the C++ standard libraries we need for the header
#include <vector>
#include <iostream>
#include <map>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

// include the unit with Cartesian 3-vectors
#include "Cartesian3.h"
// the render parameters
#include "RenderParameters.h"

// define a macro for "not used" flag
#define NO_SUCH_ELEMENT -1

// use macros for the "previous" and "next" IDs
#define PREVIOUS_EDGE(x) ((x) % 3) ? ((x) - 1) : ((x) + 2)
#define NEXT_EDGE(x) (((x) % 3) == 2) ? ((x) - 2) : ((x) + 1)

//Set first directed edge structure
struct Vertex
{
    long FirstEdgeIndex;
    Cartesian3 vertices;

    bool isBoundary;
};

struct Edge
{
    int firstIndex;

    int secondeIndex;

    int hasEdge = -1;

};

class AttributedObject
    { // class AttributedObject
    public:
    // vector of vertices
    std::vector<Cartesian3> vertices;

    // vector of colours stored as cartesian triples in float
    std::vector<Cartesian3> colours;
    
    // vector of normals
    std::vector<Cartesian3> normals;
    
    // vector of texture coordinates (stored as triple to simplify code)
    std::vector<Cartesian3> textureCoords;

    // vector for the "first" directed edge for each vertex
    std::vector<long> firstDirectedEdge;

    // vector of faces - doubles as the "to" array for edges
    std::vector<unsigned int> faceVertices;

    // vector of the "other half" of the directed edges
    std::vector<long> otherHalf;

    //TODO:: vector of the directed vertices
    std::vector<Vertex> TODO_vertices;

    //TODO:: vector of the directed edges
    std::vector<Edge> TODO_edges;

    std::vector<Edge> boundary_edges;

    // centre of gravity - computed after reading
    Cartesian3 centreOfGravity;

    // size of object - i.e. radius of circumscribing sphere centred at centre of gravity
    float objectSize;

    // constructor will initialise to safe values
    AttributedObject();

   //use floater's method
    void Floaters();

    //all of the vertices except boundary points
    int interior_vertices;

    // read routine returns true on success, failure otherwise
    bool ReadObjectStream(std::istream &geometryStream);

    // write routine
    void WriteObjectStream(std::ostream &geometryStream);

    void DirectedAttribute();

    //Caculate boundary length
    float Distance(Cartesian3 a, Cartesian3 b);

    float Angle(Cartesian3 a,Cartesian3 b,Cartesian3 c);

    float Area(Cartesian3 a,Cartesian3 b,Cartesian3 c);
    bool isInTriangle(Cartesian3 a, Cartesian3 b, Cartesian3 c, Cartesian3 p);
    bool isBoundary(int index);
    void TextureMapping();
    void control(RenderParameters *renderParameters);
    // routine to render
    void Render(RenderParameters *renderParameters);
    }; // class AttributedObject

// end of include guard for AttributedObject
#endif
