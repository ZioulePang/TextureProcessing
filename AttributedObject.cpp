///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  AttributedObject.cpp
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

// include the header file
#include "AttributedObject.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include <map>
#include <Eigen/Dense>

// include the Cartesian 3- vector class
#include "Cartesian3.h"
#include <set>

#define MAXIMUM_LINE_LENGTH 1024
#define REMAP_TO_UNIT_INTERVAL(x) (0.5 + (0.5*(x)))
#define REMAP_FROM_UNIT_INTERVAL(x) (-1.0 + (2.0*(x)))
#define M_PI 3.14159265358979323846
#define N_ITERATIONS 100000

// constructor will initialise to safe values
AttributedObject::AttributedObject()
    : centreOfGravity(0.0,0.0,0.0)
    { // AttributedObject()
    // force arrays to size 0
    vertices.resize(0);
    colours.resize(0);
    normals.resize(0);
    textureCoords.resize(0);
    firstDirectedEdge.resize(0);
    faceVertices.resize(0);
	otherHalf.resize(0);
    } // AttributedObject()

// read routine returns true on success, failure otherwise
bool AttributedObject::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()
    
    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];
    
    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
        // character to read
        char firstChar = geometryStream.get();
        
//         std::cout << "Read: " << firstChar << std::endl;
        
        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the character we read
        switch (firstChar)
            { // switch on first character
            case '#':       // comment line
                // read and discard the line
                geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
                break;
                
            case 'v':       // vertex data of some type
                { // some sort of vertex data
                // retrieve another character
                char secondChar = geometryStream.get();
                
                // bail if we ran out of file
                if (geometryStream.eof())
                    break;

                // now use the second character to choose branch
                switch (secondChar)
                    { // switch on second character
                    case ' ':       // space - indicates a vertex
                        { // vertex read
                        Cartesian3 vertex;
                        geometryStream >> vertex;
                        vertices.push_back(vertex);
//                         std::cout << "Vertex " << vertex << std::endl;
                        break;
                        } // vertex read
                    case 'c':       // c indicates colour
                        { // normal read
                        Cartesian3 colour;
                        geometryStream >> colour;
                        colours.push_back(colour);
//                         std::cout << "Colour " << colour << std::endl;
                        break;
                        } // normal read
                    case 'n':       // n indicates normal vector
                        { // normal read
                        Cartesian3 normal;
                        geometryStream >> normal;
                        normals.push_back(normal);
//                         std::cout << "Normal " << normal << std::endl;
                        break;
                        } // normal read
                    case 't':       // t indicates texture coords
                        { // tex coord
                        Cartesian3 texCoord;
                        geometryStream >> texCoord;
                        textureCoords.push_back(texCoord);
//                         std::cout << "Tex Coords " << texCoord << std::endl;
                        break;                  
                        } // tex coord
                    default:
                        break;
                    } // switch on second character 
                break;
                } // some sort of vertex data
                
            case 'f':       // face data
                { // face
				// make a hard assumption that we have a single triangle per line
                unsigned int vertexID;
                
                // read in three vertices
				for (unsigned int vertex = 0; vertex < 3; vertex++)
					{ // per vertex
					// read a vertex ID
					geometryStream >> vertexID;

					// subtract one and store them (OBJ uses 1-based numbering)
					faceVertices.push_back(vertexID-1);
					} // per vertex
				break;
                } // face
                
            // default processing: do nothing
            default:
                break;

            } // switch on first character

        } // not eof

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (vertices.size() != 0)
        { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];
        
        // and divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();         
            
            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;
                
            } // per vertex
        } // non-empty vertex set

// 	std::cout << "Centre of Gravity: " << centreOfGravity << std::endl;
// 	std::cout << "Object Size:       " << objectSize << std::endl;

    // return a success code
    return true;
	} // ReadObjectStream()

// write routine
void AttributedObject::WriteObjectStream(std::ostream &geometryStream)
    { // WriteObjectStream()
    geometryStream << "# " << (faceVertices.size()/3) << " triangles" << std::endl;
    geometryStream << std::endl;

    // output the vertex coordinates
    geometryStream << "# " << vertices.size() << " vertices" << std::endl;
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "v  " << std::fixed << vertices[vertex] << std::endl;

    // output the vertex colours
    geometryStream << "# " << colours.size() << " vertex colours" << std::endl;
    for (unsigned int vertex = 0; vertex < colours.size(); vertex++)
        geometryStream << "vc " << std::fixed << colours[vertex] << std::endl;

    // output the vertex normals
    geometryStream << "# " << normals.size() << " vertex normals" << std::endl;
    for (unsigned int vertex = 0; vertex < normals.size(); vertex++)
        geometryStream << "vn " << std::fixed << normals[vertex] << std::endl;

    // output the vertex coords
    geometryStream << "# " << textureCoords.size() << " vertex tex coords" << std::endl;
    for (unsigned int vertex = 0; vertex < textureCoords.size(); vertex++)
        geometryStream << "vt " << std::fixed << textureCoords[vertex] << std::endl;

    // and the faces
    for (unsigned int face = 0; face < faceVertices.size(); face+=3)
        { // per face
        geometryStream << "f";
        
        // loop through # of vertices
        for (unsigned int vertex = 0; vertex < 3; vertex++)
			{ // per vertex
            geometryStream << " ";
            geometryStream << faceVertices[face+vertex] + 1;
			} // per vertex
		// end the line
        geometryStream << std::endl;
        } // per face
    
    } // WriteObjectStream()

//caculate Attributes
void AttributedObject::DirectedAttribute()
{
    int i,j = 0;

    //Store vertices
    for (i = 0;i<vertices.size();i++)
    {
        Vertex temp_vertex;
        temp_vertex.vertices = vertices[i];
        temp_vertex.FirstEdgeIndex = 0;
        temp_vertex.isBoundary =false;
        TODO_vertices.push_back(temp_vertex);
    }

    //Store directed edges
    for (i = 0;i < faceVertices.size() / 3;i++)
    {
        Edge temp_edge1,temp_edge2,temp_edge3;
        temp_edge1.firstIndex = faceVertices[3 * i + 2];
        temp_edge1.secondeIndex = faceVertices[3 * i];

        temp_edge2.firstIndex = faceVertices[3 * i];
        temp_edge2.secondeIndex = faceVertices[3 * i + 1];

        temp_edge3.firstIndex = faceVertices[3 * i + 1];
        temp_edge3.secondeIndex = faceVertices[3 * i + 2];

        TODO_edges.push_back(temp_edge1);
        TODO_edges.push_back(temp_edge2);
        TODO_edges.push_back(temp_edge3);
    }

    //Find first directed edge
    for (i = 0;i<vertices.size();i++)
    {
        int time = 0;
        for (j = 0;j<TODO_edges.size();j++)
        {
            if(TODO_edges[j].firstIndex == i && time == 0)
            {
                TODO_vertices[i].FirstEdgeIndex = TODO_edges[j].secondeIndex;
                firstDirectedEdge.push_back(TODO_edges[j].secondeIndex);
                time++;
            }
        }
    }
    //insert otherhalf
    for (i = 0;i < TODO_edges.size();i++)
    {
        for (j = 0;j < TODO_edges.size();j++)
        {
            if((TODO_edges[j].secondeIndex == TODO_edges[i].firstIndex && TODO_edges[j].firstIndex == TODO_edges[i].secondeIndex) && i != j)
            {
                TODO_edges[i].hasEdge = 1;
                otherHalf.push_back(j);
            }
        }
    }

    //find boundary points
    for (i = 0;i<TODO_edges.size();i++)
    {
        if(TODO_edges[i].hasEdge == -1)
        {
            TODO_vertices[TODO_edges[i].firstIndex].isBoundary = true;
            TODO_vertices[TODO_edges[i].secondeIndex].isBoundary = true;
            boundary_edges.push_back(TODO_edges[i]);
        }
    }


}

//distinguish wether it is boundary
bool AttributedObject::isBoundary(int index)
{
    if(TODO_edges[index].hasEdge == -1)
    {
        return true;
    }
    return false;
}

float AttributedObject::Distance(Cartesian3 a, Cartesian3 b){
    double x=a.x-b.x;
    double y=a.y-b.y;
    double z=a.z-b.z;

    return sqrt(x*x+y*y+z*z);
}

void AttributedObject::Floaters()
{
    int i,j;
    int* total_Nei=new int[vertices.size()];
    double** weight=new double*[vertices.size()];
    for (i=0;i<vertices.size();i++)
    {
        weight[i]= new double[vertices.size()];
        for (j=0; j<vertices.size();j++)
           weight[i][j]=0;

    }
    //caculate how many boundary points we have
    std::vector<int> boundary_vertices;
    std::vector<int> interiorPoints;
    std::vector<int> exteriorPoints;
    interior_vertices = vertices.size() - boundary_vertices.size();
    for (int i = 0;i<TODO_vertices.size();i++)
    {
        if(TODO_vertices[i].isBoundary)
        {
            boundary_vertices.push_back(i);
        }else
        {
            interiorPoints.push_back(i);
        }
    }

    //Sort boundary edges
    for (i=0;i<boundary_edges.size();i++)
    {
        int count = i + 1;
        for (j=i+1;j<boundary_edges.size();j++) {
            if(boundary_edges[i].secondeIndex == boundary_edges[j].firstIndex && i != j)
            {
                Edge temp = boundary_edges[j];
                boundary_edges[j] = boundary_edges[count];
                boundary_edges[count] = temp;
                count ++;
                break;
            }

        }

    }

    std::vector<Cartesian3> temp_Vertex = vertices;
    std::vector<int> sorted_points;
    for (i =0;i<boundary_edges.size();i++)
    {
        exteriorPoints.push_back(boundary_edges[i].firstIndex);
    }
    for (i = 0;i<vertices.size();i++)
    {
        if(i<exteriorPoints.size())
        {
            sorted_points.push_back(exteriorPoints[i]);
        }
        else {
            sorted_points.push_back(interiorPoints[i - exteriorPoints.size()]);
        }
    }

    //Caculate boundary perimeter
    float *boundary_length_arr = new float[boundary_vertices.size()];
    float sum=0;
    for (i = 0;i < boundary_vertices.size();i++)
    {
        boundary_length_arr[i]=Distance(vertices[boundary_edges[i].firstIndex],vertices[boundary_edges[i].secondeIndex]);
        sum+=boundary_length_arr[i];
    }

    //parameterize to square
    Cartesian3* parameterised_arr = new Cartesian3[vertices.size()];
    float temp_distance=0;
    for (i=0;i<boundary_vertices.size();i++)
    {
        if (temp_distance<2)
        {
           parameterised_arr[i].x = 1.0f;
           parameterised_arr[i].y = temp_distance - 1.0f;
        }
        else if (temp_distance<4)
        {
            parameterised_arr[i].x = 3.0f - temp_distance;
            parameterised_arr[i].y = 1.0f;
        }
        else if (temp_distance<6)
        {
            parameterised_arr[i].x = -1.0f;
            parameterised_arr[i].y = 5.0f - temp_distance;
        }
        else
        {
            parameterised_arr[i].x = temp_distance - 7.0f;
            parameterised_arr[i].y = -1.0f;
        }
        temp_distance+= boundary_length_arr[i]*8/sum;
    }

    std::vector<float> angles;
    for (i = boundary_vertices.size();i < vertices.size();i++)
    {
        float sumAngle = 0;
        int counter = 0,first_nei = 0;
        for (j = 0;j<faceVertices.size();j++)
        {
            if(faceVertices[j] == sorted_points[i])
            {
                if(j%3 == 0)
                {
                    total_Nei[faceVertices[j + 1]] = faceVertices[j + 2];
                    counter++;
                    if(counter == 1) first_nei = faceVertices[j + 1];
                }
                else if (j%3 == 1)
                {
                   total_Nei[faceVertices[j + 1]] = faceVertices[j - 1];
                   counter++;
                   if(counter == 1) first_nei = faceVertices[j + 1];
                }
                else if (j%3 == 2)
                {
                    total_Nei[faceVertices[j - 2]] = faceVertices[j - 1];
                    counter++;
                    if(counter == 1) first_nei = faceVertices[j - 2];
                }

            }
        }

        int* nei = new int[counter];
        nei[0] = first_nei;
        for ( j=1;j<counter;j++)
        {
            nei[j]=total_Nei[nei[j-1]];
        }
        //Caculate interior sum angle



        if(counter!=0)
        {
            for (j=0;j<counter;j++)
            {
               sumAngle += Angle(vertices[nei[j]],vertices[sorted_points[i]],vertices[nei[(j+1)%counter]]);
            }
            Cartesian3* tempPosition= new Cartesian3[counter];
           tempPosition[0].x=(vertices[sorted_points[i]]-vertices[nei[0]]).length();

           float tSumAng=0;
           for (j=1;j<counter;j++)
           {
             float target_Angle=2*M_PI*Angle(vertices[nei[j - 1]],vertices[sorted_points[i]],vertices[nei[j]])/sumAngle;
             tSumAng += target_Angle;
             float length_01 =  (vertices[sorted_points[i]] - vertices[nei[j - 1]]).length();
             float length_02 = (vertices[sorted_points[i]] - vertices[nei[j]]).length();
             float phy = length_02 / length_01;
             tempPosition[j].x=length_02*cos(tSumAng);
             tempPosition[j].y=length_02*sin(tSumAng);
             //tempPosition[j].x = phy*cos(tSumAng)*tempPosition[j - 1].x - phy*sin(tSumAng)*tempPosition[j-1].y;
             //tempPosition[j].y = phy*sin(tSumAng)*tempPosition[j - 1].x + phy*cos(tSumAng)*tempPosition[j-1].y;
           }
           //caculate weight
           int time = 0;
           for (j = 0;j<counter- 2;j++)
           {
              for (int k = j + 1;k<counter - 1;k++)
              {
                for (int m = k + 1;m <counter;m++)
                {
                   if(isInTriangle(tempPosition[j],tempPosition[k] ,tempPosition[m], Cartesian3(0,0,0)))
                   {
                      float total = Area(tempPosition[j],tempPosition[k],tempPosition[m]);
                      float w1 = Area(Cartesian3(0,0,0),tempPosition[k],tempPosition[m]);
                      float w2 = Area(Cartesian3(0,0,0),tempPosition[j],tempPosition[m]);
                      float w3 = Area(Cartesian3(0,0,0),tempPosition[j],tempPosition[k]);
                      if(w1/total<1 && w2/total < 1 && w3/total<1)
                      {
                          weight[sorted_points[i]][nei[j]]+=w1/total;
                          weight[sorted_points[i]][nei[k]]+=w2/total;
                          weight[sorted_points[i]][nei[m]]+=w3/total;
                          time++;
                      }

                   }

                }
              }
           }
           for (j = 0;j<counter;j++)
          {
              if(weight[sorted_points[i]][nei[j]] > 0)
              {
                 weight[sorted_points[i]][nei[j]] /= time;
              }
          }

        }



     }
      Eigen::MatrixXf A = Eigen::MatrixXf::Identity(vertices.size(),vertices.size());
      for (i = exteriorPoints.size();i < vertices.size();i++)
      {
          for (j = 0;j<vertices.size();j++)
          {
             if(i!=j)
             {
                 A(i,j) = -1 * weight[sorted_points[i]][sorted_points[j]];
             }
          }
      }


      //caculate B
      Eigen::MatrixXf B = Eigen::MatrixXf::Zero(vertices.size(),2);
      for (i = 0;i<boundary_vertices.size();i++)
      {
          B(i,0) = parameterised_arr[i].x;
          B(i,1) = parameterised_arr[i].y;
      }

      Eigen::MatrixXf c =  A.inverse() * B;
      for (i = 0;i<vertices.size();i++)
      {
          vertices[sorted_points[i]] = Cartesian3(c(i,0),c(i,1),0);
      }

      std::cout <<textureCoords.size()<< " "<< colours.size()<<std::endl;

}

void AttributedObject::TextureMapping()
{
    int i,j;
    for (i = 0;i<vertices.size();i++)
    {
        colours[i] = (vertices[i] + Cartesian3(1,1,0))/2.0f;
    }


}

float AttributedObject::Angle(Cartesian3 a,Cartesian3 b,Cartesian3 c)
{
    Cartesian3 v1(a.x-b.x,a.y-b.y,a.z-b.z);
    Cartesian3 v2(c.x-b.x,c.y-b.y,c.z-b.z);
    float d1=Distance(a,b);
    float d2=Distance(c,b);
    float cross = v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
    return acos(cross/(d1 * d2));

}

float AttributedObject::Area(Cartesian3 a,Cartesian3 b,Cartesian3 c)
{
    Cartesian3 v1 = a - b;
    Cartesian3 v2 = c - b;
    float cosB = v1.dot(v2)/(v1.length()*v2.length());

    float area = (v1.length()*v2.length()) * sqrt(1 - cosB*cosB) / 2.0f;
    return area;

}

bool AttributedObject::isInTriangle(Cartesian3 a, Cartesian3 b, Cartesian3 c, Cartesian3 p)
{
    Cartesian3 v0 = a - b;
    Cartesian3 v1 = c - a;
    Cartesian3 v2 = b - c;

    Cartesian3 p0 = a - p;
    Cartesian3 p1 = c - p;
    Cartesian3 p2 = b - p;

    Cartesian3 c1 = v0.cross(p0);
    Cartesian3 c2 = v1.cross(p1);
    Cartesian3 c3 = v2.cross(p2);

    if((c1.dot(c2) > 0 && c2.dot(c3) > 0 && c1.dot(c3)>0) ||(c1.dot(c2) < 0 && c2.dot(c3) < 0 && c1.dot(c3)<0))
    {
        return  true;
    }
    else
    {
        return false;
    }
}

void AttributedObject::control(RenderParameters *renderParameters)
{
    if(renderParameters->useTexCoords)
    {
        DirectedAttribute();
        Floaters();
        if(renderParameters->useTexCoords)
        {
            TextureMapping();
        }
    }


}
// routine to render
void AttributedObject::Render(RenderParameters *renderParameters)
    { // Render()
	// make sure that textures are disabled
	glDisable(GL_TEXTURE_2D);

	float scale = renderParameters->zoomScale;
	scale /= objectSize;
	// Scale defaults to the zoom setting
	glTranslatef(-centreOfGravity.x * scale, -centreOfGravity.y * scale, -centreOfGravity.z * scale);
		
	if (renderParameters->useWireframe)
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
	else
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    // start rendering
    glBegin(GL_TRIANGLES);
	
    // loop through the faces: note that they may not be triangles, which complicates life

    for (unsigned int face = 0; face < faceVertices.size(); face+=3)
        { // per face

        // now do a loop over three vertices
        for (unsigned int vertex = 0; vertex < 3; vertex++)
            { // per vertex
            // set colour using vertex ID

            glColor3f
                (

                colours[faceVertices[face+vertex]].x,
                colours[faceVertices[face+vertex]].y,
                colours[faceVertices[face+vertex]].z
                );

            // use scaled xyz for vertex position
            glVertex3f
                (
                scale * vertices[faceVertices[face+vertex]].x,
                scale * vertices[faceVertices[face+vertex]].y,
                scale * vertices[faceVertices[face+vertex]].z
                );
            } // per vertex
        } // per face

    // close off the triangles
    glEnd();

    // revert render mode  
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    } // Render()


