//
//  Faces.h
//  LPM
//
//  Created by Peter Bosler on 11/5/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

/** @file Faces.h 
	@brief Faces class header file
	@author Peter Bosler, Sandia National Laboratories, Multiphysics Applications
*/

#ifndef __LPM__Faces__
#define __LPM__Faces__

#include <iostream>
#include <vector>
#include "OutputMessage.h"
#include "Logger.h"
#include "Particles.h"
#include "Edges.h"
#include "xyzVector.h"

/** @class Faces
	@brief Structure of arrays, where each index corresponds to one polygonal face which connects Edges and Particles.  
	
	Faces have vertices (indices of Particles) and edges (indices of Edges), and one center index also to Particles.
	Faces are responsible for setting the area represented by their center particle and maintaining counter-clockwise 
	orientation of the vertices and edges.
	
	Faces are maintained in a quad-tree structure; divided faces create 4 children.
	
	Faces objects with constant type (e.g., all quadrilateral or all triangular) take advantage of that known structure
	by enforcing a constant counter-clockwise ordering of edges and vertices relative to each individual face:
	
	@image html TriFace.png "Vertex and edge ordering relative to an individual triangular face."
	@image latex TriFace.eps "Vertex and edge ordering relative to an individual triangular face." width=5in
	
	@image html QuadFace.png "Vertex and edge ordering relative to an individual quadrilateral face."
	@image latex QuadFace.eps "Vertex and edge ordering relative to an individual quadrilateral face." width=5in
	
	@todo Add Voronoi faces (variable types of polygons)
*/
class Faces
{
	public:
		enum facesKind { triangularFaces, ///< each face is a triangle
						quadrilateralFaces ///< each face is a quadrilateral
					   };
	
		/** @brief main constructor.  Allocates memory but does not initialize geometry.
			@param nMax
			@param kind
			@param procRank
			@param numProcs
			@param gkin
		*/
		Faces( const int nMax = 6, const facesKind kind=triangularFaces, const int procRank = 0, const int numProcs = 1,
			   const xyzVector::geometryKind gkin = xyzVector::EuclideanGeometry);	   
	
		/// Returns the current number of Faces
		int N() const { return _N; }
		/// Returns the current number of undivided Faces
		int nActive() const { return _nActive; }
		/** @brief Sets the current number of undivided Faces.
			@warning For use ONLY during mesh initialization, for example, in @ref PolyMesh2d::initializeFromSeed.
			Otherwise, _nActive is set in member functions Faces::divideTriFace and Faces::divideQuadFace.
		*/
		void setNActive( const int nA ){ _nActive = nA; }
		
		/// Returns the maximum number of Faces allowed in memory
		int nMax() const { return _nMax; }
	
		/// Returns the area of a face
		double area( const int index, const Particles& particles ) const { return particles.area[ _centerParticle[index] ]; }
		
		/// Returns the index of a Particles member associated with a Faces center.
		int centerParticle( const int index ) const { return _centerParticle[index]; }
		
		/** @brief Sets the area represnted by face._centerParticle[index]

			Computes the sum of a face's planar subtriangles.  
	
			@param index
			@param particles
			@param edges
		*/
		void setArea( const int index, Particles& particles, const Edges& edges ) const;
		/// Sets the area of a face to zero, sets the value of area carried by that face's center particle to zero.
		void setAreaToZero( const int index, Particles& particles ) const;
		
		/** @brief Returns a counter-clockwise ordered list of edges around a face
			@param index
		*/		
		std::vector<int> edges( const int index) const { return _edges[index]; }
		/** @brief Returns a counter-clockwise ordered list of vertices around a face
			@param index
		*/
		std::vector<int> vertices( const int index) const{ return _vertices[index]; }
		int nVerts( const int index ) const { return _vertices[index].size(); }
		
		std::vector<int> children( const int index) const{ return _children[index]; }
		
		int parent( const int index ) const { return _parent[index]; }
		
		/** Sets the edges around a face, in counter-clockwise order
			@param index
			@param newEdges CCW-ordered list of indices to @ref Edges
		*/
		void setEdges( const int index, const std::vector<int> newEdges ){ _edges[index] = newEdges; }
		
		/** Sets the vertices around a face, in counter-clockwise order
			@param index
			@param newVerts CCW-order list of indices to Particles
		*/
		void setVertices( const int index, const std::vector<int> newVerts ){ _vertices[index] = newVerts; }
		
		/** @brief Inserts a face to the end of the Faces data structure.
			@param cntrParticle = index to particle at center of new face (index to a member of @ref Particles )
			@param vertices = CCW list of vertices around new face (indices to member of @ref Particles)
			@param edges = CCW list of edges around new face (indices to members of @ref Edges )
		*/
		void insertFace( const int cntrParticle, std::vector<int> vertices, std::vector<int> edges );
		
		/** @brief Divides a triangular face to create 4 child Faces
			
			Automatically adds @ref Edges and @ref Particles as required.
			
			@image html DivideTriFace.png  "Division of a parent face to create four children."
			@image latex DivideTriFace.eps "Division of a parent face to create four children." width=5in
		*/
		void divideTriFace( const int index, Particles& particles, Edges& edges );
		
		/** @brief Divides a quadrilateral face to create 4 child Faces
			
			Automatically adds @ref Edges and @ref Particles as required.
			
			@image html DivideQuadFace.png  "Division of a parent face to create four children."
			@image latex DivideQuadFace.eps "Division of a parent face to create four children." width=5in
		*/
		void divideQuadFace( const int index, Particles& particles, Edges& edges );
		
		/** @brief Returns true if Edge[edgeIndex] is oriented positively with respect to Face[faceIndex].
			@param faceIndex
			@param edgeIndex
			@param edges
		*/
		bool positiveEdge( const int faceIndex, const int edgeIndex, const Edges& edges ) const;
		
		/** @brief Returns true if Face has not been subdivided.
			@param faceIndex
		*/
		bool isActive( const int faceIndex ) const { return !_hasChildren[faceIndex]; }
		bool hasChildren(const int index) const { return _hasChildren[index]; }
		
		/// @brief Returns the faceKind attribute of a Faces data object.
		facesKind kind() const { return _kind; }
		
		/** @brief Reduces the memory used by a Faces data object to only the number of initialized faces.  
			@warning This function should only be called when no additional faces will be added, i.e., when the
			connectivity of Faces, Edges, and Particles will be constant for the rest of their scope lifetimes.
		*/
		void shrinkMemory();
		
		/** @brief Writes the state of a Faces data object to an output stream in a format readable by Matlab.
			@param fs = ostream object already opened
		*/
		void writeVariablesToMatlab( std::ostream& fs ) const;
		
		void writePolygonsToVTK( std::ostream& fs ) const;
		
		int countDivided() const;
		
		/// returns general info about a Faces object for logging and console output.
		std::vector< std::string > getInfo() const;
		
	protected:
		std::vector<int> _centerParticle; ///< indices to members of a @ref Particles data object
		std::vector< std::vector<int> > _vertices; ///< indices to members of a @ref Particles data object
		std::vector< std::vector<int> > _edges; ///< indices to members of an @ref Edges data object
		std::vector<bool> _hasChildren; ///< _hasChildren[i] = true if Face [i] has been divided
		std::vector< std::vector<int> > _children; ///< indices to members of a @ref Faces data object, quadtree organization
		std::vector<int> _parent; ///< indices to members of a @ref Faces data object, quadtree organization
		facesKind _kind; ///< faceKind identifier
		int _N; ///< Current number of initialized faces in memory
		int _nMax; ///< Maximum number of faces allowed in memory
		int _nActive; ///< Current number of initialized and undivided faces
		xyzVector::geometryKind gk;
		OutputMessage::priority logLevel; 
		Logger* log;
};

//std::ostream& operator << ( std::ostream& os, const Faces& faces){};

#endif