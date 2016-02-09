//
//  Edges.h
//  LPM
//
//  Created by Peter Bosler on 11/5/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

/** @file Edges.h 
	@brief Edges class header file
	@author Peter Bosler, Sandia National Laboratories, Multiphysics Applications
*/

#ifndef __LPM__Edges__
#define __LPM__Edges__

#include <iostream>
#include <vector>
#include "OutputMessage.h"
#include "Logger.h"
#include "xyzVector.h"

class Particles;

/** @class Edges
	@brief Directed edges to connect Particles, with dual mesh information available.
	This class is designed to be used by PolyMesh2d objects, not by itself; hence, many of its attributes are protected,
	unlike the Particles class, which can be used alone.
	
	Edge i is accessed with the same index for all arrays.
	Its origin is given by `edges.orig(i)`, which returns an integer that corresponds to the origin particles
	in a Particles object.  
	Similarly, `edges.dest(i)` returns the index of the destination Particle.
	`edges.leftFace(i)` and `edges.rightFace(i)` return indices to elements of a Faces object.
	
	Edges are maintained in a binary tree structure.  
	Divided edges create two child edges.  
	
	Edges may be subclassed to add more particles to each edge, for high-order mesh reconstructions.  	
	Note: this will affect the convexity of faces.  
	
	Dual mesh information is recorded at each vertex in the Particles class via the Edges::recordIncidentEdgeAtParticles method.  
*/
class Edges
{
	public:
		/** @brief Constructor.
		*/
		Edges( const int nMax = 12, const int procRank = 0, const int numProcs = 1, 
			   const xyzVector::geometryKind gkin = xyzVector::EuclideanGeometry);
		/// Destructor
		virtual ~Edges(){};
		
		/** @brief Returns the origin (index to a member of a Particles object) of an edge.
			@param index index of edge in Edges object
			@return index of origin particle in Particles object
		*/
		int orig( const int index ) const { return _orig[index]; }
		/** @brief Returns the destination (index to a member of a Particles object) of an edge.
			@param index index of edge in Edges object
			@return index of destination particle in Particles object
		*/		
		int dest( const int index ) const { return _dest[index]; }
		/** @brief Returns the face (index to a member of a Faces object) to the left of a directed edge.
			@param index index of edge in Edges object
			@return index of face to the left of directed edge in Faces object
		*/
		int leftFace( const int index ) const { return _leftFace[index]; }
		/** @brief Returns the face (index to a member of a Faces object) to the right of a directed edge.
			@param index index of edge in Edges object
			@return index of face to the right of directed edge in Faces object
		*/
		int rightFace( const int index) const { return _rightFace[index]; }
		/** @brief Returns the parent edge (if any) of an edge.
			@param index of edge in Edges object
			@return index of parent edge in Edges object
		*/
		int parent( const int index ) const { return _parent[index]; }
		/** @brief Returns the first child of an edge (an edge and its first child have the same origin)
			@param index index of edge in Edges object
			@param index of child edge 1 in Edges object
		*/
		int child1( const int index ) const { return _child1[index]; }
		/** @brief Returns the second child of an edge (an edge and its second child have the same destination)
			@param index index of edge in Edges object
			@param index of child edge 2 in Edges object
		*/
		int child2( const int index ) const { return _child2[index]; }
		/** @brief Returns the number of edges currently defined in an Edges object.
		*/
		int N() const { return _N; }
		/** @brief Returns the maximum number of edges allowed in memory for an Edges object.
		*/
		int nMax() const { return _nMax; }		
		
		/** @brief Divides an edge to create two child edges.
			@param index index of edge in Edges to be divided.
			@param particles Particles object associated with Edges
		*/		
		virtual void divideEdge( const int index, Particles& particles );
		/** @brief Inserts a new edge into an Edges object (inserted edges have no parents).
			@param origIndex index to origin vertex (in Particles) of new edge
			@param destIndex index to destination vertex (in Particles) of new edge
			@param left index to left face of new edge (in Faces) 
			@param right index to right face of new edge (in Faces)
			@param particles Particles object assocated with Edges
		*/
		virtual void insertEdge( const int origIndex, const int destIndex, const int left, const int right, Particles& particles );
		
		/** @brief sets the left face of an edge.
			@param index index of edge (in Edges)
			@param newLeft index of new left face (in Faces)
		*/
		void setLeftFace( const int index, const int newLeft ) { _leftFace[index] = newLeft; }
		/** @brief sets the right face of an edge.
			@param index index of edge (in Edges)
			@param newRight index of new right face (in Faces)
		*/
		void setRightFace( const int index, const int newRight) { _rightFace[index] = newRight; }
		
		/** @brief Releases unused memory, if too much was allocated.
			Call this method after an Edges object has been initialized, and only if no additional particles will be added.
		*/
		void shrinkMemory();
		
		/// True if edge(index) has been divided
		bool isDivided( const int index ) const { return _hasChildren[index]; }
		/// True if edge(index) has not been divided, which indicates it is a member of the current active mesh.
		bool isActive( const int index ) const { return !_hasChildren[index]; }
		/// True if edge(index) has been divided
		bool hasChildren( const int index ) const { return _hasChildren[index]; }
		
		/// True for planar meshes if an edge is a boundary edge of the mesh.
		bool onBoundary( const int index ) const { return ( _leftFace[index] < 0 || _rightFace[index] < 0 ); }
		
		/** @brief Records an edge at each of its vertices in a Particles object to provide dual mesh information.
			@param edgeIndex index of edge (in Edges) to be recorded at each of its vertices (in Particles)
			@param particles Particles associated with *this Edges.
		*/
		void recordIncidentEdgeAtParticles( const int edgeIndex, Particles& particles ) const;
		
		/// Writes Edges data in a Matlab-readable format.
		void writeVariablesToMatlab( std::ostream& fs ) const;
		
		/// returns general info about a Edges object for logging and console output.
		std::vector< std::string > getInfo() const;
		
	protected:
		std::vector<int> _orig;
		std::vector<int> _dest;
		std::vector<int> _rightFace;
		std::vector<int> _leftFace;
		std::vector<int> _parent;
		std::vector<bool> _hasChildren;
		std::vector<int> _child1;
		std::vector<int> _child2;
		int _nActive;
		int _N;
		int _nMax;
		OutputMessage::priority logLevel;
		Logger* log;
		xyzVector::geometryKind gk;
		
		void replaceIncidentEdgeWithChild( const int parentIndex, const int child1Index, Particles& particles ) const;
		double edgeAngleAtOrig( const int edgeIndex, const Particles& particles ) const;
		double edgeAngleAtDest( const int edgeIndex, const Particles& particles ) const;
};

std::ostream& operator << ( std::ostream& os, const Edges& edges);

#endif