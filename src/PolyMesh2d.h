//
//  PolyMesh2d.h
//  LPM
//
//  Created by Peter Bosler on 2/15/15.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

/** @file PolyMesh2d.h 
	@brief PolyMesh2d class header file
	@author Peter Bosler, Sandia National Laboratories, Multiphysics Applications
*/
#ifndef __LPM__PolyMesh2d__
#define __LPM__PolyMesh2d__

#include "Particles.h"
#include "Edges.h"
#include "Faces.h"
#include "OutputMessage.h"
#include "Logger.h"
#include "Field.h"
#include <iostream>
#include <vector>

/** @class PolyMesh2d
	@brief Manages meshes of Particles, Edges, and Faces for 2d topologies (e.g., plane and sphere).
	
	Current implementation for __planar meshes__ includes a mesh of quadrilateral faces from a square or a mesh of triangular faces from 
	a hexagon, with free boundaries.
	Users choose the initial mesh from one of the predefined "seeds," defined to approximately cover the domain
	@f$ D = \{ (x,y) : x,y\in [-1,1] \} @f$, depending on the chosen seed.  
	The seed is then multiplied by a user-supplied scalar `maxR` to cover larger subsets of the plane.  
	
	Implementations for __spherical meshes__ include a mesh of quadrilaterals from the cubed sphere or a mesh of triangles from
	an icosahedral triangulation of the sphere.
	Users choose one of the spherical seeds (which discretize a unit sphere) and may set maxR = @ref GlobalConstants::EarthRadius()
	to define an Earth-sized sphere.
	
	The mesh is recursively refined until a user-supplied initial depth is achieved in the Faces quadtree.  
	This refinement increases spatial resolution.  
	
	Adaptive refinement depends on the context of each problem and is handled by subclasses that include relevant @ref Field objects
	for storing and handling physical data.  
	
	@todo implement periodic boundary conditions
	@todo implement polar coordinate planer mesh
	@todo implement Voronoi/Delaunay meshes
	
	Currently implemented mesh seeds are illustrated below. @n 
	Particles are drawn as circles and annotated by @f$p_i@f$ where the subscript i indicates 
	the index of that particle in the Particles structure.@n
	Edges are illustrated as arrows @f$e_i@f$ and Faces by @f$F_i@f$.
	
	Seeds are recursively refined by either Faces::divideTriFace or Faces::divideQuadFace 
	until a desired spatial resolution is achieved.
	
	@image html QuadRectSeed.png "Seed for an initially square mesh of quadrilateral faces."
	@image latex QuadRectSeed.eps "Seed for an initially square mesh of quadrilateral faces." width=5in
	
	and
	
	@image html TriHexSeed.png "Seed for an initially hexagonal mesh of triangular faces." 
	@image latex TriHexSeed.eps "Seed for an initially hexagonal mesh of triangular faces." width=5in
*/
class PolyMesh2d
{
	public : 
		enum meshSeed { triHexSeed, ///< seed for a mesh of triangular faces whose union is a hexagon
					quadRectSeed, ///< seed for a mesh of quadrilateral faces whose union is a square
					polarDiscSeed, ///< seed for a mesh of polar coordinate rectangles (sectors of annuli)
					quadRectPeriodic, ///< seed for a mesh of quadrilateral faces whose union is a doubly-periodic square
					icosTriSphereSeed, ///< seed for an icosahedral triangulation of the unit sphere`
					cubedSphereSeed ///< seed for a cubed sphere discretization of the unit sphere
				 };
		/** @brief constructor
			@param initialNest number of times each face of the meshSeed will be recursively divided.
			@param meshSeed seed to define the initial mesh
			@param maxR amplification factor applied to initial meshSeed
			@param procRank
			@param numProcs
			@param gkin coordinate kind definition (e.g, Cartesian or Spherical)
		*/		 
		PolyMesh2d( const int intialNest = 0, const meshSeed seed = triHexSeed,
					const double maxR = 1.0, const int procRank = 0, const int numProcs = 1 );
		
		/// Returns the number of Faces in a PolyMesh2d.
		int nFaces() const { return faces.N(); }
		/// Returns the number of undivided Faces in a PolyMesh2d
		int nActiveFaces() const { return faces.nActive(); }
		/// Returns the number of Edges in a PolyMesh2d
		int nEdges() const {return edges.N(); }
		/// Returns the number of Particles in a PolyMesh2d
		int nParticles() const {return particles.N();}
		/// Returns the number of parent faces
		int nDividedFaces() const { return faces.countDivided(); }
		/// Returns the mesh seed used to initialize a PolyMesh2d object
		meshSeed getSeed() const { return _seed; }
		/// Returns the number of components in particles' coordinate vectors
		int nDim() const { return particles.nDim(); }
		/// Returns the faceKind 
		Faces::facesKind getFaceKind() const { return faces.kind(); }
		
		Particles* getParticles() { return &particles; }
		
		double particleArea( const int particleIndex ) const { return particles.area[particleIndex]; }
		
		/** @brief Reduces the memory used by a PlanarMesh to only the number of initialized members of its members.  
			@warning This function should only be called when no additional faces will be added, i.e., when the
			connectivity of Faces, Edges, and Particles will be constant for the rest of the PlanarMesh lifetime.
		*/
		virtual void shrinkMemory();
		
		/** @brief Returns the index (in a Particles object) of a face's center particle.
		
			@param faceIndex index of face whose center particle is needed
			@return index to Particles corresponding to the face's center particle.
		*/	
		int faceParticleIndex( const int faceIndex ) const { return faces.centerParticle(faceIndex);}
	
		/** @brief returns the physical position of the active particle interior to a face.
			@param faceIndex index to a face of Faces
		*/					
		xyzVector facePosition( const int faceIndex ) const;

		/** @brief returns a vector corresponding to the centroid of a face (defined by its vertices).
			@param faceIndex index to a face in Faces
		*/
		xyzVector faceCentroid( const int faceIndex ) const;

		/** @brief returns the phsyical position of a particle in Particles associated with this mesh.
			@param particleIndex index to a particle (either a vertex or a center) in Particles
		*/
		xyzVector particlePosition( const int particleIndex ) const { return particles.physCoord(particleIndex); }
		
		xyzVector particleLagCoord( const int particleIndex ) const { return particles.lagCoord(particleIndex); }

		/** @brief Sets the surface area of a mesh. 
			@warning This method assumes the surface area has already been computed.
		*/
		void setSurfaceAreaFromParticles();
		
		/** @brief Sets the surface area of a mesh, using the vertices of each face.
		*/
		void setSurfaceAreaByCalculation();
	
		/** @brief Returns the surface area in physical coordinates covered by a PlanarMesh.
			@return _surfaceArea
		*/
		double surfaceArea() const { return _surfaceArea; }
		
		/** @brief Output mesh data to a legacy VTK file.  Filename should be `*.vtk`
			@param filename with extension `.vtk`
		*/
		virtual void outputToVTK( const std::string filename, const std::vector<Field> fields = std::vector<Field>() ) const;
	
	    void writeActiveParticlesToCSV(const std::string filename) const;
	    
	    void writeActiveParticlesToVTK(const std::string filename, const std::string title = " ") const;
	
		/** @brief Output mesh data in a matlab-readable ascii file.  Filename should be `*.m`
			@param filename with extension `.m`
		*/
		virtual void outputToMatlab( const std::string filename ) const;
	
		/** @brief Output mesh data in NetCDF format.  Filename should be `*.nc`
			@param filename with extension `.nc`
		*/	
		virtual void outputToNetCDF( const std::string filename ) const{};
	
		/** @brief Output mesh data in Exodus format.  Filename should be `*.g`
			@param filename with extension `.g`
		*/
		virtual void outputToExodus( const std::string filename ) const{};
	
		/** @brief Returns a collection of strings containing info about the state of a PlanarMesh object.
	
			Use for building LongMessage objects to send to a Logger.
		*/
		std::vector< std::string > getInfo() const ;

		/** @brief Returns the index to the face (in Faces) that contains the query point.
			Returns -1 if the queryPt is outside the mesh.
			@param queryPt
			@return index to face closest to queryPt.  
		*/
		int locatePointInMeshFace( const xyzVector& queryPt ) const;
		
		/** @brief Returns true if a point is outside of the mesh (Planar geometries only )
			@param queryPt
		*/
		bool pointIsOutsideMesh( const xyzVector& queryPt ) const;
		
		/** @brief Returns a counter-clockwise ordered list of edges (in Edges) around a face.
			@todo Change this function to support AMR, differentiate it from Faces::edges.
			@param faceIndex index to a face (in Faces) whose edges are needed
			@return CCW-ordered list of indices to edges (in Edges)
		*/
		std::vector<int> ccwEdgesAroundFace( const int faceIndex ) const { return faces.edges(faceIndex); }
		
		/** @brief Returns a counter-clockwise ordered list of faces (in Faces) sharing an edge with the input face.
			@param faceIndex face whose neighbors are needed
			@return CCW-ordered list of adjacent faces
		*/
		std::vector<int> ccwAdjacentFaces( const int faceIndex ) const;
		
		/** @brief Returns a counter-clockwise order list of vertices (in Particles) in a face.
			@todo Change this function to support AMR, differentiate it from Faces::vertices.	
			@param faceIndex face whose vertices are needed
			@return CCW-ordered list of vertices (in Particles) around face
		*/
		std::vector<int> ccwVerticesAroundFace( const int faceIndex ) const { return faces.vertices(faceIndex); }
		
		/** @brief Returns a counter-clockwise ordered list of faces (in Faces) around a vertex (in Particles), i.e., like a dual mesh.
			@param particleIndex index to a vertex (in Particles) whose incident faces are needed
			@return ccw-ordered list of indices to faces (in Faces) that touch the input vertex.
		*/
		std::vector<int> ccwFacesAroundVertex( const int particleIndex ) const;
	
		/** @brief True if face has children */
		bool faceIsDivided( const int faceIndex ) const { return faces.hasChildren(faceIndex); }
		
	protected : 
		Particles particles;
		Edges edges;
		Faces faces;
		int _nCoordDim;
		xyzVector::geometryKind gk;
		OutputMessage::priority logLevel;
		Logger* log;
		double _surfaceArea;
		double _maxR; 
		double _meshSpacing;
		meshSeed _seed;
		Faces::facesKind _faceKind;
		int _initNest;
		
		/** @brief Initializes the root mesh from the chosen mesh seed.  Mesh must have been constructed
			prior to calling this function.
		
			Mesh seeds are defined here.
			Each mesh seed begins with maximum particle distance from the origin = 1.0.
			Meshes that need to cover more area should use an ampFactor > 1.
		
			@param ampFactor
		*/
		void initializeFromSeed(const double ampFactor = 1.0 );


		/** @brief Returns the number of vertices in a PlanarMesh whose Faces tree has depth = initNest.
			@param initNest
			@return nVertices
		*/
		int nVertices( const int initNest ) const;
	
		/** @brief Returns the number of faces in a PlanarMesh whose Faces tree has depth = initNest.
			@param initNest
			@return nFaces
		*/
		int nFaces( const int initNest ) const;
	
		/** @brief Returns the number of edges in a PlanarMesh whose Faces tree has depth = initNest.
		
			Uses Euler's formula for an open polyhedron, @f$ F + V - E = 1 @f$.  
			Therefore, the number of vertices and the number of faces in the mesh must be known prior to
			calling this function (use PlanarMesh::nVertices and PlanarMesh::nFaces).
			@param nVert = V
			@param nFace = F
			@return nEdges = E
		*/
		int nEdges( const int nVert, const int nFace ) const;
	
		/** @brief Sorts the incident edges at each vertex (locally) into counter-clockwise order.
		*/
		void sortEdgesAroundVertices();
		
		/** @brief Estimates the mesh size (units of length) using the surface area and the number of faces */
		void estimateMeshSize(){ _meshSpacing = std::sqrt( _surfaceArea / (double)faces.N() ); };

		/** @brief Used to initialize a tree search with PolyMesh2d::locatePointTreeSearch.  
			Locates the nearest root face (a face from the initial mesh seed) to a given query point.
			@param queryPt  
			@return index of nearest root face (in Faces)
		*/
		int nearestRootFace( const xyzVector& queryPt ) const;
		
		/** @brief Recursive function. 
			Performs a tree search point-query algorithm.  
			Must be initialized by PolyMesh2d::nearestRootFace
			@param queryPt
			@param startIndex
			@return closest leaf in Faces tree to queryPt
		*/
		int locatePointTreeSearch( const xyzVector& queryPt, int startIndex ) const;
		
		/** @brief Recursive function.
			Performs a walk-search point query algorithm.
			Performance is determined by the skill of the initial guess.
			@param queryPt
			@param startIndex
			@return index to face (in Faces) closest to queryPt
		*/
		int locatePointWalkSearch( const xyzVector& queryPt, const int startIndex ) const;
};

std::ostream& operator<<( std::ostream& os, const PolyMesh2d& aMesh);

#endif