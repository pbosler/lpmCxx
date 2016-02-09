//
//  Particles.h
//  LPM-Testing
//
//  Created by Bosler, Peter Andrew on 11/3/14.
//  Copyright (c) 2014 Bosler, Peter Andrew. All rights reserved.
//

/** @file Particles.h
	@brief Particles class header file.
	@author Peter Bosler <pabosle@sandia.gov>
 */

#ifndef __LPM_Testing__Particles__
#define __LPM_Testing__Particles__

#include <iostream>
#include <vector>
#include "OutputMessage.h"
#include "Logger.h"
#include "xyzVector.h"

/** @class Particles
	@brief Particles class defines a "structure of arrays" for organizing a set of particles to discretize space.
	
	Particles are identified by their index in the class, so the ith particle has x-coordinate x[i], y-coordinate y[i], 
	z-coordinate z[i], etc.  
	Most particles will have zero area and zero volume, indicating that they are vertices that do not correspond to a region of space.
	If nonzero, area and volume values are defined by another class (Faces or Cells).  
	Methods from those classes may assign a nonzero value to a particle's area or volume for use with quadrature schemes.  
	
	Each particle with zero area (all vertices) maintains a record of edges incident to that vertex.  
	These records are double-integer pairs.
	The double value is between @f$-\pi@f$ and @f$pi@f$ and represents the angle of incidence of the edge
	relative to a coordinate system centered at that particle.
	The integer value is the edge index. 
	This enables quick sorting of the edges about each vertex into ccw order.  
	
	The Field class contains scalar and vector field data defined upon a particle set.  
	Each particle has a 1-to-1 correspondence with an element of each Field.   

*/
class Particles {
    
private:
    int _N; ///< number of initialized particles
    int _nMax; ///< maximum number of particles allowed in memory
    int _nDim; ///< coordinate dimension (must be 2 or 3)
    OutputMessage::priority logLevel; ///< base priority level for Particles log output
    Logger* log; ///< Logger for Particles objects
    xyzVector::geometryKind gk; ///< coordinate type (Euclidean, spherical, etc.)
       
public:
    std::vector<double> x; ///< x-coordinates in physical space
    std::vector<double> y; ///< y-coordinates in physical space
    std::vector<double> z; ///< z-coordinates in physical space
    std::vector<double> x0;///< x-coordinates in Lagrangian space
    std::vector<double> y0;///< y-coordinates in Lagrangian space
    std::vector<double> z0;///< z-coordinates in Lagrangian space
    std::vector<double> area; ///< area represented by each Particles
    std::vector<double> volume; ///< volume represented by each particle (3d only)
    std::vector< std::vector< std::pair<double,int> > > incidentEdges; ///< incident edge record for each particle (null for non-vertices)
    
    void shrinkMemory();
    
    /** @brief Constructor.
    	@param nDim
    	@param nMax
    	@param procRank (for log initialization)
    	@param numProcs (for log initialization)
    */
    Particles( const int nDim = 2, const int nMax = 7, 
    		   const int procRank = 0, const int numProcs = 1, 
    		   const xyzVector::geometryKind gkin = xyzVector::EuclideanGeometry);
    
	/// returns the coordinate dimension
    int nDim() const { return _nDim; }
    
    /// returns the current number of particles
    int N() const { return _N; }
    
    /// returns the maximum allowable number of particles
    int nMax() const {return _nMax; }
    
    /** converts a particle's physical coordinates to an xyzVector
    	@param index
    */
    xyzVector physCoord( const int index ) const;
    
    /** converts a particle's Lagrangian coordinates to an xyzVector
    	@param index
    */
    xyzVector lagCoord( const int index ) const; 
    
    /** @brief Inserts a new particle at the end of the data structure.  
    	@param physCoords
    	@param lagCoords
    */
    void insertParticle( const xyzVector& physCoords, const xyzVector& lagCoords);
    
    /** @brief Sets the physical coordinates of the particle at given index.
    	@param index
    	@param physCoords
    */
    void setPhysicalCoords( const int index, const xyzVector& physCoords );
    
    /** @brief Sets the physical coordinates of the particle at given index.
    	@param index
    	@param nx
    	@param ny
    	@param nz
    */
    void setPhysicalCoords( const int index, const double nx, const double ny, const double nz = 0.0);
    
    /** @brief Sets the Lagrangian coordinates of a particle at the given index.
    	@param index
    	@param lagCoords
    */
    void setLagrangianCoords( const int index, const xyzVector& lagCoords);
    
	/** @brief Sets the Lagrangian coordinates of a particle at the given index.
		@param index
    	@param ax
    	@param ay
    	@param az
    */
    void setLagrangianCoords( const int index, const double ax, const double ay, const double az = 0.0);
    
    /** @brief amplifies each coordinate component by ampFactor
    	@param ampFactor
    */
	void rescale( const double ampFactor );
  
    /// Returns the minimum x-coordinate value in physical space.
    double minX() const;
    /// Returns the maximum x-coordinate value in physical space.
    double maxX() const;
    /// Returns the minimum y-coordinate value in physical space.
    double minY() const;
    /// Returns the maximum y-coordinate value in physical space.
    double maxY() const;
    /// Returns the minimum z-coordinate value in physical space.
    double minZ() const;
    /// Returns the maximum z-coordinate value in physical space.
    double maxZ() const;
    
    /// Returns the minimum x-coordinate value in Lagrangian space.
    double minX0() const;
    /// Returns the maximum x-coordinate value in Lagrangian space.
    double maxX0() const;
    /// Returns the minimum y-coordinate value in Lagrangian space.
    double minY0() const;
    /// Returns the maximum y-coordinate value in Lagrangian space.
    double maxY0() const;
    /// Returns the minimum z-coordinate value in Lagrangian space.
    double minZ0() const;
    /// Returns the maximum z-coordinate value in Lagrangian space.
    double maxZ0() const;
    /// Returns the total area of the space represented by a 2d particle set
    double totalArea() const;
    /// Returns the total volume of the space represented by a 3d particle set
    double totalVolume() const;
    
    /// returns general info about a Particles object for logging and console output.
    std::vector< std::string > getInfo() const;
    
//    void recordIncidentEdgeAtParticle(const int particleIndex, const int edgeIndex);
    
    /// Sorts the incident edges at a vertex into ccw order.
    void sortIncidentEdgesAtParticle( const int particleIndex );
    
    /// Writes Particles data to a file in a format readable by Matlab.  
    void writeVariablesToMatlab( std::ostream& fs ) const;
    
    /// Writes particle positions to a .vtk file.
    void writePointsToVTK( std::ostream& fs, std::string title = "--" ) const;
    
    void writeLagCoordToVTK( std::ostream& fs ) const;
    
    std::vector<int>  edgesAtVertex( const int particleIndex ) const;
    
    void printIncidentEdges() const;
};

/// basic console output
std::ostream& operator << ( std::ostream& os, const Particles& aParticles);


#endif /* defined(__LPM_Testing__Particles__) */
