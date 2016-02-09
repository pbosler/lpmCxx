/** @file ParticlesAndFieldTest.cpp
	@author Peter Bosler, Sandia National Laboratories
	@brief Unit test for mesh primitives, Particles and Field.
*/

#include <iostream>
#include "Particles.h"
#include "xyzVector.h"
#include "Field.h"
#include <cmath>
#include "lpmConfig.h"
#include <sstream>

int main ( int argc, char* argv[] )
{
	
	const int procRank = 0;
	const int numProcs = 1;
	const OutputMessage::priority loggingLevel = OutputMessage::debugPriority;
	
	Logger* exeLog = Logger::Instance(loggingLevel, procRank, numProcs);

	std::stringstream ss;
    std::string s;
    ss << "-- LPM Version " << LPM_VERSION_MAJOR << "." << LPM_VERSION_MINOR << " -- \n ";
    ss << "   Unit Test " << argv[0] << " \n";
    ss << "		 This test builds a Cartesian particle set, defines a scalar field and a vector field on those particles,\n"
       << " then outputs the particles and fields to a Matlab-readable .m file. \n";
    OutputMessage introMsg( ss.str(), OutputMessage::remarkPriority, "main");
    exeLog->logMessage(introMsg);

	const int nn = 101;
	const double dx = 0.1;
	Particles p( 2, nn*nn);

	std::cout << "returned from particles constructor.\n";
	
	std::cout << "particles N = " << p.N() << std::endl;
	std::cout << " size x = " << p.x.size() << std::endl;
	std::cout << p << std::endl;
	
	double xi;
	double yj;
	int k = 0;
	xyzVector phys, lag;
	for ( int i = 0; i < nn; ++i)
	{
		xi = -5.0 + i * dx;
		for ( int j = 0; j < nn; ++j)
		{
			yj = -5.0 + j * dx;	
			phys = xyzVector(xi, yj);
			lag	= phys;
			
			p.insertParticle( phys, lag);

			k++;
		}
	}
	
	LongMessage particleStatus("Particle set built", OutputMessage::tracePriority, "main", p.getInfo() );
	exeLog->logMessage(particleStatus);

	Field f("testScalar", "n/a", 1, p.N());
	for (  int  k = 0; k < p.N(); ++k)
	{
		xi = p.x[k];
		yj = p.y[k];
		f.scalar[k] = std::exp( -xi * xi - yj * yj );
	}
	
	f.outputForMatlab( "testScalarField.m", p);
	std::cout << "scalar field test complete.\n";

	Field vf("testVector", "n/a", 2, p.N());
	for ( int k = 0; k < p.N(); ++k)
	{
		xi = p.x[k];
		yj = p.y[k];
		vf.insertVectorToField( k, -yj, xi );
	}
	
	vf.outputForMatlab( "testVectorField.m", p);
	
	std::cout << "vector field test complete.\n";
	
	return 0;
};