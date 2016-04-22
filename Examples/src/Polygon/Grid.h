//----------------------------------------------------------------------
// Description: Grid class 
// Author:      Pierpaolo Necchi
// Email:       pierpaolo.necchi@gmail.com
// Date:        ven 22 apr 2016 16:19:52 CEST
//----------------------------------------------------------------------

#ifndef GRID_H
#define GRID_H

#include <ostream>
#include <vector>
#include <string>
#include <memory>
#include "Polygon.hpp"

/** \class Grid
 *	\brief Grid is a set of vertices which form a set of polygons. 
 */
class Grid
{

	friend std::ostream& operator<<(std::ostream &os, const Grid& grid);

public:
	// Constructor
	Grid()=default;
	
	/** 
	 * Constructor from file 
	 * Initializes a Grid object by reading the grid stored in a file.  
	 * @param filename the name of the input file (with relative path) 
	 */
	Grid(const std::string &filename);

	// Copy-constructor
	Grid(const Grid&)=default;

	// Assigment operator
	Grid& operator=(const Grid&)=default;
	
	// Destructor
	virtual ~Grid ()=default;

	// Size
	size_t nVertices() const { return vertices.size(); }
	size_t nPolygons() const { return polygonVec.size(); }

private:
	std::vector<Geometry::Point2D> vertices;
	std::vector<std::shared_ptr<Geometry::AbstractPolygon>> polygonVec;
};

// Friend functions declaration
std::ostream& operator<<(std::ostream &os, const Grid& grid);

#endif /* end of include guard: GRID_H */
