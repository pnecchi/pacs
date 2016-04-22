#include "Grid.h"
#include <string>   /* std::string */
#include <sstream>  /* std::ifstream */
#include <fstream>  /* std::istringstream */
#include <iostream> 

std::ostream& operator<<(std::ostream &os, const Grid& grid)
{
	// Print vertices
	os << "Vertices:" << std::endl;
	os << "X\tY" << std::endl;
	for(auto const &i : grid.vertices)
	{
		os << i.x() << "\t" << i.y() << std::endl;
	}

	// Print polygons
	os << "Polygons:" << std::endl;
	for(auto const &i : grid.polygonVec)
	{
		i->showMe();
	}	
	return os;
}

Grid::Grid(const std::string &filename)
{
	// Initialize filestream from filename
	std::ifstream ifs(filename);
	std::string line;
	
	// Read number of points and number of polygons from first line 
	size_t NPoints = 0;
	size_t NPolygons = 0;
	if (getline(ifs, line))
	{
		std::istringstream linestream(line);
		linestream >> NPoints >> NPolygons;
	}
	
	// Resize vector of vertices and of polygons
	vertices.resize(NPoints);
	polygonVec.reserve(NPolygons);
	
	// Read vertices 
	int idx = 0;
	double xCoord = 0.0;
	double yCoord = 0.0; 
	for(size_t i = 0; i < NPoints && getline(ifs, line); ++i)
	{
		std::istringstream linestream(line);
		linestream >> idx >> xCoord >> yCoord;
		vertices[idx] = Geometry::Point2D(xCoord, yCoord); 
	}
	
	// Read polygons 
	int polygonType = 0;
	int idxVertex = 0; 
	Geometry::Vertices polygonVertices; 
	for(size_t i = 0; i < NPolygons && getline(ifs, line); ++i)
	{
		polygonVertices.resize(0);  // Reset vertices to empty vector
		
		// Read polygon index and type
		std::istringstream linestream(line);
		linestream >> idx >> polygonType;
		
		// Read vertices and store them 
		while(linestream >> idxVertex)
		{
			polygonVertices.push_back(vertices[idxVertex]);
		}
	
		// Create polygon object of appropriate type and store it
		switch (polygonType) 
		{
			case 0:  // Triangle
				polygonVec.emplace_back(new Geometry::Triangle(polygonVertices)); 
				break;
			case 1:  // Square
				polygonVec.emplace_back(new Geometry::Square(polygonVertices)); 
				break;
			case 2:  // Generic polygon
				polygonVec.emplace_back(new Geometry::Polygon(polygonVertices));
				break;
		}
	}
}


double Grid::area() const
{
	double totArea = 0.0;
	for (auto const &i : polygonVec)
	{
		totArea += i->area();
	}
	return totArea;
}
