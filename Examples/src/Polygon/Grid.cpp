#include "Grid.h"
#include <string>    /* std::string */
#include <sstream>   /* std::ifstream */
#include <fstream>   /* std::istringstream */
#include <set>       /* std::set */
#include <algorithm> /* std::set_difference */
#include <iterator>  /* std::inserter */
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
	for(auto const &i : grid.polygons)
	{
		i->showMe();
	}

	// Print edges 
	os << "Edges:" << std::endl;
	for(auto const &i : grid.edges)
	{
		os << i.first << " -- " << i.second << std::endl;
	}

	// Print boundary
	os << "Boundary:" << std::endl;
	for(auto const &i : grid.boundary)
	{
		os << i.first << " -- " << i.second << std::endl;
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
	polygons.reserve(NPolygons);
	
	// Read vertices 
	unsigned int idx = 0;
	double xCoord = 0.0;
	double yCoord = 0.0; 
	for(size_t i = 0; i < NPoints && getline(ifs, line); ++i)
	{
		std::istringstream linestream(line);
		linestream >> idx >> xCoord >> yCoord;
		vertices[idx] = Geometry::Point2D(xCoord, yCoord); 
	}
	
	// Read polygons 
	unsigned int polygonType = 0;
	unsigned int idxVertex = 0;
	Geometry::Vertices polygonVertices; 
	std::vector<unsigned int> idxVerticesVec;
	std::set<Edge> setEdges;
	std::set<Edge> setNonBoundaryEdges; 
	std::pair<std::set<Edge>::iterator, bool> ret;

	for(size_t i = 0; i < NPolygons && getline(ifs, line); ++i)
	{
		// Reset vertices to empty vector
		idxVerticesVec.resize(0);   
		polygonVertices.resize(0); 		

		// Read polygon index and type
		std::istringstream linestream(line);
		linestream >> idx >> polygonType;

		// Read vertices and edges and store them 
		while(linestream >> idxVertex)
		{
			// Insert vertex
			idxVerticesVec.push_back(idxVertex);
			polygonVertices.push_back(vertices[idxVertex]);
		}

		// Insert first vertex in the back, to close polygon 
		idxVerticesVec.push_back(idxVerticesVec[0]);
	
		// Create polygon object of appropriate type and store it
		switch (polygonType) 
		{
			case 0:  // Triangle
				polygons.emplace_back(new Geometry::Triangle(polygonVertices)); 
				break;
			case 1:  // Square
				polygons.emplace_back(new Geometry::Square(polygonVertices)); 
				break;
			case 2:  // Generic polygon
				polygons.emplace_back(new Geometry::Polygon(polygonVertices));
				break;
		}

		// Create edges 
		for(size_t i = 0; i < idxVerticesVec.size() - 1; ++i)
		{
			Edge edge(idxVerticesVec[i], idxVerticesVec[i+1]);
			ret = setEdges.insert(edge);

			if (!ret.second) // Edge does not belong to the boundary
			{
				setNonBoundaryEdges.insert(edge);
			}
		}
	}

	// Compute boundary edges by set difference
	std::set<Edge> setBoundaryEdges;
	std::set_difference(setEdges.begin(), setEdges.end(),
						setNonBoundaryEdges.begin(), setNonBoundaryEdges.end(),
						std::inserter(setBoundaryEdges, setBoundaryEdges.end()));

	// Convert sets to vectors
	edges.resize(setEdges.size());
	boundary.resize(setBoundaryEdges.size());
	std::copy(setEdges.begin(), setEdges.end(), edges.begin());
	std::copy(setBoundaryEdges.begin(), setBoundaryEdges.end(), boundary.begin());					
}


double Grid::area() const
{
	double totArea = 0.0;
	for (auto const &i : polygons)
	{
		totArea += i->area();
	}
	return totArea;
}
