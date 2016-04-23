/*#include "Polygon.hpp"*/
//#include <iostream>
////! Main program
//int main()
//{
  //using namespace Geometry;

  ////! five vertices
  //Vertices v(5);

  //v[0]={0.5,0.5}; //C++11 sintax!
  //v[1]={0.5,0.8};
  //v[2]={0.0,0.0};
  //v[3]={0.0,-1.0};
  //v[4]={0.3,-0.5};

  //Polygon aPolygon(v);
  //aPolygon.showMe();
  //std::cout<<"Area: "<<aPolygon.area()<<std::endl;
  //// A triangle built with the first 3 vertices
  //Triangle aTriangle(Vertices(v.begin(),v.begin()+3));
  //AbstractPolygon * p_ab=&aTriangle;

  //p_ab->showMe();
  //std::cout<<"Area: "<<aTriangle.area()<<std::endl;
  ////! Unit Square
  ////C++11 syntax. I expoit the fact that the is an implicit conversion
  ////between initializer lists with two double and a Point2d, since the latter has a  constructor taking two doubles
  ////as argument 
  //// In C++98 I would have written
  ////Square aSquare(Point2D(0.0,0.0),1.0); 
  //Square aSquare({0.0,0.0},1.0);
  //Square s2(aSquare); 
  //AbstractPolygon & r_ab=s2;
  //r_ab.showMe();
  //std::cout<<"Area: "<<r_ab.area()<<std::endl;
/*}*/
 
#include <iostream>
#include "Grid.h"
#include "Edge.h"

int main()
{
	// Read grid from file
	std::string filename = "mesh.dat";
	std::cout << ">> Read grid from " << filename << std::endl;
	Grid grid(filename);
	
	// Output grid to terminal
	std::cout << ">> Output grid" << std::endl;
	std::cout << grid << std::endl;
	
	// Compute grid total area
	std::cout << ">> Grid area: " << grid.area() << std::endl;	
	
	// Output edges to file
	std::cout << ">> Print edges to file" << std::endl;
	grid.outputEdges();

	return 0;
}
