#include "Mesh.h"

Mesh::Mesh()
{
    //ctor
}

Mesh::Mesh(double xMin, double xMax, unsigned long NNodes_)
    : NNodes(NNodes_),
      step((xMax - xMin) / (NNodes_-1)),
      nodes(NNodes_),
      elements(NNodes_-1, std::vector<unsigned long>(2))
{
    // Initialize nodes
    for (unsigned long i = 0; i < NNodes; i++)
    {
        nodes[i] = xMin + i * step;
    }

    // Initialize element -> node mapping
    for (unsigned long i = 0; i < NNodes-1; i++)
    {
        elements[i][0] = i;
        elements[i][1] = i+1;
    }
}


Mesh::Mesh(const Mesh& Mesh_)
    : NNodes(Mesh_.NNodes),
      step(Mesh_.step),
      nodes(Mesh_.nodes),
      elements(Mesh_.elements)
{
    //
}

Mesh::~Mesh()
{
    //dtor
}

unsigned long Mesh::size() const
{
    return NNodes;
}

double Mesh::getStep() const
{
    return step;
}

double Mesh::getNode(unsigned long idx) const
{
    return nodes[idx];
}

std::vector<unsigned long> Mesh::getElementNodes(unsigned long idx) const
{
    return elements[idx];
}
