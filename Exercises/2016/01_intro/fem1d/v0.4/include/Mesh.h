#ifndef MESH_H
#define MESH_H

#include <vector>

class Mesh
{
    public:
        // Constructors
        Mesh();
        Mesh(double xMin, double xMax, unsigned long NNodes);
        Mesh(const Mesh& Mesh_);

        // Destructors
        virtual ~Mesh();

        // Other methods
        unsigned long size() const;
        double getStep() const;
        double getNode(unsigned long idx) const;
        std::vector<unsigned long> getElementNodes(unsigned long idx) const;

    private:
        unsigned long NNodes;
        double step;
        std::vector<double> nodes;
        std::vector<std::vector<unsigned long>> elements;
};

#endif // MESH_H
