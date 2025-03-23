#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include "ProblemSetup.h"
#include <fstream>
class FEMSolver
{
public:
    FEMSolver(ProblemSetup& problemSetup);
    void Solve();
    void OutputResultsToFile(const std::string& outputFile);
    void WriteMeshData(const std::string& outputFile);
    void OutputSigmaXXToFile(const std::string& outputFile);
private:
    void CalculateD();
    void CalculateBC(const Eigen::VectorXf& nodesX, const Eigen::VectorXf& nodesY, const Element& element,
        Eigen::Matrix<float, 3, 6>& B, Eigen::Matrix3f& C);
    void ApplyConstraints(Eigen::SparseMatrix<float>& K, const std::vector<Constraint>& constraints);
    void SetConstraints(Eigen::SparseMatrix<float>::InnerIterator& it, int index);

    Eigen::Matrix3f D;
    Eigen::VectorXf nodesX, nodesY, loads, displacements;
    std::vector<Element> elements;
    std::vector<Constraint> constraints;
    int nodesCount;
    float poissonRatio, youngModulus;
};

#endif // FEM_SOLVER_H
