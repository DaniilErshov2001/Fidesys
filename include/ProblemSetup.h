#ifndef PROBLEM_SETUP_H
#define PROBLEM_SETUP_H

#include <Eigen/Dense>
#include <vector>
#include <string>
#include "Element.h"
#include "Constraint.h"
#include <regex>

extern float poissonRatio_global;
extern float youngModulus_global;
extern float F_load_global;

class ProblemSetup
{
public:
    ProblemSetup(const std::string& inputFile);

    const Eigen::VectorXf& GetNodesX() const;
    const Eigen::VectorXf& GetNodesY() const;
    const std::vector<Element>& GetElements() const;
    const std::vector<Constraint>& GetConstraints() const;
    const Eigen::VectorXf& GetLoads() const;
    int GetNodesCount() const;
    float GetPoissonRatio() const;
    float GetYoungModulus() const;

private:
    void ReadInputData(const std::string& inputFile);

    Eigen::VectorXf nodesX, nodesY, loads;
    std::vector<Element> elements;
    std::vector<Constraint> constraints;
    int nodesCount;
    float poissonRatio = poissonRatio_global, youngModulus = youngModulus_global;
    float F_load = F_load_global;
};

#endif // PROBLEM_SETUP_H
