#include "ProblemSetup.h"
#include "FEMSolver.h"

#include <iostream>
float poissonRatio_global = 0.3f, youngModulus_global = 2e5f; //MPa
float F_load_global = 1.f; //MN/m

int main()
{
    std::string inputFile = "../../mesh/task_mesh_extrafine.k";
    //std::string inputFile = "C:/fydesis/fydesis/input.txt";
    std::string outputFile = "../../result/output.txt";
    std::string inputCheck = "../../result/input_check.txt";
    std::string sigmaXX = "../../result/sigmaXX.txt";

    ProblemSetup problemSetup(inputFile);
    FEMSolver solver(problemSetup);

    solver.Solve();
    solver.OutputResultsToFile(outputFile);
    solver.WriteMeshData(inputCheck);
    solver.OutputSigmaXXToFile(sigmaXX);
    return 0;
}
