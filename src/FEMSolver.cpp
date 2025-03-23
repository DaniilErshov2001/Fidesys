#include "FEMSolver.h"
#include <unordered_set>
#include <Eigen/SparseLU>

FEMSolver::FEMSolver(ProblemSetup& problemSetup)
    : nodesX(problemSetup.GetNodesX()),
    nodesY(problemSetup.GetNodesY()),
    elements(problemSetup.GetElements()),
    constraints(problemSetup.GetConstraints()),
    loads(problemSetup.GetLoads()),
    nodesCount(problemSetup.GetNodesCount()),
    poissonRatio(problemSetup.GetPoissonRatio()),
    youngModulus(problemSetup.GetYoungModulus())
{
}

void FEMSolver::Solve()
{
    // Расчет матрицы D (матрица упругости)
    CalculateD();

    // Массив для хранения всех ненулевых элементов глобальной матрицы жесткости
    std::vector<Eigen::Triplet<float>> triplets;

    for (auto& element : elements)
    {
        Eigen::Matrix<float, 3, 6> B;
        Eigen::Matrix3f C;
        CalculateBC(nodesX, nodesY, element, B, C);

        // Расчет матрицы жесткости для элемента
        Eigen::Matrix<float, 6, 6> K = B.transpose() * D * B * std::abs(C.determinant()) / 2.0;

        // Заполнение глобальной матрицы жесткости
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                triplets.push_back(Eigen::Triplet<float>(2 * element.nodesIds[i] + 0, 2 * element.nodesIds[j] + 0, K(2 * i + 0, 2 * j + 0)));
                triplets.push_back(Eigen::Triplet<float>(2 * element.nodesIds[i] + 0, 2 * element.nodesIds[j] + 1, K(2 * i + 0, 2 * j + 1)));
                triplets.push_back(Eigen::Triplet<float>(2 * element.nodesIds[i] + 1, 2 * element.nodesIds[j] + 0, K(2 * i + 1, 2 * j + 0)));
                triplets.push_back(Eigen::Triplet<float>(2 * element.nodesIds[i] + 1, 2 * element.nodesIds[j] + 1, K(2 * i + 1, 2 * j + 1)));
            }
        }
    }

    // Создание глобальной матрицы жесткости
    Eigen::SparseMatrix<float> globalK(2 * nodesCount, 2 * nodesCount);
    globalK.setFromTriplets(triplets.begin(), triplets.end());

    // Применение граничных условий
    ApplyConstraints(globalK, constraints);

    // Решение системы линейных уравнений
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver(globalK);
    Eigen::VectorXf displacements = solver.solve(loads);

    // Сохраняем результат в переменной displacements
    this->displacements = displacements;
}

void FEMSolver::OutputResultsToFile(const std::string& outputFile)
{
    std::ofstream outfile(outputFile);
    if (outfile.is_open())
    {
        // Записываем в файл результат перемещения узлов
        outfile << displacements << std::endl;
    }

    // Записываем напряжения для каждого элемента
    for (auto& element : elements)
    {
        Eigen::Matrix<float, 3, 6> B;
        Eigen::Matrix3f C;
        CalculateBC(nodesX, nodesY, element, B, C);

        // Расчет деформаций в узлах элемента
        Eigen::Matrix<float, 6, 1> delta;
        delta << displacements.segment<2>(2 * element.nodesIds[0]),
            displacements.segment<2>(2 * element.nodesIds[1]),
            displacements.segment<2>(2 * element.nodesIds[2]);

        // Расчет напряжений
        Eigen::Vector3f sigma = D * B * delta;
        float sigma_mises = sqrt(sigma[0] * sigma[0] - sigma[0] * sigma[1] + sigma[1] * sigma[1] + 3.0f * sigma[2] * sigma[2]);

        // Запись напряжений в файл
        outfile << sigma_mises << std::endl;
    }

}

void FEMSolver::OutputSigmaXXToFile(const std::string& outputFile)
{
    std::ofstream outfile(outputFile);
    if (!outfile.is_open())
    {
        return;
    }

    // Записываем напряжения sigma_xx для каждого элемента
    for (const auto& element : elements)
    {
        Eigen::Matrix<float, 3, 6> B;
        Eigen::Matrix3f C;
        CalculateBC(nodesX, nodesY, element, B, C);

        // Расчет деформаций в узлах элемента
        Eigen::Matrix<float, 6, 1> delta;
        delta << displacements.segment<2>(2 * element.nodesIds[0]),
            displacements.segment<2>(2 * element.nodesIds[1]),
            displacements.segment<2>(2 * element.nodesIds[2]);

        // Расчет напряжений
        Eigen::Vector3f sigma = D * B * delta;
        float sigma_xx = sigma[0]; // Выбираем первую компоненту напряжений (σxx)

        // Запись sigma_xx в файл
        outfile << sigma_xx << std::endl;
    }
}


void FEMSolver::WriteMeshData(const std::string& outputFile)
{
    std::ofstream outfile(outputFile);

    // Запись коэффициентов
    outfile << poissonRatio << " " << youngModulus << "\n";

    // Запись узлов
    outfile << nodesCount << "\n";
    for (int i = 0; i < nodesCount; ++i)
    {
        outfile << nodesX[i] << " " << nodesY[i] << "\n";
    }

    // Запись элементов
    outfile << elements.size() << "\n";
    for (int i = 0; i < elements.size(); ++i)
    {
        outfile << elements[i].nodesIds[0] << " " << elements[i].nodesIds[1] << " " << elements[i].nodesIds[2] << "\n";    }

    // Запись ограничений
    outfile << constraints.size() << "\n";
    for (int i = 0; i < constraints.size(); ++i)
    {
        outfile << constraints[i].node << " " << constraints[i].type << "\n";    }

    // Запись нагрузок
    int loadsCount = 0;
    for (int i = 0; i < nodesCount; ++i)
    {
        if (loads[2 * i + 0] != 0 || loads[2 * i + 1] != 0) {
            loadsCount++;
        }
    }
    outfile << loadsCount << "\n";

    for (int i = 0; i < nodesCount; ++i)
    {
        if (loads[2 * i + 0] != 0 || loads[2 * i + 1] != 0) {
            outfile << i << " "  << loads[2 * i + 0] << " " << loads[2 * i + 1] << "\n";
        }
    }
}


void FEMSolver::CalculateD()
{
    D << 1.0f, poissonRatio, 0.0f,
        poissonRatio, 1.0f, 0.0f,
        0.0f, 0.0f, (1.0f - poissonRatio) / 2.0f;
    D *= youngModulus / (1.0f - std::pow(poissonRatio, 2.0f));
}

void FEMSolver::CalculateBC(const Eigen::VectorXf& nodesX, const Eigen::VectorXf& nodesY, const Element& element,
    Eigen::Matrix<float, 3, 6>& B, Eigen::Matrix3f& C)
{
    Eigen::Vector3f x, y;
    x << nodesX[element.nodesIds[0]], nodesX[element.nodesIds[1]], nodesX[element.nodesIds[2]];
    y << nodesY[element.nodesIds[0]], nodesY[element.nodesIds[1]], nodesY[element.nodesIds[2]];

    // Заполняем матрицу C
    C << Eigen::Vector3f(1.0f, 1.0f, 1.0f), x, y;
    Eigen::Matrix3f IC = C.inverse();

    // Заполняем матрицу B
    for (int i = 0; i < 3; ++i)
    {
        B(0, 2 * i + 0) = IC(1, i);
        B(0, 2 * i + 1) = 0.0f;
        B(1, 2 * i + 0) = 0.0f;
        B(1, 2 * i + 1) = IC(2, i);
        B(2, 2 * i + 0) = IC(2, i);
        B(2, 2 * i + 1) = IC(1, i);
    }
}

void FEMSolver::ApplyConstraints(Eigen::SparseMatrix<float>& K, const std::vector<Constraint>& constraints)
{
    std::unordered_set<int> constraintIndices;

    // Добавляем индексы, которые нужно зафиксировать
    for (const auto& constraint : constraints)
    {
        if (constraint.type == Constraint::UX || constraint.type == Constraint::UXUY)
        {
            constraintIndices.insert(2 * constraint.node);
        }
        if (constraint.type == Constraint::UY || constraint.type == Constraint::UXUY)
        {
            constraintIndices.insert(2 * constraint.node + 1);
        }

    }

    // Применяем ограничения к глобальной матрице жесткости
    for (int k = 0; k < K.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<float>::InnerIterator it(K, k); it; ++it)
        {
            if (constraintIndices.count(it.row()) || constraintIndices.count(it.col()))
            {
                SetConstraints(it, it.row());
                SetConstraints(it, it.col());
            }
        }
    }
}

void FEMSolver::SetConstraints(Eigen::SparseMatrix<float>::InnerIterator& it, int index)
{
    if (it.row() == index || it.col() == index)
    {
        it.valueRef() = (it.row() == it.col()) ? 1.0f : 0.0f;
    }
}
