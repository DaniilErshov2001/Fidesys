#include "ProblemSetup.h"
#include <fstream>
#include <iostream>

ProblemSetup::ProblemSetup(const std::string& inputFile)
{
    ReadInputData(inputFile);
}




void ProblemSetup::ReadInputData(const std::string& inputFile)
{
    std::ifstream file(inputFile);
    if (!file) {
        std::cerr << "Ошибка открытия файла!\n";
        return;
    }

    std::string line;
    std::regex nodeRegex(R"(^\s*(\d+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+))");
    std::regex elementRegex(R"((\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+))");

    std::unordered_map<int, int> nodeIdMap; // Соответствие внешнего ID узла внутреннему индексу

    std::smatch match;
    bool inNodeSection = false, inElementSection = false;
    int internalNodeIndex = 0; // Внутренний счетчик узлов

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '$') continue; // Пропускаем пустые строки и комментарии

        if (line.find("*NODE") != std::string::npos) {
            inNodeSection = true;
            inElementSection = false;
            continue;
        }
        if (line.find("*ELEMENT_SHELL") != std::string::npos) {
            inElementSection = true;
            inNodeSection = false;
            continue;
        }

        // Парсинг узлов
        if (inNodeSection && std::regex_search(line, match, nodeRegex)) {
            int originalId = std::stoi(match[1].str()); // Исходный ID узла
            double x = std::stod(match[2].str());
            double y = std::stod(match[3].str());

            // Сохраняем соответствие исходного ID и внутреннего индекса
            nodeIdMap[originalId] = internalNodeIndex;

            // Добавляем координаты в массивы
            nodesX.conservativeResize(internalNodeIndex + 1);
            nodesY.conservativeResize(internalNodeIndex + 1);
            nodesX[internalNodeIndex] = x;
            nodesY[internalNodeIndex] = y;

            internalNodeIndex++; // Увеличиваем индекс
        }
        // Парсинг элементов
        else if (inElementSection && std::regex_search(line, match, elementRegex)) {
            Element elem;
            int id = std::stoi(match[1].str());
            int pid = std::stoi(match[2].str());

            // Преобразуем исходные ID узлов элемента в внутренние индексы
            elem.nodesIds[0] = nodeIdMap[std::stoi(match[3].str())];
            elem.nodesIds[1] = nodeIdMap[std::stoi(match[4].str())];
            elem.nodesIds[2] = nodeIdMap[std::stoi(match[5].str())];

            // Добавляем элемент в массив
            elements.push_back(elem);
        }
    }
    file.close();

    // Вывод узлов
    std::cout << "!!!NODES\n";
    for (int i = 0; i < nodesX.size(); i++) {
        std::cout << i << " " << nodesX[i] << " " << nodesY[i] << "\n";
    }

    // Вывод элементов
    std::cout << "\n\n\n!!!ELEMENTS\n";
    for (int i = 0; i < elements.size(); i++) {
        std::cout << i << " " << elements[i].nodesIds[0] << " "
            << elements[i].nodesIds[1] << " " << elements[i].nodesIds[2] << "\n";
    }

    // Создание ограничений
    for (int i = 0; i < nodesX.size(); i++)
    {
        if (nodesX[i] == 0.) {
            constraints.push_back({ i, static_cast<Constraint::Type>(1) });
        }
        if (nodesY[i] == 0.) {
            constraints.push_back({ i, static_cast<Constraint::Type>(2) });
        }
    }

    nodesCount = nodesX.size();
    loads.resize(2 * nodesCount);
    loads.setZero();

    float minY = nodesY.minCoeff();
    float maxY = nodesY.maxCoeff();
    float minX = nodesX.minCoeff();
    float maxX = nodesX.maxCoeff();

    for (int i = 0; i < nodesCount; ++i)
    {
        if (nodesX[i] == maxX) {
            if (nodesY[i] == minY || nodesY[i] == maxY) {
                loads[2 * i + 0] = F_load / 2.0f;
                loads[2 * i + 1] = 0;
            }
            else {
                loads[2 * i + 0] = F_load;
                loads[2 * i + 1] = 0;
            }
        }
    }
}

const Eigen::VectorXf& ProblemSetup::GetNodesX() const { return nodesX; }
const Eigen::VectorXf& ProblemSetup::GetNodesY() const { return nodesY; }
const std::vector<Element>& ProblemSetup::GetElements() const { return elements; }
const std::vector<Constraint>& ProblemSetup::GetConstraints() const { return constraints; }
const Eigen::VectorXf& ProblemSetup::GetLoads() const { return loads; }
int ProblemSetup::GetNodesCount() const { return nodesCount; }
float ProblemSetup::GetPoissonRatio() const { return poissonRatio; }
float ProblemSetup::GetYoungModulus() const { return youngModulus; }
