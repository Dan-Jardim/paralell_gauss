#include "linearSystem.h"

LinearSystem::LinearSystem(string fileName)
{
    string line;
    ifstream sysFile(fileName);
    getline(sysFile,line);
    this->size = stod(line);
    cout << "Size: "<< this->size << endl;
    for (int i = 0; i<this->size; i++)
    {
        this->matrix.push_back(vector<double>{});
        for (int j = 0; j < this->size; j++){
            getline(sysFile,line);
            this->matrix.at(i).push_back(stod(line));
        }
    }
    for (int i = 0; i < this->size; i++)
    {
        getline(sysFile,line);
        this->vecRight.push_back(stod(line));
        this->vecResult.push_back(0);
    }

    sysFile.close();
}

LinearSystem::LinearSystem(vector<vector<double>> matrix, vector<double> rightVector)
{
    this->matrix = matrix;
    this->vecRight = rightVector;
    this->size = matrix.size();
}

LinearSystem::~LinearSystem()
{

}

void LinearSystem::printSystem()
{
    for(int i = 0; i< (this->size) ; i++)
    {
        cout << "|";
        for(auto num : this->matrix.at(i))
        {
            cout << num << " ";
        }
        cout << "| |";
        cout << this->vecResult.at(i);
        cout << "| = | ";
        cout << this->vecRight.at(i) << " |" << endl;
    }
}