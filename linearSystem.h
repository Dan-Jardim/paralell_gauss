#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

class LinearSystem
{
private:
    vector<vector<double>> matrix;
    vector<double> vecRight;
    vector<double> vecResult;
    int size;

public:
    LinearSystem(string fileName);
    LinearSystem(vector<vector<double>> matrix, vector<double> rightVector);
    ~LinearSystem();

    void printSystem();
};