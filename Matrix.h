#pragma once
#include <iostream>
#include <vector>

class Matrix{
    friend class DataFrame;
private:
    std::vector<std::vector<double>> matrix;
    int numRow = 0; // 행
    int numCol = 0; // 열

public:    
    Matrix operator+(Matrix const& other_matrix);
    Matrix operator-(Matrix const& other_matrix);
    Matrix operator*(Matrix const& other_matrix);
    
    Matrix Transpose() const;
    Matrix GetSubVectorbyColumn(int column) const;

    Matrix zero_maker(const Matrix& X);

    void Print();
    int GetNumRow();
    int GetNumColumn();
    double GetVal(int x, int y);

    void Set(int i, int j, double val);
    std::vector<double> nx1_vector_converter() const;
    void SetnumRow(int x);
    void SetnumCol(int x);
    void Setmatrix(const std::vector<std::vector<double>>& matrix_);
};

int Sign(double x);
double PearsonCorrelation(const std::vector<double>& x, const std::vector<double>& y);
double KendallTau(const std::vector<double>& x, const std::vector<double>& y);

Matrix Cor(Matrix& mat, int method);
//Matrix SimpleLinearRegression(const Matrix& X, const Matrix& Y);

class DataFrame : public Matrix{
private:
    std::vector<std::vector<double>> data;
    std::vector<std::string> colName;
    bool switch_0 = false;

public:
    int ReadData(std::string FileName, char sep, char comment, bool IsHeader);
    Matrix GetMatrix(int index[], int nColumn);
};