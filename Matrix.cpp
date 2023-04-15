#include "Matrix.h"
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <exception>

using namespace std;

// sep; delimiter, comment; 주석 문자(#)
int DataFrame::ReadData(std::string FileName, char sep, char comment, bool IsHeader) {
    // constructor로 ifstream type의 'file' object 생성
    // ifstream은 stream class to read from files
    std::ifstream file(FileName);

    // 예외 처리, file open error check
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << FileName << endl;
        return -1; 
    }

    // IsHeader 값에 따라 switch 만들기
    if (IsHeader){switch_0 = true;}

    std::string line;
    while (std::getline(file, line)) { // reads a single line of 'file'(input stearm), 'line'이라는 string variable로 저장
        if (line.empty() || line[0] == comment) { // skip empty lines or comment lines ('#'이 나오면 csv 파일에 대한 설명인 경우가 대부분이므로 skip)
            continue; 
        }

        //istringstream class의 'iss' object 생성
        std::istringstream iss(line); //treat a string like a stream, and read values from it using 
                                      //the same input operations you would use on a std::istream
        std::string token;
        
        // IsHeader도 True이고, data도 비어있다면, stirng type을 원소로 갖는 colName vector에 저장
        // data는 double type을 원소로 갖는 vector를 원소로 갖는 vector, data.empty()에서 True가 retrun된다는 것은 비어있다는 것.
        if (IsHeader && data.empty() && (switch_0 == true)){
            while(std::getline(iss, token, sep)){
                colName.push_back(token);      
            }
            // column의 갯수는 colName의 갯수와 같다.
            // 왜냐하면 column data를 ',' 단위로 나눴기 때문에(token) 열의 갯수가 된다.
            numCol = colName.size(); 
            // colName 저장하고 나서 switch 끄기. 첫 번째 line만 header이기 때문에 switch를 off해주면 다음부터는 이 조건에 안 걸린다.                   
            switch_0 = false;        
        }
        else{
            vector<double> row; 
            while (std::getline(iss, token, sep)) { // 'sep'에 오는 argument로 line delimiter를 사용한다. default는 '\n'
                // Throw an exception when input string is not a valid
                try {
                    double value = std::stod(token); // converts a string to an double type value
                    row.push_back(value);            // add the double to a vector
                }
                catch (std::invalid_argument const& e) { // exception: a function receives an argument that is not valid
                    std::cerr << "Invalid argument: " << e.what() << endl;
                    return -1;
                }
            }
            // IsHeader가 false이었기 때문에 numCol 값이 아직 지정이 안 되어있다면, 즉 아직도 default 값인 0이라면
            if (!numCol){numCol=row.size();}

            //예외처리, row의 원소의 갯수가 열의 갯수와 안 맞는 row라면 오류를 발생시킨다.
            if ((row.size() != numCol) && (row.size() != 0)) { 
                std::cerr<<"Input row's size: " << row.size() << std::endl;
                std::cerr << "Row has wrong number of columns" << std::endl;
                return -1;
            }
            data.push_back(row); // row 1개가 1줄의 값을 double type으로 저장해둔 vector이기 때문에 data는 원소 1개가 1행을 의미하게 된다.
        }
        
        // header line을 넣었을 때는 data.empty()가 true로 나오고, header는 행의 갯수로 세면 안 되기 때문이다.
        // data에 1줄씩 넣어줄 때마다 행의 갯수는 증가시킨다.
        if(!data.empty()){numRow++;}
    }
    return 0;
}

Matrix DataFrame::GetMatrix(int index[], int nColumn){
    Matrix mat;
    mat.numRow = numRow;
    mat.numCol = nColumn;
    for (int i = 0; i < numRow; i++){
        vector<double> row;
        for (int j = 0; j < nColumn; j++){
            row.push_back(data[i][index[j]]);
        }
        mat.matrix.push_back(row);
    }
    return mat;
}

// Segmentation Fault를 해결하고 다른 연산을 하기 편한 방법은,  원소가 0이고 MXN인 Matrix를 만드는 것.
Matrix Matrix::zero_maker(const Matrix& X){
    Matrix zero_matrix;
    vector<vector<double>> temp(X.numRow, vector<double>(X.numCol, 0.0));
    zero_matrix.matrix = temp;
    return zero_matrix;
}

Matrix Matrix::operator+(Matrix const& other_matrix) {
    //예외 처리
    if ((this->numRow != other_matrix.numRow) || (this->numCol != other_matrix.numCol)) {
        throw std::invalid_argument("Matrix dimensions must match.");
    }

    Matrix result;
    result = zero_maker(other_matrix);
    result.numRow = this->numRow;
    result.numCol = this->numCol;
    for (int i = 0; i < numRow; i++) {
        for (int j = 0; j < numCol; j++) {
            result.matrix[i][j] = matrix[i][j] + other_matrix.matrix[i][j];
        }
    }
    return result;
}

Matrix Matrix::operator-(Matrix const& other_matrix) {
    //예외 처리
    if ((this->numRow != other_matrix.numRow) || (this->numCol != other_matrix.numCol)) {
        throw std::invalid_argument("Matrix dimensions must match.");
    }    

    Matrix result;
    result = zero_maker(*this);
    result.numRow = this->numRow;
    result.numCol = this->numCol;
    for (int i = 0; i < numRow; i++) {
        for (int j = 0; j < numCol; j++) {
            result.matrix[i][j] = matrix[i][j] - other_matrix.matrix[i][j];
        }
    }
    return result;
}

Matrix Matrix::operator*(Matrix const& other_matrix) {
    if (numCol != other_matrix.numRow) {
        throw std::invalid_argument("Matrix dimensions must match.");
    }

    Matrix zero_matrix;
    vector<vector<double>> temp(numRow, vector<double>(other_matrix.numCol, 0.0));
    zero_matrix.matrix = temp;
    zero_matrix.numRow = this->numRow;
    zero_matrix.numCol = other_matrix.numCol;

    for (int i = 0; i < numRow; i++) {
        for (int j = 0; j < other_matrix.numCol; j++) {
            double sum = 0;

            for (int k = 0; k < numCol; k++) {
                sum += matrix[i][k] * other_matrix.matrix[k][j];
            }

            zero_matrix.matrix[i][j] = sum;
        }
    }
    return zero_matrix;
}

Matrix Matrix::Transpose() const
{
    Matrix result;
    vector<vector<double>> temp(this->numCol, vector<double>(this->numRow, 0.0));
    result.matrix = temp;
    result.numRow = this->numCol;
    result.numCol = this->numRow;
    
    for (int i = 0; i < numRow; ++i)
    {
        for (int j = 0; j < numCol; ++j)
        {
            result.matrix[j][i] = matrix[i][j];
        }
    }
    return result;
}

Matrix Matrix::GetSubVectorbyColumn(int column) const
{
    //예외 처리
    if (this->numCol < column) {
        throw std::invalid_argument("Please choose a number smaller than the column number.");
    }

    Matrix result;
    vector<vector<double>> temp(this->numRow, vector<double>(1, 0.0));
    result.matrix = temp;
    result.numRow = this->numRow;
    result.numCol = 1;

    for (int i = 0; i < this->numRow; ++i)
    {
        result.matrix[i][0] = this->matrix[i][column];
    }
    return result;
}

void Matrix::Print(){
   for (const auto& row : matrix) {
        for (const auto& element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
}

int Matrix::GetNumRow(){return numRow;}

int Matrix::GetNumColumn(){return numCol;}

double Matrix::GetVal(int x, int y){
    // 예외 처리
    if (x < 0 || x >= this->numRow || y < 0 || y >= this->numCol) {
        throw std::out_of_range("Invalid indices");
    }   
    return this->matrix[x][y];
}

void Matrix::Set(int x, int y, double val){
    this->matrix[x][y] = val;
}

std::vector<double> Matrix::nx1_vector_converter() const{
    // 예외 처리, nX1 matrix만 받아야 한다.
    if (this->numCol != 1){
        cout<<"Please Input nX1 matrix."<<endl;
    }
    std::vector<double> result;
    for (int i=0; i < this->numRow; i++){
        result.push_back(this->matrix[i][0]);
    }
    return result;
}

void Matrix::SetnumRow(int x){
    this->numRow = x;
}

void Matrix::SetnumCol(int x)
{
    this->numCol = x;
}
void Matrix::Setmatrix(const std::vector<std::vector<double>>& matrix_){
    this->matrix = matrix_;
}

double Mean(const vector<double>& x){
    double sum = 0.0;
    for (double val : x){
        sum += val;
    }
    return sum / x.size();
}

double StdDev(const vector<double>& x){
    double mean = Mean(x);
    double sumSqDiff = 0.0;
    for (double val : x) {
        double diff = val - mean;
        sumSqDiff += diff * diff;
    }
    double variance = sumSqDiff / x.size();
    return std::sqrt(variance);
}

int Sign(double x) {
    if (x > 0) {return 1;} 
    else if (x < 0) {return -1;} 
    else {return 0;}
}

double PearsonCorrelation(const vector<double>& x, const vector<double>& y) {
    double meanX = Mean(x);
    double meanY = Mean(y);
    double stdX = StdDev(x);
    double stdY = StdDev(y);
    // scalar-vector substraction이 안 돼서 for 문으로 구현.. 
    // double covariance = (x - meanX) * (y - meanY);이 왜 안 될까..
    double covariance = 0.0;
    for (int i = 0; i < x.size(); i++) {
        covariance += (x[i] - meanX) * (y[i] - meanY);
    }
    double corrCoeff = covariance / (stdX * stdY);
    return corrCoeff;
}

double KendallTau(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    // binomial coefficient for the number of ways to choose two items from n items
    int numPairs = n * (n - 1) / 2;
    int numConcordant = 0;
    int numDiscordant = 0;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            int xSign = Sign(x[j] - x[i]);
            int ySign = Sign(y[j] - y[i]);
            if (xSign == ySign) {
                numConcordant++;
            } else {
                numDiscordant++;
            }
        }
    }
    double tau = (double)(numConcordant - numDiscordant) / (double)numPairs;
    return tau;
}

Matrix Cor(Matrix &mat, int method = 1) {
    // 결과값으로 출력할 m X m matrix 만들기. 원소는 대각은 1, 나머지는 0으로 초기화
    int n = mat.GetNumRow();
    int m = mat.GetNumColumn();
    Matrix corr;
        vector<vector<double>> temp(m, vector<double>(m, 0.0));
        for (int k=0; k < m; k++){
            temp[k][k] = 1;
        }
        corr.Setmatrix(temp);
        corr.SetnumCol(m);
        corr.SetnumRow(m);

    for (int i = 0; i < m; i++) {
        for (int j = i + 1; j < m; j++) {
            double corrCoeff = 0.0;
            // Pearson correlation
            if (method == 1) {
                corrCoeff = PearsonCorrelation(mat.GetSubVectorbyColumn(i).nx1_vector_converter(), 
                                            mat.GetSubVectorbyColumn(j).nx1_vector_converter());
            }
            // Kendall's tau 
            else if (method == 2) {
                corrCoeff = KendallTau(mat.GetSubVectorbyColumn(i).nx1_vector_converter(), 
                                    mat.GetSubVectorbyColumn(j).nx1_vector_converter());
            }
            // correlation matrix는 symmetric이기 때문에 이렇게 만들어줘야 한다.
            corr.Set(i, j, corrCoeff);
            corr.Set(j, i, corrCoeff);
        }
    }
    return corr;
}

// Matrix Matrix::SimpleLinearRegression(Matrix &X, Matrix &Y) {
//     int n = X.GetNumRow();
//     Matrix ones(n, 1, 1.0);
//     Matrix X_with_ones = ones.AppendMatrix(X);
//     Matrix X_transpose = X_with_ones.Transpose();
//     Matrix XtX = X_transpose.Multiply(X_with_ones);
//     Matrix XtY = X_transpose.Multiply(Y);
//     Matrix beta = XtX.Inverse().Multiply(XtY);
//     return beta;
// }