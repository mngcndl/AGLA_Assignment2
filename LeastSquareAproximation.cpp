// Liubov Smirnova, CS-05
#include <iostream>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <iomanip>
using namespace std;

class Matrix{
public:
    int rows;
    int columns;
    vector<vector<double>> matrix;
    explicit Matrix(int rows, int columns, vector<vector<double>> arrayOfInput) {
        this->rows = rows;
        this->columns = columns;
        matrix.resize(rows, vector<double>(columns));
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrix[i][j] = arrayOfInput[i][j];
            }
        }
    }
    explicit Matrix(int rows, int columns) {
        this->rows = rows;
        this->columns = columns;
        matrix.resize(rows, vector<double>(columns));
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrix[i][j] = 0;
            }
        }
    }

    Matrix() = default;
    friend istream &operator>>(istream &in, Matrix &input) {
        vector<vector<double>> matrixElements;
        for (int j = 0; j < input.rows; j++) {
            vector<double> matrixRow;
            int counter = 0;
            for (int k = 0; k < input.columns; k++){
                int value;
                in >> value;
                matrixRow.push_back(value);
                counter++;
            }
            matrixElements.push_back(matrixRow);
        }
        input.matrix = matrixElements;
        return in;
    }
    friend ostream &operator<<(ostream &out, Matrix &output) {
        for(int i = 0; i < output.rows; ++i) {
            for (int j = 0; j < output.columns - 1; ++j) {
                if(output.matrix[i][j] > -0.0000000001 && output.matrix[i][j] <= 0) {
                    output.matrix[i][j] = 0;
                }
                cout << fixed << setprecision(4) << output.matrix[i][j] << " ";
            }
            if(output.matrix[i][output.columns-1] <= 0 && output.matrix[i][output.columns-1] > -0.0000000001) {
                output.matrix[i][output.columns-1] = 0;}
            cout << fixed << setprecision(4) << output.matrix[i][output.columns-1] << endl;
        }
        return out;
    }

    void matrixToString(){
        for(int i = 0; i < this->rows; ++i) {
            for (int j = 0; j < this->columns - 1; ++j) {
                if(this->matrix[i][j] > -0.0000000001 && this->matrix[i][j] <= 0) {
                    this->matrix[i][j] = 0;
                }
                cout << fixed << setprecision(4) << this->matrix[i][j] << " ";
            }
            if(this->matrix[i][this->columns-1] <= 0 && this->matrix[i][this->columns-1] > -0.0000000001) {
                this->matrix[i][this->columns-1] = 0;}
            cout << fixed << setprecision(4) << this->matrix[i][this->columns-1] << endl;
        }
    }

    Matrix operator+(Matrix& arr) {
        Matrix resultingMatrix (arr.rows, arr.columns);
        if (arr.rows != rows || arr.columns != columns) {
            cout << "Error: the dimensional problem occurred\n";
            return resultingMatrix;
        }
        for (int i = 0; i < arr.rows; i++) {
            for (int j = 0; j < arr.columns; j++) {
                resultingMatrix.matrix[i][j] = arr.matrix[i][j] + this->matrix[i][j];
            }
        }
        return resultingMatrix;
    }

    Matrix(vector<vector<double>>newM){
        matrix = newM;
        rows = newM.size();
        columns = newM[0].size();
    }

    Matrix operator-(Matrix& arr) {
        Matrix resultingMatrix (arr.rows, arr.columns);
        if (arr.rows != rows || arr.columns != columns) {
            cout << "Error: the dimensional problem occurred\n";
            return resultingMatrix;
        }
        for (int i = 0; i < arr.rows; i++) {
            for (int j = 0; j < arr.columns; j++) {
                resultingMatrix.matrix[i][j] = this->matrix[i][j] - arr.matrix[i][j];
            }
        }
        return resultingMatrix;
    }

    Matrix operator*(Matrix& arr) {
        Matrix resultingMatrix = *new Matrix(rows, arr.columns);
        if (columns != arr.rows) {
            cout<< "Error: the dimensional problem occurred\n";
            return resultingMatrix;
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < arr.columns; j++) {
                for(int k = 0; k < columns; k++){
                    resultingMatrix.matrix[i][j] += this->matrix[i][k] * arr.matrix[k][j];
                }
            }
        }
        return resultingMatrix;
    }

    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            rows = other.rows;
            columns = other.columns;
            matrix = other.matrix;
        }
        return *this;
    }
    Matrix transposition() {
        int newRowNumber = this->columns;
        int newColumnNumber = this->rows;
        Matrix resultingMatrix(newRowNumber, newColumnNumber);
        for (int i = 0; i < newRowNumber; i++) {
            for (int j = 0; j < newColumnNumber; j++) {
                resultingMatrix.matrix[i][j] = this->matrix[j][i];
            }
        }
        return resultingMatrix;
    }
};

class SquareMatrix : public Matrix {
public:
    SquareMatrix(int rows, vector<vector<double>>& squareMatrix) : Matrix(rows, rows, squareMatrix) {}
    explicit SquareMatrix(int rows) : Matrix(rows, rows) {}
    SquareMatrix() = default;
    SquareMatrix(vector<vector<double>> m){
        matrix = m;
        rows = m.size();
        columns = m[0].size();
    }

    friend istream &operator>>(istream &in, SquareMatrix input) {
        vector<vector<double>> matrixElements;
        for (int j = 0; j < input.rows; j++) {
            vector<double> matrixRow;
            int counter = 0;
            for (int k = 0; k < input.rows; k++){
                int value;
                in >> value;
                matrixRow.push_back(value);
                counter++;
            }
            if (counter == input.rows) {
                matrixElements.push_back(matrixRow);
            }
        }
        input.matrix = matrixElements;
        return in;
    }
    friend ostream &operator<<(ostream &out, const SquareMatrix &output) {
        for (int j = 0; j < output.rows; j++) {
            for(int k = 0; k < output.rows; k++) {
                out << output.matrix[j][k] << " ";
            }
            out << endl;
        }
        return out;
    }
    void squareMatrixToString(){
        for(int i = 0; i < this->rows; ++i) {
            for (int j = 0; j < this->rows - 1; ++j) {
                if(this->matrix[i][j] <= 0 && this->matrix[i][j] > -0.0000000001) {
                    this->matrix[i][j] = 0;}
                cout << fixed << setprecision(4) << this->matrix[i][j] << " ";
            }
            if(this->matrix[i][this->columns-1] <= 0 && this->matrix[i][this->rows-1] > -0.0000000001) {
                this->matrix[i][this->columns-1] = 0;
            }
            cout << fixed << setprecision(4) << this->matrix[i][this->rows-1] << endl;
        }
    }

    SquareMatrix operator+(Matrix* arr) {
        SquareMatrix resultingMatrix (arr->rows);
        if (arr->rows != arr->columns) {
            cout << "Error: the dimensional problem occurred\n";
            return resultingMatrix;
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                resultingMatrix.matrix[i][j] = arr->matrix[i][j] + this->matrix[i][j];
            }
        }
        return resultingMatrix;
    }
    SquareMatrix operator-(Matrix* arr) {
        SquareMatrix resultingMatrix (arr->rows);
        if (arr->rows != arr->columns) {
            cout << "Error: the dimensional problem occurred\n";
            return resultingMatrix;
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                resultingMatrix.matrix[i][j] = this->matrix[i][j] - arr->matrix[i][j];
            }
        }
        return resultingMatrix;
    }
    SquareMatrix operator*(Matrix* arr) {
        SquareMatrix resultingMatrix = *new SquareMatrix(rows);
        if (arr->columns != arr->rows) {
            cout<< "Error: the dimensional problem occurred\n";
            return resultingMatrix;
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < arr->columns; j++) {
                for(int k = 0; k < rows; k++){
                    resultingMatrix.matrix[i][j] += this->matrix[i][k] * arr->matrix[k][j];
                }
            }
        }
        return resultingMatrix;
    }
    SquareMatrix operator*(SquareMatrix* arr) {
        SquareMatrix resultingMatrix = *new SquareMatrix(rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < arr->columns; j++) {
                for(int k = 0; k < rows; k++){
                    resultingMatrix.matrix[i][j] += this->matrix[i][k] * arr->matrix[k][j];
                }
            }
        }
        return resultingMatrix;
    }
    SquareMatrix& operator=(Matrix other) {
        if (this != &other) {
            rows = other.rows;
            matrix = other.matrix;
        }
        return *this;
    }
    SquareMatrix transposition() {
        SquareMatrix resultingMatrix(rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                resultingMatrix.matrix[i][j] = this->matrix[j][i];
            }
        }
        return resultingMatrix;
    }
};

class EliminationMatrix : public Matrix{
public:
    int rows;
    int columns;
    vector<vector<double>> matrix;
    EliminationMatrix(int rows, int columns, Matrix &givenMatrix, int eliminationLine, int column, int pivotLineNumber) {
        this->rows = rows;
        this->columns = columns;
        matrix.resize(rows, vector<double>(columns));
        for (int k = 0; k < rows; k++) {
            for (int m = 0; m < columns; m++) {
                if (k != m) {matrix[k][m] = 0;}
                else{matrix[k][m] = 1;}
            }
        }
        matrix[eliminationLine][column] = -1.0 * (givenMatrix.matrix[eliminationLine][column] / givenMatrix.matrix[column][column]);
    }

    Matrix operator+(Matrix& arr) {
        return this->Matrix::operator+(arr);
    }

    Matrix operator-(Matrix& arr) {
        return this->Matrix::operator-(arr);
    }

    Matrix operator*(Matrix& arr) {
        Matrix resultingMatrix = *new Matrix(rows, arr.columns);
        if (this->columns != arr.rows) {
            cout<< "Error: the dimensional problem occurred\n";
            return resultingMatrix;
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                for(int k = 0; k < rows; k++){
                    resultingMatrix.matrix[i][j] += this->matrix[i][k] * arr.matrix[k][j];
                }
            }
        }
        if (arr.columns == arr.rows * 2) {
            for (int i = 0; i < rows; i++) {
                for (int j = arr.rows; j < arr.columns; j++) {
                    for(int k = 0; k < rows; k++){
                        resultingMatrix.matrix[i][j] += this->matrix[i][k] * arr.matrix[k][j];
                    }
                }
            }
        }
        return resultingMatrix;
    }
    EliminationMatrix& operator=(Matrix other)
    {
        if (this != &other) {
            rows = other.rows;
            columns = other.columns;
            matrix = other.matrix;
        }
        return *this;
    }
};

class PermutationMatrix : public Matrix{
public:
    PermutationMatrix(int rows, int columns) {
        this->rows = rows;
        this->columns = columns;
        matrix.resize(this->rows, vector<double>(columns));
        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < this->columns; j++) {
                if (i == j) {
                    matrix[i][j] = 1;
                }
                else {matrix[i][j] = 0;}
            }
        }
    }

    void permutationMatrixToString(){
        this->Matrix::matrixToString();
    }

    void conductPermutationOfRowsIJ(int rowIIndex, int rowJIndex) {
        swap(matrix[rowIIndex], matrix[rowJIndex]);
    }

    Matrix operator+(Matrix& arr) {
        return this->Matrix::operator+(arr);
    }

    Matrix operator-(Matrix& arr) {
        return this->Matrix::operator-(arr);
    }

    Matrix operator*(Matrix& arr) {
        Matrix resultingMatrix = *new Matrix(rows, arr.columns);
        if (this->columns != arr.rows) {
            cout<< "Error: the dimensional problem occurred\n";
            return resultingMatrix;
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                for(int k = 0; k < rows; k++){
                    resultingMatrix.matrix[i][j] += this->matrix[i][k] * arr.matrix[k][j];
                }
            }
        }
        if (arr.columns == arr.rows * 2) {
            for (int i = 0; i < rows; i++) {
                for (int j = arr.rows; j < arr.columns; j++) {
                    for(int k = 0; k < rows; k++){
                        resultingMatrix.matrix[i][j] += this->matrix[i][k] * arr.matrix[k][j];
                    }
                }
            }
        }
        return resultingMatrix;
    }

    PermutationMatrix& operator=(Matrix other)
    {
        if (this != &other) {
            rows = other.rows;
            columns = other.columns;
            matrix = other.matrix;
        }
        return *this;

    }
};

class GaussElimination{
public:
    SquareMatrix sM;
    Matrix squarePlusAugmentedMatrix;

    explicit GaussElimination(SquareMatrix &squareMatrix){
        sM = squareMatrix;
        squarePlusAugmentedMatrix = Matrix(sM.rows, sM.rows * 2);
        squarePlusAugmentedMatrix.rows = squareMatrix.rows;
        squarePlusAugmentedMatrix.columns = squareMatrix.rows * 2;
        for(int i = 0; i < squareMatrix.rows; i++) {
            for (int j = 0; j < squareMatrix.rows; j++) {
                squarePlusAugmentedMatrix.matrix[i][j] = squareMatrix.matrix[i][j];
            }
        }
        for(int i = 0; i < sM.rows; i++) {
            for (int j = sM.rows; j < sM.rows * 2; j++) {
                if (i == j-sM.rows){
                    squarePlusAugmentedMatrix.matrix[i][j] = 1;
                } else {
                    squarePlusAugmentedMatrix.matrix[i][j] = 0;
                }
            }
        }
    }

    int lineIndexOfTheFirstMaximumAbsolutePivot(int pivotRowNumber, int pivotColumnNumber, Matrix& givenMatrix) {
        double maximalPivot = abs(givenMatrix.matrix[pivotRowNumber][pivotColumnNumber]);
        int maximalPivotRowNumber = pivotRowNumber;
        for (int i = pivotRowNumber + 1; i < sM.rows; i++) {
            if (abs(givenMatrix.matrix[i][pivotColumnNumber]) > maximalPivot){
                maximalPivot = abs(givenMatrix.matrix[i][pivotColumnNumber]);
                maximalPivotRowNumber = i;
            }
        }
        return maximalPivotRowNumber;
    }
    Matrix inverseMatrix(){
        int stepCounter = 1;
        for (int j = 0; j < this-> sM.rows; j++) {
            for(int i = 0; i < this->sM.rows; i++) {
                if (j == i) {
                    int maximalIthPivotLineNumber = lineIndexOfTheFirstMaximumAbsolutePivot(j, i, squarePlusAugmentedMatrix);
                    if (squarePlusAugmentedMatrix.matrix[maximalIthPivotLineNumber][i] != abs(squarePlusAugmentedMatrix.matrix[i][i])) {
                        PermutationMatrix pM(sM.rows, sM.rows);
                        pM.conductPermutationOfRowsIJ(i, maximalIthPivotLineNumber);
                        squarePlusAugmentedMatrix = pM * squarePlusAugmentedMatrix;
                        maximalIthPivotLineNumber = lineIndexOfTheFirstMaximumAbsolutePivot(j + 1, i, squarePlusAugmentedMatrix);
                    }
                }
                if (i > j && squarePlusAugmentedMatrix.matrix[i][j] != 0) {
                    EliminationMatrix eM = EliminationMatrix(sM.rows, sM.rows, squarePlusAugmentedMatrix, i, j, i);
                    squarePlusAugmentedMatrix = eM * squarePlusAugmentedMatrix;
                }

            }
        }
        for (int i =  this->sM.rows - 1; i > 0; i--) {
            int j = i;
            for (int k = i-1; k >= 0; k--) {
                EliminationMatrix eM = EliminationMatrix(sM.rows, sM.rows, squarePlusAugmentedMatrix, k, j, i);
                squarePlusAugmentedMatrix = eM * squarePlusAugmentedMatrix;
            }
        }
        for (int i = 0; i < sM.rows; i++) {
            double divider = squarePlusAugmentedMatrix.matrix[i][i];
            for (int j = i; j < squarePlusAugmentedMatrix.columns; j++) {
                squarePlusAugmentedMatrix.matrix[i][j] = squarePlusAugmentedMatrix.matrix[i][j] / divider;
            }
        }
        Matrix res(squarePlusAugmentedMatrix.rows, squarePlusAugmentedMatrix.rows);
        for (int i = 0; i < res.rows; i++) {
            for (int j = squarePlusAugmentedMatrix.rows; j < squarePlusAugmentedMatrix.columns; j++) {
                res.matrix[i][j - squarePlusAugmentedMatrix.rows] = squarePlusAugmentedMatrix.matrix[i][j];
            }
        }
        return res;
    }
};

class ColumnVector : public Matrix{
public:
    ColumnVector(int rows, Matrix &givenMatrix, int eliminationLine, int column, int pivotLineNumber) {
        this->rows = rows;
        this->columns = 1;
        matrix.resize(rows, vector<double>(columns));
        for (int k = 0; k < rows; k++) {
            this->matrix[k] = givenMatrix.matrix[k];
        }
    }
    ColumnVector(int rows){
        this->rows = rows;
        this->columns = 1;
        matrix.resize(rows, vector<double>(columns));
    }
    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            rows = other.rows;
            columns = 1;
            matrix = other.matrix;
        }
        return *this;
    }
};

#define GNUPLOT_NAME "gnuplot -persist"
int main() {
    FILE* pipe = popen(GNUPLOT_NAME, "w");
    int numberOfLines, dimensionality;
    cin >> numberOfLines;
    vector<vector<double>> matrix;
    matrix.resize(numberOfLines, vector<double>(2));
    Matrix inputMatrix(numberOfLines, 2, matrix);
    cin >> inputMatrix;
    cin >> dimensionality;
    Matrix A(inputMatrix.rows, dimensionality + 1);
    for (int i = 0; i < numberOfLines;i++) {
        for (int j = 0; j < dimensionality + 1; j++){
            A.matrix[i][j] = pow(inputMatrix.matrix[i][0], j);
        }
    }

    Matrix A_T = A.transposition();
    SquareMatrix A_TMul_A;
    A_TMul_A = A_T * A;
    GaussElimination getInverse(A_TMul_A);
    Matrix inverseMatrix = getInverse.inverseMatrix();
    ColumnVector b(numberOfLines);
    for (int i = 0; i < numberOfLines; i++) {
        b.matrix[i][0] = inputMatrix.matrix[i][1];
    }
    Matrix A_TMul_b;
    A_TMul_b = A_T * b;
    Matrix x;
    x = inverseMatrix * A_TMul_b;
    cout << "A:\n" << A;
    cout << "A_T*A:\n" << A_TMul_A;
    cout << "(A_T*A)^-1:\n" << inverseMatrix;
    cout << "A_T*b:\n" << A_TMul_b;
    cout << "x~:\n" << x;

    vector<double> xs;
    for (int i = 0; i<numberOfLines; i++) {
        xs.push_back(inputMatrix.matrix[i][0]);
    }

    vector<double> ys;
    for (int i = 0; i<numberOfLines; i++) {
        ys.push_back(inputMatrix.matrix[i][1]);
    }
    double lowX = xs[0];
    double bigX = xs[0];
    double lowY = ys[0];
    double bigY = ys[0];

    cout << lowX << "; " << bigX <<"; " << lowY << "; " << bigY << endl;
    for(int i = 1; i < xs.size(); i++)
    {
        if(xs[i] < lowX){lowX = xs[i];}
        if(xs[i] > bigX){bigX = xs[i];}
        if(ys[i] < lowY){lowY = ys[i];}
        if(ys[i] > bigY){bigY = ys[i];}
    }
    if (dimensionality == 3) {
        fprintf(pipe, "plot [%lf : %lf] [%lf : %lf] %lf*x**3 + %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n", lowX * 2, bigX * 2, lowY * 2, bigY * 2, x.matrix[3][0], x.matrix[2][0], x.matrix[1][0], x.matrix[0][0]);
        for (int i = 0; i < numberOfLines; i++) {
            fprintf(pipe, "%f\t%f\n", inputMatrix.matrix[i][0], inputMatrix.matrix[i][1]);
        }
    }
    fprintf(pipe, "e\n");
    fflush(pipe);
    pclose(pipe);
    return 0;
}
