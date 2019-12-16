#include "matrix.h"

Matrix::Matrix(unsigned int dim) : dim(dim)
{
    myMatrix = new int*[dim];
    for (unsigned int i = 0; i < dim; i++) {
        myMatrix[i] = new int[dim];
    }
    /*float matrix [dim_x] [dim_y];
    myMatrix = *matrix;*/
}

Matrix::~Matrix()
{
    for (unsigned int i = 0; i < dim; i++) {
        delete [] myMatrix[i];
    }
    delete [] myMatrix;
}

void Matrix::setTo(unsigned int x, unsigned int y, int value)
{
    myMatrix[x][y] = value;
}

int Matrix::getFrom(unsigned int x, unsigned int y)
{
    return myMatrix[x][y];
}

int **Matrix::getMyMatrix() const
{
    return myMatrix;
}

void Matrix::setMyMatrix(int **value)
{
    myMatrix = value;
}
