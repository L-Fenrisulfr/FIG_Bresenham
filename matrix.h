#ifndef MATRIX_H
#define MATRIX_H

class Matrix
{
public:
    Matrix(unsigned int dim);
    ~Matrix();

    void setTo(unsigned int x, unsigned int y, int value);
    int getFrom(unsigned int x, unsigned int y);

    int **getMyMatrix() const;
    void setMyMatrix(int **value);

private:
    unsigned int dim;
    int** myMatrix; // ids
};

#endif // MATRIX_H
