#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits.h>
#include <algorithm>
#include <cstdlib>
using namespace std;

typedef vector<unsigned long long int> ROW;
typedef vector<ROW> MATRIX;



void fillMatrix(MATRIX& v, const unsigned long long int& rc)
{
    for (unsigned long long int i = 0; i < rc; i++)
    {
        for (unsigned long long int j = 0; j < rc; j++)
        {
            unsigned long long int temp = rand() % 25 + 1;
            v[i][j] = temp;
        }
    }
}

void getQuad(const MATRIX& IN,
             MATRIX& A,
             MATRIX& B,
             MATRIX& C,
             MATRIX& D)
{
    int  sz = IN.size(),
         mid = sz/2,
         count=0;
    for (auto i : IN)
    {
        ROW temp1(i.begin(),i.begin()+mid),
            temp2(i.begin()+mid,i.end());
       if (count < mid)
       {
           A.push_back(temp1);
           B.push_back(temp2);
           ++count;
       }
       else
       {
           C.push_back(temp1);
           D.push_back(temp2);
       }

    }

}

MATRIX matrixMult(const MATRIX& A, const MATRIX& B)
{
    MATRIX C=A;
    auto  M1 = (A[0][0] + A[1][1]) * (B[0][0] + B[1][1]),
          M2 = (A[1][0] + A[1][1]) * (B[0][0]),
          M3 = (A[0][0]) * (B[0][1] - B[1][1]),
          M4 = (A[1][1]) * (B[1][0] - B[0][0]),
          M5 = (A[0][0] + A[0][1]) * (B[1][1]),
          M6 = (A[1][0] - A[0][0]) * (B[0][0] + B[0][1]),
          M7 = (A[0][1] - A[1][1]) * (B[1][0] + B[1][1]);

    C[0][0] = M1 + M4 - M5 + M7;
    C[0][1] = M3 + M5;
    C[1][0] = M2 + M4;
    C[1][1] = M1 - M2 + M3 + M6;
    return C;
}

MATRIX addMatrix(const MATRIX& A, const MATRIX& B)
{
    MATRIX C=A;
    for (unsigned long long int i =0; i < A.size(); ++i)
    {
        for(unsigned long long int j = 0; j < B.size(); ++j)
            C[i][j] = A[i][j] + B[i][j];
    }
    return C;
}

MATRIX subMatrix(const MATRIX& A, const MATRIX& B)
{
    MATRIX C = A;
    for (unsigned long long int i =0; i < A.size(); ++i)
    {
        for(unsigned long long int j = 0; j < B.size(); ++j)
            C[i][j] = A[i][j] - B[i][j];
    }
    return C;
}

MATRIX getC11(const MATRIX& M1, const MATRIX& M4,
              const MATRIX& M5, const MATRIX& M7)
{
    MATRIX temp = addMatrix(M1,M4);
    temp = subMatrix(temp,M5);
    return addMatrix(temp,M7);
}

MATRIX getC12(const MATRIX& M3, const MATRIX& M5)
{
    return addMatrix(M3,M5);
}

MATRIX getC21(const MATRIX& M2, const MATRIX& M4)
{
    return addMatrix(M2,M4);
}

MATRIX getC22(const MATRIX& M1, const MATRIX& M2,
              const MATRIX& M3, const MATRIX& M6)
{
    MATRIX temp = subMatrix(M1,M2);
    temp = addMatrix(temp,M3);
    return addMatrix(temp,M6);
}

MATRIX constructC(const MATRIX& C11, const MATRIX& C12,
                  const MATRIX& C21, const MATRIX& C22)
{
    //ROW R;
    MATRIX C;
    int size = C11.size();
    for(unsigned long long int i=0; i < size*2; i++)
    {
        ROW R;
        if (i < size)
        {
            R = C11[i];
            R.insert(R.end(),C12[i].begin(),C12[i].end());
        }
        else
        {
            unsigned long long int idx = i - size;
            R = C21[idx];
            R.insert(R.end(),C22[idx].begin(),C22[idx].end());
        }
        C.push_back(R);
    }

    return C;
}

/*
 * strassen's algorithm, taking into account M1-M7,
 * has a recurrence relation of T(n) = 7*T(n/2) + n^2
 *
 * solving using master theorem:
 *
 * a = 7, b = 2, f(n) = n^2
 *
 * compare:
 *
 * n^(log2(7)) and n^2
 * log2(7) ~= 2.81, thus n^2.81 is larger than n^2, making
 * the runtime Theta(n^log2(7))
*/
MATRIX strassen(const MATRIX& A, const MATRIX& B)
{
    if (A.size() == 2 && A[0].size() == 2)
    {
       return matrixMult(A,B);
    }
    else
    {
        /*
            ____________
            | X11 | X12 |
            ------------
            | X21 | X22 |
            ------------

        */
        MATRIX A11,A12,A21,A22;
        getQuad(A,A11,A12,A21,A22);
        MATRIX B11,B12,B21,B22;
        getQuad(B,B11,B12,B21,B22);

        MATRIX M1 = strassen(addMatrix(A11,A22), addMatrix(B11,B22)),
               M2 = strassen(addMatrix(A21,A22), B11),
               M3 = strassen(A11, subMatrix(B12,B22)),
               M4 = strassen(A22, subMatrix(B21,B11)),
               M5 = strassen(addMatrix(A11,A12),B22),
               M6 = strassen(subMatrix(A21,A11),addMatrix(B11,B12)),
               M7 = strassen(subMatrix(A12,A22),addMatrix(B21,B22));
        MATRIX C11 = getC11(M1,M4,M5,M7),
               C12 = getC12(M3,M5),
               C21 = getC21(M2,M4),
               C22 = getC22(M1,M2,M3,M6);
        return constructC(C11,C12,C21,C22);
    }
}


// included as a means to compare, runs at O(n^3)
void naive(const MATRIX& v1, const MATRIX& v2)
{
    MATRIX res;

    for (unsigned long long int i  = 0; i < v1.size(); ++i)
    {
        ROW temp1;
        for (unsigned long long int j = 0; j < v2[0].size(); ++j)
        {
            unsigned long long int temp2=0;
            for (unsigned long long int k = 0; k < v1[0].size(); ++k)
            {
                temp2 += v1[i][k]*v2[k][j];
            }
            temp1.push_back(temp2);
        }
        res.push_back(temp1);
    }
    cout << ' ' << endl;
}

int main()
{
    int n = 32;
    MATRIX X(n,vector<unsigned long long int>(32)), Y(n,vector<unsigned long long int>(32));
    fillMatrix(X,n);
    fillMatrix(Y,n);
    strassen(X,Y);
    return 1;
}
