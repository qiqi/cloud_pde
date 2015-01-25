#ifndef DIAMOND_H
#define DIAMOND_H

#include<cstring>
#include<cmath>
#include<cstdlib>
#include<cassert>
#include<mpi.h>

// ====== common interface for local schemes ======= //
//
template<size_t numVar>
class LocalInputs {
    public:
    virtual double val(int iVar) const = 0;
    virtual double nbr(int iVar, int iNbr) const = 0;
};

template<size_t numVar>
class LocalOutputs {
    public:
    virtual void val(int iVar, double newVal) = 0;
};

struct LocalMesh {
    double x, dx;
    LocalMesh(double x, double dx) : x(x), dx(dx) {}
    LocalMesh() : x(nan("")), dx(nan("")) {}
};

// ====== common interface of Discretization ======= //

class Discretization {
    public:
    template<size_t numVar>
    Discretization(int numGrids, double dx,
         void (&localOperator)(
               LocalOutputs<numVar>&, const LocalMesh&));

    template<size_t numInput, size_t numOutput>
    void applyOp(void (&localOperator)(
                 const LocalInputs<numInput>& inputs,
                 LocalOutputs<numOutput>& outputs,
                 const LocalMesh& mesh));
};

// ====== common utility for diamond_scheme and classic_scheme ====== //
//
const char* mpiRankString() {
    static char mpiRankStr[256] = "";
    if (mpiRankStr[0] == 0) {
        int iProc, numProc;
        MPI_Comm_size(MPI_COMM_WORLD, &numProc);
        MPI_Comm_rank(MPI_COMM_WORLD, &iProc);

        sprintf(mpiRankStr, "%d", numProc);
        int maxStrLen = strlen(mpiRankStr);

        sprintf(mpiRankStr, "%d", iProc);
        int strLen = strlen(mpiRankStr);
        int padLen = maxStrLen - strLen;
        memmove(mpiRankStr + padLen, mpiRankStr, strLen + 1);
        memset(mpiRankStr, '0', padLen);
    }
    return mpiRankStr;
}

template<size_t numVar>
class LocalInputs1D : public LocalInputs<numVar> {
    private:
    const double *_pVal, *_pValL, *_pValR;

    public:
    LocalInputs1D() : _pVal(0), _pValL(0), _pValR(0) {}

    LocalInputs1D(const double* pVal, const double* pValL, const double* pValR)
        : _pVal(pVal), _pValL(pValL), _pValR(pValR) {}

    virtual double val(int iVar) const {
        assert(_pVal);
        assert(iVar < numVar);
        return _pVal[iVar];
    }
    virtual double nbr(int iVar, int iNbr) const {
        assert(iVar < numVar);
        switch(iNbr) {
            case 0:
                assert(_pValL);
                return _pValL[iVar];
            case 1:
                assert(_pValR);
                return _pValR[iVar];
            default:
                assert(0);
                return nan("");
        }
    }
};

template<size_t numVar>
class LocalOutputs1D : public LocalOutputs<numVar> {
    private:
    double *_pVal;

    public:
    LocalOutputs1D() : _pVal(0) {}

    LocalOutputs1D(double* pVal) : _pVal(pVal) {}

    virtual void val(int iVar, double newVal) {
        assert(_pVal);
        assert(iVar < numVar);
        _pVal[iVar] = newVal;
    }
};

#endif
