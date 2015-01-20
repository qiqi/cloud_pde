#include<vector>
#include<cstring>
#include<array>
#include<cassert>
#include<cstdlib>
#include<cmath>
#include<mpi.h>
#include"PngWriter.hpp"

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
class LocalVariables {
    public:
    virtual double val(int iVar) const = 0;
    virtual void val(int iVar, double newVal) = 0;
    virtual double nbr(int iVar, int iNbr) const = 0;
};

struct LocalMesh {
    double x, dx;
    LocalMesh(double x, double dx) : x(x), dx(dx) {}
    LocalMesh() : x(nan("")), dx(nan("")) {}
};

// sub-step

template<size_t numInput, size_t numOutput>
class SubStep {
    public:
    typedef void (&LocalOperator)(
         const LocalVariables<numInput>& inputs,
         LocalVariables<numOutput>& output, const LocalMesh& mesh);

    SubStep (LocalOperator localOp)
    :
        _localOp(localOp)
    {
    }

    private:
    LocalOperator _localOp;
};

// ====== common interface of Discretization ======= //

class Discretization {
    public:
    template<size_t numVar>
    Discretization(int numGrids, double dx,
         void (&localOperator)(
               LocalVariables<numVar>&, const LocalMesh&));

    template<size_t numInput, size_t numOutput>
    void applyOp(void (&localOperator)(
                 const LocalVariables<numInput>& inputs,
                 LocalVariables<numOutput>& outputs,
                 const LocalMesh& mesh));
};

// ====== implementation of ClassicDiscretization1D ======= //

template<size_t numVar>
class LocalVariables1D : public LocalVariables<numVar> {
    private:
    double *_pVal, *_pValL, *_pValR;

    public:
    LocalVariables1D() : _pVal(0), _pValL(0), _pValR(0) {}

    LocalVariables1D(double* pVal) : _pVal(pVal), _pValL(0), _pValR(0) {}

    LocalVariables1D(double* pVal, double* pValL, double* pValR)
        : _pVal(pVal), _pValL(pValL), _pValR(pValR) {}

    void assignPointers(double* pVal, double* pValL, double* pValR) {
        _pVal = pVal; _pValL = pValL; _pValR = pValR;
    }

    virtual double val(int iVar) const {
        assert(_pVal);
        return _pVal[iVar];
    }
    virtual void val(int iVar, double newVal) {
        assert(_pVal);
        _pVal[iVar] = newVal;
    }
    virtual double nbr(int iVar, int iNbr) const {
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

class ClassicDiscretization1D {
    private:
    int _numGrids;
    double _dx;
    double _x0;
    int _numVariables;
    double * _variablesData;

    char* _pngFilename;
    PngWriter _png;

    void _commonInit() {
        MPI_Init(0, 0);

        int iProc;
        MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
        _x0 = _numGrids * _dx * iProc;
    }

    public:
    ClassicDiscretization1D()
    : _numGrids(100), _dx(1.), 
      _variablesData(0), _pngFilename(0), _png(0,0)
    {
        _commonInit();
    }

    virtual ~ClassicDiscretization1D() {
        if (_variablesData) {
            delete _variablesData;
        }
        MPI_Finalize();
    }

    private:
    template<size_t numVar>
    inline void _applyInitialization(
         void (&localOperator)(
               LocalVariables<numVar>&, const LocalMesh&),
         double *data, int iGrid)
    {
        double (*pData)[numVar] = (double (*)[numVar]) data;
        LocalVariables1D<numVar> localVars(pData[iGrid]);
        LocalMesh localMesh(_x0 + iGrid * _dx, _dx);
        localOperator(localVars, localMesh);
    }

    template<size_t numVar>
    class _Syncer {
        private:
        double (*_data)[numVar];
        int _numGrids;
        MPI_Request _reqs[4];
    
        public:
        _Syncer(double * data, int numGrids)
        : _data((double(*)[numVar])data), _numGrids(numGrids)
        {
            static int numProc, iProc, iProcLeft, iProcRight;
                // static vars are auto-initialized to 0
            if (numProc == 0) {
                MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
                MPI_Comm_size(MPI_COMM_WORLD, &numProc);
                iProcLeft = (iProc + numProc - 1) % numProc;
                iProcRight = (iProc + 1) % numProc;
            }
    
            const int iGridLeft = 1, iGridRight = _numGrids;
            MPI_Isend(_data[iGridLeft], numVar, MPI_DOUBLE,
                      iProcLeft, 0, MPI_COMM_WORLD, _reqs);
            MPI_Isend(_data[iGridRight], numVar, MPI_DOUBLE,
                      iProcRight, 1, MPI_COMM_WORLD, _reqs + 1);
    
            const int iGridLeftGhost = 0, iGridRightGhost = _numGrids + 1;
            MPI_Irecv(_data[iGridRightGhost], numVar, MPI_DOUBLE,
                      iProcRight, 0, MPI_COMM_WORLD, _reqs + 2);
            MPI_Irecv(_data[iGridLeftGhost], numVar, MPI_DOUBLE,
                      iProcLeft, 1, MPI_COMM_WORLD, _reqs + 3);
        }
    
        void waitTillDone() {
            MPI_Status stats[4];
            MPI_Waitall(4, _reqs, stats);
        }
    
        ~_Syncer() {
            waitTillDone();
        }
    };


    public:
    template<size_t numVar>
    ClassicDiscretization1D(int numGrids, double dx,
         void (&localOperator)(
               LocalVariables<numVar>&, const LocalMesh&))
    :
        _numGrids(numGrids), _dx(dx), _numVariables(numVar),
        _pngFilename(0), _png(0,0)
    {
        _commonInit();
        _variablesData = new double[numVar * (_numGrids + 2)];
    
        const int iGridLeft = 1, iGridRight = _numGrids;
        _applyInitialization(localOperator, _variablesData, iGridLeft);
        _applyInitialization(localOperator, _variablesData, iGridRight);

        _Syncer<numVar> sync(_variablesData, _numGrids);

        for (int iGrid = iGridLeft + 1; iGrid < iGridRight; ++iGrid) {
            _applyInitialization(localOperator, _variablesData, iGrid);
        }

        sync.waitTillDone();
    }

    private:
    template<size_t numInput, size_t numOutput>
    void _applyLocalOp(void (&localOperator)(
                     const LocalVariables<numInput>& inputs,
                     LocalVariables<numOutput>& outputs,
                     const LocalMesh& mesh),
                 double * input, double * output, int iGrid)
    {
        double (*pInput)[numInput] = (double(*)[numInput])input;
        double (*pOutput)[numOutput] = (double(*)[numOutput])output;

        LocalVariables1D<numInput> localInputs(pInput[iGrid],
                                               pInput[iGrid - 1],
                                               pInput[iGrid + 1]);
        LocalVariables1D<numOutput> localOutputs(pOutput[iGrid]);
        LocalMesh localMesh(iGrid * _dx, _dx);

        localOperator(localInputs, localOutputs, localMesh);
    }

    public:
    template<size_t numInput, size_t numOutput>
    void applyOp(void (&localOperator)(
                 const LocalVariables<numInput>& inputs,
                 LocalVariables<numOutput>& outputs,
                 const LocalMesh& mesh))
    {
        assert(numInput == _numVariables);
        double * newVariablesData = new double[numOutput * (_numGrids + 2)];

        const int iGridLeft = 1, iGridRight = _numGrids;

        _applyLocalOp(localOperator, _variablesData, newVariablesData, iGridLeft);
        _applyLocalOp(localOperator, _variablesData, newVariablesData, iGridRight);

        _Syncer<numOutput> sync(newVariablesData, _numGrids);

        for (int iGrid = iGridLeft + 1; iGrid < iGridRight; ++iGrid) {
            _applyLocalOp(localOperator, _variablesData, newVariablesData, iGrid);
        }

        delete[] _variablesData;
        _variablesData = newVariablesData;
        _numVariables = numOutput;

        sync.waitTillDone();
    }

    // ------------ write to png file ------------- //

    public:

    struct {
        class _VarColor {
            private:
            int _iVar;
            double _low, _high;

            public:
            _VarColor() : _iVar(0), _low(0.0), _high(INFINITY) {}
            void set(int iVar, double low, double high) {
                _iVar = iVar;
                _low = low;
                _high = high;
            }

            void assertIVarLessThan(int numVars) {
                assert(_iVar < numVars);
            }

            double map(double* pVal) {
                return (pVal[_iVar] - _low) / (_high - _low);
            }
        } red, green, blue;

        void assertIVarLessThan(int numVars) {
            red.assertIVarLessThan(numVars);
            green.assertIVarLessThan(numVars);
            blue.assertIVarLessThan(numVars);
        }
    } colorMap;

    void variablesToColor(int iStep)
    {
        for (int iGrid = 0; iGrid < _numGrids; ++iGrid) {
            colorMap.assertIVarLessThan(_numVariables);
            double* pGrid = _variablesData + iGrid * _numVariables;
            double r = colorMap.red.map(pGrid),
                   g = colorMap.green.map(pGrid),
                   b = colorMap.blue.map(pGrid);
            _png.set(iGrid, iStep, r, g, b);
        }
    }

    void writePng() {
        if (_pngFilename) {
            _png.write(_pngFilename);
        }
    }

    void writePng(const char* filename) {
        if (_pngFilename) {
            delete[] _pngFilename;
        }

        _pngFilename = new char[strlen(filename) + strlen(mpiRankString()) + 5];
        strcpy(_pngFilename, filename);
        strcat(_pngFilename, mpiRankString());
        strcat(_pngFilename, ".png");
        writePng();
    }

};


