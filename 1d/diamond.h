#include<vector>
#include<cstring>
#include<array>
#include<cassert>
#include<cmath>
#include"PngWriter.hpp"

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
    int _numVariables;
    double * _variablesData;

    char* _pngFilename;
    PngWriter _png;

    public:
    ClassicDiscretization1D()
        : _variablesData(0), _pngFilename(0), _png(0,0) {}

    virtual ~ClassicDiscretization1D() {
        if (_variablesData) {
            delete _variablesData;
        }
    }

    template<size_t numVar>
    ClassicDiscretization1D(int numGrids, double dx,
         void (&localOperator)(
               LocalVariables<numVar>&, const LocalMesh&))
    :
        _numGrids(numGrids), _dx(dx), _numVariables(numVar),
        _pngFilename(0), _png(0,0)
    {
        _variablesData = new double[numVar * _numGrids];
    
        LocalVariables1D<numVar> localVars;
        for (int iGrid = 0; iGrid < _numGrids; ++iGrid) {
            localVars.assignPointers(_variablesData + iGrid * numVar, 0, 0);
            LocalMesh localMesh(iGrid * _dx, _dx);
            localOperator(localVars, localMesh);
        }
    }

    template<size_t numInput, size_t numOutput>
    void applyOp(void (&localOperator)(
                 const LocalVariables<numInput>& inputs,
                 LocalVariables<numOutput>& outputs,
                 const LocalMesh& mesh))
    {
        assert(numInput == _numVariables);
        double * newVariablesData = new double[numOutput * _numGrids];

        LocalVariables1D<numInput> localInputs;
        LocalVariables1D<numOutput> localOutputs;

        for (int iGrid = 0; iGrid < _numGrids; ++iGrid) {
            double * pInput = _variablesData + iGrid * numInput;

            int iGridLeft = (iGrid + _numGrids - 1) % _numGrids;
            double * pInputL = _variablesData + iGridLeft * numInput;

            int iGridRight = (iGrid + 1) % _numGrids;
            double * pInputR = _variablesData + iGridRight * numInput;

            localInputs.assignPointers(pInput, pInputL, pInputR);

            double * pOutput = newVariablesData + iGrid * numOutput;
            localOutputs.assignPointers(pOutput, 0,0);

            LocalMesh localMesh(iGrid * _dx, _dx);
            localOperator(localInputs, localOutputs, localMesh);
        }

        delete[] _variablesData;
        _variablesData = newVariablesData;
        _numVariables = numOutput;
    }

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
        _pngFilename = new char[strlen(filename) + 1];
        strcpy(_pngFilename, filename);
        writePng();
    }

};


