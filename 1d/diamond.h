#include<vector>
#include<array>
#include<cassert>
#include<cmath>

class LocalVariable {
    public:
    virtual double val() const = 0;
    virtual void val(double newVal) = 0;
    virtual double nbr(int i) const = 0;
};

class LocalMesh {
    public:
    double dx;
};

// sub-step

template<size_t numInput, size_t numOutput>
class SubStep {
    public:
    typedef void (&LocalOperator)(
         std::array<const LocalVariable*, numInput> inputs,
         LocalVariable* output, LocalMesh& mesh);

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
               std::array<LocalVariable*, numVar>, LocalMesh&));

    template<size_t numInput, size_t numOutput>
    void applyOp(void (&localOperator)(
                 std::array<const LocalVariable*, numInput> inputs,
                 std::array<LocalVariable*, numOutput> outputs,
                 LocalMesh& mesh));
};

// ====== implementation of ClassicDiscretization1D ======= //

class LocalVariable1D : public LocalVariable {
    private:
    double *_pVal, *_pValL, *_pValR;

    public:
    LocalVariable1D() : _pVal(0), _pValL(0), _pValR(0) {}

    void assignPointers(double* pVal, double* pValL, double* pValR) {
        _pVal = pVal; _pValL = pValL; _pValR = pValR;
    }

    virtual double val() const {
        assert(_pVal);
        return *_pVal;
    }
    virtual void val(double newVal) {
        assert(_pVal);
        *_pVal = newVal;
    }
    virtual double nbr(int i) const {
        switch(i) {
            case 0:
                assert(_pValL);
                return *_pValL;
            case 1:
                assert(_pValR);
                return *_pValR;
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

    public:
    ClassicDiscretization1D() : _variablesData(0) {}

    virtual ~ClassicDiscretization1D() {
        if (_variablesData) {
            delete _variablesData;
        }
    }

    template<size_t numVar>
    ClassicDiscretization1D(int numGrids, double dx,
         void (&localOperator)(
               std::array<LocalVariable*, numVar>, LocalMesh&))
    :
        _numGrids(numGrids), _dx(dx)
    {
        _numVariables = numVar;
        _variablesData = new double[numVar * _numGrids];
    
        LocalMesh localMesh;
        std::array<LocalVariable1D, numVar> localVars;
        std::array<LocalVariable*, numVar> pLocalVars;
        for (int iVar = 0; iVar < numVar; ++iVar) {
            pLocalVars[iVar] = &localVars[iVar];
        }
        for (int iGrid = 0; iGrid < _numGrids; ++iGrid) {
            double *pVal = _variablesData + iGrid * numVar;
            for (int iVar = 0; iVar < numVar; ++iVar) {
                localVars[iVar].assignPointers(pVal + iVar, 0, 0);
            }
            localOperator(pLocalVars, localMesh);
        }
    }

    template<size_t numInput, size_t numOutput>
    void applyOp(void (&localOperator)(
                 std::array<const LocalVariable*, numInput> inputs,
                 std::array<LocalVariable*, numOutput> outputs,
                 LocalMesh& mesh))
    {
    }
};


