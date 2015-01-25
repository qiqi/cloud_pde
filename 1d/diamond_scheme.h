#ifndef DIAMOND_SCHEME_H
#define DIAMOND_SCHEME_H

#include<vector>
#include<deque>
#include<cstring>
#include<array>
#include<cassert>
#include<cstdlib>
#include<cmath>
#include<mpi.h>
#include"diamond.h"
#include"PngWriter.hpp"

int iProc() {
    int iPr;
    MPI_Comm_rank(MPI_COMM_WORLD, &iPr);
    return iPr;
}

int iProcLeft() {
    int iProc, numProc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    return (iProc + numProc - 1) % numProc;
}

int iProcRight() {
    int iProc, numProc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    return (iProc + 1) % numProc;
}

//       . .
//     . . . .
//   . . . . . .  
// + + . . . . # #
// + + + . . # # # #   
// + + + + # # # # # # 
// + + + ` ` # # # # . .
// + + ` ` ` ` # # . . . .
// . ` ` ` ` ` ` . . . . . .
// _ _ _ _ _ _ _ _ _ _ _ _ _   initial data

// ====== implementation of DiamondDiscretization1D ======= //
//
class LocalOperatorBase {
    // I define this base class, so that a list of LocalOperatrs
    // with different numbers of inputs and outputs can be
    // stored in a queue, and called sequentially using the same
    // virtual function "apply".
    // The first two arguments of "apply" should be pointers to
    // LocalInputs<numInput> and LocalOutputs<numOutput>
    // instances.
    public:
    virtual void apply(const double* pInputs, const double* pInputsL,
                       const double* pInputsR, double* pOutputs,
                       const LocalMesh& mesh) const = 0;
    virtual size_t numInputs() const = 0;
    virtual size_t numOutputs() const = 0;
};

template<size_t numInput, size_t numOutput>
class LocalOperator : public LocalOperatorBase {
    public:
    virtual size_t numInputs() const { return numInput; }
    virtual size_t numOutputs() const { return numOutput; }
    private:
    void (*_operator)(
         const LocalInputs<numInput>& inputs,
         LocalOutputs<numOutput>& outputs,
         const LocalMesh& mesh);

    public:
    LocalOperator(void (&localOperator)(
                  const LocalInputs<numInput>& inputs,
                  LocalOutputs<numOutput>& outputs,
                  const LocalMesh& mesh))
    : _operator(&localOperator) {}

    virtual void apply(const double* pInputs, const double* pInputsL,
                       const double* pInputsR, double* pOutputs,
                       const LocalMesh& mesh) const
    {
        LocalInputs1D<numInput> inputs(pInputs, pInputsL, pInputsR);
        LocalOutputs1D<numOutput> outputs(pOutputs);
        _operator(inputs, outputs, mesh);
    }
};

class LocalVariablesQueue {
    private:
    size_t _maxBytes;
    char* _pData;
    size_t _enqueueByte, _dequeueByte;
    MPI_Request _req;
    int _hasSendOrRecvBeenCalled;

    public:
    LocalVariablesQueue(size_t maxBytes)
    : _maxBytes(maxBytes), _enqueueByte(0), _dequeueByte(0),
      _hasSendOrRecvBeenCalled(0)
    {
        _pData = (char*)malloc(maxBytes);
    }

    virtual ~LocalVariablesQueue() {
        free(_pData);
    }

    void Irecv(int iProc, int tag) {
        assert(!_hasSendOrRecvBeenCalled);
        _hasSendOrRecvBeenCalled = 1;
        MPI_Irecv(_pData, _maxBytes, MPI_BYTE, iProc, tag, MPI_COMM_WORLD, &_req);
        // std::cout << "Irecv()..." << _req << " iProc:" << iProc << ",Tag:" << tag << std::endl;
    }

    void Isend(int iProc, int tag) {
        assert(!_hasSendOrRecvBeenCalled);
        _hasSendOrRecvBeenCalled = 1;
        MPI_Isend(_pData, _maxBytes, MPI_BYTE, iProc, tag, MPI_COMM_WORLD, &_req);
        // std::cout << "Isend()..." << _req << " iProc:" << iProc << ",Tag:" << tag << std::endl;
    }

    void waitForSendOrRecv() {
        assert(_hasSendOrRecvBeenCalled);
        // std::cout << "waitForSendOrRecv()..." << _req << std::endl;
        MPI_Wait(&_req, MPI_STATUS_IGNORE);
        // std::cout << "waitForSendOrRecv() done." << std::endl;
    }

    int isSendOrRecvComplete() {
        assert(_hasSendOrRecvBeenCalled);
        int testResult;
        MPI_Test(&_req, &testResult, MPI_STATUS_IGNORE);
        return testResult;
    }

    void enqueue(size_t numVar, const double* pData) {
        *(size_t*)(_pData + _enqueueByte) = numVar;
        memcpy(_pData + _enqueueByte + sizeof(size_t), pData,
               sizeof(double) * numVar);
        _enqueueByte += sizeof(size_t) + sizeof(double) * numVar;
        assert(_enqueueByte <= _maxBytes);
    }

    size_t dequeue(double* pData) {
        size_t numVar = *(size_t*)(_pData + _dequeueByte);
        memcpy(pData, _pData + _dequeueByte + sizeof(size_t),
               sizeof(double) * numVar);
        _dequeueByte += sizeof(size_t) + sizeof(double) * numVar;
        assert(_dequeueByte <= _maxBytes);
        return numVar;
    }
};

class ClassicSyncer1D {
    // Only used for syncing initial data
    private:
    double *_data;
    size_t _numGrids, _numVar;
    MPI_Request _reqs[4];

    double * _pGrid(int iGrid) {
        return _data + iGrid * _numVar;
    }

    public:
    ClassicSyncer1D(double * data, int numGrids, int numVar)
    : _data(data), _numGrids(numGrids), _numVar(numVar)
    {
        const int iGridLeft = 1, iGridRight = _numGrids;
        MPI_Isend(_pGrid(iGridLeft), numVar, MPI_DOUBLE,
                  iProcLeft(), 0, MPI_COMM_WORLD, _reqs);
        MPI_Isend(_pGrid(iGridRight), numVar, MPI_DOUBLE,
                  iProcRight(), 1, MPI_COMM_WORLD, _reqs + 1);

        const int iGridLeftGhost = 0, iGridRightGhost = _numGrids + 1;
        MPI_Irecv(_pGrid(iGridRightGhost), numVar, MPI_DOUBLE,
                  iProcRight(), 0, MPI_COMM_WORLD, _reqs + 2);
        MPI_Irecv(_pGrid(iGridLeftGhost), numVar, MPI_DOUBLE,
                  iProcLeft(), 1, MPI_COMM_WORLD, _reqs + 3);
    }

    void waitTillDone() {
        MPI_Status stats[4];
        MPI_Waitall(4, _reqs, stats);
    }

    ~ClassicSyncer1D() {
        waitTillDone();
    }
};

class Diamond1D {
    // An instance of this class computes one space-time diamond
    // When constructing the instance, one feeds in a series of LocalMesh
    // instances, spaning the spatial domain of this diamond, as well
    // as a series of LocalOperator instances, spanning the time domain
    // of this diamond. From these inputs, we know how much input data to
    // expect from either the initial condition or the "foundation" of
    // the diamonds.  Computation of the whole diamond is performed
    // upon its construction. Afterwards, call "getRoof" function to obtain
    // the outputs.
    private:
    std::vector<LocalMesh> _localMeshes;
    std::vector<const LocalOperatorBase*> _localOperators;

    // This is where the actual computation is done
    int _numVariables;
    double * _variablesData;

    inline double * _varAtGrid(size_t iGrid) {
        return _variablesData + iGrid * _numVariables;
    }

    LocalVariablesQueue * _pLeftFoundation, * _pRightFoundation;
    LocalVariablesQueue * _pLeftRoof, * _pRightRoof;

    void _step(const double* pInput, double* pOutput,
               size_t iFirst, size_t iLast,
               const LocalOperatorBase* pOp)
    {
        size_t numInput = pOp->numInputs();
        size_t numOutput = pOp->numOutputs();
        for (size_t iGrid = iFirst; iGrid <= iLast; ++iGrid)
        {
            const double* input = pInput + iGrid * numInput;
            const double* inputL = pInput + (iGrid - 1) * numInput;
            const double* inputR = pInput + (iGrid + 1) * numInput;
            double* output = pOutput + iGrid * numOutput;
            pOp->apply(input, inputL, inputR, output, _localMeshes[iGrid]);
        }
    }

    void _dequeueTwoGrids(LocalVariablesQueue* pFoundation, double * pData)
    {
        size_t numVar = pFoundation->dequeue(pData);
        assert(_numVariables == numVar);
        numVar = pFoundation->dequeue(pData + numVar);
        assert(_numVariables == numVar);
    }

    bool _isHalfDiamond() {
        if (_localMeshes.size() == _localOperators.size() + 1) {
            return false;
        } else {
            assert(_localMeshes.size() == _localOperators.size() * 2);
            return true;
        }
    }

    void _computeFirstHalf() 
    {
        assert(!_isHalfDiamond());

        const size_t nGrid = _localMeshes.size();
        _numVariables = _localOperators[0]->numInputs();
        _variablesData = new double[_numVariables * (nGrid + 2)];
        _dequeueTwoGrids(_pLeftFoundation, _varAtGrid(nGrid / 2 - 1));
        _dequeueTwoGrids(_pRightFoundation, _varAtGrid(nGrid / 2 + 1));

        size_t lastStep = nGrid/2 - 2;
        for (size_t iStep = 0; iStep <= lastStep; ++ iStep)
        {
            assert(_numVariables == _localOperators[iStep]->numInputs());
            size_t numOutput = _localOperators[iStep]->numOutputs();
            double * pOutput = new double[numOutput * (nGrid + 2)];

            const size_t iFirst = nGrid/2 - iStep, iLast = nGrid/2 + iStep + 1;
            _step(_variablesData, pOutput, iFirst, iLast, _localOperators[iStep]);

            delete[] _variablesData;
            _variablesData = pOutput;
            _numVariables = numOutput;

            _dequeueTwoGrids(_pLeftFoundation, _varAtGrid(iFirst - 2));
            _dequeueTwoGrids(_pRightFoundation, _varAtGrid(iLast + 1));
        }
    }

    void enqueueTwoGrids(LocalVariablesQueue* pRoof, double * pData)
    {
        pRoof->enqueue(_numVariables, pData);
        pRoof->enqueue(_numVariables, pData + _numVariables);
    }

    void _computeSecondHalf()
    {
        const size_t nGrid = _localMeshes.size();

        size_t firstStep = nGrid/2 - 1, lastStep = nGrid - 2;
        if (_isHalfDiamond()) {
            lastStep -= firstStep;
            firstStep = 0;
        }

        for (size_t iStep = firstStep; iStep <= lastStep; ++ iStep)
        {
            assert(_numVariables == _localOperators[iStep]->numInputs());
            size_t numOutput = _localOperators[iStep]->numOutputs();
            double * pOutput = new double[numOutput * (nGrid + 2)];

            const size_t iFirst = nGrid/2 - lastStep + iStep,
                         iLast = nGrid/2 + lastStep - iStep + 1;
            _step(_variablesData, pOutput, iFirst, iLast,
                  _localOperators[iStep]);

            delete[] _variablesData;
            _variablesData = pOutput;
            _numVariables = numOutput;

            enqueueTwoGrids(_pLeftRoof, _varAtGrid(iFirst));
            enqueueTwoGrids(_pRightRoof, _varAtGrid(iLast - 1));
        }
        delete[] _variablesData;
    }

    public:
    // constructor from initial condition
    Diamond1D(std::vector<LocalMesh>::const_iterator firstGrid,
              std::vector<LocalMesh>::const_iterator lastGrid,
              std::deque<const LocalOperatorBase*>::const_iterator firstOp,
              std::deque<const LocalOperatorBase*>::const_iterator lastOp,
              double* initialData, int numVar)
    :
        _localMeshes(firstGrid, lastGrid),
        _localOperators(firstOp, lastOp)
    {
        // std::cout << "Creating Diamond with initial data" << std::endl;

        assert(_localMeshes.size() == _localOperators.size() * 2);
        _numVariables = numVar;
        _variablesData = new double[numVar * (_localMeshes.size() + 2)];

        size_t bytesToCopy = sizeof(double) * numVar * _localMeshes.size();
        memcpy(_variablesData + numVar, initialData, bytesToCopy);

        ClassicSyncer1D sync(_variablesData, _localMeshes.size(), numVar);
        sync.waitTillDone();

        // This tells the "compute" function to only process the top half
        // of the diamond
        _pLeftFoundation = _pRightFoundation = 0;
    }

    private:
    size_t _foundationBytes()
    {
        if (_isHalfDiamond()) {
            return 0;
        } else {
            size_t foundationBytes = 0;
            for (size_t iOp = 0; iOp < _localMeshes.size() / 2; ++iOp) {
                size_t localVarBytes = sizeof(size_t)
                    + sizeof(double) * _localOperators[iOp]->numInputs();
                foundationBytes += localVarBytes * 2;
            }
            return foundationBytes;
        }
    }

    size_t _roofBytes()
    {
        size_t firstOp = _localMeshes.size() / 2 - 1;
        size_t lastOp = _localMeshes.size() - 2;

        if (_isHalfDiamond()) {
            lastOp -= firstOp;
            firstOp = 0;
        }

        size_t roofBytes = 0;
        for (size_t iOp = firstOp; iOp <= lastOp; ++iOp) {
            size_t localVarBytes = sizeof(size_t)
                + sizeof(double) * _localOperators[iOp]->numOutputs();
            roofBytes += localVarBytes * 2;
        }
        return roofBytes;
    }

    public:
    // constructor from other diamonds
    Diamond1D(std::vector<LocalMesh>::const_iterator firstGrid,
              std::vector<LocalMesh>::const_iterator lastGrid,
              std::deque<const LocalOperatorBase*>::const_iterator firstOp,
              std::deque<const LocalOperatorBase*>::const_iterator lastOp,
              int iProcLeftFoundationIsFrom, int tagLeftFoundationIsFrom,
              int iProcRightFoundationIsFrom, int tagRightFoundationIsFrom)
    :
        _localMeshes(firstGrid, lastGrid),
        _localOperators(firstOp, lastOp)
    {
        // std::cout << "Creating full Diamond" << std::endl;

        assert(_localMeshes.size() % 2 == 0);
        assert(_localMeshes.size() == _localOperators.size() + 1);

        _pLeftFoundation = new LocalVariablesQueue(_foundationBytes());
        _pLeftFoundation->Irecv(iProcLeftFoundationIsFrom,
                                tagLeftFoundationIsFrom);
        _pRightFoundation = new LocalVariablesQueue(_foundationBytes());
        _pRightFoundation->Irecv(iProcRightFoundationIsFrom,
                                 tagRightFoundationIsFrom);
    }

    void compute(int iProcLeftRoofGoesTo, int tagLeftRoofGoesTo,
                 int iProcRightRoofGoesTo, int tagRightRoofGoesTo)
    {
        // std::cout << "Computing Diamond" << std::endl;

        if (_pLeftFoundation && _pRightFoundation) {
            _pLeftFoundation->waitForSendOrRecv();
            _pRightFoundation->waitForSendOrRecv();
            _computeFirstHalf();
            delete _pLeftFoundation;
            delete _pRightFoundation;
        }

        _pLeftRoof = new LocalVariablesQueue(_roofBytes());
        _pRightRoof = new LocalVariablesQueue(_roofBytes());
        _computeSecondHalf();
        _pLeftRoof->Isend(iProcLeftRoofGoesTo, tagLeftRoofGoesTo);
        _pRightRoof->Isend(iProcRightRoofGoesTo, tagRightRoofGoesTo);

        // std::cout << "Finish computing Diamond" << std::endl;
    }

    bool isSendComplete() {
        return _pLeftRoof->isSendOrRecvComplete()
            && _pRightRoof->isSendOrRecvComplete();
    }

    ~Diamond1D() {
        delete _pLeftRoof;
        delete _pRightRoof;
    }
};

class DiamondDiscretization1D {
    private:
    std::deque<const LocalOperatorBase*> _localOperators;
    std::vector<LocalMesh> _localMeshes;
    std::deque<Diamond1D> _diamonds;
    size_t _initialNumVar;
    double* _initialData;

    public:
    ~DiamondDiscretization1D() {
        MPI_Finalize();
    }

    template<size_t numVar>
    DiamondDiscretization1D(int numGrids, double dx,
         void (&localOperator)(
               LocalOutputs<numVar>&, const LocalMesh&))
    {
        MPI_Init(0, 0);

        double x0 = numGrids * dx * iProc();

        assert(numGrids % 2 == 0);
        for (int iGrid = 0; iGrid < numGrids * 3 / 2; ++iGrid) {
            _localMeshes.emplace_back(x0 + iGrid * dx, dx);
        }

        _initialNumVar = numVar;
        _initialData = new double[numGrids * numVar];
        for (int iGrid = 0; iGrid < numGrids; ++iGrid) {
            LocalOutputs1D<numVar> localVar(_initialData + iGrid * numVar);
            localOperator(localVar, _localMeshes[iGrid]);
        }
    }

    private:
    size_t _numGrids() {
        assert(_localMeshes.size() % 3 == 0);
        return _localMeshes.size() / 3 * 2;
    }

    bool _lastDiamondOnTheLeft;

    void _initialDiamond() {
        assert(_diamonds.empty());
        _diamonds.emplace_back(_localMeshes.begin(),
                               _localMeshes.begin() + _numGrids(),
                               _localOperators.begin(),
                               _localOperators.end(),
                               _initialData, _initialNumVar);

        const int tagLeftwards = 1, tagToSelf = 0;
        _diamonds.front().compute(iProcLeft(), tagLeftwards,
                                   iProc(), tagToSelf);

        _lastDiamondOnTheLeft = true;
    }

    void _newDiamond() {
        assert(!_diamonds.empty());

        const int tagLeftwards = 1, tagRightwards = 2, tagToSelf = 0;
        if (_lastDiamondOnTheLeft) {
            _diamonds.emplace_back(_localMeshes.begin(),
                                   _localMeshes.begin() + _numGrids(),
                                   _localOperators.begin(),
                                   _localOperators.end(),
                                   iProc(), tagToSelf,
                                   iProcRight(), tagLeftwards);
            _lastDiamondOnTheLeft = false;
            _diamonds.back().compute(iProc(), tagToSelf,
                                   iProcRight(), tagRightwards);
        }
        else {
            _diamonds.emplace_back(_localMeshes.begin(),
                                   _localMeshes.begin() + _numGrids(),
                                   _localOperators.begin(),
                                   _localOperators.end(),
                                   iProcLeft(), tagRightwards,
                                   iProc(), tagToSelf);
            _lastDiamondOnTheLeft = true;
            _diamonds.back().compute(iProcLeft(), tagLeftwards,
                                     iProc(), tagToSelf);
        }
    }

    public:
    template<size_t numInput, size_t numOutput>
    void applyOp(void (&localOperator)(
                 const LocalInputs<numInput>& inputs,
                 LocalOutputs<numOutput>& outputs,
                 const LocalMesh& mesh))
    {
        _localOperators.push_back(new LocalOperator<numInput, numOutput>(localOperator));

        if (_initialData && _localOperators.size() == _numGrids() / 2)
        {
            _initialDiamond();

            _localOperators.pop_front();
            delete[] _initialData;
            _initialData = 0;
        }

        if (_localOperators.size() == _numGrids() - 1)
        {
            _newDiamond();

            for (size_t iStep = 0; iStep < _numGrids() / 2; ++iStep) {
                _localOperators.pop_front();
            }
        }
    }
};

#endif
