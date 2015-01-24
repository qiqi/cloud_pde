#ifndef DIAMOND_SCHEME_H
#define DIAMOND_SCHEME_H

#include<vector>
#include<cstring>
#include<array>
#include<cassert>
#include<cstdlib>
#include<cmath>
#include<mpi.h>
#include"diamond.h"
#include"PngWriter.hpp"

//       . .
//     . . . .
//   . . . . . .  
// + + . . . . # #
// + + + . . # # # #   
// + + + + # # # # # # 
// + + + ` ` # # # # . .
// + + ` ` ` ` # # . . . .
// + ` ` ` ` ` ` . . . . . .
// _ _ _ _ _ _ _ _ _ _ _ _ _   initial data

// ====== implementation of DiamondDiscretization1D ======= //
//
class LocalOperatorBase {
    // I define this base class, so that a list of LocalOperatrs
    // with different numbers of inputs and outputs can be
    // stored in a queue, and called sequentially using the same
    // virtual function "apply".
    // The first two arguments of "apply" should be pointers to
    // LocalVariables<numInput> and LocalVariables<numOutput>
    // instances.
    public:
    virtual void apply(const void* pInputs, void* pOutputs,
                       const LocalMesh& mesh) const = 0;
    virtual size_t numInputs() const = 0;
    virtual size_t numOutputs() const = 0;
};

template<size_t numInput, size_t numOutput>
class LocalOperator : LocalOperatorBase {
    public:
    virtual size_t numInputs() { return numInput; }
    virtual size_t numOutputs() { return numOutput; }
    private:
    void (*_operator)(
         const LocalVariables<numInput>& inputs,
         LocalVariables<numOutput>& outputs,
         const LocalMesh& mesh);

    public:
    LocalOperator(void (&localOperator)(
                  const LocalVariables<numInput>& inputs,
                  LocalVariables<numOutput>& outputs,
                  const LocalMesh& mesh))
    : _operator(&localOperator) {}

    virtual void apply(const void* pInputs, void* pOutputs,
                       const LocalMesh& mesh) const
    {
        const LocalVariables<numInput>& inputs
            = *(const LocalVariables<numInput>*)pInputs;
        LocalVariables<numOutput>& outputs
            = *(const LocalVariables<numOutput>*)pOutputs;
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
    }

    void Isend(int iProc, int tag) {
        assert(!_hasSendOrRecvBeenCalled);
        _hasSendOrRecvBeenCalled = 1;
        MPI_Isend(_pData, _maxBytes, MPI_BYTE, iProc, tag, MPI_COMM_WORLD, &_req);
    }

    void waitForSendOrRecv() {
        assert(_hasSendOrRecvBeenCalled);
        MPI_Wait(&_req, MPI_STATUS_IGNORE);
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

template<size_t numVar>
class ClassicSyncer1D {
    // Only used for syncing initial data
    private:
    double (*_data)[numVar];
    int _numGrids;
    MPI_Request _reqs[4];

    public:
    ClassicSyncer1D(double * data, int numGrids)
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

    // apply initial condition
    template<size_t numVar>
    inline void _applyInitialization(
         void (&localOperator)(
               LocalVariables<numVar>&, const LocalMesh&),
         double *data, int iGrid)
    {
        double (*pData)[numVar] = (double (*)[numVar]) data;
        LocalVariables1D<numVar> localVars(pData[iGrid]);
        // pData includes the ghost grids,
        // but _localMeshes only include the interior grids
        assert(iGrid > 0 && iGrid <= _localMeshes.size());
        localOperator(localVars, _localMeshes[iGrid - 1]);
    }

    // This is where the actual computation is done
    int _numVariables;
    double * _variablesData;

    LocalVariablesQueue * _pLeftFoundation, * _pRightFoundation;
    LocalVariablesQueue * _pLeftRoof, * _pRightRoof;

    void _step(const double* pInput, double* pOutput, size_t numOutput,
               const LocalOperatorBase* pOp)
    {
    }

    void dequeueTwoGrids(LocalVariablesQueue* pFoundation, double * pData)
    {
        size_t numVar = pFoundation->dequeue(pData);
        assert(_numVariables == numVar);
        numVar = pFoundation->dequeue(pData + numVar);
        assert(_numVariables == numVar);
    }

    void _computeFirstHalf() {
        _numVariables = _localOperators[0]->numInputs();
        _variablesData = new double[_numVariables * 4];
        dequeueTwoGrids(_pLeftFoundation, _variablesData);
        dequeueTwoGrids(_pRightFoundation, _variablesData + 2 * _numVariables);

        size_t lastStep = _localMeshes.size() / 2 - 2;
        for (size_t iStep = 0; iStep <= lastStep; ++ iStep)
        {
            size_t numOutput = _localOperators[iStep]->numOutputs();
            const size_t nGridOutput = iStep * 2 + 2;
            const size_t nGridOutputIncGhost = nGridOutput + 4;
            double * pOutput = new double[numOutput * nGridOutputIncGhost];

            _step(_variablesData, pOutput + numOutput * 2, nGridOutput,
                  _localOperators[iStep]);

            delete[] _variablesData;
            _variablesData = pOutput;
            _numVariables = numOutput;

            dequeueTwoGrids(_pLeftFoundation, _variablesData);
            dequeueTwoGrids(_pRightFoundation, _variablesData
                            + (nGridOutput + 2) * _numVariables);
        }
    }

    void enqueueTwoGrids(LocalVariablesQueue* pRoof, double * pData)
    {
        pRoof->enqueue(_numVariables, pData);
        pRoof->enqueue(_numVariables, pData + _numVariables);
    }

    void _computeSecondHalf()
    {
        size_t firstStep = _localMeshes.size() / 2 - 1;
        size_t lastStep = _localMeshes.size() - 2;
        for (size_t iStep = firstStep; iStep <= lastStep; ++ iStep)
        {
            size_t numOutput = _localOperators[iStep]->numOutputs();
            const size_t nGridOutput = (lastStep - iStep) * 2 + 2;
            double * pOutput = new double[numOutput * nGridOutput];

            _step(_variablesData, pOutput, nGridOutput,
                  _localOperators[iStep]);

            delete[] _variablesData;
            _variablesData = pOutput;
            _numVariables = numOutput;

            enqueueTwoGrids(_pLeftRoof, _variablesData);
            enqueueTwoGrids(_pRightRoof, _variablesData
                            + (nGridOutput - 2) * _numVariables);
        }
        delete[] _variablesData;
    }

    public:
    // constructor from initial condition
    template<size_t numVar>
    Diamond1D(std::vector<LocalMesh>::const_iterator firstGrid,
              std::vector<LocalMesh>::const_iterator lastGrid,
              std::vector<const LocalOperatorBase*>::const_iterator firstOp,
              std::vector<const LocalOperatorBase*>::const_iterator lastOp,
              void (&localOperator)(
                    LocalVariables<numVar>&, const LocalMesh&))
    :
        _localMeshes(firstGrid, lastGrid),
        _localOperators(firstOp, lastOp)
    {
        assert(_localMeshes.size() == _localOperators.size() * 2);
        _numVariables = numVar;
        _variablesData = new double[numVar * (_localMeshes.size() + 2)];

        const int iGridLeft = 1, iGridRight = _localMeshes.size();
        _applyInitialization(localOperator, _variablesData, iGridLeft);
        _applyInitialization(localOperator, _variablesData, iGridRight);

        ClassicSyncer1D<numVar> sync(_variablesData, _localMeshes.size());

        for (int iGrid = iGridLeft + 1; iGrid < iGridRight; ++iGrid) {
            _applyInitialization(localOperator, _variablesData, iGrid);
        }
        sync.waitTillDone();

        // This tells the "compute" function to only process the top half
        // of the diamond
        _pLeftFoundation = _pRightFoundation = 0;
    }

    private:
    size_t _foundationBytes()
    {
        size_t foundationBytes = 0;
        for (size_t iOp = 0; iOp < _localMeshes.size() / 2; ++iOp) {
            size_t localVarBytes = 
                + sizeof(double) * _localOperators[iOp]->numInputs();
            foundationBytes += localVarBytes * 2;
        }
        return foundationBytes;
    }

    size_t _roofBytes()
    {
        size_t roofBytes = 0;
        for (size_t iOp = _localMeshes.size() / 2 - 1;
                iOp < _localMeshes.size() - 1; ++iOp) {
            size_t localVarBytes = 
                + sizeof(double) * _localOperators[iOp]->numOutputs();
            roofBytes += localVarBytes * 2;
        }
        return roofBytes;
    }

    public:
    // constructor from other diamonds
    Diamond1D(std::vector<LocalMesh>::const_iterator firstGrid,
              std::vector<LocalMesh>::const_iterator lastGrid,
              std::vector<const LocalOperatorBase*>::const_iterator firstOp,
              std::vector<const LocalOperatorBase*>::const_iterator lastOp,
              int iProcLeftFoundationIsFrom, int tagLeftFoundationIsFrom,
              int iProcRightFoundationIsFrom, int tagRightFoundationIsFrom)
    :
        _localMeshes(firstGrid, lastGrid),
        _localOperators(firstOp, lastOp)
    {
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

};

#endif
