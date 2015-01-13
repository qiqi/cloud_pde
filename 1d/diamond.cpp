#include"diamond.h"

template<size_t numInput, size_t numOutput>
SubStep<numInput, numOutput>::SubStep (
        SubStep<numInput, numOutput>::LocalOperator localOp)
:
    _localOp(localOp)
{

}

