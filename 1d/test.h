#include<vector>
#include<array>

class Variable {
    private:
    std::vector<VariableState> states;

    public:
};

class VariableState {
    private:
    bool hasNeighborValue;
    std::vector<Variable *> dependees;

    public:
    VariableState() : hasNeighborValue(true) {}
};

class LocalVariable {
    public:
    virtual double val() const = 0;
    virtual void val(double newVal) = 0;
    virtual double nbr(int i) const = 0;
};

class LocalVariableStencilTester : public LocalVariable {
    private:
    mutable bool _isLocal; // set to false once const function nbr is called

    public:
    LocalVariableStencilTester() : _isLocal(true) {}
    bool isLocal() {
        return _isLocal;
    }

    virtual double val() const {
        return 0.0;
    }
    virtual void val(double) {
    }
    virtual double nbr(int) const {
        _isLocal = false;
        return 0.0;
    }
};

class LocalMesh {
    public:
    double dx;
};

// single-input single-output operator

class Operator {
    public:
    typedef void (*LocalOperator)(
         const LocalVariable& input, LocalVariable& output, LocalMesh& mesh);

    private:
    LocalOperator _localOp;
    Variable _input;
    Variable _output;

    bool isLocalOperation() {
        LocalVariableStencilTester dummyLocalIn, dummyLocalOut;
        LocalMesh localMesh;
        _localOp(dummyLocalIn, dummyLocalOut, localMesh);
        return dummyLocalIn.isLocal();
    }

    friend Operator operator >> (Operator op, Variable v);
    friend Operator operator << (Operator op, Variable v);

    public:
    Operator (LocalOperator localOp)
    : _localOp(localOp) {
    }
};

Operator operator >> (Operator op, const Variable v) {
    op._input = v;
    return op;
}

Operator operator << (Operator op, Variable v) {
    op._output = v;
    return op;
}

// multi-input single-output operator

template<int numInput>
class nInputOperator {
    public:
    typedef void (*LocalOperator)(
         std::array<const LocalVariable&, numInput> input,
         LocalVariable& output, LocalMesh& mesh);

    private:
    LocalOperator _localOp;
    std::array<Variable, numInput> _inputs;
    Variable _output;

    friend Operator operator >> (Operator op, Variable v);
    friend Operator operator << (Operator op, Variable v);

    public:
    nInputOperator (LocalOperator localOp)
    : _localOp(localOp) {
    }
};

HERE==================

Operator operator >> (Operator op, const Variable v) {
    op._input = v;
    return op;
}

Operator operator << (Operator op, Variable v) {
    op._output = v;
    return op;
}

// pre-defined operators

void _predefinedUnitaryMinus(const LocalVariable& input,
                             LocalVariable& output,
                             LocalMesh&) {
    output.val(-input.val());
}

Variable operator - (Variable x) {
    static Operator op(_predefinedUnitaryMinus);
    Variable minusX;
    op >> x << minusX;
    return minusX;
}

Variable operator + (Variable x, Variable y) {
    return x;
}

Variable operator - (Variable x, Variable y) {
    return x;
}

Variable operator += (Variable x, Variable y) {
    return x;
}

Variable operator * (Variable x, Variable y) {
    return y;
}

Variable operator * (double a, Variable x) {
    return x;
}

