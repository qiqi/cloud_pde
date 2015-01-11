#include<vector>

class Variable {
    private:

    bool synced;
    double * data;

    public:
};

Variable operator - (Variable x) {
    return x;
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

class LocalVariable {
    private:
    double* _val;
    std::vector<double*> _nbr;

    public:
    double val() const {
        return *_val;
    }
    void val(double newVal) {
        *_val = newVal;
    }
    double nbr(int i) const {
        return *_nbr.at(i);
    }
};

class LocalMesh {
    public:
    double dx;
};

typedef void (*LocalOperator)(std::vector<const LocalVariable> inputs,
             std::vector<LocalVariable> outputs,
             LocalMesh& mesh);

class Operator {
    public:
    Operator (LocalOperator localOp) {
    }
};

Operator operator >> (Operator op, Variable v) {
    return op;
}

Operator operator << (Operator op, Variable v) {
    return op;
}

void laplacian(std::vector<const LocalVariable> inputs,
             std::vector<LocalVariable> outputs,
             LocalMesh& mesh) {
    double u = inputs[0].val();
    double uL = inputs[0].nbr(0);
    double uR = inputs[0].nbr(1);
    outputs[0].val((uL + uR - 2 * u) / (mesh.dx * mesh.dx));
}

void derivative(std::vector<const LocalVariable> inputs,
             std::vector<LocalVariable> outputs,
             LocalMesh& mesh) {
    double uL = inputs[0].nbr(0);
    double uR = inputs[0].nbr(1);
    outputs[0].val((uR - uL) / (2 * mesh.dx));
}

Variable calcDudt(Variable u) {
    Variable uxx, diff, conv;

    Operator(laplacian) >> u << uxx;
    Operator(laplacian) >> (u + uxx) << diff;
    Operator(derivative) >> (0.5 * u * u) << conv;
    return -conv - diff;
}

int main()
{
    double dt = 0.01;
    Variable u;
    for (int iSteps = 0; iSteps < 100; ++ iSteps) {
        Variable u2 = u + 0.5* dt * calcDudt(u);
        u += dt * calcDudt(u2);
    }
}

