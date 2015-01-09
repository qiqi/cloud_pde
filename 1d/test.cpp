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

class Operator {
};

Operator& operator >> (Operator& op, Variable v) {
    return op;
}

Operator& operator << (Operator& op, Variable v) {
    return op;
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

class Laplacian : public Operator {
    public:
    void compute(std::vector<const LocalVariable> inputs,
                 std::vector<LocalVariable> outputs,
                 LocalMesh& mesh) {
        double u = inputs[0].val();
        double uL = inputs[0].nbr(0);
        double uR = inputs[0].nbr(1);
        outputs[0].val((uL + uR - 2 * u) / (mesh.dx * mesh.dx));
    }
};
class Derivative : public Operator {
};

Variable calcDudt(Variable u) {
    Laplacian laplacian;
    Derivative derivative;

    Variable uxx, diff, conv;
    laplacian >> u << uxx;
    laplacian >> (u + uxx) << diff;
    derivative >> (0.5 * u * u) << conv;
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

