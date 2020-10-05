#ifndef _FACTOR_HPP
#define _FACTOR_HPP

#include <cassert>
#include "node.hpp"
#include "distribution.hpp"

using namespace std;

class Node;

// An abstract base class for a factor in a factor graph
class Factor {
    Factor(const Factor &) = delete;
protected:
    vector<Node*> _nodes;           // List of all associated variable nodes

    // The actual value of the factor for an assignment of all discrete variables
    virtual double potential(const vector<unsigned int>&) = 0;
public:
    // Constructor from a list of variable node pointers (copying the list and registerig the factor with each variable)
    Factor(const vector<Node*>& nodes) : _nodes(nodes) { for(auto n: _nodes) { n->add_factor(this); } }
    // Since nothing gets allocated, the destructor is standard
    virtual ~Factor() {}

    // Computes the message to the variable node using the belief propagation algorithm 
    virtual void message_to_variable(Node*, MessagePtr);
    // Returns the pointer to the i-th variable node
    Node* get_node(size_t i) { assert(i >= 0 && i < _nodes.size()); return(_nodes[i]); }
};

// A TableFactor is a factor where the factor value is given by a static (flattened) table
class TableFactor: public Factor {
    const vector<double> _tab;
protected:
    // Retrieves the value of function from a flattened table
    virtual double potential(const vector<unsigned int>&);
public:
    // Standard constructor
    TableFactor(vector<Node*>& nodes, vector<double>& tab) : Factor(nodes), _tab(tab) {}
};


// A class for the main infection spread propagation factor of the CRISP++ model (see (10) in the CRISP++ paper)
class SEIRFactor : public Factor {
    const vector<double>    _piE;               // Distribution of patients being n days in E state (see (6) and (7))
    const vector<double>    _piI;               // Distribution of patients being n days in I state (see (6) and (8))
    double                  _p0;                // Exogeneous propability of infection spread per day (see (5))
    double                  _p1;                // Probability of infection spread per contact (see (5))
    const SEIRStateSpace&   _states;            // Complete set of detailed SEIR states

    // Computes the probablity of being n days in a state (E or I) from the logP and lopPTail
    static const vector<double> init_pi(const Distribution&);

    // Computes the complete set of neighboring factors for the main SEIR infection spread factor
    static vector<Node *> init_helper(SEIRNode&, SEIRNode&, vector<SEIRNode*>);

    // Computes the forward message (see Subsection 3.2.1 in the CRISP++ paper)
    void message_forward(MessagePtr);
    // Computes the backward message (see Subsection 3.2.2 in the CRISP++ paper)
    void message_backward(MessagePtr);
    // Computes the message to a contact (see Subsection 3.2.3 in the CRISP++ paper)
    void message_vertical(Node*, MessagePtr);
public:
    // Standard constructor
    SEIRFactor(const Distribution&, const Distribution&,
               double, double, 
               SEIRNode&, SEIRNode&, 
               vector<SEIRNode*> = vector<SEIRNode*>());
    virtual ~SEIRFactor() {}

    // Computes the message to the variable node using the belief propagation algorithm in Section 3.2 of the CRISP++ paper
    void message_to_variable(Node*, MessagePtr);
    // An SEIRFactor implements its own message_to_variable() method
    virtual double potential(const vector<unsigned int>&) { assert(false); return 0.0; }
};


// A class for the infection state prior factor of the CRISP++ model (f_{u,0} in Section 3.2 of the CRISP++ paper)
class SEIRInitFactor: public Factor {
    bool                    _patient_zero;          // Is this a patient who is in state E_1 with probability 1.0?
    double                  _p0;                    // What is the exogeneous probability that infects peope a-priori
    const SEIRStateSpace&   _states;                // Complete set of detailed SEIR states
public:
    // Standard constructor
    SEIRInitFactor(SEIRNode&, bool = false, double = 0);
    // Empty destructor
    virtual ~SEIRInitFactor() {}

    // Computes the prior message taking into account the person being a patient-zero
    void message_to_variable(Node*, MessagePtr);
    // An SEIRInitFactor implements its own message_to_variable() method
    virtual double potential(const vector<unsigned int>&) { assert(false); return 0.0; }
};


// A class for the test outcome factor of the CRISP++ model (equation (19) and(12) in the CRISP++ paper)
class SEIRTestFactor: public Factor {
    bool                    _positive;          // Was this a positive test outcome?
    double                  _alpha;             // The false-negative rate of the test
    double                  _beta;              // The false-positive rate of the test
    const SEIRStateSpace&   _states;            // Complete set of detailed SEIR states
public:
    // Standard constructor
    SEIRTestFactor(SEIRNode&, bool, double, double);
    // Empty destructor
    virtual ~SEIRTestFactor() {}
    // The value of the factor for all detailed SEIR states
    virtual double potential(const vector<unsigned int>&);
};

#endif
