#ifndef F_GRAPH_HPP
#define F_GRAPH_HPP

#include "node.hpp"

using namespace std;

class Node;


class Factor {
public:
    vector<Node*> _nodes;

    virtual double potential( const vector<unsigned int> &) = 0;

public:
    Factor( const vector<Node*> &nodes);
    Factor( const Factor &other) :_nodes(other._nodes) {}
    virtual ~Factor() {}

    virtual MessagePtr message_to( Node *);
};


class TableFactor: public Factor {
    const Message _tab;
protected:
    virtual double potential( const vector<unsigned int> &);

public:
    TableFactor( vector<Node*> &nodes, Message tab);

};



class SEIRFactor : public Factor {

    //const Message& _qE; 
    //const Message& _qI;
    const Message _piE; 
    const Message _piI;
    double _p0;
    double _p1;

    static const vector<double> init_pi( const vector<double> q);

    static vector<Node *> init_helper(SEIRNode &in, SEIRNode &out, vector<SEIRNode *> contacts );

    MessagePtr message_horizontally( bool forward);
    MessagePtr message_vertically( Node *n);

    const SEIRStateSpace &_states;

public:
    SEIRFactor( const Message &qE, const Message &qI,
                double p0, double p1, 
                SEIRNode &in, SEIRNode &out, 
                vector<SEIRNode *> contacts = vector<SEIRNode*>());
    SEIRFactor( const SEIRFactor &other) : Factor(other), _piE( other._piE), _piI(other._piI), _p0(other._p0), _p1(other._p1), _states(other._states) {}
    virtual ~SEIRFactor(){}
    MessagePtr message_to( Node *);

    virtual double potential( const vector<unsigned int> &);

};


class SEIRInitFactor: public Factor {
    bool _patient_zero;
    const SEIRStateSpace &_states;

public:
    SEIRInitFactor( SEIRNode &in, bool patient_zero = false);
    virtual ~SEIRInitFactor(){}
    virtual double potential( const vector<unsigned int> &);
};


class SEIRTestFactor: public Factor {
    bool _positive;
    double _alpha;
    double _beta;
    const SEIRStateSpace &_states;

public:
    SEIRTestFactor( SEIRNode &in, bool positive, double alpha, double beta);
    virtual ~SEIRTestFactor(){}
    virtual double potential( const vector<unsigned int> &);
};

#endif