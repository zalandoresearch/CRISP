#ifndef F_GRAPH_HPP
#define F_GRAPH_HPP

#include "node.hpp"

#include <vector>
#include <any>
#include <iostream>
using namespace std;

class Node;


class Factor {
protected:
    vector<Node*> _nodes;

    virtual double potential( vector<vector<any>::const_iterator>) = 0;

public:
    vector<double> message_to( Node *);
};



class SEIRFactor : public Factor {

    const vector<double>& _qE; 
    const vector<double>& _qI;
    const vector<double> _piE; 
    const vector<double> _piI;
    double _p0;
    double _p1;

    static const vector<double> init_pi( const vector<double> q);

public:
    SEIRFactor( const vector<double> &qE, const vector<double> &qI,
                double p0, double p1, 
                SEIRNode &in, SEIRNode &out, VirusLoadNode &load);

    virtual double potential( vector<vector<any>::const_iterator>);

};

#endif