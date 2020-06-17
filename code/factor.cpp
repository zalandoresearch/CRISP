#include "factor.hpp"
#include "node.hpp"
#include "seir.hpp"
#include <algorithm>
#include <cmath>


#include<iostream>


vector<double> Factor::message_to( Node* n) {

    size_t i_node = 0; // ToDo

    vector<vector<double>::const_iterator> messages(_nodes.size());
    vector< vector<any>::const_iterator> states(_nodes.size());

    for( size_t i=0; i<_nodes.size(); i++) {
        messages[i] = _nodes[i]->_message.begin();
        states[i] = _nodes[i]->_state.begin();
    }
    vector<double> message( n->_message.size(), 0.0);
    vector<double>::iterator msg_it = message.begin();

    while( true) {

        double p = 0.0;
        for( int i=0; i<_nodes.size(); i++) {
            if( i!=i_node) {
                p += *(messages[i]);
            }
        }
        *(msg_it) += p * potential( states);

        for( int i=_nodes.size()-1; i>=0; i--) {
            messages[i]++;
            states[i]++;
            if(_nodes[i] == n) msg_it++;
            
            if( messages[i] == _nodes[i]->_message.end()) {
                messages[i] = _nodes[i]->_message.begin();
                states[i] = _nodes[i]->_state.begin();  
                if(_nodes[i] == n) msg_it = message.begin(); 

                if( i==0) goto done;           
            }          
            else break;
        }
    } 
    done:

    return message;
}



SEIRFactor::SEIRFactor( const vector<double> &qE, const vector<double> &qI, 
                        double p0, double p1, 
                        SEIRNode &in, SEIRNode &out, VirusLoadNode &load) : 
_qE(qE), _qI(qI), _piE(init_pi(qE)), _piI(init_pi(qI)), _p0(p0), _p1(p1) {
    _nodes.resize(3);
    _nodes[0] = &in;
    _nodes[1] = &out;
    _nodes[2] = &load;
}


const vector<double> SEIRFactor::init_pi( const vector<double> q) {

    auto res = vector<double>(q.size());
    double scale = 0;
    for( int i=q.size()-1; i>=0; i--) 
        scale += q[i];

    double sum = 0;
    for( int i=q.size()-1; i>=0; i--) {
        sum += q[i]/scale;
        res[i] = q[i]/scale/sum;
    }
    return res;
}

double SEIRFactor::potential( vector<vector<any>::const_iterator> states) {

    SEIRState in = std::any_cast<SEIRState>(*(states[0]));
    SEIRState out = std::any_cast<SEIRState>(*(states[1]));
    double load = std::any_cast<double>(*(states[2]));

    cout << in << " " << out << " " << load << endl;
    switch( in.phase()) {
        case SEIRState::S: {
            if( in.next(/*change = */ true)  == out) return 1.0 - (1.0-_p0)*pow(1.0-_p1, load);
            if( in.next(/*change = */ false) == out) return (1.0-_p0)*pow(1.0-_p1, load);
        } break;
        case SEIRState::E: {
            if( in.next(/*change = */ true)  == out) return _piE[in.days()];
            if( in.next(/*change = */ false) == out) return 1.0-_piE[in.days()];

        } break;
        case SEIRState::I: {
            if( in.next(/*change = */ true)  == out) return _piI[in.days()];
            if( in.next(/*change = */ false) == out) return 1.0-_piI[in.days()];
        } break;
        case SEIRState::R: {
            return out.phase()==SEIRState::R ? 1.0 : 0.0;
        } break;

    }
    return 0.0;
}