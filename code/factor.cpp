
#include "factor.hpp"
#include "node.hpp"
#include "seir.hpp"

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>

using namespace std;


Factor::Factor( const vector<Node*> &nodes, const vector<Node*> &child_nodes) : 
    _nodes(nodes) 
{
    if( child_nodes.size()>0) {
        for( auto n: child_nodes) {
            n->_factors.push_back(this);
            n->_messages.emplace_back( new Message(n->_N, 1.0) );
        }
    }
    else {
        for( auto n: _nodes) {
            n->_factors.push_back(this);
            n->_messages.emplace_back( new Message(n->_N, 1.0) );
        }
    }
}


MessagePtr Factor::message_to( Node* n) {

    size_t i_node = 0; 
    for( ; i_node<_nodes.size() && _nodes[i_node]!=n; i_node++){
    }
    assert( i_node<_nodes.size()); // n needs to be in _nodes


    vector< MessagePtr>  incoming_messages;
    vector< Message::const_iterator>   msg_its;

    vector< unsigned int> state_its(_nodes.size());

    for( unsigned int i=0; i!=_nodes.size(); i++ ) {
        incoming_messages.push_back(_nodes[i]->message_to(this));
        msg_its.push_back(incoming_messages[i]->begin());
    }

    MessagePtr outgoing_message( new Message(n->_N, 0.0));
    Message::iterator out_msg_it = outgoing_message->begin();

    while( true) {

        double p = 1.0;
        for(  size_t i=0; i<_nodes.size(); i++) {
            if( i!=i_node) {
                p *= *(msg_its[i]);
            }
        }
        double upd = p * potential( state_its);
        // cout << *(out_msg_it) << " " << upd << " "<< p << " " << endl;
        *(out_msg_it) += upd;

        // generate the next state configuration
        for( int i=_nodes.size()-1; i>=0; i--) {
            msg_its[i]++;
            state_its[i]++;
            if(_nodes[i] == n) out_msg_it++;
            
            if( msg_its[i] == incoming_messages[i]->end()) {
                // if they overrun, wrap them over to the beginning 
                msg_its[i]   = incoming_messages[i]->begin();
                state_its[i] = 0;  
                if(_nodes[i] == n) 
                    out_msg_it = outgoing_message->begin(); 

                // if the first iterator wrapped around we are done
                if( i==0) goto done;           
            }          
            else break;
        }
    } 
    done:

    return outgoing_message;
}

TableFactor::TableFactor( vector<Node*> &nodes, Message tab) :
    Factor(nodes), _tab(tab)
{
}

double TableFactor::potential( const vector<unsigned int> & state_its) {
    assert( state_its.size()== _nodes.size());

    int idx = 0;
    for( unsigned int i=0; i<_nodes.size(); i++) {
        if( i>0) idx *= _nodes[i-1]->size();
        idx += state_its[i];
        //cout << _nodes[i]->size() << " " <<  std::any_cast<int>(*(state_its[i])) << endl;
    }
    //cout << "ïdx: " << idx << endl;
    return _tab[idx];
}


SEIRFactor::SEIRFactor( const Distribution &qE, const Distribution &qI, 
                        double p0, double p1, 
                        SEIRNode &in, SEIRNode &out, vector<SEIRNode *> contacts,
                        const vector<Node*> &child_nodes
                        ) : 
    Factor(init_helper(in, out, contacts), child_nodes), _piE(init_pi(qE)), _piI(init_pi(qI)), _p0(p0), _p1(p1), _states(in.states()) 
{
}


vector<Node *> SEIRFactor::init_helper(SEIRNode &in, SEIRNode &out, vector<SEIRNode *> contacts ) {
    vector<Node *> nodes({ &in, &out});
    nodes.insert( nodes.end(), contacts.begin(), contacts.end());
    return nodes;
}


const vector<double> SEIRFactor::init_pi( const Distribution &q) {

    auto res = Message(q.getMaxOutcomeValue()+1);
    for( int i=0; i<res.size(); i++)
        res[i] = exp( q.getLogP(i)-q.getLogPTail(i)); 
    return res;
}


MessagePtr SEIRFactor::message_horizontally( bool forward) {

    VirusLoad load;
    for( size_t i=2; i<_nodes.size(); i++) {
        MessagePtr  contact_message = _nodes[i]->message_to(this);
        double p = 0, p_tot = 0;
        for( size_t i=0; i< contact_message->size(); i++) {
            p_tot += (*contact_message)[i];
            if( _states[i].phase() == SEIRState::I) 
                p += (*contact_message)[i];
        }
        load.add_source( p/p_tot, 1.0);
       }

    MessagePtr outgoing_message( new Message(_nodes[ forward ? 1 : 0]->size(), 0.0));

    MessagePtr input_message;
    MessagePtr output_message;
    if( forward) {
        input_message = _nodes[0]->message_to(this);
        output_message = outgoing_message;
    }
    else {
        input_message = outgoing_message;
        output_message = _nodes[1]->message_to(this);
    }

    Message::const_iterator it_input_message;
    Message::const_iterator it_output_message;

    it_input_message = input_message->begin();
    for( auto it_input_states = _states.cbegin(); it_input_states != _states.cend(); ++it_input_states) 
    {
        auto in = *it_input_states;

        it_output_message = output_message->begin();
        for( auto it_output_states = _states.cbegin(); it_output_states != _states.cend(); ++it_output_states) 
        {
            auto out = *it_output_states;

            const double& in_value = *(forward ? it_input_message : it_output_message);
            double &out_value = const_cast<double &>(*(forward ? it_output_message : it_input_message));
    
            switch( in.phase() ) {
                case SEIRState::S: {
                    double p_keep = 0.0;
                    for( auto it_load = load.cbegin(); it_load != load.cend(); ++it_load) {
                        p_keep += (it_load->second) * pow(1.0-_p1, it_load->first);
                    }
                    p_keep *= (1.0-_p0);
                    if( in.next(/*change = */ true)  == out) {
                        out_value += in_value * (1.0 - p_keep);
                        //cerr << in << "," << out << "=" << in_value <<", "<<out_value<<endl;
                        }
                    if( in.next(/*change = */ false) == out) {
                        out_value += in_value * p_keep;
                        //cerr << in << "," << out << "=" << in_value <<", "<<out_value<<endl;
                        }
                } break;
                case SEIRState::E: {
                    if( in.next(/*change = */ true)  == out) {
                        out_value += in_value * _piE[in.days()];
                        //cerr << in << "," << out << "=" << in_value <<", "<<out_value<<endl;
                        }
                    if( in.next(/*change = */ false) == out) {
                        out_value += in_value * (1.0-_piE[in.days()]);
                        //cerr << in << "," << out << "=" << in_value <<", "<<out_value<<endl;
                        }

                } break;
                case SEIRState::I: {
                    if( in.next(/*change = */ true)  == out) {
                        out_value += in_value * _piI[in.days()];
                        //cerr << in << "," << out << "=" << in_value <<", "<<out_value<<endl;
                        }
                    if( in.next(/*change = */ false) == out) {
                        out_value += in_value * (1.0-_piI[in.days()]);
                        //cerr << in << "," << out << "=" << in_value <<", "<<out_value<<endl;
                        }
                } break;
                case SEIRState::R: {
                    out_value += in_value * (out.phase()==SEIRState::R ? 1.0 : 0.0);
                } break;
            }
            ++it_output_message;
        }
        ++it_input_message;
    }
    // cerr << "input_message" << &(*input_message) << endl;
    //cerr << "outgoing_message" << *outgoing_message << endl;
    return outgoing_message;
}


MessagePtr SEIRFactor::message_vertically( Node *n) {

    VirusLoad load;
    for( size_t i=2; i<_nodes.size(); i++) {
        if( _nodes[i] != n) {

            MessagePtr  contact_message = _nodes[i]->message_to(this);
            double p = 0, p_tot = 0;
            for( size_t i=0; i< contact_message->size(); i++) {
                p_tot += (*contact_message)[i];
                if( _states[i].phase() == SEIRState::I) 
                    p += (*contact_message)[i];
            }
            load.add_source( p/p_tot, 1.0);
        }
    }
    // if(_nodes.size()>2)
    //     cerr <<"load: " << load << endl;

    auto input_message = _nodes[0]->message_to(this);
    auto output_message = _nodes[1]->message_to(this);


    Message::const_iterator it_input_message;
    Message::const_iterator it_output_message;


    Message p_outgoing(2, 0.0);
    for(int j=0; j<2; j++) {
        it_input_message = input_message->begin();
        for( auto it_input_states = _states.cbegin(); it_input_states != _states.cend(); ++it_input_states) 
        {
            auto in = *it_input_states;

            it_output_message = output_message->begin();
            for( auto it_output_states = _states.cbegin(); it_output_states != _states.cend(); ++it_output_states) 
            {
                auto out = *it_output_states;

                const double& in_value = *it_input_message;
                const double &out_value = *it_output_message;
        
                switch( in.phase() ) {
                    case SEIRState::S: {
                        double p_keep = 0.0;
                        for( auto it_load = load.cbegin(); it_load != load.cend(); ++it_load) {
                            p_keep += (it_load->second) * pow(1.0-_p1, it_load->first + (double) j); 
                                                                                            // ^^^
                                                                                            // here is the load increase by the contact node
                        }
                        p_keep *= (1.0-_p0);
                        if( in.next(/*change = */ true)  == out) {
                            p_outgoing[j] += out_value * in_value * (1.0 - p_keep);
                            //cerr << in << "," << out << "=" << in_value <<", "<<out_value<<endl;
                            }
                        if( in.next(/*change = */ false) == out) {
                            p_outgoing[j] += out_value * in_value * p_keep;
                            //cerr << in << "," << out << "=" << in_value <<", "<<out_value<<endl;
                            }
                    } break;
                    case SEIRState::E: {
                        if( in.next(/*change = */ true)  == out) {
                            p_outgoing[j] += out_value * in_value * _piE[in.days()];
                            //cerr << in << "," << out << "=" << in_value <<", "<<out_value<<endl;
                            }
                        if( in.next(/*change = */ false) == out) {
                            p_outgoing[j] += out_value * in_value * (1.0-_piE[in.days()]);
                            //cerr << in << "," << out << "=" << in_value <<", "<<out_value<<endl;
                            }

                    } break;
                    case SEIRState::I: {
                        if( in.next(/*change = */ true)  == out) {
                            p_outgoing[j] += out_value * in_value * _piI[in.days()];
                            //cerr << in << "," << out << "=" << in_value <<", "<<out_value<<endl;
                            }
                        if( in.next(/*change = */ false) == out) {
                            p_outgoing[j] += out_value * in_value * (1.0-_piI[in.days()]);
                            //cerr << in << "," << out << "=" << in_value <<", "<<out_value<<endl;
                            }
                    } break;
                    case SEIRState::R: {
                        p_outgoing[j] += out_value * in_value * (out.phase()==SEIRState::R ? 1.0 : 0.0);
                    } break;
                }
                ++it_output_message;
            }
            ++it_input_message;
        }
    //cerr << p_outgoing << endl;
    } 

    p_outgoing = normalize(p_outgoing);
    MessagePtr outgoing_message( new Message( n->_N));
    for(int i=0; i<n->_N; i++) {
        (*outgoing_message)[i] = (_states[i].phase()==SEIRState::I ? p_outgoing[1] : p_outgoing[0]);
    }

    return outgoing_message;

}


MessagePtr SEIRFactor::message_to( Node *n) {

    if( n==_nodes[0]) return message_horizontally( /*forward =*/ false);
    if( n==_nodes[1]) return message_horizontally( /*forward =*/ true);
    return message_vertically( n);
}


double SEIRFactor::potential( const vector<unsigned int> & ) {

    assert(false); // SEIRFactor implements its own message_to() method
}


SEIRInitFactor::SEIRInitFactor( SEIRNode &out, bool patient_zero) :
    Factor({&out}), _patient_zero(patient_zero), _states(out.states())
{
}

double SEIRInitFactor::potential( const vector<unsigned int> &state_its) {

    assert(state_its.size()==1);

    SEIRState out = _states[state_its[0]];

    if( _patient_zero) 
        return out == SEIRState(SEIRState::E,1) ? 1.0 : 0.0;
    else 
        return out == SEIRState(SEIRState::S) ? 1.0 : 0.0;
}


SEIRTestFactor::SEIRTestFactor( SEIRNode &out, bool positive, double alpha, double beta) :
    Factor({&out}), _states(out.states())
{
    _positive = positive;
    _alpha = alpha;
    _beta = beta;
}

double SEIRTestFactor::potential( const vector<unsigned int> &state_its) {

    SEIRState out = _states[state_its[0]];

    if( _positive) 
        return out.phase() == SEIRState::I ? 1.0-_alpha : _beta;
    else 
        return out.phase() == SEIRState::I ? _alpha : 1.0-_beta;
}