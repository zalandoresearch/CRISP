
#include "factor.hpp"
#include "node.hpp"
#include "seir.hpp"

#include <vector>
#include <iostream>
#include <cmath>
#include <any>
#include <cstdint>
#include <functional>
#include <memory>

using namespace std;


Factor::Factor( const vector<Node*> &nodes) : 
    _nodes(nodes) 
{
    for( auto n: _nodes) {
        n->_factors.push_back(this);
        n->_messages.emplace_back( new Message(n->_N, 1.0) );
    }
}

MessagePtr Factor::message_to( Node* n) {

    size_t i_node = 0; 
    for( ; i_node<_nodes.size() && _nodes[i_node]!=n; i_node++){
    }
    assert( i_node<_nodes.size()); // n needs to be in _nodes


    vector< MessagePtr>  incoming_messages;
    vector< vector<std::any>> input_states;
    vector< Message::const_iterator>   msg_its;
    vector< vector<std::any>::const_iterator> state_its;

    for( unsigned int i=0; i!=_nodes.size(); i++ ) {
        incoming_messages.push_back(_nodes[i]->message_to(this));
        msg_its.push_back(incoming_messages[i]->begin());

        input_states.push_back(_nodes[i]->_states);
        state_its.push_back( input_states[i].begin());
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
                state_its[i] = input_states[i].begin();  
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

double TableFactor::potential( vector<vector<any>::const_iterator> state_its) {
    assert( state_its.size()== _nodes.size());

    int idx = 0;
    for( unsigned int i=0; i<_nodes.size(); i++) {
        if( i>0) idx *= _nodes[i-1]->size();
        idx += std::any_cast<int>(*(state_its[i]));
        //cout << _nodes[i]->size() << " " <<  std::any_cast<int>(*(state_its[i])) << endl;
    }
    //cout << "ïdx: " << idx << endl;
    return _tab[idx];
}


SEIRFactor::SEIRFactor( const Message &qE, const Message &qI, 
                        double p0, double p1, 
                        SEIRNode &in, SEIRNode &out, vector<SEIRNode *> contacts) : 
    Factor(init_helper(in, out, contacts)), /* _qE(qE), _qI(qI), */ _piE(init_pi(qE)), _piI(init_pi(qI)), _p0(p0), _p1(p1) 
{
}


vector<Node *> SEIRFactor::init_helper(SEIRNode &in, SEIRNode &out, vector<SEIRNode *> contacts ) {
    vector<Node *> nodes({ &in, &out});
    nodes.insert( nodes.end(), contacts.begin(), contacts.end());
    return nodes;
}


const vector<double> SEIRFactor::init_pi( const vector<double> q) {

    auto res = Message(q.size());
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


MessagePtr SEIRFactor::message_horizontally( bool forward) {

    const vector<std::any>& input_states = _nodes[0]->states();
    const vector<std::any> output_states = _nodes[1]->states();

    VirusLoad load;
    for( size_t i=2; i<_nodes.size(); i++) {
        MessagePtr  contact_message = _nodes[i]->message_to(this);
        const vector<std::any>& contact_states = _nodes[i]->states();
            double p = 0, p_tot = 0;
            for( size_t i=0; i< contact_message->size(); i++) {
                p_tot += (*contact_message)[i];
                if( std::any_cast<SEIRState>(contact_states[i]).phase() == SEIRState::I) 
                    p += (*contact_message)[i];
            }
            load.add_source( p/p_tot, 1.0);
       }
    //if(_nodes.size()>2)
    //    cerr <<"load: " << load << endl;

    MessagePtr outgoing_message( new Message(_nodes[ forward ? 1 : 0]->size(), 0.0));

    MessagePtr input_message;
    MessagePtr output_message;
    if( forward) {
        //cerr << "forward\n";
        input_message = _nodes[0]->message_to(this);
        output_message = outgoing_message;
    }
    else {
        //cerr << "backward\n";
        input_message = outgoing_message;
        output_message = _nodes[1]->message_to(this);
    }

    // cerr << "input_message" << input_message << endl;
    // cerr << "outgoing_message" << outgoing_message << endl;

    // cerr << "input_message" << *input_message << endl;
    // cerr << "output_message" << *output_message << endl;
    // cerr << "forward="<< forward << endl;

    Message::const_iterator it_input_message;
    Message::const_iterator it_output_message;

    it_input_message = input_message->begin();
    for( auto it_input_states = input_states.cbegin(); it_input_states != input_states.cend(); ++it_input_states) 
    {
        auto in = std::any_cast<SEIRState>(*it_input_states);

        it_output_message = output_message->begin();
        for( auto it_output_states = output_states.cbegin(); it_output_states != output_states.cend(); ++it_output_states) 
        {
            auto out = std::any_cast<SEIRState>(*it_output_states);

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
            const vector<std::any>& contact_states = _nodes[i]->states();
            double p = 0, p_tot = 0;
            for( size_t i=0; i< contact_message->size(); i++) {
                p_tot += (*contact_message)[i];
                if( std::any_cast<SEIRState>(contact_states[i]).phase() == SEIRState::I) 
                    p += (*contact_message)[i];
            }
            load.add_source( p/p_tot, 1.0);
        }
    }
    // if(_nodes.size()>2)
    //     cerr <<"load: " << load << endl;

    auto input_message = _nodes[0]->message_to(this);
    auto output_message = _nodes[1]->message_to(this);

    const vector<any>& input_states = _nodes[0]->_states;
    const vector<any>& output_states = _nodes[1]->_states;


    Message::const_iterator it_input_message;
    Message::const_iterator it_output_message;


    Message p_outgoing(2, 0.0);
    for(int j=0; j<2; j++) {
        it_input_message = input_message->begin();
        for( auto it_input_states = input_states.cbegin(); it_input_states != input_states.cend(); ++it_input_states) 
        {
            auto in = std::any_cast<SEIRState>(*it_input_states);

            it_output_message = output_message->begin();
            for( auto it_output_states = output_states.cbegin(); it_output_states != output_states.cend(); ++it_output_states) 
            {
                auto out = std::any_cast<SEIRState>(*it_output_states);

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
    cerr << p_outgoing << endl;
    } 

    p_outgoing = normalize(p_outgoing);
    MessagePtr outgoing_message( new Message( n->_N));
    for(int i=0; i<n->_N; i++) {
        (*outgoing_message)[i] = (any_cast<SEIRState>(_nodes[0]->_states[i]).phase()==SEIRState::I ? p_outgoing[1] : p_outgoing[0]);
    }

    return outgoing_message;

}


MessagePtr SEIRFactor::message_to( Node *n) {

    if( n==_nodes[0]) return message_horizontally( /*forward =*/ false);
    if( n==_nodes[1]) return message_horizontally( /*forward =*/ true);
    //return MessagePtr( new Message(53,1.0)); 
    return message_vertically( n);
}


double SEIRFactor::potential( vector<vector<any>::const_iterator> ) {

    assert(false); // SEIRFactor implements its own message_to() method
}


SEIRInitFactor::SEIRInitFactor( SEIRNode &out, bool patient_zero) :
    Factor({&out}), _patient_zero(patient_zero)
{
}

double SEIRInitFactor::potential( vector<vector<any>::const_iterator> state_its) {

    assert(state_its.size()==1);
    SEIRState out = std::any_cast<SEIRState>(*(state_its[0]));

    if( _patient_zero) 
        return out == SEIRState(SEIRState::E,1) ? 1.0 : 0.0;
    else 
        return out == SEIRState(SEIRState::S) ? 1.0 : 0.0;
}