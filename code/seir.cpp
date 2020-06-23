#include"seir.hpp"

#include <vector>
#include <iostream>
#include <map>
using namespace std;



SEIRState::SEIRState( Phase p, int d) : _p(p), _d(d) {};
SEIRState::SEIRState( const SEIRState &other) : _p(other._p), _d(other._d) {};

SEIRState SEIRState::next( bool change) const {
    if(change) {
        switch(_p) {
            case S: return SEIRState(E, 1);
            case E: return SEIRState(I, 1);
            case I: return SEIRState(R);
            case R: throw std::invalid_argument( "Pahse R cannot be changed from!" );
        }
    }
    else {
        switch(_p) {
            case S: return SEIRState(S);
            case E: return SEIRState(E, _d+1);
            case I: return SEIRState(I, _d+1);
            case R: return SEIRState(R);
        }
    }
}


bool SEIRState::operator == (const SEIRState &other) const {
    if(_p==S && other._p==S) return true;
    if(_p==R && other._p==R) return true;
    if(_d != other._d) return false;
    if(_p==E && other._p==E) return true;
    if(_p==I && other._p==I) return true;

    return false;
}

bool SEIRState::operator == (SEIRState::Phase other) const {
    return _p == other;
}

SEIRStateSpace SEIRState::all_states( int dE, int dI) {

    SEIRStateSpace res;

    res.push_back( SEIRState(S));
    for( int d=1; d<=dE; res.push_back( SEIRState(E,d++)));
    for( int d=1; d<=dI; res.push_back( SEIRState(I,d++)));
    res.push_back( SEIRState(R));

    return res;
}

ostream& operator<<(ostream& os, const SEIRState &s) {
    switch(s.phase()) {
        case SEIRState::S: os << "S"; break;
        case SEIRState::E: os << "E(" << s.days() << ")"; break;
        case SEIRState::I: os << "I(" << s.days() << ")"; break;
        case SEIRState::R: os << "R"; break;
    }
    return os;
}


Message basic_states( const Message &message, const SEIRStateSpace &states) {

    assert(message.size() == states.size());

    Message out(4);

    for( size_t i=0; i<message.size(); i++) {
        out[states[i].phase()] += message[i]; 
    }
    return out;
}



VirusLoad::VirusLoad() {
    _Px[0] = 1.0;
}

void VirusLoad::add_source( double p, double x) {
    
    auto Px(_Px);

    for( auto it = _Px.begin(); it!=_Px.end(); ++it)
        it->second *= 1.0-p;

    for(auto it = Px.begin(); it!=Px.end(); ++it )
        _Px[it->first+x] += p * it->second;

}

ostream& operator<<(ostream& os, const VirusLoad &l) {
    double p=0;
    os << "[";
    for( auto it=l.cbegin(); it!=l.cend(); ++it) {
        os << "(" <<it->first <<": "<< it->second << ") ";
        p += it->second;
    }
    os << "], p = " << p; 
    return os;
}
