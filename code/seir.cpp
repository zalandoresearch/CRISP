#include"seir.hpp"

#include <vector>
#include <iostream>
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

vector<SEIRState> SEIRState::all_states( int dE, int dI) {

     vector<SEIRState> res;

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
