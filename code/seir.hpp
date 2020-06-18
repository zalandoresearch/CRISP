#ifndef _SEIR_HPP
#define _SEIR_HPP

#include <vector>
#include <map>
using namespace std;

class SEIRState {
public:
    enum Phase { S, E, I, R};

    SEIRState( Phase p, int d=0);
    SEIRState( const SEIRState &other);

    SEIRState next( bool change=false) const;
    bool operator== (const SEIRState &other) const;

    Phase phase() const { return _p;}
    int days() const { return _d;}

    static vector<SEIRState> all_states( int dE, int dI);

private:
    Phase _p;
    int _d;

};

ostream& operator<<(ostream& os, const SEIRState &s);




class VirusLoad {

    map<double, double> _Px;

public:
    VirusLoad(); 
    void add_source( double p, double x);

    map<double, double>::const_iterator cbegin() const {return _Px.cbegin();};
    map<double, double>::const_iterator cend() const {return _Px.cend();};

};

ostream& operator<<(ostream& os, const VirusLoad &l);


#endif