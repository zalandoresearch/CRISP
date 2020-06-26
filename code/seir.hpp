#ifndef _SEIR_HPP
#define _SEIR_HPP


#include <vector>
#include <map>
using namespace std;




class SEIRState {
public:
    enum Phase { S=0, E=1, I=2, R=3};

    SEIRState( Phase p, int d=0);
    SEIRState( const SEIRState &other);

    SEIRState next( bool change=false) const;
    bool operator== (const SEIRState &other) const;
    bool operator== (Phase other) const;

    Phase phase() const { return _p;}
    int days() const { return _d;}

private:
    Phase _p;
    int _d;

};

ostream& operator<<(ostream& os, const SEIRState &s);

class SEIRStateSpace: public vector<SEIRState> {
    SEIRStateSpace( const SEIRStateSpace&) = delete; //avoid copying this object, one instance should always be sufficient

    vector<SEIRState> _states;
public:
    const int dEMax;
    const int dIMax;

    SEIRStateSpace( int dE, int dI);
    int operator[] (const SEIRState &s) const;
    const SEIRState &operator[] (int i) const {return _states[i];};
};







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