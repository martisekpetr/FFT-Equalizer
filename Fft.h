//------------------------------------
// Petr Martisek
// FFT ekvalizer
// ADS II, cvicici Martin Mares
// 2012/13
//------------------------------------

#ifndef FFT_H
#define FFT_H

#include <vector>
#include <complex>
#include <tuple>

typedef std::vector<std::vector<std::complex<float>>> cplx2DVector;

class Fft
{
public:
    Fft (unsigned long numSamples);
    ~Fft (){}
    void    Transform ();
    void    Equalize(std::vector<std::tuple<int,int,float>> eqParam, int SampleRate);
    void    InverseTransform();
    void    CopyIn (const std::vector<long>& in);
    const std::vector<std::complex<float>> & GetResult(){return _inPlaceBuffer; }
private:
    unsigned long _numData; //aktualni velikost vstupnich dat
    unsigned long _numSamples; // nejblizsi vyssi mocnina dvojky - velikost bufferu, v nemz bude probihat vypocet
    unsigned long _numLevels; // pocet urovni FFT (dvojkovy logaritmus _numSamples)
    std::vector<unsigned long> _bitReverseMap; // mapa pro bitove reverzni setrideni
    std::vector<std::complex<float>> _inPlaceBuffer;  // in-place vektor, v nemz probiha transformace
    cplx2DVector _omega;             // vektor predpocitanych hodnot odmocnin z jednicky, tj. komplexnich exponencial pro vsechny urovne
};

#endif
