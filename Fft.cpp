//------------------------------------
// Petr Martisek
// FFT ekvalizer
// ADS II, cvicici Martin Mares
// 2012/13
//------------------------------------

#define _USE_MATH_DEFINES
#include "Fft.h"
#include <math.h>
#define PI (2.0 * asin(1.0))

Fft::Fft (unsigned long numSamples) : _numSamples (numSamples) // numSamples musi byt mocnina dvojky
{    
    _numData = 0;
    _inPlaceBuffer = std::vector<std::complex<float>>(_numSamples,std::complex<float>(0,0));

    // vypocet dvojkoveho logaritmu = pocet urovni FFT
    _numLevels = 0;
    int n = _numSamples;
    while (n != 1)
    {
        n /= 2;
        _numLevels++;
    }

    // priprava mapy pro bitove reverzni setrideni vzorku
    // napr. 000, 100, 010, 110, 001, 101, 011, 111 je bitove reverzne setridena posl.
    // postupne vhodne prepiname bity pomoci masky

    _bitReverseMap = std::vector<unsigned long>(_numSamples);
    int aktual = 0;
    int halfnumSamples = _numSamples/2;
    _bitReverseMap.at(0) = 0;

    for (unsigned long i = 1; i < _numSamples; i++)
    {
        int mask = halfnumSamples;
        // zakladni maska - 10000...
        while (aktual >= mask)
        {
            aktual -= mask; // vypneme bit dany maskou
            mask >>= 1;     // a posuneme masku (1000.. -> 0100...)
        }
        aktual += mask;
        _bitReverseMap.at(i) = aktual;
    }

    // predpocitame si komplexni exponencialy, tj. mocniny primitivni odmocniny z jednicky
    // pro inverzi pouziju tytez hodnoty, jen komplexne sdruzene
    // v praxi to znamena pridat pred sin znamenko minus

    _omega = cplx2DVector(_numLevels + 1);
    int denominator = 2;
    for (unsigned long level = 1; level <= _numLevels; level++)
    {
        _omega.at(level) = std::vector<std::complex<float>>(_numSamples);

        for ( unsigned long i = 0; i < _numSamples; i++ )
        {
            float re =  (float)cos (2 * M_PI * i / denominator);
            float im = (float)sin (2 * M_PI * i / denominator);
            _omega.at(level).at(i) = std::complex<float>(re, im);
        }
        denominator *= 2;
    }
}

void Fft::CopyIn (const std::vector<long>& in)
{
    //nakopirovani dat do bufferu, avsak uz vhodne (tj. reverzne bitove) setridenych
    auto iter = in.begin();
    _numData = (unsigned long)in.size();
    if (_numData > _numSamples)
        return;
    for (unsigned long i = 0; i < _numData; i++, iter++)
        _inPlaceBuffer[_bitReverseMap[i]] = std::complex<float>(*iter);
}

void Fft::Transform ()
{
    // samotna transformace pomoci nerekurzivniho in-place algoritmu
    // tj. pocitame hodnoty v hladinach FFT site postupne zleva doprava

    int step = 1; //sirka "butterfly" vypoctu, roztec mezi dvema prvky, jejichz linearni kombinace davaji hodnotu na dalsi urovni
    for (unsigned long level = 1; level <= _numLevels; level++)
    {
        int increm = step * 2;  // velikost "bloku" - graficky je to cast dane vrstvy, kde se "butterflies" prekryvaji
                                // po provedeni vsech vypoctu v jednom bloku je treba skocit prave o tuto hodnotu na dalsi blok
        for (int j = 0; j < step; j++)
        {
            std::complex<float> U = _omega.at(level).at(j);
            for (unsigned long i = j; i < _numSamples; i += increm)
            {
                std::complex<float> T = U;

                //in-place prepocet dalsi vrstvy
                T *= _inPlaceBuffer.at(i+step);
                _inPlaceBuffer.at(i+step) = _inPlaceBuffer.at(i) - T;
                _inPlaceBuffer.at(i) = _inPlaceBuffer.at(i) + T;

            }
        }
        step *= 2; //zdvojnasobeni sirky jednoho "butterfly" vypoctu
    }
}
void Fft::InverseTransform(){

    // data v _inPlaceBuffer je treba pred inverzni transformaci opet reverzne bitove setridit
    // pouzijeme pomocny vektor _tmpBuffer a abychom redukovali pocet kopirovani, bude transformace probihat v nem

    std::vector<std::complex<float>> _tmpBuffer(_numSamples, std::complex<float>(0,0));
    for (unsigned long i = 0; i < _numSamples; i++)
        _tmpBuffer.at(_bitReverseMap.at(i)) = _inPlaceBuffer.at(i);

    //jinak je princip transformace prakticky totozny, az na pouziti kompl. sdruzene odm. z jednicky a zaverecnou normalizaci
    int step = 1;
    for (unsigned long level = 1; level <= _numLevels; level++)
    {
        int increm = step * 2;
        for (int j = 0; j < step; j++)
        {
            std::complex<float> U =  std::complex<float>(_omega.at(level).at(j).real(), -_omega.at(level).at(j).imag()); //kompl. sdruzene cislo
            for (unsigned long i = j; i < _numSamples; i += increm)
            {
                std::complex<float> T = U;
                T *= _tmpBuffer.at(i+step);
                _tmpBuffer.at(i+step) = _tmpBuffer.at(i) - T;
                _tmpBuffer.at(i) = _tmpBuffer.at(i) + T;

            }
        }
        step *= 2;
    }

    //normalizace - vydeleni vsech hodnot _numSamples
    for(unsigned long j = 0; j < _numSamples; j++){
        _tmpBuffer.at(j) /= _numSamples;
        //soucasne hodnoty prekopirujeme do puvodniho bufferu
        _inPlaceBuffer.at(j) = _tmpBuffer.at(j);
    }
}

void Fft::Equalize(std::vector<std::tuple<int,int,float>> eqParam, int SampleRate){

    // druha polovina frekvencniho spektra je zrcadlovy obraz prvni, abychom zabranili sumum ve vysledku,
    // vynulujeme druhou polovinu spektra a prvky prvni vynasobime dvema (abychom vykompenzovali ubytek intenzity)
    for(unsigned long i = 1; i<_numSamples/2; i++)
        _inPlaceBuffer.at(i) *= 2;
    for(unsigned long i = _numSamples/2 + 1; i<_numSamples; i++){
        _inPlaceBuffer.at(i) = 0;
    }

    for(auto it = eqParam.begin(); it != eqParam.end(); it++){
        // prevod frekvencnich mezi na spravny index
        // spodni mez je zaokrouhlena nahoru
        unsigned long lower_bound = (unsigned long)(0.5 + (std::get<0>(*it)/(float)SampleRate *_numSamples));
        unsigned long upper_bound = (unsigned long)(std::get<1>(*it)/(float)SampleRate *_numSamples);

        //kontrola velikosti mezi, pripadna oprava na prijatelnou hodnotu
        if (lower_bound <= 0)
            lower_bound = 1;
        if (lower_bound > _numSamples/2)
            lower_bound = _numSamples/2;
        if (upper_bound <= 0)
            upper_bound = 1;
        if (upper_bound > _numSamples/2)
            upper_bound = _numSamples/2;

        //samotna ekvalizace
        for(unsigned long i = lower_bound; i <= upper_bound; i++){
            _inPlaceBuffer.at(i) *= std::get<2>(*it);
        }
    }
}
