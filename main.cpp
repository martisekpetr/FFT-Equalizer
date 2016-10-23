//------------------------------------
// Petr Martisek
// FFT ekvalizer
// ADS II, cvicici Martin Mares
// 2012/13
//------------------------------------

#include <iostream>
#include "Fft.h"
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

using namespace std;

//prevod char* -> int
int si(char* s){
    stringstream ss(s);
    int n;
    ss >> n;
    return n;
}

//prevod char* -> float
float sd(char* s){
    stringstream ss(s);
    float n;
    ss >> n;
    return n;
}

void DisplayHelp(){
    cout << "Usage: " << endl
         << "ffteq [-i INPUT_FILE] [-o OUTPUT_FILE] [-e LOWER_BOUND UPPER_BOUND QUOTIENT] [-e LOWER_BOUND2 UPPER_BOUND2 QUOTIENT2] ..." << endl << endl
         << "where INPUT_FILE is *.wav file in PCM encoding, OUTPUT_FILE will be in the same format, "
         << "LOWER_BOUND and UPPER_BOUND are frequencies in Hertzs which define the frequency range on which the equalization will be done, "
         << "and QUOTIENT defines the amount of change, where 1.0 means no change, 0.5 means lower the presence of given frequencies to 50 % etc." <<endl;
}

int main(int, char** argv)
{
    string iFile = "";
    string oFile = "";
    vector<tuple<int,int,float>> eqParam;

    /* parametry vstupu:
     *  -o    vystupni soubor
     *  -h    napoveda
     *  -e    parametry ekvalizace ve formatu DOLNI_MEZ_FREKVENCE_V_Hz, HORNI_MEZ_FREKVENCE_V_Hz, KVOCIENT_ZMENY
     *
     * napr:
     * vstup: fft -i test.wav -o test_copy.wav -e 100 500 1.2 -e 2000 3000 0.2
     */

    //============================================
    // Parsovani vstupnich parametru
    //============================================

    if(!*++argv){
        cerr << "Input file not specified! See help (ffteq -h) or program documentation."  ;
        return 1;
    }
    else if(strcmp(*argv, "-h") == 0){
        DisplayHelp();
        return 0;
    }
    else
        iFile = *argv;

    while(*++argv && **argv == '-')
    {
        switch (argv[0][1]){
        case 'h':
            DisplayHelp();
            return 0;
        case 'o':
            if (*++argv)
            {oFile = (*argv);}
            else
                argv--;
            break;
        case 'e':
            if (*++argv && *++argv && *++argv)
            {
                int a = si(*(--(--argv)));
                int b = si(*++argv);
                float c = sd(*++argv);
                eqParam.push_back(make_tuple(a,b,c));
            }
            else
                --(--(--argv));
        default:
            break;
        }
    }

    //pojmenovani vystupniho souboru neni-li definovano uzivatelsky
    if(oFile == ""){
        size_t pos = iFile.find_last_of('.');
        oFile = iFile.substr(0,pos) + "_eq.wav";
    }

    //=======================================================
    // Precteni a naparsovani vstupniho *.wav souboru
    //=======================================================

    streampos fileSize;
    // je treba cist a zapisovat v binarnim modu, jinak muze dojit nechtene modifikaci dat operacnim systemem
    // napr. automaticke nahrazovani LF za LF+CR
    ifstream file(iFile, ios::binary);
    if(!file)
    {
        cerr << "Cannot open file.";
        return 1;
    }

    // ziskani velikosti souboru
    file.seekg(0, file.end);
    fileSize = file.tellg();
    file.seekg(0, file.beg);

    //alokujeme buffer, do nejz nakopirujeme data
    char * buffer = new char[fileSize];

    //precteni dat
    file.read(buffer, fileSize);

    //kontrola, zda se precetlo vse
    if(!file){
        cerr << "Error while reading file.";
        return 1;
    }
    file.close();

    //ziskavani informaci z hlavicky
    cout << "Chunk ID: ";
    char* part = new char;
    strncpy(part, buffer, 4);
    part[4] = '\0';
    if (strcmp(part, "RIFF") != 0) {
        cerr << "Not a *.wav file, cannot apply FFT.";
        return 1;
    };
    cout << part << endl;

    unsigned int ChunkSize =  *reinterpret_cast<unsigned int*>(buffer + 4) ;
    cout << "ChunkSize: " << ChunkSize << endl;
    cout << "File size: " << (float)(ChunkSize + 8) / (1024 * 1024) << " MB" << endl;
    if (ChunkSize+8 != fileSize) {
        cerr << "Corrupted file, size of the file does not correspond with the information in header.";
        return 1;
    }

    cout << "Format: ";
    strncpy(part, buffer + 8, 4);
    part[4] = '\0';
    if (strcmp(part, "WAVE") != 0) {
        cerr << "Not a *.wav file, cannot apply FFT.";
        return 1;
    };
    cout << part << endl;

    cout << "SubChunk1ID: ";
    strncpy(part, buffer + 12, 4);
    part[4] = '\0';
    cout << part << endl;

    unsigned long Subchunk1Size = *reinterpret_cast<unsigned long*>( buffer + 16 );
    cout << "Subchunk1Size: " << Subchunk1Size << endl;

    unsigned short AudioFormat = *reinterpret_cast<unsigned short*>( buffer + 20 );
    cout << "AudioFormat: " << AudioFormat << endl;
    if (AudioFormat != 1) {
        cerr << "Compressed file, cannot apply FFT." << endl;
        return 1;
    }

    unsigned short NumChannels = *reinterpret_cast<unsigned short*>( buffer + 22 );
    cout << "NumChannels: " << NumChannels << endl;

    unsigned long SampleRate = *reinterpret_cast<unsigned long*>( buffer + 24 );
    cout << "SampleRate: " << SampleRate << endl;

    unsigned long ByteRate = *reinterpret_cast<unsigned long*>( buffer + 28 );
    cout << "ByteRate: " << ByteRate << endl;

    unsigned short BlockAlign = *reinterpret_cast<unsigned short*>( buffer + 32 );
    cout << "BlockAlign: " << BlockAlign << endl;

    unsigned short BitsPerSample = *reinterpret_cast<unsigned short*>( buffer + 34 );
    cout << "BitsPerSample: " << BitsPerSample << endl;

    cout << "SubChunk2ID: ";
    strncpy(part, buffer + 20 + Subchunk1Size, 4);
    part[4] = '\0';
    cout << part << endl;

    unsigned long Subchunk2Size= *reinterpret_cast<unsigned long*>( buffer + 24 + Subchunk1Size );
    cout << "Subchunk2Size: " << Subchunk2Size << endl;

    vector<vector<signed long>> channels(NumChannels);

    for (unsigned long i = Subchunk1Size + 28; i < (ChunkSize + 8); i += BlockAlign)
    {
        unsigned int j = 0;

        for(auto it = channels.begin(); it!=channels.end(); it++,j += BitsPerSample/8){
            if(BitsPerSample <= 8)
                it->push_back(*reinterpret_cast<unsigned char*>(buffer + i + j));
            else if(BitsPerSample <= 16)
                it->push_back(*reinterpret_cast<signed short*>(buffer + i + j));
            else if(BitsPerSample <= 32)
                it->push_back(*reinterpret_cast<signed long*>(buffer + i + j));
            else{
                cerr << "Bit depth too large.";
                return 1;
            }
        }
    }


    ofstream out;
    out.open(oFile, ios::binary);
    if (out.is_open())
    {
        for (unsigned short i = 0; i < Subchunk1Size + 28; i++)
        {
            //prekopirovani hlavicky souboru
            out << buffer[i];
        }
        for (unsigned long j = 0; j < Subchunk2Size / BlockAlign; j++)
        {
            unsigned short i=0;
            for(auto it = channels.begin(); it != channels.end(); it++, i++)
            {
                channels.at(i).at(j) = (signed long)(*it).at(j).real();

                //pretypovani vysledku podle hodnoty BitsPerSample a pripadna zmena kodovani (Big -> Little Endian)
                if(BitsPerSample <= 8)
                    out << (char)(channels.at(i).at(j));
                else if(BitsPerSample <= 16){
                    signed short pom = (signed short)(channels.at(i).at(j));
                    out << (char)(pom & 0xFF) << (char)((pom & 0xFF00)>>8);
                }
                else if(BitsPerSample <= 32){
                    signed long pom = (signed long)(channels.at(i).at(j));
                    out << (char)(pom & 0xFF) << (char)((pom & 0xFF00)>>8) << (char)((pom & 0xFF0000) >> 16) << (char)((pom & 0xFF000000) >> 24);
                }
            }
        }
      }
    cout << "Equalization finished." << endl;
    out.close();
    delete [] buffer;
    return 0;
}

