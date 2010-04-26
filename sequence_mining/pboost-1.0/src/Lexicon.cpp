
#ifndef	LEXICON_H
#define	LEXICON_H


#include <iostream>
#include <map>
#include <sstream>
#include <fstream>

using namespace std;

class Lexicon {
    
public:
    Lexicon::Lexicon(const char* filename, const char* replacements_filename="replacements.txt") {
        ifstream inFile (filename);
        cerr << "Opening keyfile " << filename << "\n";
        if (!inFile) {
            cerr << "Cannot open keyfile.  Aborting" << endl;
            return;
        }
        

        
        
        while (!inFile.eof()) {
            for (string line; getline(inFile,line); ) {
                int id = -1;
                if (line.size() == 0) { continue; }
                
                istringstream ss(line);
                for (string field; getline(ss,field,'\t'); ) {
                    if (id == -1) {
                        id = atoi(field.c_str());
                        //cout << field << endl;
                        if ((unsigned int)id > max_id) { max_id = id; }
                    } else {
                        //cout << field << endl;
                        lexicon[(unsigned int)id] = field;
                    }
                }
            }
        }
        // read in a map of the elements to replace
        ifstream rFile (replacements_filename);
        cerr << "Opening replacement file " << replacements_filename << endl;
        if (!rFile) {
            cerr << "Cannot open replacements file.  Aborting" << endl;
            return;
        }
    

        while (!rFile.eof()) {
            for (string line; getline(rFile,line); ) {
                if (line.size() < 2) { continue; }
                line.erase(line.length()-1); // remove trailing char
                istringstream ss(line);
                bool first = true;
                unsigned int id = 0;
                for (string field; getline(ss,field,'\t'); ) {
                    if (field[0] == '#') {
                        continue;
                    } else {
                        
                        if (first) {
                            first = false;
                            id = atoi(field.c_str());
                        } else {
                            replacements[id] = atoi(field.c_str());
                        }
                    }
                    
                }
            }
        }
        // store max data value
        data_max = max_id;
        inFile.close();
    }
    
    void Lexicon::print() {
        for (unsigned int i=0; i < max_id; i++) {
            cout << i << "    " << lexicon[i] << endl;
        }
    }    
    
    string Lexicon::getString(unsigned int i) {
         map<unsigned int, string>::iterator it = lexicon.find(i);
         if (it == lexicon.end()) {
             return string("not found");
         } else {
             return it->second;
         } 
    }
    
    
    unsigned int Lexicon::getId( unsigned int inID) {
        //cerr << "Getting ID" << inID;
        map<unsigned int, unsigned int>::iterator it = replacements.find(inID);
        if (it == replacements.end()) 
            return inID;
        else 
            return it->second;
    }
    
    private:
        map<unsigned int,string> lexicon;
        map<unsigned int, unsigned int> replacements;
        unsigned int max_id;
        unsigned int data_max;
    
};

#endif

