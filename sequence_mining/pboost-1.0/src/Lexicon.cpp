

#include <iostream>
#include <map>
#include <sstream>
class Lexicon {
    
public:
    Lexicon::Lexicon(const char* filename) {
        std::ifstream inFile (filename);
          
        while (!inFile.eof()) {
            for (std::string line; std::getline(inFile,line); ) {
                int id = -1;
                if (line.size() == 0) { continue; }
                
                std::istringstream ss(line);
                for (std::string field; std::getline(ss,field,'\t'); ) {
                    if (id == -1) {
                        id = atoi(field.c_str());
                        std::cout << field << std::endl;
                        if (id > max_id) { max_id = id; }
                    } else {
                        std::cout << field << std::endl;
                        field.erase(field.length()-1); // remove trailing char
                        lexicon[id] = field;
                    }
                }
            }
            
        }
        
        inFile.close();
        for (int i=0; i < max_id; i++) {
            std::cout << i << "    " << lexicon[i] << std::endl;
        }
    }
    
    private:
        std::map<int,std::string> lexicon;  
        unsigned int max_id;
    
};



