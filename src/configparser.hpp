#ifndef CONFIGPARSER_HPP_

#include <vector>
#include <string>
#include <sstream>
#include <map>

class ConfigParser {
    public:
        ConfigParser(std::string file, std::string delim=" ");
        
        bool exists(std::string key);
        int getNumParams(std::string key);
        
        int getInt(std::string key, int index, int defaultVal);
        double getDouble(std::string key, int index, double defaultVal);
        char getChar(std::string key, int index, char defaultVal);
        std::string getString(std::string key, int index, std::string defaultVal);
    
    private:
        int getIndex(std::string key);
        
        std::map<std::string, std::vector<std::string> > configs;
};

#endif
