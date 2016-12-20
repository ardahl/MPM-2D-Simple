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
        int getNumOptions();
        
        int getInt(std::string key, int index=0, int defaultVal=0);
        double getDouble(std::string key, int index=0, double defaultVal=0.0);
        char getChar(std::string key, int index=0, char defaultVal=0);
        std::string getString(std::string key, int index=0, std::string defaultVal="");
    
    private:
        std::map<std::string, std::vector<std::string> > configs;
};

#endif
