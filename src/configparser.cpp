#include "configparser.hpp"
#include <fstream>
#include <limits>
#include <cstring>

ConfigParser::ConfigParser(std::string file, std::string delim) {
    std::ifstream scn;
    scn.open(file);
    if(scn.fail()) {
        printf("File could not be opened: %s\n", file.c_str());
        exit(0);
    }
    
    std::string str, param, token;
    while(std::getline(scn, str)) {
        //Ignore empty lines
        if(str.empty()) {
            continue;
        }
        //Ignore Comments
        else if(str[0] == '#') {
            continue;
        }

        //Tokenize
        std::vector<std::string> options;
        size_t pos = 0;
        //Option name
        pos = str.find(delim);
        param = str.substr(0, pos);
        str.erase(0, pos + delim.length());
        //Option values
        while ((pos = str.find(delim)) != std::string::npos) {
            token = str.substr(0, pos);
            if(!token.empty()) {
                options.push_back(token);
            }
            str.erase(0, pos + delim.length());
        }
        //Loop skips last value. If string isn't empty there's something left
        if(!str.empty()) {
            options.push_back(str);
        }
        //Put in map
        configs[param] = options;
    }
}

bool ConfigParser::exists(std::string key) {
    return configs.find(key) != configs.end();
}

int ConfigParser::getNumParams(std::string key) {
    std::map<std::string, std::vector<std::string> >::iterator it = configs.find(key);
    if(it == configs.end()) {
        return -1;
    }
    return it->second.size();
}

int ConfigParser::getNumOptions() {
    return configs.size();
}
        
int ConfigParser::getInt(std::string key, int index, int defaultVal) {
    std::map<std::string, std::vector<std::string> >::iterator it = configs.find(key);
    if(it == configs.end()) {
        return defaultVal;
    }
    std::vector<std::string> vals = it->second;
    return std::stoi(vals[index]);
}

double ConfigParser::getDouble(std::string key, int index, double defaultVal) {
    std::map<std::string, std::vector<std::string> >::iterator it = configs.find(key);
    if(it == configs.end()) {
        return defaultVal;
    }
    std::vector<std::string> vals = it->second;
    return std::stod(vals[index]);
}

char ConfigParser::getChar(std::string key, int index, char defaultVal) {
    std::map<std::string, std::vector<std::string> >::iterator it = configs.find(key);
    if(it == configs.end()) {
        return defaultVal;
    }
    std::vector<std::string> vals = it->second;
    return vals[index][0];
}

std::string ConfigParser::getString(std::string key, int index, std::string defaultVal) {
    std::map<std::string, std::vector<std::string> >::iterator it = configs.find(key);
    if(it == configs.end()) {
        return defaultVal;
    }
    std::vector<std::string> vals = it->second;
    return vals[index];
}
