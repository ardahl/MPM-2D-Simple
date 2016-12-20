#include "configparser.hpp"

ConfigParser::ConfigParser(std::string file, char delim) {
    std::ifstream scn;
    scn.open(config);
    if(scn.fail()) {
        printf("File could not be opened: %s\n", config.c_str());
        exit(0);
    }
    
    std::string str, tok;
    while(std:getline(scn, str)) {
        //Ignore empty lines
        if(str.empty()) {
            continue;
        }
        //Ignore Comments
        else if(str[0] == '#') {
            scn.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }
        
        //Tokenize
        std::vector<std::string> options;
        char *p = strtok(str.c_str(), delim);
        //param
        std::string param = std::string(p);
        p = strtok(NULL, delim);
        //values
        while (p) {
            std::string v(p);
            if(!v.empty()) {
                options.push_back(std::string(p));
            }
            p = strtok(NULL, delim);
        }
        configs[param] = options;
    }
}

bool ConfigParser::exists(std::string key) {
    return configs.find(key) != configs.end();
}

int getNumParams(std::string key) {
    std::vector<std::string>::iterator it = configs.find(key);
    if(it == configs.end()) {
        return -1;
    }
    return it->second.size();
}
        
int ConfigParser::getInt(std::string key, int index, int defaultVal) {
    std::vector<std::string>::iterator it = configs.find(key);
    if(it == configs.end()) {
        return defaultVal;
    }
    std::vector<std::string> vals = it->second;
    return std::stoi(vals[index]);
}

double ConfigParser::getDouble(std::string key, int index, double defaultVal) {
    std::vector<std::string>::iterator it = configs.find(key);
    if(it == configs.end()) {
        return defaultVal;
    }
    std::vector<std::string> vals = it->second;
    return std::stod(vals[index]);
}

char ConfigParser::getChar(std::string key, int index, char defaultVal) {
    std::vector<std::string>::iterator it = configs.find(key);
    if(it == configs.end()) {
        return defaultVal;
    }
    std::vector<std::string> vals = it->second;
    return vals[index][0];
}

std::string ConfigParser::getString(std::string key, int index, std::string defaultVal) {
    std::vector<std::string>::iterator it = configs.find(key);
    if(it == configs.end()) {
        return defaultVal;
    }
    std::vector<std::string> vals = it->second;
    return vals[index];
}
