#ifndef STOCHASTIC_ANALYTIC_CONTINUATION_STRINGS_H
#define STOCHASTIC_ANALYTIC_CONTINUATION_STRINGS_H
#pragma once

#include <vector>

void split(const std::string& s, std::vector<std::string>& tokens, const std::string& delimiters = " ") {
    std::string::size_type lastPos = s.find_first_not_of(delimiters, 0);
    std::string::size_type pos = s.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos) {
        tokens.emplace_back(s.substr(lastPos, pos - lastPos));
        lastPos = s.find_first_not_of(delimiters, pos);
        pos = s.find_first_of(delimiters, lastPos);
    }
}

double str2double(const std::string& str) {
    if (!str.empty()) {
        double output;
        std::stringstream sstream;
        sstream << str;
        sstream >> output;
        return output;
    }
    else { return 0;}
}

int str2int(const std::string& str) {
    if (!str.empty()) {
        int output;
        std::stringstream sstream;
        sstream << str;
        sstream >> output;
        return output;
    }
    else { return 0;}
}

#endif //STOCHASTIC_ANALYTIC_CONTINUATION_STRINGS_H