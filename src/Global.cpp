
#include "Global.h"

std::string intToString(const int& number, const int& digits) {
    std::string format;
    if(digits <= 0){
        format = "%d";
    }else{
        char b[10];
        sprintf(b,"%d",digits);
        format = b;
        format = "%0"+format+"d";
    }
    char buffer[100];
    sprintf(buffer, format.c_str(), number);
    return buffer;
}

std::string doubleToString(const double& number, const int& after) {
    char buffer[100];
    switch (after) {
        case 0: sprintf(buffer, "%f", number);
            break;
        case 1: sprintf(buffer, "%.1f", number);
            break;
        case 2: sprintf(buffer, "%.2f", number);
            break;
        case 3: sprintf(buffer, "%.3f", number);
            break;
        case 4: sprintf(buffer, "%.4f", number);
            break;
        case 5: sprintf(buffer, "%.5f", number);
            break;
        case 6: sprintf(buffer, "%.6f", number);
            break;
        case 7: sprintf(buffer, "%.7f", number);
            break;
        case 8: sprintf(buffer, "%.8f", number);
            break;
        default: sprintf(buffer, "%.9f", number);
            break;
    }
    return buffer;
}

std::string MTToString(MonomerType a){
    std::string result = "";
    switch(a){
        case W: result = "W"; break;
        case H: result = "H"; break;
        case P: result = "P"; break;
    }
    return result;
}

void getOutputStream(std::ofstream& fout, std::string filename, bool app) {
    filename = OUTPUT_DIR + filename;
    if(app){
        fout.open(filename.c_str(), std::ios::app);
    }else{
        fout.open(filename.c_str());
    }
    if (!fout.is_open()) {
        std::cout << "Output file: " << filename << " open failed!" << std::endl;
        exit(1);
    }
}

void getInputStream(std::ifstream& fin, std::string filename) {
    filename = INPUT_DIR + filename;
    fin.open(filename.c_str());
    if (!fin.is_open()) {
        std::cout << "Input file: " << filename << " open failed!" << std::endl;
        exit(1);
    }
}


