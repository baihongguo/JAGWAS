#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator> 
#include <numeric>
#include <string> 
#include <armadillo>
#include <cstdlib> 


#include </HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/pgamma.c>
#include </HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/lgamma.c>
#include </HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/fmax2.c>
#include </HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/mlutils.c>
//#include </HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/dpois.c>
#include </HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/dnorm.c>
#include </HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/pnorm.c>
#include </HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/ftrunc.c>
#include </HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/gamma.c>
#include </HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/lgammacor.c>
#include </HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/stirlerr.c>
#include </HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/bd0.c>
#include </HGCNT95FS/ADDLIE/work/Test/JAGWAS_C++/chebyshev.c>


void processFiles(const std::string& outputPath, const std::string& Cor_M, int nrow , double MAF, bool score_test, bool Beta_se , bool logP, const std::vector<std::string>& fileNames) {
 
    int nline = nrow;
    double AF;
    double AF1;
    arma::mat cor;
    std::string value;
    std::string Value;
    std::string line;
    std::string Line;
    std::string data;
    std::size_t dim = fileNames.size();
    std::vector<std::ifstream> file_handles;
    std::ofstream results(outputPath);
    std::ostringstream oss;  //buffer
    if (logP == true) {
    oss << "CHR" << "\t" << "SNP" << "\t" << "POS" << "\t" << "A1" << "\t" << "A2" << "\t" << "N" << "\t" << "AF1" << "\t" << "log10(P)"<< "\n";
    }
    else
    {
    oss << "CHR" << "\t" << "SNP" << "\t" << "POS" << "\t" << "A1" << "\t" << "A2" << "\t" << "N" << "\t" << "AF1" << "\t" << "P"<< "\n";
    }
    
    
    for (size_t f = 0; f < fileNames.size(); f++) {
        // std::ifstream currentFile = openFile(fileNames[f]);
        file_handles.emplace_back(fileNames[f]);
    }
   
    if (!cor.load(Cor_M, arma::auto_detect)) {
        std::cerr << "Error loading the matrix from file." << std::endl;
    }
    arma::mat cor_i = inv(cor);

    int p = 0;

//To skip the header
std::getline(file_handles[0], line);
for (size_t f = 1; f < fileNames.size(); f++) {
        std::getline(file_handles[f], Line);
    }



if (score_test == true) {

while (std::getline(file_handles[0], line)) {
  
         
        std::istringstream iss(line);
        std::vector<std::string> data;
        while (std::getline(iss, value, '\t')) data.push_back(value);
        AF = std::stod(data[4]);
    if (AF>MAF && AF<(1-MAF))
    {
        std::vector<double> zscores(dim);
        oss << "NA" << "\t" << data[0] << "\t" << "NA" << "\t" << data[1] << "\t" << data[2] << "\t" << data[3] << "\t" << data[4] << "\t";
        double score = std::stod(data[5]);
        double var = std::stod(data[6]);
        zscores[0] = (score/var)/sqrt(1/var);
// chunk into different fields, when reading, save fields to string,  

       for (size_t f = 1; f < fileNames.size(); f++) {
        std::string Value;
        std::string Line;
        std::getline(file_handles[f], Line);
        std::istringstream iss(Line);
        std::vector<std::string> data;
        while (std::getline(iss, Value, '\t')) data.push_back(Value);

        AF1 = std::stod(data[4]);
        if (AF1 != AF) {
        std::cerr << "Error: inconsistent variants!" << std::endl;
        }
       double score = std::stod(data[5]);
        double var = std::stod(data[6]);
        zscores[f] = (score/var)/sqrt(1/var);
    } 
   arma::vec Zscores = zscores;
   arma::mat z = arma::reshape(Zscores, 1, Zscores.size());
   //arma::mat zt = arma::reshape(Zscores, Zscores.size(), 1);
   arma::mat zt = trans(z);
   arma::mat chisq = z * cor_i * zt;
   double chisqs = chisq(0,0);
   double pValue = pgamma(chisqs/4, dim/2, 0.5, 0, logP);
   if (logP == true) {
   pValue = pValue*(1/log(10));      
   }
   oss << pValue << "\n";

   p++;
    
        if (p % nline == 0) {
            results << oss.str(); // dump oss to results, every n line 
            oss.str(std::string());
            oss.clear();
        }	
    }

    else {
        for (size_t f = 1; f < fileNames.size(); f++) {
        std::string Line;
        std::getline(file_handles[f], Line);
    }
    }
    // do the matrix production, get the chi-sq test statistics and output the results  
   

        }    
    results << oss.str();    //for the remainder
    oss.str(std::string());
    oss.clear();
    results.close();

}







else {

    if (Beta_se == false) {
    while (std::getline(file_handles[0], line)) {
        std::istringstream iss(line);
        std::vector<std::string> data;
        while (std::getline(iss, value, '\t')) data.push_back(value);
        AF = std::stod(data[6]);
    if (AF>MAF && AF<(1-MAF))
    {
        std::vector<double> zscores(dim);
        zscores[0] = std::stod(data[7]);
        oss << data[0] << "\t" << data[1] << "\t" << data[2] << "\t" << data[3] << "\t" << data[4] << "\t" << data[5] << "\t" << data[6] << "\t";
    
// chunk into different fields, when reading, save fields to string,  

       for (size_t f = 1; f < fileNames.size(); f++) {
        std::getline(file_handles[f], Line);
        std::istringstream iss(Line);
        std::vector<std::string> data;
        while (std::getline(iss, Value, '\t')) data.push_back(Value);
        zscores[f] = std::stod(data[7]); 
    }
  
    // do the matrix production, get the chi-sq test statistics and output the results  
   arma::vec Zscores = zscores;
   arma::mat z = arma::reshape(Zscores, 1, Zscores.size());
   //arma::mat zt = arma::reshape(Zscores, Zscores.size(), 1);
   arma::mat zt = trans(z);
   arma::mat chisq = z * cor * zt;
   double chisqs = chisq(0,0);
   double pValue = pgamma(chisqs/4, dim/2, 0.5, 0, logP);
   if (logP == true) {
   pValue = pValue*(1/log(10));      
   }
   oss << pValue << "\n";

   p++;
    
        if (p % nline == 0) {
            results << oss.str(); // dump oss to results, every nline 
            oss.str(std::string());
            oss.clear();
        }	
    }
    else {
        for (size_t f = 1; f < fileNames.size(); f++) {
        std::string Line;
        std::getline(file_handles[f], Line);
    }
         }
    }
    }
       
    




if (Beta_se == true) {
    while (std::getline(file_handles[0], line)) {
  
         
        std::istringstream iss(line);
        std::vector<std::string> data;
        while (std::getline(iss, value, '\t')) data.push_back(value);
        AF = std::stod(data[6]);
    if (AF>MAF && AF<(1-MAF))
    {
        std::vector<double> zscores(dim);
        oss << data[0] << "\t" << data[1] << "\t" << data[2] << "\t" << data[3] << "\t" << data[4] << "\t" << data[5] << "\t" << data[6] << "\t";
        double beta = std::stod(data[7]);
        double se = std::stod(data[8]);
        zscores[0] = beta/se;
// chunk into different fields, when reading, save fields to string,  

       for (size_t f = 1; f < fileNames.size(); f++) {
        std::string Value;
        std::string Line;
        std::getline(file_handles[f], Line);
        std::istringstream iss(Line);
        std::vector<std::string> data;
        while (std::getline(iss, Value, '\t')) data.push_back(Value);

        AF1 = std::stod(data[6]);
        if (AF1 != AF) {
        std::cerr << "Error: inconsistent variants!" << std::endl;
        }
        double beta = std::stod(data[7]);
        double se = std::stod(data[8]);
        zscores[f] = beta/se; 
    } 
   arma::vec Zscores = zscores;
   arma::mat z = arma::reshape(Zscores, 1, Zscores.size());
   //arma::mat zt = arma::reshape(Zscores, Zscores.size(), 1);
   arma::mat zt = trans(z);
   arma::mat chisq = z * cor_i * zt;
   double chisqs = chisq(0,0);
   double pValue = pgamma(chisqs/4, dim/2, 0.5, 0, logP);
   if (logP == true) {
   pValue = pValue*(1/log(10));      
   }
   oss << pValue << "\n";

   p++;
    
        if (p % nline == 0) {
            results << oss.str(); // dump oss to results, every n line 
            oss.str(std::string());
            oss.clear();
        }	
    }

    else {
        for (size_t f = 1; f < fileNames.size(); f++) {
        std::string Line;
        std::getline(file_handles[f], Line);
    }
    }
    // do the matrix production, get the chi-sq test statistics and output the results  
   

        }    

    }

    results << oss.str();    //for the remainder
    oss.str(std::string());
    oss.clear();
    results.close();
}
}

        



int main(int argc, char* argv[]) {
    
    std::string outputFilePath(argv[1]);
    std::string cor_m(argv[2]);
    int nrow = std::atoi(argv[3]);
    double MAF = std::stod(argv[4]);
    bool score_test = std::stoi(argv[5]) != 0;
    bool beta_se = std::stoi(argv[6]) != 0;
    bool logP = std::stoi(argv[7]) != 0;
    
    std::vector<std::string> fileNames;
    for (int i = 8; i < argc; ++i) {
        fileNames.push_back(argv[i]);
    }

    processFiles(outputFilePath, cor_m, nrow,  MAF, score_test, beta_se, logP, fileNames);
    return 0;

}


