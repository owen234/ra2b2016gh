#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

bool create_model_pars_data3()
{

    std::string filein_name = "outputfiles/data-chi2-fit-model-pars.txt";
    std::string filein2_name = "outputfiles/model-pars-qcdmc3.txt";
    std::string fileout_name = "outputfiles/model-pars-data3.txt";

    
    ifstream filein(filein_name);
    ofstream fileout(fileout_name);
    
    if ( !filein ) { cout << "Warning: file " << filein_name <<" doesn't exist" << endl; return 0; }
    if ( !fileout ) { cout << "Warning: cannot open/create file" << fileout_name << endl; return 0; }
    
    std::string line;
    while( std::getline(filein,line) )
    {
        
        std::size_t index = line.find("+/-");
        if (index == std::string::npos) { cout << "Warning: The structure of file " << filein_name << "is not as expected";filein.close();fileout.close(); return 0; }
        line.replace( index, 3, "   ");
        
        index = line.find("(");
        if (index == std::string::npos) { cout << "Warning: The structure of file " << filein_name << "is not as expected";filein.close();fileout.close(); return 0; }
        line.replace(line.find('('),line.find(')')-line.find('(')+1,"0.00");
        
        string val_str, name_str;       
        stringstream convert_temp(line);
        convert_temp >> name_str;
        convert_temp >> val_str;
        stringstream convert(val_str);
        double val;
        if ( !(convert >> val) )  { cout << val << "Warning: The structure of file " << filein_name << "is not as expected";filein.close();fileout.close(); return 0; }
        fileout << line << std::endl;
    }
   filein.close();


    ifstream filein2( filein2_name );
    
    if ( !filein2 ) { cout << "Warning: file " << filein2_name << " doesn't exist" << endl; return 0; }
    
    while( std::getline(filein2,line) )
    {
        
        std::string val_str, name_str;       
        stringstream convert_temp(line);
        convert_temp >> name_str;
        convert_temp >> val_str;
        stringstream convert(val_str);
        double val;
        if ( !(convert >> val) )  { cout << val << "Warning: The structure of file " << filein_name << "is not as expected";filein.close();fileout.close(); return 0; }
        if ( name_str.compare(5,3,"mht") == 0 || name_str.compare(5,2,"nb") == 0 )  
           fileout << line << std::endl;
    }
   filein.close();





   return 1;

}
