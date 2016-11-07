#ifndef num_to_str_h
#define num_to_str_h

TString num_to_str(int value)
{
 
   TString str;
   str.Form("%d",value);
   return str;
   
}

TString num_to_str(double value, int digis_no)
{

   TString str, digis_no_str;

   char a[10] = "%";
   digis_no_str.Form("%s.%df", a, digis_no);
   std::cout << "amin: " << digis_no_str << " " << value << std::endl;

   str.Form(digis_no_str,value);
   std::cout << str << std::endl;
   return str;
		          
}

TString num_to_str(float value, int digis_no)
{

   TString str, digis_no_str;

   char a[10] = "%";
   digis_no_str.Form("%s.%df", a, digis_no);
   std::cout << "amin: " << digis_no_str << std::endl;

   str.Form(digis_no_str,value);
   return str;
   
}

#endif
