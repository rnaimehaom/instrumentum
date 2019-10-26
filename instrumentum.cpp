#include "molecular_assembler.h"

std::map<int,std::string> element_table;

int main(int argc,char** argv)
{
  if (argc != 2) {
    std::cout << "Usage: ./instrumentum parameter_file" << std::endl;
    return 0;
  }

  // Element names...
  element_table[1] =   "H";
  element_table[6] =   "C";
  element_table[7] =   "N";
  element_table[8] =   "O";
  element_table[9] =   "F";
  element_table[15] =  "P";
  element_table[16] =  "S";
  element_table[17] = "Cl";
  element_table[35] = "Br";
  element_table[47] = "Ag";
  element_table[53] =  "I";

  std::string filename(argv[1]);
  Molecular_Assembler masm(filename);
  masm.run();

  return 0;
}

