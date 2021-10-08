#include "molecular_assembler.h"

int main(int argc,char** argv)
{
  if (argc != 2) {
    std::cerr << "Usage: ./instrumentum <parameter file>" << std::endl;
    return 1;
  }

  std::string filename(argv[1]);
  Molecular_Assembler masm(filename);
  masm.assemble();

  return 0;
}

