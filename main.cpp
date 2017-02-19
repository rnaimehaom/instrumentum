#include "molecular_assembler.h"

int main(int argc,char** argv)
{
  if (argc != 2) {
    std::cout << "Usage: ./instrumentum parameter_file" << std::endl;
    return 0;
  }
  boost::timer::auto_cpu_timer timer(3);

  Molecular_Assembler masm(argv[1]);
  masm.run();

  return 0;
}

