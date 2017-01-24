#include "masm.h"

int main(int argc,char** argv)
{
  time_t t1,t2;
  MASM* m;
#if DBASE
  m = new MASM;
  m->set_jobid(boost::lexical_cast<int>(argv[1]));
  m->retrieve_db();
#else
  m = new MASM(argv[1]);
#endif

  t1 = std::time(NULL);
  m->run();
  t2 = std::time(NULL);
  std::cout << "The total runtime for this version of Werkzeug is " << t2 - t1 << " seconds." << std::endl;

  delete m;

  return 0;
}

