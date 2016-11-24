#include <THtml.h>
#include <TSystem.h>

int main()
{
  gSystem->Load("../rootclasses/.libs/libStrangecalc");

  THtml h;
  h.SetInputDir("../rootclasses");
  h.SetOutputDir("../htmldoc/");
  h.SetEtcDir("/home/tom/Install/root/etc/html");

  h.SetProductName("strangecalc-wrapper");

  h.MakeIndex();
  h.MakeAll();

}
