#include <memory>
#include <vector>

#include <Pipeline.h>



void Pipeline::setup(std::vector<dict> dicts)
{
  int i = 0;
  
  for (auto& params : dicts)
  {
    if (i==0)
    {
      this->pipes[i]->setup(params);
    }
    else
    {
      this->pipes[i]->setup(params, this->pipes[i-1]);
    }
    i++;
  }
}

void Pipeline::start(void)
{
  int iter = 0;
  for (auto& pipe : this->pipes)
  {
    pipe->run();
  }
}