#include <memory>
#include <vector>

#include <Pipeline.h>

void Pipeline::setup(std::vector<std::shared_ptr<dict>> dicts)
{
  int i = 0;
  for (auto& params : dicts)
  {
    this->pipes[i]->setup(*params);
    i++;
  }
}

void Pipeline::start(void)
{
  for (auto& pipe : this->pipes)
  {
    pipe->run();
  }
}