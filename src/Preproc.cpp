#include <Preproc.h>

// ------------------------------------------------------------ //
// --------------------------- Preproc ------------------------ //
// ------------------------------------------------------------ //
template <typename MeshClass, typename... Args>
Preproc<MeshClass, Args...>::Preproc(dict& params)
{
  this->params = params;
}

template <typename MeshClass, typename... Args>
Preproc<MeshClass, Args...>::~Preproc()
{

}

template <typename MeshClass>
void Preproc<MeshClass>::setup(void)
{
  this->mesh = MeshClass {this->params}; 
}

template <typename MeshClass, typename... Args>
void Preproc<MeshClass, Args...>::setup(void)
{
  this->mesh = MeshClass<Args..> {this->params}; 
}

template <typename MeshClass, typename... Args>
void Preproc<MeshClass, Args...>::run(void)
{
  // pass
}
// ------------------------------------------------------------ //
