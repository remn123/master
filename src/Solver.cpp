#include <Dict.h>
#include <Solver.h>

// ------------------------------------------------------------ //
// ---------------------------- Solver ------------------------ //
// ------------------------------------------------------------ //
template <template <typename, typename...> class SpaceClass,
          template <typename> class TimeClass,
          typename SpaceType, typename TimeType, typename... Args>
Solver<SpaceClass, TimeClass, SpaceType, TimeType, Args...>::Solver(dict& params)
{
  this->space_params = split_params(params, "s_");
  this->time_params  = split_params(params, "t_");
  this->params       = split_params(params, "d_");
}

template <template <typename, typename...> class SpaceClass,
          template <typename> class TimeClass,
          typename SpaceType, typename TimeType, typename... Args>
Solver<SpaceClass, TimeClass, SpaceType, TimeType, Args...>::~Solver()
{

}

template <template <typename, typename...> class SpaceClass,
          template <typename> class TimeClass,
          typename SpaceType, typename TimeType, typename... Args>
void Solver<SpaceClass, TimeClass, SpaceType, TimeType, Args...>::setup(void)
{
  this->space = SpaceClass<SpaceType, Args..> {this->space_params}; 
  this->time  = TimeClass<TimeType> {this->time_params};
}

template <template <typename, typename...> class SpaceClass,
          template <typename> class TimeClass,
          typename SpaceType, typename TimeType, typename... Args>
void Solver<SpaceClass, TimeClass, SpaceType, TimeType, Args...>::run(void)
{
  // pass
}
// ------------------------------------------------------------ //
