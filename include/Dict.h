#pragma once

#include <map>
#include <string>

// class IValue
// {
// private:

// public:
//   IValue();
//   virtual ~IValue();

// };


// template <typename T>
// class Value : public IValue
// {
// private:
//   T value;
// public:
//   Value(T val) : value(val){};
//   ~Value();

//   T get (void)
//   {
//     return this->value;
//   }
// };


typedef std::map<const std::string, const std::string> dict;


dict split_params(const dict& params, const std::string& pattern)
{
  dict splitted_dict;
  std::size_t found;

  std::string key, val;
  for( auto const& p : params )
  {
    key = p.first;
    val = p.second;
    found = key.find(pattern);
    if (found!=std::string::npos)
    {
      splitted_dict.emplace(key.substr(found+pattern.length()), val);
    }
  }
  return splitted_dict;
}