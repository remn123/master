#include <vector>
#include <string>
#include <map>

typedef std::map<const std::string, const void> dict;

/* PIPE */
class Pipe
{
public:
    dict params;
    std::string message;
public:
    Pipe(/* args */);
    virtual ~Pipe();

    void set(dict params);
    void get(void);
    void run(void);
};

/* PREPROC */
class Preproc : public Pipe
{
private:
    
public:
    Preproc();
    ~Preproc();
};

/* SOLVER */
template<typename T>
class Solver : public Pipe, public Method<T>
{
private:
    /* data */
public:
    Solver();
    ~Solver();
};

/* POSPROC */
class Posproc : public Pipe
{
private:
    /* data */
public:
    Posproc();
    ~Posproc();
};

/* PIPELINE */
class Pipeline
{
private:
    std::vector<Pipe> pipes;
    int verbose;
    std::string message;
public:
    Pipeline(std::vector<Pipe>);
    ~Pipeline();

    void set(void);
    void start(void);
    void log(const std::string&);
    void print(const std::string&);
};

void Pipeline::set(void)
{
    for (auto pipe : this->pipes)
    {
        if (pipe.params["verbose"]==1)
        {
            this->print(pipe.message);
        }
        else if(pipe.params["verbose"]==2)
        {
            this->log(pipe.message);
        }
        pipe.set();
    }
}

void Pipeline::start(void)
{
    for (auto pipe : this->pipes)
    {
        if (pipe.params["verbose"]==1)
        {
            this->print(pipe.message);
        }
        else if(pipe.params["verbose"]==2)
        {
            this->log(pipe.message);
        }
        pipe.run();
    }
}

void Pipeline::log(const std::string& message)
{
    //std::cout << ""
}
