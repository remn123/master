/* 
  ******************************** SPECTRAL DIFFERENCE **********************************
  ***************************************************************************************

   * 1) Apply initial condition (Q)
   * 2) Find out flux (FP) and solution (SP) points on numerical space
   * 3) Interpolate solution (Q) from SP to FP
   * 4) Solve Riemann-Problem at each edge/face (Rusanov, Roe, etc) ... flux calculation
   * 5) Interpolate corrected flux solution from FP to SP
   * 6) Calculate Residue
   * 7) Iterate in time (explicit/implicit)

  ***************************************************************************************
*/

template<typename T>
concept Solverkind = requires(T a){
  {std::function_traits<a>::}
};


template<typename T>
concept Solverkind = requires(T a){
  {std::function_traits<a>::}
};

template<Solverkind T, Solverkind U>
class Solver : public Pipe, public Method<T>
{
public:

  Solver;
  ~Solver;

private:

}


/* 
 * Method<SD, ExplicitRungeKutta>
 * Solver<SD, ExplicitEuler>
 * */

class Solver : public Pipe
{
public:
  SpacialDiscretization spacial;
  TimeDiscretization time;
private:
  
public:
  Solver();
  ~Solver();

}


/* TimeDiscretization */
class TimeDiscretization
{
public:
  double delta_t;
  double CFL;

public:
  TimeDiscretization();
  virtual ~TimeDiscretization();

  virtual void setup(void);
  virtual void step(void);
}

class ExplicitEuler : public TimeDiscretization
{
public:

public:
  ExplicitEuler();
  ~ExplicitEuler();

  void setup(void);
  void step(void);
}



/* SpacialDiscretization */
class SpacialDiscretization
{
public:

public:
  SpacialDiscretization();
  virtual ~SpacialDiscretization();

  virtual void setup(void);
  virtual void residue(void);
}

class SD : public SpacialDiscretization
{
public:
  std::vector<Node> fnodes; // flux points (FP)
  std::vector<Node> snodes; // solution points (SP)

private:

public:
  SD();
  ~SD();

  void setup(void);
  void initial_condition(Element&);
  void calc_high_order_nodes(Element&);
  void interpolate_sp2fp(Element&);
  void riemann_solver(Element&);
  void interpolate_fp2sp(Element&);
  void residue(Element&);

}
