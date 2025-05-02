function mykspsetup(ksp)
  pc       = Ref{GridapPETSc.PETSC.PC}()
  mumpsmat = Ref{GridapPETSc.PETSC.Mat}()
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPPREONLY)
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCLU)
  @check_error_code GridapPETSc.PETSC.PCFactorSetMatSolverType(pc[],GridapPETSc.PETSC.MATSOLVERMUMPS)
  @check_error_code GridapPETSc.PETSC.PCFactorSetUpMatSolverType(pc[])
  @check_error_code GridapPETSc.PETSC.PCFactorGetMatrix(pc[],mumpsmat)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  4, 0)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 28, 2)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 29, 2)
  @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
end

function check_solver_accuracy(A,x,b)
  nr = norm(A*x - b)
  nb = norm(b)
  nx = norm(x)
  # @show nr, nr/nb, nr/nx
  tol_warn = 1.0e-12
  if nr > tol_warn && nr/nb > tol_warn && nr/nx > tol_warn
    @warn "Solver not accurate"
    @show nr, nr/nb, nr/nx
  end
end

function _solve(op,ps)
  A = get_matrix(op)
  b = get_vector(op)
  x = similar(b)
  ss = symbolic_setup(ps,A)
  ns = numerical_setup(ss,A)
  solve!(x,ns,b)
  check_solver_accuracy(A,x,b)
  x
end

function solve_stokes(op,X,ps)
  x = _solve(op,ps)
  u¹,p¹,_ = FEFunction(X,x)
  u¹,p¹
end

function _assemble_problem(a,b,assem,X,Y,A)
  dv = get_fe_basis(Y)
  du = get_trial_fe_basis(X)
  if A === nothing
    assemble_matrix_and_vector(a(du,dv),b(dv),assem,X,Y)
  else
    A,assemble_vector(b(dv),assem,Y)
  end
end

function _solve_problem(A,b,Xᵘ,ps)
  x = similar(b)
  ss = symbolic_setup(ps, A)
  ns = numerical_setup(ss, A)
  solve!(x,ns,b)
  check_solver_accuracy(A,x,b)
  FEFunction(Xᵘ,x)
end