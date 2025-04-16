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
  @check_error_code GridapPETSc.PETSC.MatMumpsSetCntl(mumpsmat[], 3, 1.0e-6)
  @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
end

check_solver_accuracy(A,x,b) = _check_solver_accuracy(A,x,b)

function check_solver_accuracy(A::PSparseMatrix,x,b)
  x = GridapDistributed.change_ghost(x,axes(A,2))
  _check_solver_accuracy(A,x,b)
end

function _check_solver_accuracy(A,x,b)
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

function solve_surface(op,X,ps)
  x = _solve(op,ps)
  u¹,_ = FEFunction(X,x)
  u¹
end

function solve_stokes(op,X,ps)
  x = _solve(op,ps)
  u¹,p¹,_ = FEFunction(X,x)
  u¹,p¹
end

function solve_stokes_stokes(op,X,ps)
  x = _solve(op,ps)
  u¹,p¹,_,u²,p²,_ = FEFunction(X,x)
  u¹,p¹,u²,p²
end

function solve_stokes_with_rigid_body(op,X,ps)
  x = _solve(op,ps)
  u¹,p¹,_,u² = FEFunction(X,x)
  u¹,p¹,u²
end