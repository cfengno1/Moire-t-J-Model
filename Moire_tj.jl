using ITensors
using LinearAlgebra
BLAS.set_num_threads(1)
using ITensors.HDF5
ITensors.Strided.disable_threads()
ITensors.enable_threaded_blocksparse()
using StatsBase: sample

Base.Sys.set_process_title(string("#CORES=",Threads.nthreads()))
M=4000
Ny = 4
Nx = 36
phi= pi*7/12
inv_delta=12
sys=string("N",Ny,"x",Nx,"_delta",inv_delta,"_7pi12_M",M)
text=string(sys,"_node04")

delta=1/inv_delta
J1=1
t1=3
N = Nx*Ny
Ne=convert(Int64,N*(1-delta))
println()
@show sys
@show text
@show Threads.nthreads()

let
  sites = siteinds("tJ", N;
                   conserve_qns = true)

  lattice = triangular_lattice_Moire(Nx, Ny; yperiodic = true)

  # Define the Heisenberg spin Hamiltonian on this lattice
  ampo = OpSum()
  for b in lattice
    if b.type=="1" || b.type=="2"
      #x- or y-bond

      ampo .+= J1/2*exp(-2*im*phi), "S+", b.s1, "S-", b.s2
      ampo .+= J1/2*exp(2*im*phi), "S-", b.s1, "S+", b.s2
      ampo .+= J1,   "Sz", b.s1, "Sz", b.s2
      ampo .+= -0.25*J1, "Ntot", b.s1, "Ntot", b.s2
      ampo .+= -t1*exp(-im*phi), "Cdagup", b.s1, "Cup", b.s2
      ampo .+= -t1*exp(im*phi), "Cdagup", b.s2, "Cup", b.s1
      ampo .+= -t1*exp(im*phi), "Cdagdn", b.s1, "Cdn", b.s2
      ampo .+= -t1*exp(-im*phi), "Cdagdn", b.s2, "Cdn", b.s1
    else
      #diagonal bond

      ampo .+= J1/2*exp(2*im*phi), "S+", b.s1, "S-", b.s2
      ampo .+= J1/2*exp(-2*im*phi), "S-", b.s1, "S+", b.s2
      ampo .+= J1,   "Sz", b.s1, "Sz", b.s2
      ampo .+= -0.25*J1, "Ntot", b.s1, "Ntot", b.s2
      ampo .+= -t1*exp(im*phi), "Cdagup", b.s1, "Cup", b.s2
      ampo .+= -t1*exp(-im*phi), "Cdagup", b.s2, "Cup", b.s1
      ampo .+= -t1*exp(-im*phi), "Cdagdn", b.s1, "Cdn", b.s2
      ampo .+= -t1*exp(im*phi), "Cdagdn", b.s2, "Cdn", b.s1
    end
  end

  H = MPO(ampo,sites)
  H = splitblocks(linkinds, H)

  f=h5open("H_MPO.h5","w")
  write(f,"H",H)
  close(f)

  ### Choose a initial state with randomly distributed electron

  sites_emp=inv_delta:inv_delta:N
	sites_Ne=setdiff(1:N,sites_emp)
	sites_up=sites_Ne[1:2:end]
	sites_dn=sites_Ne[2:2:end]

  state=String[]
  for i in 1:N
  	if i in sites_up
			state = vcat(state,"Up")
		elseif i in sites_dn
			state = vcat(state,"Dn")
		else
			state = vcat(state,"Emp")
		end
	end

  psi0 = randomMPS(sites,state,600)

  sweeps = Sweeps(36)
  maxdim!(sweeps,800, 800, 800, 800,1000,1000,1000,1000,1000,1000,1000,1000,1000,2000,2000,2000,2000,2000,2000,2000,2000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,4000,4000,4000,4000,4000,4000,4000,4000)
  cutoff!(sweeps,1E-8)
  noise!(sweeps,1E-5,1E-5,1E-5,1E-5,1E-5,1E-5,1E-5,1E-5,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-9,1E-9,1E-9,1E-9,0)

  @show sweeps

  file=string("psi_",sys)
  energy,psi = dmrg(H,psi0,sweeps; write_when_maxdim_exceeds=4000, write_step=file)
  mv(string(file,"_tmp.h5"),string(file,".h5");force=true)

  include("send_email.jl")

  using Plots
  NN = expect(psi,"Ntot";sites = 1:Ny:N)
  open(string("Ne_",sys,".dat"), "w") do io
      for i in 1:1
      	for j in 1:Nx
          println(io, i,"\t",j,"\t",real(NN[j]))
        end
      end
  end
  FIG=plot(1:Nx,NN);
  savefig(FIG,string("Ne_",sys,".pdf"))

  return
end
