using Plots

function Measure_func(Nx::Int, Ny::Int, sys::String ; kwargs...)

	N=Nx*Ny
	ini=convert(Int64,8)
	shift = get(kwargs, :shift, 1)
  Occu = get(kwargs, :Occu, false)
  Sz = get(kwargs, :Sz, false)
  Gr = get(kwargs, :Gr, false)
  Szz = get(kwargs, :Szz, false)
  Sxy = get(kwargs, :Sxy, false)
  Sxy_M = get(kwargs, :Sxy_M, false)
  Szz_M = get(kwargs, :Szz_M, false)
  Dr = get(kwargs, :Dr, false)
  Paa = get(kwargs, :Paa, false)
  Pbb = get(kwargs, :Pbb, false)
  Pcc = get(kwargs, :Pcc, false)
  Pbb_shift = get(kwargs, :Pbb_shift, false)
  Pba = get(kwargs, :Pba, false)
  Pbc = get(kwargs, :Pbc, false)
  Ptud_bb = get(kwargs, :Ptud_bb, false)
  Ptuu_bb = get(kwargs, :Ptuu_bb, false)
  Ptuu_bb_shift = get(kwargs, :Ptuu_bb_shift, false)
  Ptuu_bb_y = get(kwargs, :Ptuu_bb_y, false)
  Ptdd_bb = get(kwargs, :Ptdd_bb, false)
  Ptdd_bb_y = get(kwargs, :Ptdd_bb_y, false)
  Ptdd_bb_shift = get(kwargs, :Ptdd_bb_shift, false)
  Ptud_aa = get(kwargs, :Ptud_aa, false)
  Ptuu_aa = get(kwargs, :Ptuu_aa, false)
  Ptdd_aa = get(kwargs, :Ptdd_aa, false)
  Ptud_cc = get(kwargs, :Ptud_cc, false)
  Ptuu_cc = get(kwargs, :Ptuu_cc, false)
  Ptud_ba = get(kwargs, :Ptud_ba, false)
  Ptuu_ba = get(kwargs, :Ptuu_ba, false)
  Ptud_bc = get(kwargs, :Ptud_bc, false)
  Ptuu_bc = get(kwargs, :Ptuu_bc, false)
  Ptud_ac = get(kwargs, :Ptud_ac, false)
  Ptuu_ac = get(kwargs, :Ptuu_ac, false)

  f = h5open(string("psi_",sys,".h5"),"r")
  psi=read(f,"psi",MPS)
  close(f)

  """
  ##############################################################################################################################
  Occu
  """
  if Occu
		Ne = expect(psi,"Ntot";sites = 1:Ny:N)
		open(string("Ne_",sys,".dat"), "w") do io
		    for i in 1:1
		    	for j in 1:Nx
		        println(io, i,"\t",j,"\t",real(Ne[j]))
		      end
		    end
		end
		FIG=plot(1:Nx,Ne);
		savefig(FIG,string("Ne_",sys,".pdf"))
  end
  	
  """
  ##############################################################################################################################
  Sz
  """
	if Sz
  	Sz = expect(psi,"Sz";sites = 1:Ny:N)
  	open(string("Sz_",sys,".dat"), "w") do io
  	    for i in 1:1
  	    	for j in 1:Nx
  	        println(io, i,"\t",j,"\t",real(Sz[j]))
  	      end
  	    end
  	end
	end
  
  """
  ##############################################################################################################################
  Gr
  """
	if Gr
  	open(string("Gr_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  	    for j=1:1
  	        range_x= ini*Ny+j
  	        range_y= collect((ini+1):(Nx-ini)).*Ny .+ j
  	        Gr = correlation_vector_1N(psi,"Cdagup","Cup";site_x = range_x, sites_y = range_y)
  	        Gr = Gr + correlation_vector_1N(psi,"Cdagdn","Cdn";site_x = range_x, sites_y = range_y)
  	
  	        for i in 1:length(range_y)
  	            println(io, range_x,"\t",range_y[i],"\t",real(Gr[i]),"\t",imag(Gr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	        end
  	    end
  	end
	end
  	
  """
  ##############################################################################################################################
  Szz
  """
	if Szz
  	open(string("Szz_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  	    for j=1:1
  	        range_x= ini*Ny+j
  	        range_y= collect((ini+1):(Nx-ini)).*Ny .+ j
  	        Sr = correlation_vector_1N(psi,"Sz","Sz";site_x = range_x, sites_y = range_y)
  	
  	        for i in 1:length(range_y)
  	            println(io, range_x,"\t",range_y[i],"\t",real(Sr[i]),"\t",imag(Sr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	        end
  	    end
  	end
	end
  	
  """
  ##############################################################################################################################
  Szz_M
  """
	if Szz_M
  	range_x= collect(1:N)
		#range_x= collect((ini*Ny+1):(Nx-ini)*Ny)
  	Sr = correlation_matrix(psi,"Sz","Sz";site_range = range_x)
  	open(string("Szz_M_",sys,".dat"), "w") do io
  	  println(io, "r1","\t","r2","\t","Real","\t","Imag")
  		for i in 1:size(range_x,1)
  		  for j in 1:size(range_x,1)
  		    println(io, range_x[i],"\t",range_x[j],"\t",real(Sr[i,j]))
  		  end
  		end
  	end
	end

  """
  ##############################################################################################################################
  Sxy_M
  """
	if Sxy_M
  	range_x= collect(1:N)
		#range_x= collect((ini*Ny+1):(Nx-ini)*Ny)
  	Sr = correlation_matrix(psi,"S+","S-";site_range = range_x)
  	open(string("Sxy_M_",sys,".dat"), "w") do io
  	  println(io, "r1","\t","r2","\t","Real","\t","Imag")
  		for i in 1:size(range_x,1)
  		  for j in 1:size(range_x,1)
  		    println(io, range_x[i],"\t",range_x[j],"\t",real(Sr[i,j]))
  		  end
  		end
  	end
	end

  """
  ##############################################################################################################################
  Sxy
  """
	if Sxy
  	open(string("Sxy_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  	    for j=1:1
  	        range_x= ini*Ny+j
  	        range_y= collect((ini+1):(Nx-ini)).*Ny .+ j
  	        Sr = correlation_vector_1N(psi,"S+","S-";site_x = range_x, sites_y = range_y)
  	
  	        for i in 1:length(range_y)
  	            println(io, range_x,"\t",range_y[i],"\t",real(Sr[i]),"\t",imag(Sr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	        end
  	    end
  	end
	end
  	
  """
  ##############################################################################################################################
  Dr
  """
	if Dr
  	open(string("Dr_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Real","\t","Imag")
  	    for j=1:1
  	        range_x= ini*Ny+j
  	        range_y= collect((ini+1):(Nx-ini)).*Ny .+ j
  	        Dr = correlation_vector_1N(psi,"Ntot","Ntot";site_x = range_x, sites_y = range_y)
  	
  	        for i in 1:length(range_y)
  	            println(io, range_x,"\t",range_y[i],"\t",real(Dr[i]),"\t",imag(Dr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	        end
  	    end
  	end
	end

  """
  ##############################################################################################################################
  Paa
  """
	if Paa
  	open(string("Paa_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for j=1:1
  	        range_x= ini*Ny+j
  	        range_y= collect((ini+2):(Nx-ini)).*Ny .+ j
  	        bl=Ny
  	        Pr1 = correlation_4p(psi,"Cdagup","Cdagdn","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr2 = correlation_4p(psi,"Cdagup","Cdagdn","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr3 = correlation_4p(psi,"Cdagdn","Cdagup","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr4 = correlation_4p(psi,"Cdagdn","Cdagup","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr = Pr1 - Pr2 + Pr3 - Pr4
  	    
  	       for i in 1:length(range_y)
  	           println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
  	       end
  	    end
  	end
	end
  	
  """
  ##############################################################################################################################
  Pbb_shift
  """
	if Pbb_shift
  	open(string("Pbb_",sys,"_shift.dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for j=1:1
  	        range_x= ini*Ny+j
  	        range_y= collect((ini+1):(Nx-ini)).*Ny .+ j .+1
  	        bl=1
  	        Pr1 = correlation_4p(psi,"Cdagup","Cdagdn","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr2 = correlation_4p(psi,"Cdagup","Cdagdn","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr3 = correlation_4p(psi,"Cdagdn","Cdagup","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr4 = correlation_4p(psi,"Cdagdn","Cdagup","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr = Pr1 - Pr2 + Pr3 - Pr4
  	        
  	        for i in 1:length(range_y)
  	            println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
  	        end
  	    end
  	end
	end
  	
  """
  ##############################################################################################################################
  Pbb
  """
	if Pbb
  	open(string("Pbb_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for j=1:1
  	        range_x= ini*Ny+j
  	        range_y= collect((ini+1):(Nx-ini)).*Ny .+ j
  	        bl=1
  	        Pr1 = correlation_4p(psi,"Cdagup","Cdagdn","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr2 = correlation_4p(psi,"Cdagup","Cdagdn","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr3 = correlation_4p(psi,"Cdagdn","Cdagup","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr4 = correlation_4p(psi,"Cdagdn","Cdagup","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	        Pr = Pr1 - Pr2 + Pr3 - Pr4
  	        
  	        for i in 1:length(range_y)
  	            println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
  	        end
  	    end
  	end
	end
  	
  """
  ##############################################################################################################################
  Pcc
  """
	if Pcc
		open(string("Pcc_",sys,".dat"), "w") do io
		    println(io, "r0","\t","r1","\t","Norm")
		    for j=1:1
		        range_x= ini*Ny+j
		        range_y= collect((ini+2):(Nx-ini)).*Ny .+ j
		        bl1=1+Ny
		        bl2=1+Ny
		        Pr1 = correlation_4p(psi,"Cdagup","Cdagdn","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		        Pr2 = correlation_4p(psi,"Cdagup","Cdagdn","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		        Pr3 = correlation_4p(psi,"Cdagdn","Cdagup","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		        Pr4 = correlation_4p(psi,"Cdagdn","Cdagup","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		        Pr = Pr1 - Pr2 + Pr3 - Pr4
		        
		        for i in 1:length(range_y)
		            println(io, "(",range_x,",",range_x+bl1,")","\t","(",range_y[i],",",range_y[i]+bl2,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
		        end
		    end
		end
	end

  """
  ##############################################################################################################################
  Pba
  """
	if Pba
  	open(string("Pba_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for j=1:1
  	        range_x= ini*Ny+j
  	        range_y= collect((ini+1):(Nx-ini)).*Ny .+ j
  	        bl1=1
  	        bl2=Ny
  	        Pr1 = correlation_4p(psi,"Cdagup","Cdagdn","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
  	        Pr2 = correlation_4p(psi,"Cdagup","Cdagdn","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
  	        Pr3 = correlation_4p(psi,"Cdagdn","Cdagup","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
  	        Pr4 = correlation_4p(psi,"Cdagdn","Cdagup","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
  	        Pr = Pr1 - Pr2 + Pr3 - Pr4
  	        
  	        for i in 1:length(range_y)
  	            println(io, "(",range_x,",",range_x+bl1,")","\t","(",range_y[i],",",range_y[i]+bl2,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
  	        end
  	    end
  	end
	end
  	
  """
  ##############################################################################################################################
  Pbc
  """
	if Pbc
		open(string("Pbc_",sys,".dat"), "w") do io
		    println(io, "r0","\t","r1","\t","Norm")
		    for j=1:1
		        range_x= ini*Ny+j
		        range_y= collect((ini+1):(Nx-ini)).*Ny .+ j
		        bl1=1
		        bl2=1+Ny
		        Pr1 = correlation_4p(psi,"Cdagup","Cdagdn","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		        Pr2 = correlation_4p(psi,"Cdagup","Cdagdn","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		        Pr3 = correlation_4p(psi,"Cdagdn","Cdagup","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		        Pr4 = correlation_4p(psi,"Cdagdn","Cdagup","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		        Pr = Pr1 - Pr2 + Pr3 - Pr4
		        
		        for i in 1:length(range_y)
		            println(io, "(",range_x,",",range_x+bl1,")","\t","(",range_y[i],",",range_y[i]+bl2,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
		        end
		    end
		end
	end

  """
  ##############################################################################################################################
  Ptud_bb
  """
  
	if Ptud_bb 
  	range_x= ini*Ny+1
  	range_y= collect((ini+1):(Nx-ini)).*Ny .+ 1
  	bl=1
  	Pr1 = -correlation_4p(psi,"Cdagup","Cdagdn","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	Pr2 = -correlation_4p(psi,"Cdagup","Cdagdn","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	Pr3 = -correlation_4p(psi,"Cdagdn","Cdagup","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	Pr4 = -correlation_4p(psi,"Cdagdn","Cdagup","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	Pr = Pr1 + Pr2 + Pr3 + Pr4
  	
  	open(string("Ptud_bb_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for i in 1:length(range_y)
  	        println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
  	    end
  	end
	end
  
  """
  ##############################################################################################################################
  Ptuu_bb_shift
  """

	if Ptuu_bb_shift 
  	bl1=1
  	range_x= ini*Ny+1
		if shift < Ny - 1
  		range_y= collect((ini+1):(Nx-ini)).*Ny .+ 1 .+shift
  		bl2=1
  		Pr = -correlation_4p(psi,"Cdagup","Cdagup","Cup","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		elseif shift == Ny - 1
  		range_y= collect((ini+1):(Nx-ini)).*Ny .+ 1
  		bl2=Ny-1
  		Pr = correlation_4p(psi,"Cdagup","Cdagup","Cup","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		end
  	
  	open(string("Ptuu_bb_shift",shift,"_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for i in 1:length(range_y)
  	        println(io, "(",range_x,",",range_x+bl1,")","\t","(",range_y[i],",",range_y[i]+bl2,")","\t",real(Pr[i]),"\t",imag(Pr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	    end
  	end
	end
  
  """
  ##############################################################################################################################
  Ptuu_bb
  """

	if Ptuu_bb 
  	range_x= ini*Ny+1
  	range_y= collect((ini+1):(Nx-ini)).*Ny .+ 1
  	bl=1
  	Pr = -correlation_4p(psi,"Cdagup","Cdagup","Cup","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	
  	open(string("Ptuu_bb_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for i in 1:length(range_y)
  	        println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i]),"\t",imag(Pr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	    end
  	end
	end
  
  """
  ##############################################################################################################################
  Ptdd_bb_shift
  """

	if Ptdd_bb_shift 
  	bl1=1
  	range_x= ini*Ny+1
		if shift < Ny - 1
  		range_y= collect((ini+1):(Nx-ini)).*Ny .+ 1 .+shift
  		bl2=1
  		Pr = -correlation_4p(psi,"Cdagdn","Cdagdn","Cdn","Cdn";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		elseif shift == Ny - 1
  		range_y= collect((ini+1):(Nx-ini)).*Ny .+ 1
  		bl2=Ny-1
  		Pr = correlation_4p(psi,"Cdagdn","Cdagdn","Cdn","Cdn";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		end
  	
  	open(string("Ptdd_bb_shift",shift,"_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for i in 1:length(range_y)
  	        println(io, "(",range_x,",",range_x+bl1,")","\t","(",range_y[i],",",range_y[i]+bl2,")","\t",real(Pr[i]),"\t",imag(Pr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	    end
  	end
	end
  """
  ##############################################################################################################################
  Ptdd_bb
  """

	if Ptdd_bb 
  	range_x= ini*Ny+1
  	range_y= collect((ini+1):(Nx-ini)).*Ny .+ 1
  	bl=1
  	Pr = -correlation_4p(psi,"Cdagdn","Cdagdn","Cdn","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	
  	open(string("Ptdd_bb_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for i in 1:length(range_y)
  	        println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i]),"\t",imag(Pr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	    end
  	end
	end
  
  """
  ##############################################################################################################################
  Ptdd_bb_y
  """

	if Ptdd_bb_y 
  	range_x= ini*Ny+1
  	range_y= collect((ini+8)*Ny .+ (1:Ny-1))
  	bl=1
  	Pr = -correlation_4p(psi,"Cdagdn","Cdagdn","Cdn","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	
  	open(string("Ptdd_bb_y_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for i in 1:length(range_y)
  	        println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i]),"\t",imag(Pr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	    end
  	end
	end

  """
  ##############################################################################################################################
  Ptuu_bb_y
  """

	if Ptuu_bb_y 
  	range_x= ini*Ny+1
  	range_y= collect((ini+8)*Ny .+ (1:Ny-1))
  	bl=1
  	Pr = -correlation_4p(psi,"Cdagup","Cdagup","Cup","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	
  	open(string("Ptuu_bb_y_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for i in 1:length(range_y)
  	        println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i]),"\t",imag(Pr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	    end
  	end
	end
  
	"""
	##############################################################################################################################
	Ptud_ba
	"""
	
	if Ptud_ba
		range_x= ini*Ny+1
		range_y= collect((ini+1):(Nx-ini)).*Ny .+ 1
		bl1=1
		bl2=Ny
		Pr1 = -correlation_4p(psi,"Cdagup","Cdagdn","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		Pr2 = -correlation_4p(psi,"Cdagup","Cdagdn","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		Pr3 = -correlation_4p(psi,"Cdagdn","Cdagup","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		Pr4 = -correlation_4p(psi,"Cdagdn","Cdagup","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		Pr = Pr1 + Pr2 + Pr3 + Pr4
		
		open(string("Ptud_ba_",sys,".dat"), "w") do io
		    println(io, "r0","\t","r1","\t","Norm")
		    for i in 1:length(range_y)
		        println(io, "(",range_x,",",range_x+bl1,")","\t","(",range_y[i],",",range_y[i]+bl2,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
		    end
		end
	end

	"""
	##############################################################################################################################
	Ptuu_ba
	"""
	
	if Ptuu_ba
		range_x= ini*Ny+1
		range_y= collect((ini+1):(Nx-ini)).*Ny .+ 1
		bl1=1
		bl2=Ny
		Pr = -correlation_4p(psi,"Cdagup","Cdagup","Cup","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		
		open(string("Ptuu_ba_",sys,".dat"), "w") do io
		    println(io, "r0","\t","r1","\t","Norm")
		    for i in 1:length(range_y)
		        println(io, "(",range_x,",",range_x+bl1,")","\t","(",range_y[i],",",range_y[i]+bl2,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
		    end
		end
	end

	"""
	##############################################################################################################################
	Ptuu_ac
	"""
	
	if Ptuu_ac
		range_x= ini*Ny+1
		range_y= collect((ini+2):(Nx-ini)).*Ny .+ 1
		bl1=Ny
		bl2=1+Ny
		Pr = -correlation_4p(psi,"Cdagup","Cdagup","Cup","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		
		open(string("Ptuu_ac_",sys,".dat"), "w") do io
		    println(io, "r0","\t","r1","\t","Norm")
		    for i in 1:length(range_y)
		        println(io, "(",range_x,",",range_x+bl1,")","\t","(",range_y[i],",",range_y[i]+bl2,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
		    end
		end
	end

	"""
	##############################################################################################################################
	Ptuu_bc
	"""
	
	if Ptuu_bc
		range_x= ini*Ny+1
		range_y= collect((ini+1):(Nx-ini)).*Ny .+ 1
		bl1=1
		bl2=1+Ny
		Pr = -correlation_4p(psi,"Cdagup","Cdagup","Cup","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		
		open(string("Ptuu_bc_",sys,".dat"), "w") do io
		    println(io, "r0","\t","r1","\t","Norm")
		    for i in 1:length(range_y)
		        println(io, "(",range_x,",",range_x+bl1,")","\t","(",range_y[i],",",range_y[i]+bl2,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
		    end
		end
	end

	"""
	##############################################################################################################################
	Ptuu_bc
	"""
#	
#	if Ptuu_bc
#		range_x= ini*Ny+1
#		range_y= collect((ini+2):(Nx-ini)).*Ny  .- Ny
#		bl1=1
#		bl2=1+Ny
#		Pr = correlation_4p(psi,"Cdagup","Cdagup","Cup","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
#		
#		open(string("Ptuu_bc_",sys,".dat"), "w") do io
#		    println(io, "r0","\t","r1","\t","Norm")
#		    for i in 1:length(range_y)
#		        println(io, "(",range_x,",",range_x+bl1,")","\t","(",range_y[i],",",range_y[i]+bl2,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
#		    end
#		end
#	end

	"""
	##############################################################################################################################
	Ptud_ac
	"""
	
	if Ptud_ac
		range_x= ini*Ny+1
		range_y= collect((ini+2):(Nx-ini)).*Ny .+ 1
		bl1=Ny
		bl2=1+Ny
		Pr1 = -correlation_4p(psi,"Cdagup","Cdagdn","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		Pr2 = -correlation_4p(psi,"Cdagup","Cdagdn","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		Pr3 = -correlation_4p(psi,"Cdagdn","Cdagup","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		Pr4 = -correlation_4p(psi,"Cdagdn","Cdagup","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		Pr = Pr1 + Pr2 + Pr3 + Pr4
		
		open(string("Ptud_ac_",sys,".dat"), "w") do io
		    println(io, "r0","\t","r1","\t","Norm")
		    for i in 1:length(range_y)
		        println(io, "(",range_x,",",range_x+bl1,")","\t","(",range_y[i],",",range_y[i]+bl2,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
		    end
		end
	end

	"""
	##############################################################################################################################
	Ptud_bc
	"""
	
	if Ptud_bc
		range_x= ini*Ny+1
		range_y= collect((ini+1):(Nx-ini)).*Ny .+ 1
		bl1=1
		bl2=1+Ny
		Pr1 = -correlation_4p(psi,"Cdagup","Cdagdn","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		Pr2 = -correlation_4p(psi,"Cdagup","Cdagdn","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		Pr3 = -correlation_4p(psi,"Cdagdn","Cdagup","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		Pr4 = -correlation_4p(psi,"Cdagdn","Cdagup","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl1, bond2=bl2)
		Pr = Pr1 + Pr2 + Pr3 + Pr4
		
		open(string("Ptud_bc_",sys,".dat"), "w") do io
		    println(io, "r0","\t","r1","\t","Norm")
		    for i in 1:length(range_y)
		        println(io, "(",range_x,",",range_x+bl1,")","\t","(",range_y[i],",",range_y[i]+bl2,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
		    end
		end
	end

  """
  ##############################################################################################################################
  Ptud_cc
  """

	if Ptud_cc 
  	range_x= ini*Ny+1
  	range_y= collect((ini+2):(Nx-ini)).*Ny .+ 1
  	bl=Ny+1
  	Pr1 = -correlation_4p(psi,"Cdagup","Cdagdn","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	Pr2 = -correlation_4p(psi,"Cdagup","Cdagdn","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	Pr3 = -correlation_4p(psi,"Cdagdn","Cdagup","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	Pr4 = -correlation_4p(psi,"Cdagdn","Cdagup","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	Pr = Pr1 + Pr2 + Pr3 + Pr4
  	
  	open(string("Ptud_cc_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for i in 1:length(range_y)
  	        println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
  	    end
  	end
	end
  
  """
  ##############################################################################################################################
  Ptud_aa
  """

	if Ptud_aa 
  	range_x= ini*Ny+1
  	range_y= collect((ini+2):(Nx-ini)).*Ny .+ 1
  	bl=Ny
  	Pr1 = -correlation_4p(psi,"Cdagup","Cdagdn","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	Pr2 = -correlation_4p(psi,"Cdagup","Cdagdn","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	Pr3 = -correlation_4p(psi,"Cdagdn","Cdagup","Cup","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	Pr4 = -correlation_4p(psi,"Cdagdn","Cdagup","Cdn","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	Pr = Pr1 + Pr2 + Pr3 + Pr4
  	
  	open(string("Ptud_aa_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for i in 1:length(range_y)
  	        println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i])/2,"\t",imag(Pr[i])/2,"\t",round((range_y[i]-range_x)/Ny))
  	    end
  	end
	end
  
  """
  ##############################################################################################################################
  Ptuu_cc
  """
	if Ptuu_cc 
  	range_x= ini*Ny+1
  	range_y= collect((ini+2):(Nx-ini)).*Ny .+ 1
  	bl=Ny+1
  	Pr = -correlation_4p(psi,"Cdagup","Cdagup","Cup","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	
  	open(string("Ptuu_cc_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for i in 1:length(range_y)
  	        println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i]),"\t",imag(Pr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	    end
  	end
	end
  """
  ##############################################################################################################################
  Ptuu_aa
  """
	if Ptuu_aa 
  	range_x= ini*Ny+1
  	range_y= collect((ini+2):(Nx-ini)).*Ny .+ 1
  	bl=Ny
  	Pr = -correlation_4p(psi,"Cdagup","Cdagup","Cup","Cup";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	
  	open(string("Ptuu_aa_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for i in 1:length(range_y)
  	        println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i]),"\t",imag(Pr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	    end
  	end
	end
  
  """
  ###############################################################################################################################
  #Ptdd_aa
  """
	if Ptdd_aa 
  	range_x= ini*Ny+1
  	range_y= collect((ini+2):(Nx-ini)).*Ny .+ 1
  	bl=Ny
  	Pr = -correlation_4p(psi,"Cdagdn","Cdagdn","Cdn","Cdn";site_x = range_x, sites_y = range_y, bond1=bl, bond2=bl)
  	
  	open(string("Ptdd_aa_",sys,".dat"), "w") do io
  	    println(io, "r0","\t","r1","\t","Norm")
  	    for i in 1:length(range_y)
  	        println(io, "(",range_x,",",range_x+bl,")","\t","(",range_y[i],",",range_y[i]+bl,")","\t",real(Pr[i]),"\t",imag(Pr[i]),"\t",round((range_y[i]-range_x)/Ny))
  	    end
  	end
	end
  
## Other measurements
  #Sz = expect(psi,"Sz";sites = 1:N)
  #open(string("Sz_",sys,".dat"), "w") do io
  #    for i in 1:Ny
  #    	for j in 1:Nx
  #        println(io, i,"\t",j,"\t",real(Sz[(j-1)*Ny+i]))
  #      end
  #    end
  #end
  
  #Cc = correlation_matrix(psi,"Ntot","Ntot")
  #open(string("Cc_",sys,".dat"), "w") do io
  #    for i in 1:N
  #    	for j in 1:N
  #        println(io, i,"\t",j,"\t",real(Cc[i,j]))
  #      end
  #    end
  #end
  
end
  
